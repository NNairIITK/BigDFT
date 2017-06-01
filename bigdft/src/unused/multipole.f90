
!!!!    !> Calculates the multipole moments for each atom 
!!!!    subroutine multipoles_from_density(iproc, nproc, at, lzd, smats, smatl, orbs, &
!!!!               npsidim, lphi, norbsPerType, collcom, collcom_sr, orthpar, &
!!!!               ovrlp, kernel, meth_overlap, multipoles_out)
!!!!      use module_base
!!!!      use module_types
!!!!      use sparsematrix_base, only: sparsematrix_malloc0, SPARSE_FULL, assignment(=)
!!!!      use sparsematrix_init, only: matrixindex_in_compressed
!!!!      use orthonormalization, only: orthonormalizeLocalized
!!!!      use yaml_output
!!!!      use locreg_operations
!!!!      use bounds, only: geocode_buffers
!!!!      implicit none
!!!!
!!!!      ! Calling arguments
!!!!      integer,intent(in) :: iproc, nproc, npsidim, meth_overlap
!!!!      type(atoms_data),intent(in) :: at
!!!!      type(local_zone_descriptors),intent(in) :: lzd
!!!!      type(sparse_matrix),intent(inout) :: smats, smatl
!!!!      type(orbitals_data),intent(in) :: orbs
!!!!      real(kind=8),dimension(npsidim),intent(in) :: lphi
!!!!      integer,dimension(at%astruct%ntypes),intent(in) :: norbsPerType
!!!!      type(comms_linear),intent(in) :: collcom
!!!!      type(comms_linear),intent(in) :: collcom_sr
!!!!      type(orthon_data),intent(in) :: orthpar
!!!!      type(matrices),intent(in) :: kernel
!!!!      type(matrices),intent(inout) :: ovrlp !in principle also intent(in)
!!!!      real(kind=8),dimension(-lmax:lmax,0:lmax,1:at%astruct%nat),intent(out),optional :: multipoles_out
!!!!
!!!!      ! Local variables
!!!!      integer :: ist, istr, iorb, iiorb, ilr, ii, natp, isat, nr, jproc, iat, n, norb_get, istr_get
!!!!      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
!!!!      integer :: ii1, ii2, ii3, jjorb, i, itype
!!!!      real(kind=8),dimension(:),allocatable :: psir
!!!!      real(kind=8),dimension(:),pointer :: phit_c, phit_f
!!!!      type(workarr_sumrho) :: w
!!!!      integer,dimension(:),allocatable :: nat_par, norb_list
!!!!      real(kind=8),dimension(:),allocatable :: psir_get, locrad, rmax
!!!!      real(kind=8),dimension(:,:),allocatable :: locregcenter, psir_get_fake
!!!!      real(kind=8),dimension(:,:,:,:),allocatable :: phi
!!!!      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
!!!!      integer,dimension(:,:),allocatable :: comms
!!!!      logical :: can_use_transposed, arr_allocated
!!!!      real(kind=8) :: ddot, x, y, z, tt, rnorm, factor, max_error, q!, get_normalization, get_test_factor
!!!!      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
!!!!      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
!!!!      real(kind=8),dimension(:,:),allocatable :: weight_centers
!!!!      integer,dimension(:),allocatable :: n1i, n2i, n3i, ns1i, ns2i, ns3i
!!!!      !real(kind=8),dimension(-lmax:lmax,0:lmax,at%astruct%nat) :: multipoles
!!!!      real(kind=8),dimension(:,:,:),allocatable :: multipoles
!!!!      real(kind=8) :: factor_normalization, hxh, hyh, hzh, weight
!!!!      character(len=20) :: atomname
!!!!      real(kind=8),dimension(-lmax:lmax,0:lmax) :: norm
!!!!      integer :: nl1, nl2, nl3
!!!!      !real(kind=8),parameter :: rmax=5.d0
!!!!      !testarr = 0.d0
!!!!
!!!!      call f_routine(id='multipoles_from_density')
!!!!
!!!!      if (iproc==0) call yaml_comment('Atomic multipole analysis, old approach',hfill='=')
!!!!
!!!!      multipoles = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.at%astruct%nat/),id='multipoles')
!!!!
!!!!
!!!!      call unitary_test()
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_mapping_open('Multipole analysis')
!!!!      end if
!!!!
!!!!      ! Orthogonalize the support functions
!!!!      can_use_transposed = .false.
!!!!      methTransformOverlap = 1020
!!!!      phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
!!!!      phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
!!!!      phi_ortho = f_malloc(npsidim,id='phi_ortho')
!!!!      call vcopy(npsidim, lphi(1), 1, phi_ortho(1), 1)
!!!!      if (iproc==0) then
!!!!          call yaml_map('Orthonormalizing support functions',.true.)
!!!!      end if
!!!!      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
!!!!           1.d-8, npsidim, orbs, lzd, &
!!!!           smats, smatl, collcom, orthpar, &
!!!!           phi_ortho, phit_c, phit_f, &
!!!!           can_use_transposed)
!!!!      call f_free_ptr(phit_c)
!!!!      call f_free_ptr(phit_f)
!!!!
!!!!
!!!!      ! Transform the support functions to real space
!!!!      psir = f_malloc(max(collcom_sr%ndimpsi_c,1),id='psir')
!!!!      ist=1
!!!!      istr=1
!!!!      do iorb=1,orbs%norbp
!!!!          iiorb=orbs%isorb+iorb
!!!!          ilr=orbs%inwhichlocreg(iiorb)
!!!!          call initialize_work_arrays_sumrho(1,[lzd%Llr(ilr)],.true.,w)
!!!!          call daub_to_isf(lzd%Llr(ilr), w, phi_ortho(ist), psir(istr))
!!!!          call deallocate_work_arrays_sumrho(w)
!!!!          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
!!!!          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!      end do
!!!!      do i=1,ist-1
!!!!          write(445,*) i, phi_ortho(i)
!!!!      end do
!!!!      do i=1,istr-1
!!!!          write(444,*) i, psir(i)
!!!!      end do
!!!!      call f_free(phi_ortho)
!!!!      if(istr/=collcom_sr%ndimpsi_c+1) then
!!!!          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
!!!!          stop
!!!!      end if
!!!!
!!!!
!!!!      !! NEW: CALCULATE THE WEIGHT CENTER OF EACH SUPPORT FUNCTION ############################
!!!!      !weight_centers = f_malloc0((/3,orbs%norb/),id='weight_centers')
!!!!      !hxh = 0.5d0*lzd%hgrids(1)
!!!!      !hyh = 0.5d0*lzd%hgrids(2)
!!!!      !hzh = 0.5d0*lzd%hgrids(3)
!!!!
!!!!      !istr = 1
!!!!      !do iorb=1,orbs%norbp
!!!!      !    iiorb=orbs%isorb+iorb
!!!!      !    ilr=orbs%inwhichlocreg(iiorb)
!!!!      !    call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
!!!!      !    !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
!!!!      !    !com(1:3,iorb) = 0.d0
!!!!      !    weight = 0.d0
!!!!      !    do i3=1,lzd%llr(ilr)%d%n3i
!!!!      !        ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
!!!!      !        z = ii3*hzh
!!!!      !        do i2=1,lzd%llr(ilr)%d%n2i
!!!!      !            ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
!!!!      !            y = ii2*hyh
!!!!      !            do i1=1,lzd%llr(ilr)%d%n1i
!!!!      !                ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
!!!!      !                x = ii1*hxh
!!!!      !                tt = psir(istr)**2
!!!!      !                weight_centers(1,iiorb) = weight_centers(1,iiorb) + x*tt
!!!!      !                weight_centers(2,iiorb) = weight_centers(2,iiorb) + y*tt
!!!!      !                weight_centers(3,iiorb) = weight_centers(3,iiorb) + z*tt
!!!!      !                weight = weight + tt
!!!!      !                istr = istr + 1
!!!!      !            end do
!!!!      !        end do
!!!!      !    end do
!!!!      !    !call yaml_map('weight',weight)
!!!!      !    weight_centers(1:3,iorb) = weight_centers(1:3,iorb)/weight
!!!!      !end do
!!!!      !call mpiallred(weight_centers, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!!!      !if (iproc==0) then
!!!!      !    do iorb=1,orbs%norb
!!!!      !        write(*,'(a,i4,3es13.4)') 'iorb, weight_centers(1:3,iorb)', iorb, weight_centers(1:3,iorb)
!!!!      !    end do
!!!!      !end if
!!!!      !call f_free(weight_centers)
!!!!      !! ######################################################################################
!!!!
!!!!
!!!!      ! Switch to a partition over the atoms
!!!!      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
!!!!      ii = at%astruct%nat/nproc
!!!!      nat_par(0:nproc-1) = ii
!!!!      ii = at%astruct%nat-nproc*ii
!!!!      nat_par(0:ii-1) = nat_par(0:ii-1) + 1
!!!!      if (sum(nat_par(:))/=at%astruct%nat) then
!!!!          call f_err_throw('wrong partition of the atoms',err_name='BIGDFT_RUNTIME_ERROR')
!!!!      end if
!!!!      natp = nat_par(iproc)
!!!!      if (iproc==0) then
!!!!          isat = 0
!!!!      else
!!!!          isat = sum(nat_par(0:iproc-1))
!!!!      end if
!!!!      call f_free(nat_par)
!!!!      comms = f_malloc((/4,orbs%norb/),id='comms')
!!!!      !write(*,'(a,i5,3x,13i4)') 'iproc, orbs%onwhichatom', iproc, orbs%onwhichatom
!!!!      nr = 0
!!!!      norb_get = 0
!!!!      istr = 1
!!!!      istr_get = 1
!!!!      do iat=isat+1,isat+natp
!!!!          do jproc=0,nproc-1
!!!!              istr = 0
!!!!              do iorb=1,orbs%norb_par(jproc,0)
!!!!                  iiorb = iorb + orbs%isorb_par(jproc)
!!!!                  ilr = orbs%inwhichlocreg(iiorb)
!!!!                  iiat = orbs%onwhichatom(iiorb)
!!!!                  n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!                  !write(*,'(a,5i8)') 'iproc, jproc, iiorb, iat', iproc, jproc, iiorb, iat
!!!!                  !if (iat>=isat+1 .and. iat<=isat+natp) then
!!!!                  if (iiat==iat) then
!!!!                      norb_get = norb_get + 1
!!!!                      comms(1,norb_get) = jproc
!!!!                      comms(2,norb_get) = n
!!!!                      comms(3,norb_get) = istr
!!!!                      comms(4,norb_get) = istr_get
!!!!                      nr = nr + n
!!!!                      !write(*,'(a,5i8)') 'iproc, jproc, n, istr, istr_get', iproc, jproc, n, istr, istr_get
!!!!                      istr = istr + n
!!!!                      istr_get = istr_get + n
!!!!                  else
!!!!                      istr = istr + n
!!!!                  end if
!!!!              end do
!!!!          end do
!!!!      end do
!!!!      !write(*,*) 'iproc, nr', iproc, nr
!!!!      !do iorb=1,norb_get
!!!!      !    write(*,'(a,2i7,4i9)') 'iproc, iorb, comm',iproc, iorb, comms(:,iorb)
!!!!      !end do
!!!!      psir_get = f_malloc(nr,id='psir_get')
!!!!      if (nproc>1) then
!!!!          window = mpiwindow(size(psir), psir(1), bigdft_mpi%mpi_comm)
!!!!          do iorb=1,norb_get
!!!!              jproc = comms(1,iorb)
!!!!              n = comms(2,iorb)
!!!!              ioffset = comms(3,iorb)
!!!!              istr = comms(4,iorb)
!!!!              !write(*,'(5(a,i0))') 'task ',iproc,' gets ',n,' elements at position ', &
!!!!              !                     istr,' from position ',ioffset+1,' on task ',jproc
!!!!              call mpiget(psir_get(istr), n, jproc, int(ioffset,kind=mpi_address_kind), window)
!!!!          end do
!!!!          call mpi_fenceandfree(window)
!!!!      else
!!!!          do iorb=1,norb_get
!!!!              n = comms(2,iorb)
!!!!              ioffset = comms(3,iorb)
!!!!              istr = comms(4,iorb)
!!!!              call vcopy(n, psir(ioffset+1), 1, psir_get(istr), 1)
!!!!          end do
!!!!      end if
!!!!      call f_free(psir)
!!!!      call f_free(comms)
!!!!      istr = 1
!!!!      !write(*,*) 'iproc, isat, natp', iproc, isat, natp
!!!!      do iorb=1,orbs%norb
!!!!          ilr = orbs%inwhichlocreg(iorb)
!!!!          iat = orbs%onwhichatom(iorb)
!!!!          if (iat>=isat+1 .and. iat<=isat+natp) then
!!!!              n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!              !write(*,'(a,4i8,es16.6)') 'iproc, iorb, n, istr, ddot', &
!!!!              !    iproc, iorb, n, istr, ddot(n, psir_get(istr), 1, psir_get(istr), 1)
!!!!              !testarr(2,iorb) = ddot(n, psir_get(istr), 1, psir_get(istr), 1)
!!!!              istr = istr + n
!!!!          end if
!!!!      end do
!!!!      !do iorb=1,size(psir_get)
!!!!      !    write(300+iproc,'(a,i7,es16.7)') 'i, psir_get(i)', iorb, psir_get(iorb)
!!!!      !end do
!!!!      !call mpiallred(testarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!!!      !if (iproc==0) then
!!!!      !    do iorb=1,orbs%norb
!!!!      !        write(*,*) 'DIFF, iorb, val', iorb, abs(testarr(1,iorb)-testarr(2,iorb))
!!!!      !    end do
!!!!      !end if
!!!!
!!!!
!!!!      ! Calculate the kernel for orthormal support functions
!!!!      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='kernel_ortho')
!!!!      call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, meth_overlap, smats, smatl, &
!!!!           ovrlp, kernel, kernel_ortho)
!!!!
!!!!
!!!!      ! Apply the spherical harmonics to the suport functions.
!!!!      ! For the support functions on atom A we only need to apply 
!!!!      ! the spherical harmonics centered as well on atom A.
!!!!
!!!!
!!!!
!!!!      !lmax = 1
!!!!      n1i = f_malloc(orbs%norb,id='n1i')
!!!!      n2i = f_malloc(orbs%norb,id='n2i')
!!!!      n3i = f_malloc(orbs%norb,id='n3i')
!!!!      ns1i = f_malloc(orbs%norb,id='ns1i')
!!!!      ns2i = f_malloc(orbs%norb,id='ns2i')
!!!!      ns3i = f_malloc(orbs%norb,id='ns3i')
!!!!      locrad = f_malloc(orbs%norb,id='locrad')
!!!!      locregcenter = f_malloc((/3,orbs%norb/),id='locregcenter')
!!!!      do iorb=1,orbs%norb
!!!!          n1i(iorb) = lzd%llr(iorb)%d%n1i
!!!!          n2i(iorb) = lzd%llr(iorb)%d%n2i
!!!!          n3i(iorb) = lzd%llr(iorb)%d%n3i
!!!!          ns1i(iorb) = lzd%llr(iorb)%nsi1
!!!!          ns2i(iorb) = lzd%llr(iorb)%nsi2
!!!!          ns3i(iorb) = lzd%llr(iorb)%nsi3
!!!!          locrad(iorb) = lzd%llr(iorb)%locrad
!!!!          locregcenter(1:3,iorb) = lzd%llr(iorb)%locregcenter(1:3)
!!!!      end do
!!!!      rmax = f_malloc(at%astruct%nat,id='rmax')
!!!!      call multipole_analysis_core(iproc, nproc, natp, isat, at%astruct%nat, at%astruct%ntypes, orbs%norb, &
!!!!           0, at%astruct%iatype, norbsPerType, orbs%inwhichlocreg, orbs%onwhichatom, &
!!!!           n1i, n2i, n3i, ns1i, ns2i, ns3i, locrad, lzd%hgrids, locregcenter, &
!!!!           nr, psir_get, psir_get, &
!!!!           smatl%nvctr, kernel_ortho, 1, multipoles, rmax, 101, smatl)!, matrixindex)
!!!!      call f_free(psir_get_fake)
!!!!      call f_free(n1i)
!!!!      call f_free(n2i)
!!!!      call f_free(n3i)
!!!!      call f_free(ns1i)
!!!!      call f_free(ns2i)
!!!!      call f_free(ns3i)
!!!!      call f_free(locrad)
!!!!      call f_free(locregcenter)
!!!!
!!!!      ! The monopole term should be the net charge, i.e. add the positive atomic charges
!!!!      do iat=1,at%astruct%nat
!!!!          itype = at%astruct%iatype(iat)
!!!!          q = real(at%nelpsp(itype),kind=8)
!!!!          multipoles(0,0,iat) = multipoles(0,0,iat) + q
!!!!      end do
!!!!
!!!!      if (iproc==0) then
!!!!          !!call write_multipoles(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%atomnames, &
!!!!          !!     multipoles, rmax, lzd%hgrids, without_normalization=.false.)
!!!!          call write_multipoles_new(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%atomnames, &
!!!!               at%astruct%rxyz, at%astruct%units, multipoles)
!!!!      end if
!!!!      call f_free(rmax)
!!!!      call f_free(psir_get)
!!!!
!!!!      call f_free(kernel_ortho)
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_mapping_close()
!!!!      end if
!!!!
!!!!
!!!!      if (present(multipoles_out)) then
!!!!          call f_memcpy(src=multipoles,dest=multipoles_out)
!!!!      end if
!!!!      call f_free(multipoles)
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_comment('Atomic multipole analysis done',hfill='=')
!!!!      end if
!!!!
!!!!      call f_release_routine()
!!!!
!!!!
!!!!
!!!!      contains
!!!!
!!!!
!!!!        subroutine unitary_test()
!!!!          implicit none
!!!!          real(kind=8),dimension(1) :: rmax
!!!!          integer,parameter :: n1i=101, n2i=81, n3i=91
!!!!          integer,parameter :: nsi1=0, nsi2=10, nsi3=20
!!!!          !integer,parameter :: n1i=100, n2i=100, n3i=100
!!!!          !integer,parameter :: nsi1=0, nsi2=0, nsi3=0
!!!!          real(kind=8),dimension(3) :: locregcenter
!!!!          integer :: nr
!!!!          real(kind=8) :: factor_normalization, r2, r
!!!!
!!!!          if (iproc==0) then
!!!!              call yaml_mapping_open('Unitary test for multipoles')
!!!!          end if
!!!!
!!!!          locregcenter(1) = (ceiling(real(n1i,kind=8)/2.d0)+nsi1-14-1)*0.5d0*lzd%hgrids(1)
!!!!          locregcenter(2) = (ceiling(real(n2i,kind=8)/2.d0)+nsi2-14-1)*0.5d0*lzd%hgrids(2)
!!!!          locregcenter(3) = (ceiling(real(n3i,kind=8)/2.d0)+nsi3-14-1)*0.5d0*lzd%hgrids(3)
!!!!
!!!!          !psir_get_fake = f_malloc0((/lzd%llr(1)%d%n1i*lzd%llr(1)%d%n2i*lzd%llr(1)%d%n3i,2/),id='psir_get_fake')
!!!!          nr = n1i*n2i*n3i
!!!!          psir_get_fake = f_malloc0((/nr,2/),id='psir_get_fake')
!!!!          psir_get_fake(:,2) = 1.d0
!!!!          rmax(1) = min(n1i*0.25d0*lzd%hgrids(1), &
!!!!                        n2i*0.25d0*lzd%hgrids(2), &
!!!!                        n3i*0.25d0*lzd%hgrids(3))
!!!!          ! Normalized to within a sphere of radius rmax
!!!!          factor_normalization = 3.d0/(4.d0*pi*rmax(1)**3)*0.5d0*lzd%hgrids(1)*0.5d0*lzd%hgrids(2)*0.5d0*lzd%hgrids(3)
!!!!          do i3=1,n3i!lzd%llr(1)%d%n3i
!!!!              ii3 = nsi3 + i3 - 14 - 1
!!!!              z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3)
!!!!              do i2=1,n2i!lzd%llr(1)%d%n2i
!!!!                  ii2 = nsi2 + i2 - 14 - 1
!!!!                  y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2)
!!!!                  do i1=1,n1i!lzd%llr(1)%d%n1i
!!!!                      ii1 = nsi1 + i1 - 14 - 1
!!!!                      x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1)
!!!!                      r2 = x**2+y**2+z**2
!!!!                      if (r2>rmax(1)**2) cycle
!!!!                      r = sqrt(r2)
!!!!                      r = max(0.5d0,r)
!!!!                      ind = (i3-1)*n2i*n1i + (i2-1)*n1i + i1
!!!!                      do l=0,lmax
!!!!                          do m=-l,l
!!!!                              factor = get_test_factor(l,m)*factor_normalization*sqrt(4.d0*pi*real(2*l+1,kind=8))
!!!!                              ! The minus sign is necessary to compensate the fact that the multipole calculation uses itself a minus sign (since it
!!!!                              ! assumes that the charge density (which is a negative quantity) is positive)
!!!!                              psir_get_fake(ind,1) = psir_get_fake(ind,1) - factor*solid_harmonic(-2, 1.0d0, l, m , x, y, z)
!!!!                          end do
!!!!                      end do
!!!!                  end do
!!!!              end do
!!!!          end do
!!!!          ! Only do it for one MPI task
!!!!          if (iproc==0) then
!!!!              call multipole_analysis_core(0, 1, 1, 0, 1, 1, 1, &
!!!!                   0, (/1/), (/1/), (/1/), (/1/), &
!!!!                   (/n1i/), (/n2i/), (/n3i/), &
!!!!                   (/nsi1/), (/nsi2/), (/nsi3/), rmax, &
!!!!                   lzd%hgrids, locregcenter, &
!!!!                   nr, psir_get_fake(:,1), psir_get_fake(:,2), &
!!!!                   1, (/1.d0/), 1, multipoles, rmax, 102, matrixindex=(/1/))
!!!!          end if
!!!!
!!!!          if (iproc==0) then
!!!!              call write_multipoles_new(1, 1, (/1/), (/'testatom'/), (/0.d0,0.d0,0.d0/), 'fake', &
!!!!                   multipoles, check_values_=.true.)
!!!!              call yaml_mapping_close()
!!!!          end if
!!!!
!!!!        end subroutine unitary_test
!!!!
!!!!
!!!!    end subroutine multipoles_from_density


    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, is1, ie1, is2, ie2, is3, ie3, iis3, iie3, hx, hy, hz, shift, pot)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3 !< parallelized box bounds
      integer,intent(in) :: iis3, iie3 !< non-parallelized box bounds (z direction)
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp
      real(dp),dimension(3) :: r
      real(kind=8),parameter :: buffer = 1.d0
      real(kind=8) :: xxs, xxe, yys, yye, zzs, zze

      !stop 'deprecated'

      xxs = real(is1-15,kind=8)*hx + shift(1) - buffer
      xxe = real(ie1-15,kind=8)*hx + shift(1) + buffer
      yys = real(is2-15,kind=8)*hy + shift(2) - buffer
      yye = real(ie2-15,kind=8)*hy + shift(2) + buffer
      zzs = real(iis3-15,kind=8)*hz + shift(3) - buffer
      zze = real(iie3-15,kind=8)*hz + shift(3) + buffer

      !!$omp parallel &
      !!$omp default(none) &
      !!$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !!$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !!$omp do
      do impl=1,ep%nmpl
          ! Only take this atom into account if it lies outside of the simulation box (plus some buffer).
          ! Otherwise it has already been taken into account by the potential generated via the Poisson Solver
          if ( ep%mpl(impl)%rxyz(1) <= xxs .or. ep%mpl(impl)%rxyz(1) >= xxe .or. &
               ep%mpl(impl)%rxyz(2) <= yys .or. ep%mpl(impl)%rxyz(2) >= yye .or. &
               ep%mpl(impl)%rxyz(3) <= zzs .or. ep%mpl(impl)%rxyz(3) >= zze ) then
              !write(*,*) ep%mpl(impl)%rxyz(1) <= xxs
              !write(*,*) ep%mpl(impl)%rxyz(1) >= xxe
              !write(*,*) ep%mpl(impl)%rxyz(2) <= yys
              !write(*,*) ep%mpl(impl)%rxyz(2) >= yye
              !write(*,*) ep%mpl(impl)%rxyz(3) <= zzs
              !write(*,*) ep%mpl(impl)%rxyz(3) >= zze
              !write(*,*) 'ok for impl',impl, ep%mpl(impl)%rxyz(1:3), xxs, xxe, yys, yye, zzs, zze
              do i3=is3,ie3
                  ii3 = i3 - 15
                  z = real(ii3,kind=8)*hz + shift(3)
                  !write(*,'(a,i7,2es16.7)') 'i3, z, ep%mpl(1)%rxyz(3)', i3, z, ep%mpl(1)%rxyz(3)
                  do i2=is2,ie2
                      ii2 = i2 - 15
                      y = real(ii2,kind=8)*hy + shift(2)
                      do i1=is1,ie1
                          ii1 = i1 - 15
                          x = real(ii1,kind=8)*hx + shift(1)
                          r(1) = ep%mpl(impl)%rxyz(1) - x
                          r(2) = ep%mpl(impl)%rxyz(2) - y
                          r(3) = ep%mpl(impl)%rxyz(3) - z 
                          rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                          ! To avoid floating point exception
                          rnrm2=max(rnrm2,1.d-6)
                          rnrm1 = sqrt(rnrm2)
                          rnrm3 = rnrm1*rnrm2
                          rnrm5 = rnrm3*rnrm2
                          mp = 0.0_dp
                          do l=0,lmax
                              if (associated(ep%mpl(impl)%qlm(l)%q)) then
                                  select case(l)
                                  case (0)
                                      mp = mp + calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                      !write(*,'(a,3es12.4,es16.8)') 'x, y, z, monopole', x, y, z, calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                  case (1)
                                      mp = mp + calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                      !write(*,*) 'dipole', calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                  case (2)
                                      mp = mp + calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                      !write(*,*) 'quadrupole', calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                  case (3)
                                      call f_err_throw('octupole not yet implemented', err_name='BIGDFT_RUNTIME_ERROR')
                                      !multipole_terms(l) = calc_octopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                  case default
                                      call f_err_throw('Wrong value of l', err_name='BIGDFT_RUNTIME_ERROR')
                                  end select
                              end if
                          end do
                          pot(i1,i2,i3) = pot(i1,i2,i3) + mp
                          !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
                      end do
                  end do
              end do
              !!$omp end do
              !!$omp end parallel
          end if
      end do

    end subroutine potential_from_multipoles



    subroutine multipole_analysis_core(iproc, nproc, natp, isat, nat, ntypes, norb, &
               r_exponent, iatype, norbsPerType, inwhichlocreg, onwhichatom, &
               n1i, n2i, n3i, nsi1, nsi2, nsi3, locrad, hgrids, locregcenter, &
               nr, psir1_get, psir2_get, &
               nvctr, matrix_compr, verbosity, &
               multipoles, rmax, get_index, smatl, matrixindex)
      use locreg_operations, only: workarr_sumrho      
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      !real(kind=8),intent(in) :: rmax
      integer,intent(in) :: natp, isat, nat, ntypes, norb, nr, nvctr, get_index, r_exponent
      integer,dimension(nat),intent(in) :: iatype
      integer,dimension(ntypes),intent(in) :: norbsPerType
      integer,dimension(norb),intent(in) :: inwhichlocreg, onwhichatom
      !logical,intent(in) :: with_rl !if true, the spherical harmonics are multiplied by r^l, otherwise not
      integer,dimension(norb),intent(in) :: n1i, n2i, n3i, nsi1, nsi2, nsi3
      real(kind=8),dimension(norb),intent(in) :: locrad
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(3,norb) :: locregcenter
      real(kind=8),dimension(nr),intent(in) :: psir1_get, psir2_get
      real(kind=8),dimension(nvctr),intent(in) :: matrix_compr
      integer,intent(in) ::verbosity
      real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(out) :: multipoles
      real(kind=8),dimension(nat),intent(out) :: rmax
      type(sparse_matrix),intent(in),optional :: smatl
      integer,dimension(norb,norb),intent(in),optional :: matrixindex

      ! Local variables
      integer,parameter :: INDEX_AUTOMATIC=101
      integer,parameter :: INDEX_MANUALLY=102
      integer :: ist, istr, iorb, iiorb, ilr, ii, jproc, iat, n, norb_get, istr_get, i
      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
      integer :: ii1, ii2, ii3, jjorb, iilr
      real(kind=8),dimension(:),allocatable :: psir
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      type(workarr_sumrho) :: w
      integer,dimension(:),allocatable :: nat_par, norb_list
      !real(kind=8),dimension(:),allocatable :: psir_get
      real(kind=8),dimension(:,:,:,:),allocatable :: phi1, phi2
      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
      integer,dimension(:,:),allocatable :: comms
      logical :: can_use_transposed, arr_allocated
      real(kind=8) :: ddot, x, y, z, tt, rnorm, rnorm_maxdev, tt2
      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
      real(kind=8) :: factor_normalization
      character(len=20) :: atomname
      real(kind=8),dimension(-lmax:lmax,0:lmax) :: norm
      !real(kind=8),dimension(3) :: com
      real(kind=8) :: dnrm2
      !real(kind=8) :: rmax
      real(kind=8),dimension(7,7,4) :: khack

 khack(1:7,1,1) = (/ 2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,1) = (/-1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,1) = (/-3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,1) = (/-5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,1) = (/-5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,1) = (/ 7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,1) = (/-6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,2) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,2) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,2) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,2) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,2) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,2) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,2) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,3) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,3) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,3) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,3) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,3) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,3) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,3) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,4) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,4) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,4) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,4) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,4) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,4) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,4) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)


      !! Check that rmax remains within the box.
      !if (rmax>=0.5d0*(0.5d0*hgrids(1)*maxval(n1i)+0.5d0*hgrids(2)*maxval(n2i)+0.5d0*hgrids(3)*maxval(n3i))) then
      !    call f_err_throw('The radius for the multipole analysis is too small', err_name='BIGDFT_RUNTIME_ERROR')
      !end if



      rmax = 0.d0
      norb_list = f_malloc(maxval(norbsPerType(:)),id='norb_list')
      call f_zero(multipoles)
      rnorm_maxdev = 0.d0
      ist = 0
      !call mpi_barrier(mpi_comm_world,iat)
      do iat=1,natp
          iiat = iat + isat
          ityp=iatype(iiat)
          norb_per_atom = norbsPerType(ityp)
          arr_allocated = .false.
          iiorb = 0
          do iorb=1,norb
              ilr = inwhichlocreg(iorb)
              jat = onwhichatom(iorb)
              !rmax = locrad(iorb)
              if (jat==iiat) then
                  iilr = ilr
                  if (.not.arr_allocated) then
                      rmax(iiat) = min(n1i(ilr)*0.25d0*hgrids(1),n2i(ilr)*0.25d0*hgrids(2),n3i(ilr)*0.25d0*hgrids(3))
                      !rmax(iiat) = 7.d0
                      !write(*,*) 'iiat, rmax', iiat, rmax(iiat)
                      phi1 = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                       1.to.norb_per_atom/),id='ph1i')
                      phi2 = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                       1.to.norb_per_atom/),id='phi2')
                      sphi = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                        -lmax.to.lmax,0.to.lmax,1.to.norb_per_atom/),id='sphi')
                      !write(*,*) 'sphi allocated, natp, size',natp, size(sphi)
                      !write(*,'(a,4i9)') 'phi2 allocated, natp, n1, n2, n3',natp, n1i(ilr), n2i(ilr), n3i(ilr)
                      arr_allocated = .true.
                  else
                      ! Check the dimensions
                      if (n1i(ilr)/=size(phi1,1)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n2i(ilr)/=size(phi1,2)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n3i(ilr)/=size(phi1,3)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n1i(ilr)/=size(phi2,1)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n2i(ilr)/=size(phi2,2)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n3i(ilr)/=size(phi2,3)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                  end if
                  iiorb = iiorb + 1
                  norb_list(iiorb) = iorb
                  ! Apply the spherical harmonic
                  rnorm = 0.d0
                  norm = 0.d0
                  factor_normalization = sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3))
                  do i3=1,n3i(ilr)
                      ii3 = nsi3(ilr) + i3 - 14 - 1
                      z = ii3*0.5d0*hgrids(3) - locregcenter(3,ilr)
                      do i2=1,n2i(ilr)
                          ii2 = nsi2(ilr) + i2 - 14 - 1
                          y = ii2*0.5d0*hgrids(2) - locregcenter(2,ilr)
                          do i1=1,n1i(ilr)
                              ii1 = nsi1(ilr) + i1 - 14 - 1
                              x = ii1*0.5d0*hgrids(1) - locregcenter(1,ilr)
                              ind = (i3-1)*n2i(ilr)*n1i(ilr) + (i2-1)*n1i(ilr) + i1
                              phi1(i1,i2,i3,iiorb) = psir1_get(ist+ind)
                              phi2(i1,i2,i3,iiorb) = psir2_get(ist+ind)
                              if (x**2+y**2+z**2>rmax(iiat)**2) cycle
                              do l=0,lmax
                                  do m=-l,l
                                      tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                                      tt = tt*sqrt(4.d0*pi/real(2*l+1,kind=8))
                                      sphi(i1,i2,i3,m,l,iiorb) = tt*phi1(i1,i2,i3,iiorb)
                                      !write(*,*) 'i1, i1, i2, tt, phi1', i1, i2, i3, tt, phi1(i1,i2,i3,iiorb)
                                      norm(m,l) = norm(m,l) + (tt*factor_normalization)**2*&
                                          real((2*l+3)*(2*l+1),kind=8)/(4.d0*pi*rmax(iiat)**(2*l+3)) !normalization of a solid harmonic within a sphere of radius rmax... hopefully correct
                                      !write(*,*) 'iorb, i1, i2, i3, tt, phi', iorb, i1, i2, i3, tt, phi1(i1,i2,i3,iiorb)
                                  end do
                              end do
                          end do
                      end do
                  end do
                  do l=0,lmax
                      do m=-l,l
                          rnorm_maxdev = max(rnorm_maxdev,abs(1.d0-norm(m,l)))
                      end do
                  end do
                  ist = ist + ind
              end if
          end do
          ! Calculate the scalar product
          do l=0,lmax
              do m=-l,l
                  do iorb=1,norb_per_atom
                      iiorb = norb_list(iorb)
                      do jorb=1,norb_per_atom
                          jjorb = norb_list(jorb)
                          if (get_index==INDEX_AUTOMATIC) then
                              ind = matrixindex_in_compressed(smatl, iiorb, jjorb)
                          else if (get_index==INDEX_MANUALLY) then
                              ind = matrixindex(iiorb, jjorb)
                          end if
                          ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
                          tt = -ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), sphi(1,1,1,m,l,iorb), 1, phi2(1,1,1,jorb), 1)
                          !i = 0
                          !do i3=1,n3i(iilr)
                          !    do i2=1,n2i(iilr)
                          !        do i1=1,n1i(iilr)
                          !            i = i + 1
                          !            write(*,*) i, phi2(i1,i2,i3,jorb), sphi(i1,i2,i3,m,l,iorb)
                          !        end do
                          !    end do
                          !end do
                          ii = 0
                          do i3=1,n3i(iilr)
                              ii3 = nsi3(iilr) + i3 - 14 - 1
                              z = ii3*0.5d0*hgrids(3) - locregcenter(3,iilr)
                              do i2=1,n2i(iilr)
                                  ii2 = nsi2(iilr) + i2 - 14 - 1
                                  y = ii2*0.5d0*hgrids(2) - locregcenter(2,iilr)
                                  do i1=1,n1i(iilr)
                                      ii1 = nsi1(iilr) + i1 - 14 - 1
                                      x = ii1*0.5d0*hgrids(1) - locregcenter(1,iilr)
                                      if (x**2+y**2+z**2>rmax(iiat)**2) cycle
                                      !tt = tt + sphi(i1,i2,i3,m,l,iorb)*phi2(i1,i2,i3,jorb)
                                      ii = ii + 1
                                  end do
                              end do
                          end do
                          tt = tt
                          multipoles(m,l,iiat) = multipoles(m,l,iiat) + matrix_compr(ind)*tt
                          tt2 = ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), phi1(1,1,1,iorb), 1, phi2(1,1,1,jorb), 1)
                      end do
                  end do

              end do
          end do
          if (arr_allocated) then
              call f_free(phi1)
              call f_free(phi2)
              call f_free(sphi)
          end if
      end do
      if (nproc>1) then
          call mpiallred(multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(rnorm_maxdev, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          call mpiallred(rmax, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (verbosity>0 .and. iproc==0) then
          call yaml_mapping_open('Calculation of solid harmonics')
          call yaml_map('Radius of integration sphere',(/minval(rmax),maxval(rmax)/),fmt='(es8.2)')
          call yaml_map('Maximal deviation from normalization',rnorm_maxdev,fmt='(es8.2)')
          call yaml_mapping_close()
      end if


      call f_free(norb_list)
    end subroutine multipole_analysis_core

