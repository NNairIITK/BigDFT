module multipole
  use module_base
  use multipole_base, only: external_potential_descriptors,lmax
  implicit none

  private

  !> Public routines
  public :: potential_from_multipoles
  public :: multipoles_from_density

  contains

    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, pot)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(dp),intent(in) :: hx, hy, hz
      real(dp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, impl, l
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp
      real(dp),dimension(3) :: r

      !$omp parallel &
      !$omp default(none) &
      !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !$omp private(i1, i2, i3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !$omp do
      do i3=is3,ie3
          z = real(i3,kind=8)*hz
          do i2=is2,ie2
              y = real(i2,kind=8)*hy
              do i1=is1,ie1
                  x = real(i1,kind=8)*hx
                  do impl=1,ep%nmpl
                      r(1) = ep%mpl(impl)%rxyz(1) - x
                      r(2) = ep%mpl(impl)%rxyz(2) - y
                      r(3) = ep%mpl(impl)%rxyz(3) - z 
                      rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
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
                  end do
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel


      contains


        function calc_monopole(q, rnrm1) result(mpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(1),intent(in) :: q
          real(dp),intent(in) :: rnrm1
          real(dp) :: mpm

          mpm = -q(1)/rnrm1

        end function calc_monopole


        function calc_dipole(q, r, rnrm3) result(dpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(3),intent(in) :: q
          real(dp),intent(in) :: rnrm3
          real(dp),dimension(3),intent(in) :: r
          real(dp) :: dpm

          dpm = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
          dpm = -dpm/rnrm3

        end function calc_dipole


        function calc_quadropole(q, r, rnrm5) result(qpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(5),intent(in) :: q
          real(dp),intent(in) :: rnrm5
          real(dp),dimension(3),intent(in) :: r
          real(dp) :: qpm
          ! Local variables
          real(dp),dimension(3,3) :: qq

          qq(1,1) = q(1)
          qq(2,1) = q(2)
          qq(3,1) = q(3)
          qq(1,2) = qq(2,1)
          qq(2,2) = q(4)
          qq(3,2) = q(5)
          qq(1,3) = qq(3,1)
          qq(2,3) = qq(3,2)
          qq(3,3) = 1.0_dp-qq(1,1)-qq(2,2)

          qpm = qq(1,1)*r(1)*r(1) + &
               qq(2,1)*r(2)*r(1) + &
               qq(3,1)*r(3)*r(1) + &
               qq(1,2)*r(1)*r(2) + &
               qq(2,2)*r(2)*r(2) + &
               qq(3,2)*r(3)*r(2) + &
               qq(1,3)*r(1)*r(3) + &
               qq(2,3)*r(2)*r(3) + &
               qq(3,3)*r(3)*r(3)
          qpm = -0.5_dp*qpm/rnrm5

        end function calc_quadropole


        !function calc_octopole()
        !end function calc_octopole

    end subroutine potential_from_multipoles



    subroutine multipoles_from_density(iproc, nproc, at, lzd, smats, smatl, orbs, &
               npsidim, lphi, norbsPerType, collcom, collcom_sr, orthpar, &
               ovrlp, kernel, meth_overlap)
      use module_base
      use module_types
      use module_interfaces
      use sparsematrix_base, only: sparsematrix_malloc0, SPARSE_FULL, assignment(=)
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim, meth_overlap
      type(atoms_data),intent(in) :: at
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(orbitals_data),intent(in) :: orbs
      real(kind=8),dimension(npsidim),intent(in) :: lphi
      integer,dimension(at%astruct%ntypes),intent(in) :: norbsPerType
      type(comms_linear),intent(in) :: collcom
      type(comms_linear),intent(in) :: collcom_sr
      type(orthon_data),intent(in) :: orthpar
      type(matrices),intent(in) :: kernel
      type(matrices),intent(inout) :: ovrlp !in principle also intent(in)

      ! Local variables
      integer :: ist, istr, iorb, iiorb, ilr, ii, natp, isat, nr, jproc, iat, n, norb_get, istr_get
      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
      integer :: ii1, ii2, ii3, jjorb
      real(kind=8),dimension(:),allocatable :: psir
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      type(workarr_sumrho) :: w
      integer,dimension(:),allocatable :: nat_par, norb_list
      real(kind=8),dimension(:),allocatable :: psir_get
      real(kind=8),dimension(:,:,:,:),allocatable :: phi
      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
      integer,dimension(:,:),allocatable :: comms
      logical :: can_use_transposed, arr_allocated
      real(kind=8) :: ddot, x, y, z, tt
      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
      real(kind=8),dimension(-lmax:lmax,0:lmax,at%astruct%nat) :: multipoles
      real(kind=8) :: factor_normalization
      character(len=20) :: atomname
      !testarr = 0.d0

      call f_routine(id='multipoles_from_density')

      ! Orthogonalize the support functions
      can_use_transposed = .false.
      methTransformOverlap = 1020
      phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
      phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
      phi_ortho = f_malloc(npsidim,id='phi_ortho')
      call vcopy(npsidim, lphi(1), 1, phi_ortho(1), 1)
      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
           1.d-8, npsidim, orbs, lzd, &
           smats, smatl, collcom, orthpar, &
           phi_ortho, phit_c, phit_f, &
           can_use_transposed)
      call f_free_ptr(phit_c)
      call f_free_ptr(phit_f)


      ! Transform the support functions to real space
      psir = f_malloc(collcom_sr%ndimpsi_c,id='psir')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,lzd%Llr(ilr),.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi_ortho(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      call f_free(phi_ortho)
      if(istr/=collcom_sr%ndimpsi_c+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
          stop
      end if


      ! Switch to a partition over the atoms
      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
      ii = at%astruct%nat/nproc
      nat_par(0:nproc-1) = ii
      ii = at%astruct%nat-nproc*ii
      nat_par(0:ii-1) = nat_par(0:ii-1) + 1
      if (sum(nat_par(:))/=at%astruct%nat) then
          call f_err_throw('wrong partition of the atoms',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      natp = nat_par(iproc)
      if (iproc==0) then
          isat = 0
      else
          isat = sum(nat_par(0:iproc-1))
      end if
      call f_free(nat_par)
      comms = f_malloc((/4,orbs%norb/),id='comms')
      nr = 0
      norb_get = 0
      istr = 1
      istr_get = 1
      do jproc=0,nproc-1
          istr = 0
          do iorb=1,orbs%norb_par(jproc,0)
              iiorb = iorb + orbs%isorb_par(jproc)
              ilr = orbs%inwhichlocreg(iiorb)
              iat = orbs%onwhichatom(iiorb)
              n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
              !write(*,'(a,5i8)') 'iproc, jproc, iiorb, iat', iproc, jproc, iiorb, iat
              if (iat>=isat+1 .and. iat<=isat+natp) then
                  norb_get = norb_get + 1
                  comms(1,norb_get) = jproc
                  comms(2,norb_get) = n
                  comms(3,norb_get) = istr
                  comms(4,norb_get) = istr_get
                  nr = nr + n
                  !write(*,'(a,5i8)') 'iproc, jproc, n, istr, istr_get', iproc, jproc, n, istr, istr_get
                  istr = istr + n
                  istr_get = istr_get + n
              else
                  istr = istr + n
              end if
          end do
      end do
      !write(*,*) 'iproc, nr', iproc, nr
      !do iorb=1,norb_get
      !    write(*,'(a,2i7,4i9)') 'iproc, iorb, comm',iproc, iorb, comms(:,iorb)
      !end do
      psir_get = f_malloc(nr,id='psir_get')
      window = mpiwindow(size(psir), psir(1), bigdft_mpi%mpi_comm)
      do iorb=1,norb_get
          jproc = comms(1,iorb)
          n = comms(2,iorb)
          ioffset = comms(3,iorb)
          istr = comms(4,iorb)
          !write(*,'(5(a,i0))') 'task ',iproc,' gets ',n,' elements at position ', &
          !                     istr,' from position ',ioffset+1,' on task ',jproc
          call mpiget(psir_get(istr), n, jproc, int(ioffset,kind=mpi_address_kind), window)
      end do
      call mpi_fenceandfree(window)
      call f_free(psir)
      call f_free(comms)
      istr = 1
      !write(*,*) 'iproc, isat, natp', iproc, isat, natp
      do iorb=1,orbs%norb
          ilr = orbs%inwhichlocreg(iorb)
          iat = orbs%onwhichatom(iorb)
          if (iat>=isat+1 .and. iat<=isat+natp) then
              n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
              !write(*,'(a,4i8,es16.6)') 'iproc, iorb, n, istr, ddot', &
              !    iproc, iorb, n, istr, ddot(n, psir_get(istr), 1, psir_get(istr), 1)
              !testarr(2,iorb) = ddot(n, psir_get(istr), 1, psir_get(istr), 1)
              istr = istr + n
          end if
      end do
      !call mpiallred(testarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      !if (iproc==0) then
      !    do iorb=1,orbs%norb
      !        write(*,*) 'DIFF, iorb, val', iorb, abs(testarr(1,iorb)-testarr(2,iorb))
      !    end do
      !end if


      ! Calculate the kernel for orthormal support functions
      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='kernel_ortho')
      call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, meth_overlap, smats, smatl, &
           ovrlp, kernel, kernel_ortho)


      ! Apply the spherical harmonics to the suport functions.
      ! For the support functions on atom A we only need to apply 
      ! the spherical harmonics centered as well on atom A.

      !lmax = 1
      norb_list = f_malloc(maxval(norbsPerType(:)),id='norb_list')
          call f_zero(multipoles)
      factor_normalization = sqrt(4.d0*pi_param)
      do iat=1,natp
          iiat = iat + isat
          ityp=at%astruct%iatype(iiat)
          norb_per_atom = norbsPerType(ityp)
          arr_allocated = .false.
          iiorb = 0
          ist = 0
          do iorb=1,orbs%norb
              ilr = orbs%inwhichlocreg(iorb)
              jat = orbs%onwhichatom(iorb)
              if (jat==iiat) then
                  if (.not.arr_allocated) then
                      phi = f_malloc((/1.to.lzd%Llr(ilr)%d%n1i,1.to.lzd%Llr(ilr)%d%n2i,1.to.lzd%Llr(ilr)%d%n3i, &
                                       1.to.norb_per_atom/),id='sphi')
                      sphi = f_malloc((/1.to.lzd%Llr(ilr)%d%n1i,1.to.lzd%Llr(ilr)%d%n2i,1.to.lzd%Llr(ilr)%d%n3i, &
                                        1.to.norb_per_atom,-lmax.to.lmax,0.to.lmax/),id='sphi')
                      arr_allocated = .true.
                  end if
                  iiorb = iiorb + 1
                  norb_list(iiorb) = iorb
                  ! Apply the spherical harmonic
                  do i3=1,lzd%llr(ilr)%d%n3i
                      ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
                      z = ii3*0.5d0*lzd%hgrids(3) - lzd%llr(ilr)%locregcenter(3)
                      do i2=1,lzd%llr(ilr)%d%n2i
                          ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
                          y = ii2*0.5d0*lzd%hgrids(2) - lzd%llr(ilr)%locregcenter(2)
                          do i1=1,lzd%llr(ilr)%d%n1i
                              ii1 = lzd%llr(ilr)%nsi1 + i1 -14 -1
                              x = ii1*0.5d0*lzd%hgrids(1) - lzd%llr(ilr)%locregcenter(1)
                              ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                              phi(i1,i2,i3,iiorb) = psir_get(ist+ind)
                              do l=0,lmax
                                  do m=-l,l
                                      tt = spherical_harmonic(l, m, x, y, z)
                                      !if (iproc==0) write(*,'(a,3i4,3f9.3,3x,3f9.3,es20.12)') 'i1, i2, i3, xx, yy, zz, locregcenter, tt', &
                                      !    i1, i2, i3, ii3*0.5d0*lzd%hgrids(3), ii2*0.5d0*lzd%hgrids(2), ii1*0.5d0*lzd%hgrids(1), lzd%llr(ilr)%locregcenter, tt
                                      sphi(i1,i2,i3,iiorb,m,l) = factor_normalization*tt*phi(i1,i2,i3,iiorb)
                                      !sphi(i1,i2,i3,iiorb,m,l) = phi(i1,i2,i3,iiorb)
                                  end do
                              end do
                          end do
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
                          ind = matrixindex_in_compressed(smatl, iiorb, jjorb)
                          tt = ddot(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, sphi(1,1,1,iorb,m,l), 1, phi(1,1,1,jorb), 1)
                          !write(*,'(a,6i7,2es18.8)') 'iiat, m, l, iiorb, jjorb, ind, tt, K', iiat, m, l, iiorb, jjorb, ind, tt, kernel%matrix_compr(ind)
                          multipoles(m,l,iiat) = multipoles(m,l,iiat) + kernel%matrix_compr(ind)*tt
                          !write(*,'(a,5i8,es16.8)') 'iproc, l, m, iorb, jorb, ddot', iproc, l, m, iorb, jorb, tt
                      end do
                  end do

              end do
          end do
          if (arr_allocated) then
              call f_free(phi)
              call f_free(sphi)
          end if
      end do
      call f_free(psir_get)

      call f_free(norb_list)

      call mpiallred(multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
      if (iproc==0) then
          call yaml_sequence_open('Multipole analysis')
          do iat=1,at%astruct%nat
              call yaml_sequence(advance='no')
              atomname=at%astruct%atomnames(at%astruct%iatype(iat))
              call yaml_sequence_open(trim(atomname))
              do l=0,lmax
                  call yaml_sequence(advance='no')
                  call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat),fmt='(1es16.8)')
                  call yaml_newline()
              end do
              !call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
              call yaml_sequence_close()
          end do
          call yaml_sequence_close()
      end if

      !if (iproc==0) then
      !    !write(*,*) 'SUM(KERNEL)',sum(kernel_ortho), sum(kernel%matrix_compr)
      !    tt = 0.d0
      !    do iorb=1,orbs%norb
      !        do jorb=1,orbs%norb
      !             ind = matrixindex_in_compressed(smatl, iorb, jorb)
      !             write(*,*) 'iorb, jorb, val', iorb, jorb, kernel_ortho(ind)
      !             if (iorb==jorb) then
      !                 tt = tt + kernel_ortho(ind)
      !             end if
      !        end do
      !    end do
      !    !write(*,*) 'TRACE',tt
      !end if

      call f_free(kernel_ortho)

      call f_release_routine()


      contains

          !> Calculates the real spherical harmonic for given values of l, m, x, y, z
          function spherical_harmonic(l, m, x, y, z) result(sh)
            use module_base, only: pi => pi_param
            implicit none
            ! Calling arguments
            integer,intent(in) :: l, m
            real(kind=8),intent(in) :: x, y, z
            real(kind=8) :: sh

            ! Local variables
            integer,parameter :: l_max=2
            real(kind=8) :: r, r2

            if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
            if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
                err_name='BIGDFT_RUNTIME_ERROR')
            if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')

            select case (l)
            case (0)
                sh = 0.5d0*sqrt(1/pi)
            case (1)
                r = sqrt(x**2+y**2+z**2)
                select case (m)
                case (-1)
                    sh = sqrt(3.d0/(4.d0*pi))*y/r
                case (0)
                    sh = sqrt(3.d0/(4.d0*pi))*z/r
                case (1)
                    sh = sqrt(3.d0/(4.d0*pi))*x/r
                end select
            case (2)
                r2 = x**2+y**2+z**2
                select case (m)
                case (-2)
                    sh = 0.5d0*sqrt(15.d0/pi)*x*y/r2
                case (-1)
                    sh = 0.5d0*sqrt(15.d0/pi)*y*z/r2
                case (0)
                    sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
                case (1)
                    sh = 0.5d0*sqrt(15.d0/pi)*z*x/r2
                case (2)
                    sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
                end select
            end select

          end function spherical_harmonic


    end subroutine multipoles_from_density



    !> Calculate S^1/2 * K * S^1/2, which is the kernel corresponding to a
    !! orthonormal set of supoprt functions.
    subroutine kernel_for_orthonormal_basis(iproc, nproc, norbp, meth_overlap, smats, smatl, &
               ovrlp, kernel, weight_matrix_compr)
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, SPARSE_TASKGROUP, &
                                   matrices_null, assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   deallocate_matrices
      use matrix_operations, only: overlapPowerGeneral
      use sparsematrix, only: matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups
      implicit none

      ! Calling arguments
      integer :: iproc, nproc, norbp, meth_overlap
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(matrices),intent(in) :: kernel
      type(matrices),intent(inout) :: ovrlp
      real(kind=8),dimension(smatl%nvctr*smatl%nspin),intent(out) :: weight_matrix_compr

      ! Local variables
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:),allocatable :: weight_matrix_compr_tg, proj_ovrlp_half_compr
      real(kind=8) :: max_error, mean_error

      inv_ovrlp(1) = matrices_null()
      inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='inv_ovrlp(1)%matrix_compr')

      call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), -1, &
           imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
           ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
           max_error=max_error, mean_error=mean_error)
      !call f_free_ptr(ovrlp%matrix)

      proj_ovrlp_half_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='proj_mat_compr')
      if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              kernel%matrix_compr, inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr)
      end if
      weight_matrix_compr_tg = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='weight_matrix_compr_tg')
      if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr, weight_matrix_compr_tg)
      end if
      call f_free(proj_ovrlp_half_compr)

      call deallocate_matrices(inv_ovrlp(1))

      ! Maybe this can be improved... not really necessary to gather the entire matrix
      !weight_matrix_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='weight_matrix_compr')
      call gather_matrix_from_taskgroups(iproc, nproc, smatl, weight_matrix_compr_tg, weight_matrix_compr)

      call f_free(weight_matrix_compr_tg)

    end subroutine kernel_for_orthonormal_basis

end module multipole
