module postprocessing_linear
  use public_enums
  implicit none

  private

  !> Public routines
  public :: loewdin_charge_analysis
  public :: loewdin_charge_analysis_core
  public :: support_function_multipoles
  public :: build_ks_orbitals

  contains

    subroutine loewdin_charge_analysis(iproc,tmb,atoms,denspot,&
               calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   DENSE_FULL, assignment(=), &
                                   matrices_null, allocate_matrices, deallocate_matrices
      use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                              uncompress_matrix2
      use transposed_operations, only: calculate_overlap_transposed
      use matrix_operations, only: overlapPowerGeneral
      use yaml_output
      implicit none
      integer,intent(in) :: iproc
      type(dft_wavefunction),intent(inout) :: tmb
      type(atoms_data),intent(in) :: atoms
      type(DFT_local_fields), intent(inout) :: denspot
      logical,intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
      integer,intent(in) :: meth_overlap
    
      !local variables
      !integer :: ifrag,ifrag_ref,isforb,jorb
      integer :: iorb,ierr
      real(kind=gp), allocatable, dimension(:,:,:) :: proj_mat
      real(kind=gp), allocatable, dimension(:,:) :: proj_ovrlp_half, weight_matrixp
      character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
      real(kind=gp) :: max_error, mean_error
      type(matrices),dimension(1) :: inv_ovrlp
    
      ! new variables
      integer :: iat
      real(kind=8),dimension(:,:),allocatable :: weight_matrix
      !real(kind=gp),dimension(:,:),pointer :: ovrlp
      real(kind=8) :: total_charge, total_net_charge
      real(kind=8),dimension(:),allocatable :: charge_per_atom
      !logical :: psit_c_associated, psit_f_associated
    
    
      ! needs parallelizing/converting to sparse
      ! re-use overlap matrix if possible either before or after
    
      call f_routine(id='loewdin_charge_analysis')
    
      !inv_ovrlp(1) = matrices_null()
      !call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))
    
    
    
      if (calculate_overlap_matrix) then
         if(.not.tmb%can_use_transposed) then
             !!if(.not.associated(tmb%psit_c)) then
             !!    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
             !!    psit_c_associated=.false.
             !!else
             !!    psit_c_associated=.true.
             !!end if
             !!if(.not.associated(tmb%psit_f)) then
             !!    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
             !!    psit_f_associated=.false.
             !!else
             !!    psit_f_associated=.true.
             !!end if
             call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                  TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
             tmb%can_use_transposed=.true.
         end if
    
         call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
              tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
         !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
         ! This can then be deleted if the transition to the new type has been completed.
         !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr
    
    
         !!if (.not.psit_c_associated) then
         !!   call f_free_ptr(tmb%psit_c)
         !!   tmb%can_use_transposed=.false.
         !!end if
         !!if (.not.psit_f_associated) then
         !!   call f_free_ptr(tmb%psit_f)
         !!   tmb%can_use_transposed=.false.
         !!end if
      end if

      call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, &
               tmb%orbs%norb_par, tmb%orbs%isorb_par, meth_overlap, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_)
    
!!!!      if (calculate_ovrlp_half) then
!!!!         tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
!!!!         call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, &
!!!!              tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
!!!!         call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/2/), &
!!!!              tmb%orthpar%blocksize_pdsyev, &
!!!!              imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
!!!!              ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
!!!!              max_error=max_error, mean_error=mean_error)
!!!!         !!ovrlp_half=tmb%linmat%ovrlp%matrix
!!!!         call f_free_ptr(tmb%linmat%ovrlp_%matrix)
!!!!      end if
!!!!    
!!!!      ! optimize this to just change the matrix multiplication?
!!!!      proj_mat = sparsematrix_malloc0(tmb%linmat%l,iaction=DENSE_FULL,id='proj_mat')
!!!!    
!!!!      call uncompress_matrix2(iproc, bigdft_mpi%nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, proj_mat)
!!!!      !!isforb=0
!!!!      !!do ifrag=1,input%frag%nfrag
!!!!      !!   ifrag_ref=input%frag%frag_index(ifrag)
!!!!      !!   if (ifrag==ifrag_charged(1)) then
!!!!      !!      do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
!!!!      !!         proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
!!!!      !!      end do
!!!!      !!   end if
!!!!      !!   !!if (nfrag_charged==2) then
!!!!      !!   !!   if (ifrag==ifrag_charged(2)) then
!!!!      !!   !!      do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
!!!!      !!   !!         proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
!!!!      !!   !!      end do
!!!!      !!   !!   end if
!!!!      !!   !!end if
!!!!      !!   isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
!!!!      !!end do
!!!!    
!!!!      proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
!!!!      if (tmb%orbs%norbp>0) then
!!!!         call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
!!!!                tmb%orbs%norb, 1.d0, &
!!!!                proj_mat(1,1,1), tmb%orbs%norb, &
!!!!                inv_ovrlp(1)%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, 0.d0, &
!!!!                proj_ovrlp_half(1,1), tmb%orbs%norb)
!!!!      end if
!!!!      call f_free(proj_mat)
!!!!      weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
!!!!      if (tmb%orbs%norbp>0) then
!!!!         call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
!!!!              tmb%orbs%norb, 1.d0, &
!!!!              inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, &
!!!!              proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
!!!!              weight_matrixp(1,1), tmb%orbs%norb)
!!!!      end if
!!!!      !call f_free_ptr(ovrlp_half)
!!!!      call f_free(proj_ovrlp_half)
!!!!      weight_matrix=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='weight_matrix')
!!!!      if (bigdft_mpi%nproc>1) then
!!!!         call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix, &
!!!!              tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
!!!!              mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!!!      else
!!!!         call vcopy(tmb%orbs%norb*tmb%orbs%norb,weight_matrixp(1,1),1,weight_matrix(1,1),1)
!!!!      end if
!!!!      call f_free(weight_matrixp)
!!!!      !call compress_matrix(bigdft_mpi%iproc,weight_matrix)
!!!!    
!!!!      charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')
!!!!      !!do iorb=1,tmb%orbs%norb
!!!!      !!    do jorb=1,tmb%orbs%norb
!!!!      !!        if (iproc==0) write(*,'(a,2i7,es16.7)') 'iorb,jorb,weight_matrix(jorb,iorb)', iorb,jorb,weight_matrix(jorb,iorb)
!!!!      !!        if (iorb==jorb) then
!!!!      !!            total_charge = total_charge + weight_matrix(jorb,iorb)
!!!!      !!            iat=tmb%orbs%onwhichatom(iorb)
!!!!      !!            charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(jorb,iorb)
!!!!      !!        end if
!!!!      !!    end do
!!!!      !!end do
!!!!      !!if (iproc==0) then
!!!!      !!    do iat=1,atoms%astruct%nat
!!!!      !!        write(*,*) 'iat, partial total_charge', iat, charge_per_atom(iat)
!!!!      !!    end do
!!!!      !!    write(*,*) 'total total_charge',total_charge
!!!!      !!    if (iproc==0) call write_partial_charges()
!!!!      !!end if
!!!!    
!!!!      do iorb=1,tmb%orbs%norb
!!!!          iat=tmb%orbs%onwhichatom(iorb)
!!!!          charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(iorb,iorb)
!!!!      end do
!!!!      if (iproc==0) then
!!!!          !call write_partial_charges()
!!!!          call write_partial_charges(atoms, charge_per_atom)
!!!!          call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
!!!!          call calculate_dipole(iproc, atoms, charge_per_atom)
!!!!          call calculate_quadropole(iproc, atoms, charge_per_atom)
!!!!          call yaml_sequence_close()
!!!!      end if
!!!!      !!call support_function_multipoles()
!!!!    
!!!!      call deallocate_matrices(inv_ovrlp(1))
!!!!    
!!!!      call f_free(charge_per_atom)
!!!!      call f_free(weight_matrix)
      call f_release_routine()
    
    
    end subroutine loewdin_charge_analysis




    subroutine loewdin_charge_analysis_core(iproc, nproc, norb, norbp, isorb, &
               norb_par, isorb_par, meth_overlap, smats, smatl, atoms, kernel, ovrlp)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, DENSE_FULL, &
                                   deallocate_matrices, matrices_null, allocate_matrices
      use sparsematrix, only: uncompress_matrix2
      use matrix_operations, only: overlapPowerGeneral
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norb, norbp, isorb, meth_overlap
      integer,dimension(0:nproc-1),intent(in) :: norb_par, isorb_par
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(atoms_data),intent(in) :: atoms
      type(matrices),intent(in) :: kernel
      type(matrices),intent(inout) :: ovrlp

      ! Local variables
      integer :: ierr, iorb, iat
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:,:,:),allocatable :: proj_mat
      real(kind=8),dimension(:,:),allocatable :: weight_matrix, weight_matrixp, proj_ovrlp_half
      real(kind=8),dimension(:),allocatable :: charge_per_atom
      real(kind=8) :: mean_error, max_error


      inv_ovrlp(1) = matrices_null()
      call allocate_matrices(smatl, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

      ovrlp%matrix = sparsematrix_malloc_ptr(smats, iaction=DENSE_FULL, id='ovrlp%matrix')
      call uncompress_matrix2(iproc, nproc, smats, &
           ovrlp%matrix_compr, ovrlp%matrix)
      call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), -1, &
           imode=2, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
           ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
           max_error=max_error, mean_error=mean_error)
      call f_free_ptr(ovrlp%matrix)
    
      ! optimize this to just change the matrix multiplication?
      proj_mat = sparsematrix_malloc0(smatl,iaction=DENSE_FULL,id='proj_mat')
    
      call uncompress_matrix2(iproc, nproc, smatl, kernel%matrix_compr, proj_mat)

      proj_ovrlp_half=f_malloc((/norb,norbp/),id='proj_ovrlp_half')
      if (norbp>0) then
         call dgemm('n', 'n', norb, norbp, &
                norb, 1.d0, &
                proj_mat(1,1,1), norb, &
                inv_ovrlp(1)%matrix(1,isorb+1,1), norb, 0.d0, &
                proj_ovrlp_half(1,1), norb)
      end if
      call f_free(proj_mat)
      weight_matrixp=f_malloc((/norb,norbp/), id='weight_matrixp')
      if (norbp>0) then
         call dgemm('n', 'n', norb, norbp, &
              norb, 1.d0, &
              inv_ovrlp(1)%matrix(1,1,1), norb, &
              proj_ovrlp_half(1,1), norb, 0.d0, &
              weight_matrixp(1,1), norb)
      end if
      call f_free(proj_ovrlp_half)
      weight_matrix=f_malloc((/norb,norb/), id='weight_matrix')
      if (nproc>1) then
         call mpi_allgatherv(weight_matrixp, norb*norbp, mpi_double_precision, weight_matrix, &
              norb*norb_par(:), norb*isorb_par, &
              mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call vcopy(norb*norb,weight_matrixp(1,1),1,weight_matrix(1,1),1)
      end if
      call f_free(weight_matrixp)
    
      charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')
    
      do iorb=1,norb
          iat=smats%on_which_atom(iorb)
          charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(iorb,iorb)
      end do
      if (iproc==0) then
          !call write_partial_charges()
          call write_partial_charges(atoms, charge_per_atom, .true.)
          call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
          call calculate_dipole(iproc, atoms, charge_per_atom)
          call calculate_quadropole(iproc, atoms, charge_per_atom)
          call yaml_sequence_close()
      end if
      !!call support_function_multipoles()
    
      call deallocate_matrices(inv_ovrlp(1))
      call f_free(charge_per_atom)
      call f_free(weight_matrix)

    end subroutine loewdin_charge_analysis_core


    subroutine write_partial_charges(atoms, charge_per_atom, write_gnuplot)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      logical,intent(in) :: write_gnuplot
      ! Local variables
      integer :: iat, itypes, iitype, nntype, intype, iunit
      real(kind=8) :: total_charge, total_net_charge, frac_charge, range_min, range_max
      character(len=20) :: atomname, colorname
      real(kind=8),dimension(2) :: charges
      character(len=128) :: output
      character(len=2) :: backslash
      integer,parameter :: ncolors = 7
      character(len=20),dimension(ncolors),parameter :: colors=(/'violet', &
                                                                 'blue  ', &
                                                                 'cyan  ', &
                                                                 'green ', &
                                                                 'yellow', &
                                                                 'orange', &
                                                                 'red   '/)

      call yaml_sequence_open('Loewdin charge analysis (charge / net charge)')
      total_charge=0.d0
      total_net_charge=0.d0
      do iat=1,atoms%astruct%nat
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          atomname=atoms%astruct%atomnames(atoms%astruct%iatype(iat))
          charges(1)=-charge_per_atom(iat)
          charges(2)=-(charge_per_atom(iat)-real(atoms%nelpsp(atoms%astruct%iatype(iat)),kind=8))
          total_charge = total_charge + charges(1)
          total_net_charge = total_net_charge + charges(2)
          call yaml_map(trim(atomname),charges,fmt='(1es20.12)')
          call yaml_mapping_close(advance='no')
          call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
      end do
      call yaml_sequence(advance='no')
      call yaml_map('total charge',total_charge,fmt='(es16.8)')
      call yaml_sequence(advance='no')
      call yaml_map('total net charge',total_net_charge,fmt='(es16.8)')
      call yaml_sequence_close()

      if (write_gnuplot) then
          output='chargeanalysis.gp'
          call yaml_map('output file',trim(output))
          call f_open_file(iunit, file=trim(output), binary=.false.)
          write(iunit,'(a)') '# plot the fractional charge as a normalized sum of Gaussians'
          write(iunit,'(a)') 'set samples 500000'
          range_min = minval(-(charge_per_atom(:)-real(atoms%nelpsp(atoms%astruct%iatype(:)),kind=8))) - 0.1d0
          range_max = maxval(-(charge_per_atom(:)-real(atoms%nelpsp(atoms%astruct%iatype(:)),kind=8))) + 0.1d0
          write(iunit,'(a,2(es12.5,a))') 'set xrange[',range_min,':',range_max,']'
          write(iunit,'(a)') 'sigma=0.005'
          write(backslash,'(a)') '\ '
          do itypes=1,atoms%astruct%ntypes
              nntype = 0
              do iat=1,atoms%astruct%nat
                  iitype = (atoms%astruct%iatype(iat))
                  if (iitype==itypes) then
                      nntype = nntype + 1
                  end if
              end do
              write(iunit,'(a,i0,a,i0,2a)') 'f',itypes,'(x) = 1/',nntype,'.0*( '//trim(backslash)
              intype = 0
              do iat=1,atoms%astruct%nat
                  iitype = (atoms%astruct%iatype(iat))
                  if (iitype==itypes) then
                      intype = intype + 1
                      frac_charge = -(charge_per_atom(iat)-real(atoms%nelpsp(atoms%astruct%iatype(iat)),kind=8))
                      if (intype<nntype) then
                          write(iunit,'(a,es16.9,a)') '  1.0*exp(-(x-',frac_charge,')**2/(2*sigma**2)) + '//trim(backslash)
                      else
                          write(iunit,'(a,es16.9,a)') '  1.0*exp(-(x-',frac_charge,')**2/(2*sigma**2)))'
                      end if
                  end if
              end do
              atomname=atoms%astruct%atomnames(itypes)
              if (itypes<ncolors) then
                  colorname = colors(itypes)
              else
                  colorname = 'color'
              end if
              if (itypes==1) then
                  write(iunit,'(a,i0,5a)') "plot f",itypes,"(x) lc rgb '",trim(colorname), &
                      "' lt 1 lw 2 w l title '",trim(atomname),"'"
              else
                  write(iunit,'(a,i0,5a)') "replot f",itypes,"(x) lc rgb '",trim(colorname), &
                      "' lt 1 lw 2 w l title '",trim(atomname),"'"
              end if
          end do
      end if
    end subroutine write_partial_charges


    subroutine calculate_dipole(iproc, atoms, charge_per_atom)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      integer,intent(in) :: iproc
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      ! Local variables
      integer :: iat
      real(kind=8),dimension(3) :: dipole_elec, dipole_cores, dipole_net
    
      dipole_cores(1:3)=0._gp
      do iat=1,atoms%astruct%nat
         dipole_cores(1:3)=dipole_cores(1:3)+atoms%nelpsp(atoms%astruct%iatype(iat))*atoms%astruct%rxyz(1:3,iat)
      end do
    
      dipole_elec=0.d0
      do iat=1,atoms%astruct%nat
          dipole_elec(1:3) = dipole_elec(1:3) -charge_per_atom(iat)*atoms%astruct%rxyz(1:3,iat)
      end do
    
      dipole_net=dipole_cores+dipole_elec
    
      if (iproc==0) then
          !!call yaml_map('core dipole', dipole_cores)
          !!call yaml_map('electronic dipole', dipole_elec)
          call yaml_map('net dipole', dipole_net,fmt='(es12.5)')
      end if
    
    end subroutine calculate_dipole


    subroutine calculate_quadropole(iproc, atoms, charge_per_atom)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      integer,intent(in) :: iproc
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      ! Local variables
      real(kind=8),dimension(3,3) :: quadropole_elec, quadropole_cores, quadropole_net
      real(kind=8),dimension(3) :: charge_center_cores, charge_center_charge
      integer :: iat, i, j
      real(kind=8) :: delta_term, rj, ri, q, qtot
    
    
      ! charge center of the cores
      charge_center_cores(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=atoms%nelpsp(atoms%astruct%iatype(iat))
          charge_center_cores(1:3) = charge_center_cores(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_cores=charge_center_cores/qtot
    
    
      ! charge center of the charge
      charge_center_charge(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=-charge_per_atom(iat)
          charge_center_charge(1:3) = charge_center_charge(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_charge=charge_center_charge/qtot
    
    
      quadropole_cores(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=atoms%nelpsp(atoms%astruct%iatype(iat))
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = atoms%astruct%rxyz(1,iat)**2 + atoms%astruct%rxyz(2,iat)**2 + atoms%astruct%rxyz(3,iat)**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)
                 ri=atoms%astruct%rxyz(i,iat)
                 quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_cores(j,i) = quadropole_cores(j,i) + &
                 !!                        atoms%nelpsp(atoms%astruct%iatype(iat))* &
                 !!                          (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do
    
    
      quadropole_elec(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=-charge_per_atom(iat)
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = (atoms%astruct%rxyz(1,iat)+(charge_center_cores(1)-charge_center_charge(1)))**2 + &
                                  (atoms%astruct%rxyz(2,iat)+(charge_center_cores(2)-charge_center_charge(2)))**2 + &
                                  (atoms%astruct%rxyz(3,iat)+(charge_center_cores(3)-charge_center_charge(3)))**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)+(charge_center_cores(j)-charge_center_charge(j))
                 ri=atoms%astruct%rxyz(i,iat)+(charge_center_cores(i)-charge_center_charge(i))
                 quadropole_elec(j,i) = quadropole_elec(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_elec(j,i) = quadropole_elec(j,i) + &
                 !!                       -charge_per_atom(iat)* &
                 !!                         (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do
    
      quadropole_net=quadropole_cores+quadropole_elec
    
      if (iproc==0) then
          !!call yaml_sequence_open('core quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_cores(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()
    
          !!call yaml_sequence_open('electronic quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_elec(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()
    
          call yaml_sequence_open('net quadropole')
          do i=1,3
             call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es12.5)')))
          end do
          call yaml_sequence(advance='no')
          call yaml_map('trace of quadropole matrix',&
               quadropole_net(1,1)+quadropole_net(2,2)+quadropole_net(3,3),fmt='(es12.2)')
          call yaml_sequence_close()
      end if
    
    end subroutine calculate_quadropole


    subroutine support_function_multipoles(iproc, tmb, atoms, denspot)
      use module_base
      use module_types
      use yaml_output
      
      ! Calling arguments
      integer,intent(in) :: iproc
      type(DFT_wavefunction),intent(in) :: tmb
      type(atoms_data),intent(in) :: atoms
      type(DFT_local_fields), intent(inout) :: denspot
    
      integer :: ist, istr, iorb, iiorb, ilr, i, iat
      real(kind=8),dimension(3) :: charge_center_elec
      real(kind=8),dimension(:),allocatable :: phir
      type(workarr_sumrho) :: w
      character(len=20) :: atomname
      real(kind=8),dimension(:,:),allocatable :: dipole_net
      real(kind=8),dimension(:,:,:),allocatable :: quadropole_net
    
      call f_routine(id='support_function_multipoles')
    
      phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
      dipole_net = f_malloc0((/3,tmb%orbs%norb/),id='dipole_net')
      quadropole_net = f_malloc0((/3,3,tmb%orbs%norb/),id='quadropole_net')
    
      !call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
      !call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))
    
      ist=1
      istr=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          iat=tmb%orbs%onwhichatom(iiorb)
          call initialize_work_arrays_sumrho(1,tmb%lzd%Llr(ilr),.true.,w)
          ! Transform the support function to real space
          call daub_to_isf(tmb%lzd%llr(ilr), w, tmb%psi(ist), phir(istr))
          call deallocate_work_arrays_sumrho(w)
          ! Calculate the charge center
          call charge_center(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
               denspot%dpbox%hgrids, phir(istr), charge_center_elec)
          !write(*,*) 'ilr, tmb%lzd%llr(ilr)%locregcenter', iat, tmb%lzd%llr(ilr)%locregcenter
          atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
          call calculate_multipoles(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
               denspot%dpbox%hgrids, phir(istr), charge_center_elec, tmb%lzd%llr(ilr)%locregcenter, &
               dipole_net(:,iiorb), quadropole_net(:,:,iiorb))
          !write(*,*) 'charge_center', charge_center_elec
          ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
          stop
      end if
    
    
      if (bigdft_mpi%nproc>1) then
          call mpiallred(dipole_net, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(quadropole_net, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      if (iproc==0) then
          call yaml_sequence_open('Support functions moments')
          do iorb=1,tmb%orbs%norb
              iat=tmb%orbs%onwhichatom(iorb)
              atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
              call yaml_sequence_open('number'//trim(yaml_toa(iorb))// &
                   ' (atom number ='//trim(yaml_toa(iat))//', type = '//trim(atomname)//')')
              call yaml_sequence(advance='no')
              call yaml_map('net dipole',dipole_net(:,iorb),fmt='(es18.10)')
              call yaml_sequence(advance='no')
              call yaml_sequence_open('net quadropole')
              do i=1,3
                 call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3,iorb),fmt='(es15.8)')))
              end do
              call yaml_sequence_close()
              call yaml_sequence_close()
          end do
          call yaml_sequence_close()
      end if
    
      call f_free(phir)
      call f_free(dipole_net)
      call f_free(quadropole_net)
      call f_release_routine()
    
    end subroutine support_function_multipoles


    subroutine calculate_multipoles(n1i, n2i, n3i, hgrids, phir, charge_center_elec, rxyz_center, &
               dipole_net, quadropole_net)
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: n1i, n2i, n3i
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(n1i*n2i*n3i),intent(in) :: phir
      real(kind=8),dimension(3),intent(in) :: charge_center_elec, rxyz_center
      real(kind=8),dimension(3),intent(out) :: dipole_net
      real(kind=8),dimension(3,3),intent(out) :: quadropole_net
    
      integer :: i1, i2, i3, jj, iz, iy, ix, ii, i, j
      real(kind=8) :: q, x, y, z, qtot, ri, rj, delta_term
      real(kind=8),dimension(3) :: dipole_center, dipole_el
      real(kind=8),dimension(3,3) :: quadropole_center, quadropole_el
    
    
      !!call yaml_map('rxyz_center',rxyz_center,fmt='(es16.6)')
      !!call yaml_map('charge_center_elec',charge_center_elec,fmt='(es16.6)')
      !!call yaml_map('sum phir',sum(phir),fmt='(es16.6)')
    
      ! Dipole and quadropole of the support function
      dipole_el=0.d0
      quadropole_el=0.d0
      qtot=0.d0
      jj=0
      do i3=1,n3i
          do i2=1,n2i
              do i1=1,n1i
                  jj=jj+1
                  ! z component of point jj
                  iz=jj/(n2i*n1i)
                  ! Subtract the 'lower' xy layers
                  ii=jj-iz*(n2i*n1i)
                  ! y component of point jj
                  iy=ii/n1i
                  ! Subtract the 'lower' y rows
                  ii=ii-iy*n1i
                  ! x component
                  ix=ii
    
                  ! Shift the values due to the convolutions bounds
                  ix=ix-14
                  iy=iy-14
                  iz=iz-14
    
                  q = phir(jj)**2 * product(hgrids)
                  x = ix*hgrids(1) + (rxyz_center(1)-charge_center_elec(1))
                  y = iy*hgrids(2) + (rxyz_center(2)-charge_center_elec(2))
                  z = iz*hgrids(3) + (rxyz_center(3)-charge_center_elec(3))
    
                  ! Dipole part
                  dipole_el(1) = dipole_el(1) + q*x
                  dipole_el(2) = dipole_el(2) + q*y
                  dipole_el(3) = dipole_el(3) + q*z
                  qtot=qtot+q
    
                  ! Quadrupole part
                  do i=1,3
                      ri=get_r(i, x, y, z)
                      do j=1,3
                          rj=get_r(j, x, y, z)
                          if (i==j) then
                              delta_term = x**2 + y**2 + z**2
                          else
                              delta_term=0.d0
                          end if
                          quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
                      end do
                  end do
              end do
          end do
      end do
    
      ! Dipole of the center
      dipole_center(1) = -qtot*rxyz_center(1)
      dipole_center(2) = -qtot*rxyz_center(2)
      dipole_center(3) = -qtot*rxyz_center(3)
    
      ! Quadropole of the center
      quadropole_center=0.d0
      do i=1,3
          ri=rxyz_center(i)
          do j=1,3
              rj=rxyz_center(j)
              if (i==j) then
                  delta_term = rxyz_center(1)**2 + rxyz_center(2)**2 + rxyz_center(3)**2
              else
                  delta_term=0.d0
              end if
              quadropole_center(j,i) = quadropole_center(j,i) -qtot*(3.d0*rj*ri-delta_term)
          end do
      end do
    
      ! Net dipole and quadropole
      dipole_net = dipole_el + dipole_center
      quadropole_net = quadropole_el + quadropole_center
    
    
    !  call yaml_sequence_open(trim(yaml_toa(it))//'('//trim(atomname)//')')
    !  !call yaml_map('qtot',qtot)
    !  call yaml_sequence(advance='no')
    !  !call yaml_map('center dipole',dipole_center,fmt='(es16.6)')
    !  !call yaml_map('electronic dipole',dipole_el,fmt='(es18.10)')
    !  call yaml_map('net dipole',dipole_net,fmt='(es18.10)')
    !  call yaml_sequence(advance='no')
    !  !call yaml_sequence_open('center quadropole')
    !  !do i=1,3
    !  !   call yaml_sequence(trim(yaml_toa(quadropole_center(i,1:3),fmt='(es15.8)')))
    !  !end do
    !  !call yaml_sequence_close()
    !  !call yaml_sequence_open('electronic quadropole')
    !  !do i=1,3
    !  !   call yaml_sequence(trim(yaml_toa(quadropole_el(i,1:3),fmt='(es15.8)')))
    !  !end do
    !  !call yaml_sequence_close()
    !  call yaml_sequence_open('net quadropole')
    !  do i=1,3
    !     call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es15.8)')))
    !  end do
    !  call yaml_sequence_close()
    !  call yaml_sequence_close()
    
      contains
    
        function get_r(i, x, y, z)
          integer,intent(in) :: i
          real(kind=8),intent(in) :: x, y, z
          real(kind=8) :: get_r
    
          select case (i)
          case (1)
              get_r=x
          case (2)
              get_r=y
          case (3)
              get_r=z
          case default
              stop 'wrong value of i'
          end select
        end function get_r
    
    end subroutine calculate_multipoles


    subroutine build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, order_taylor, &
               energy, energyDiff, energyold)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), SPARSE_FULL
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use yaml_output
      use rhopotential, only: updatePotential, sumrho_for_TMBs, corrections_for_negative_charge
      implicit none
      
      ! Calling arguments
      integer:: iproc, nproc
      type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers), intent(inout) :: GPU
      type(energy_terms),intent(inout) :: energs
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(input_variables),intent(in) :: input
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(out) :: energy, energyDiff
      real(kind=8), intent(inout) :: energyold
    
      ! Local variables
      type(orbitals_data) :: orbs
      type(comms_cubic) :: comms
      real(gp) :: fnrm
      logical :: rho_negative
      integer :: infoCoeff, nvctrp, npsidim_global
      real(kind=8),dimension(:),pointer :: phi_global, phiwork_global
      real(kind=8),dimension(:),allocatable :: tmparr
      character(len=*),parameter :: subname='build_ks_orbitals'
      real(wp), dimension(:,:,:), pointer :: mom_vec_fake
      type(work_mpiaccumulate) :: energs_work
      integer,dimension(:,:),allocatable :: ioffset_isf
    
    
      nullify(mom_vec_fake)
    
      energs_work = work_mpiaccumulate_null()
      energs_work%ncount = 4
      call allocate_work_mpiaccumulate(energs_work)
    
    
      !debug
      !integer :: iorb, jorb, ist, jst, ierr, i
      !real(kind=8) :: ddot, tt
    
    
      ! Get the expansion coefficients
      ! Start with a "clean" density, i.e. without legacy from previous mixing steps
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
    
      tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      call f_free(tmparr)
    
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
    
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
    
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor,input%lin%max_inversion_error,input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap, energs_work)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
    
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
      end if
    
    
      !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
      !     max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      !tmb%can_use_transposed=.false.
      !call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
      !energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      !energyDiff=energy-energyold
      !energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      ! Create communication arrays for support functions in the global box
      
      call nullify_orbitals_data(orbs)
      call copy_orbitals_data(tmb%orbs, orbs, subname)
      call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)
    
    
      ! Transform the support functions to the global box
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
                         tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
      phi_global = f_malloc_ptr(npsidim_global,id='phi_global')
      phiwork_global = f_malloc_ptr(npsidim_global,id='phiwork_global')
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
           tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
           KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
      call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global(1), phiwork_global(1))
    
    
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
      call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
                 tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
      
      call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global(1), phi_global(1))  
    
      call f_free_ptr(phi_global)
    
      !!ist=1
      !!do iorb=1,KSwfn%orbs%norbp
      !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!        write(800+iproc,*) iorb, i, phiwork_global(ist)
      !!        ist=ist+1
      !!    end do
      !!end do
    
    
      !!ierr=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!do i=1,KSwfn%orbs%norb*ierr
      !!    write(401,*) i, phiwork_global(i)
      !!end do
      !!write(*,*) 'GLOBAL DDOT',ddot(KSwfn%orbs%norb*ierr, phi_global, 1, phi_global, 1)
    
      !!do i=1,KSwfn%orbs%norb*ierr
      !!     tmb%psi(i)=phi_global(i)
      !!end do
      !!call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, tmb%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !!     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
    
      !!do i=1,KSwfn%orbs%norb*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)
      !!    write(600,'(i10,es16.7)') i, tmb%psi(i)
      !!end do
    
    
      !!write(*,*) 'iproc, input%output_wf_format',iproc, WF_FORMAT_PLAIN
      call writemywaves(iproc,trim(input%dir_output)//"wavefunction", WF_FORMAT_PLAIN, &
           KSwfn%orbs, KSwfn%Lzd%Glr%d%n1, KSwfn%Lzd%Glr%d%n2, KSwfn%Lzd%Glr%d%n3, &
           KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           at, rxyz, KSwfn%Lzd%Glr%wfd, phiwork_global)
    
      if (input%write_orbitals==2) then
          call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KSDens', &
               KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz)
      end if
    
      if (input%wf_extent_analysis) then
          ioffset_isf = f_malloc0((/3,KSwfn%orbs%norbp/),id='ioffset_isf')
          call analyze_wavefunctions('Kohn Sham orbitals extent analysis', 'global', &
               KSwfn%lzd, KSwfn%orbs, KSwfn%orbs%npsidim_orbs, phiwork_global, ioffset_isf)
          call f_free(ioffset_isf)
      end if
    
    
    
       call f_free_ptr(phiwork_global)
       call deallocate_orbitals_data(orbs)
       call deallocate_comms_cubic(comms)
    
      ! To get consistent values of the energy and the Kohn-Sham residue with those
      ! which will be calculated by the cubic restart.
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      call f_free(tmparr)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor, input%lin%max_inversion_error, input%purification_quickreturn, &
           input%calculate_KS_residue, input%calculate_gap, energs_work, updatekernel=.false.)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      energyDiff=energy-energyold
      energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      call deallocate_work_mpiaccumulate(energs_work)
    
    end subroutine build_ks_orbitals
    

    !TEMPORARY, to be cleaned/removed
    subroutine build_ks_orbitals_laura_tmp(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, order_taylor, &
               energy, energyDiff, energyold, npsidim_global, phiwork_global)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix
      use yaml_output
      use rhopotential, only: updatepotential, sumrho_for_TMBs, corrections_for_negative_charge
      implicit none
      
      ! Calling arguments
      integer:: iproc, nproc
      type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers), intent(inout) :: GPU
      type(energy_terms),intent(inout) :: energs
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(input_variables),intent(in) :: input
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(out) :: energy, energyDiff
      real(kind=8), intent(inout) :: energyold
    integer, intent(in) :: npsidim_global
    real(kind=8),dimension(:),pointer :: phiwork_global
    
      ! Local variables
      type(orbitals_data) :: orbs
      type(comms_cubic) :: comms
      real(gp) :: fnrm
      logical :: rho_negative
      integer :: infoCoeff, nvctrp
      real(kind=8),dimension(:),pointer :: phi_global
      real(kind=8),dimension(:,:),allocatable :: coeffs_occs
      character(len=*),parameter :: subname='build_ks_orbitals'
      real(wp), dimension(:,:,:), pointer :: mom_vec_fake
      integer :: iorb, itmb
      type(work_mpiaccumulate) :: energs_work
    
      nullify(mom_vec_fake)
    
      energs_work = work_mpiaccumulate_null()
      energs_work%ncount = 4
      call allocate_work_mpiaccumulate(energs_work)
    
      !debug
      !integer :: iorb, jorb, ist, jst, ierr, i
      !real(kind=8) :: ddot, tt
    
    
      ! Get the expansion coefficients
      ! Start with a "clean" density, i.e. without legacy from previous mixing steps
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
    
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
    
      tmb%can_use_transposed=.false.
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor,input%lin%max_inversion_error,input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap,energs_work)
    
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
      end if
    
    
      !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
      !     max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      !tmb%can_use_transposed=.false.
      !call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
      !energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      !energyDiff=energy-energyold
      !energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      ! Create communication arrays for support functions in the global box
      
      call nullify_orbitals_data(orbs)
      call copy_orbitals_data(tmb%orbs, orbs, subname)
      call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)
    
    
      ! Transform the support functions to the global box
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      !npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
      !                   tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
      phi_global = f_malloc_ptr(npsidim_global,id='phi_global')
      !phiwork_global = f_malloc_ptr(npsidim_global,id='phiwork_global')
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
           tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
           KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
      call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global(1), phiwork_global(1))
    
    !NOT PRINTING, 2xcorrect charge?!
    !print*,'ntmb,ntmbp,norb,norbp',tmb%orbs%norb,tmb%orbs%norbp,KSwfn%orbs%norb,KSwfn%orbs%norbp
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      coeffs_occs=f_malloc0((/tmb%orbs%norb,KSwfn%orbs%norb/),id='coeffs_occs')
      !call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, KSwfn%orbs%occup(1), nvctrp, tmb%coeff(1,1), &
      !     tmb%orbs%norb, 0.d0, coeffs_occs, nvctrp)
                 do iorb=1,KSwfn%orbs%norbp
                    do itmb=1,tmb%orbs%norb
                        coeffs_occs(itmb,iorb) = sqrt(KSwfn%orbs%occup(KSwfn%orbs%isorb+iorb))*tmb%coeff(itmb,KSwfn%orbs%isorb+iorb)
    !print*,KSwfn%orbs%occup(KSwfn%orbs%isorb+iorb)
                    end do
                 end do
              call mpiallred(coeffs_occs, mpi_sum, comm=bigdft_mpi%mpi_comm)
    
      nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
      !call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
      call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, coeffs_occs(1,1), &
                 tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
      call f_free(coeffs_occs)  
    
      call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global(1), phi_global(1))  
    
      call f_free_ptr(phi_global)
    
      !!ist=1
      !!do iorb=1,KSwfn%orbs%norbp
      !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!        write(800+iproc,*) iorb, i, phiwork_global(ist)
      !!        ist=ist+1
      !!    end do
      !!end do
    
    
      !!ierr=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!do i=1,KSwfn%orbs%norb*ierr
      !!    write(401,*) i, phiwork_global(i)
      !!end do
      !!write(*,*) 'GLOBAL DDOT',ddot(KSwfn%orbs%norb*ierr, phi_global, 1, phi_global, 1)
    
      !!do i=1,KSwfn%orbs%norb*ierr
      !!     tmb%psi(i)=phi_global(i)
      !!end do
      !!call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, tmb%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !!     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
    
      !!do i=1,KSwfn%orbs%norb*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)
      !!    write(600,'(i10,es16.7)') i, tmb%psi(i)
      !!end do
    
    
      !!write(*,*) 'iproc, input%output_wf_format',iproc, WF_FORMAT_PLAIN
      !call writemywaves(iproc,trim(input%dir_output)//"wavefunction", WF_FORMAT_PLAIN, &
      !     KSwfn%orbs, KSwfn%Lzd%Glr%d%n1, KSwfn%Lzd%Glr%d%n2, KSwfn%Lzd%Glr%d%n3, &
      !     KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     at, rxyz, KSwfn%Lzd%Glr%wfd, phiwork_global)
    
       !call f_free_ptr(phiwork_global)
       call deallocate_orbitals_data(orbs)
       call deallocate_comms_cubic(comms)
    
    if (.false.) then
      ! To get consistent values of the energy and the Kohn-Sham residue with those
      ! which will be calculated by the cubic restart.
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
      tmb%can_use_transposed=.false.
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor, input%lin%max_inversion_error, input%purification_quickreturn, &
           input%calculate_KS_residue, input%calculate_gap, energs_work, updatekernel=.false.)
      energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      energyDiff=energy-energyold
      energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    end if
      call allocate_work_mpiaccumulate(energs_work)
    end subroutine build_ks_orbitals_laura_tmp

end module postprocessing_linear
