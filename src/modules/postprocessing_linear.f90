module postprocessing_linear
  use public_enums
  implicit none

  private

  !> Public routines
  public :: loewdin_charge_analysis
  public :: loewdin_charge_analysis_core
  public :: support_function_multipoles
  public :: support_function_gross_multipoles
  public :: build_ks_orbitals
  public :: calculate_theta
  !public :: supportfunction_centers
!  public :: projector_for_charge_analysis
!  public :: extract_matrix !better private

  !> Public constants
  integer,parameter,public :: CHARGE_ANALYSIS_LOEWDIN   = 501
  integer,parameter,public :: CHARGE_ANALYSIS_MULLIKEN  = 502
  integer,parameter,public :: CHARGE_ANALYSIS_PROJECTOR = 503

  contains

    subroutine loewdin_charge_analysis(iproc,tmb,atoms,denspot, &
               calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,blocksize,&
               ntheta, istheta, theta)
      use module_base
      use module_types
      !use module_interfaces
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
      integer,intent(in) :: iproc, blocksize
      type(dft_wavefunction),intent(inout) :: tmb
      type(atoms_data),intent(in) :: atoms
      type(DFT_local_fields), intent(inout) :: denspot
      logical,intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
      integer,intent(in) :: meth_overlap
      integer,intent(in),optional :: ntheta, istheta
      real(kind=8),dimension(:,:),intent(in),optional :: theta ! must have dimension (atoms%astruct%nat,ndim_theta)
    
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
      logical :: optionals_present
    
    
      ! needs parallelizing/converting to sparse
      ! re-use overlap matrix if possible either before or after
    
      call f_routine(id='loewdin_charge_analysis')


      if (present(theta)) then
          if (.not.present(ntheta)) then
              call f_err_throw('ntheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else if (.not.present(istheta)) then
              call f_err_throw('istheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else
              optionals_present = .true.
          end if
          if (size(theta,1)/=atoms%astruct%nat) then
              call f_err_throw('wrong first dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(theta,2)/=ntheta) then
              call f_err_throw('wrong second dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      else
          optionals_present = .false.
      end if

    
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

      !call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, &
      !         tmb%orbs%norb_par, tmb%orbs%isorb_par, meth_overlap, tmb%linmat%s, tmb%linmat%l, atoms, &
      !         tmb%linmat%kernel_, tmb%linmat%ovrlp_)
      if (optionals_present) then
          call loewdin_charge_analysis_core(CHARGE_ANALYSIS_LOEWDIN, bigdft_mpi%iproc, bigdft_mpi%nproc, &
               tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, tmb%linmat%s%isfvctr, &
               tmb%linmat%s%nfvctr_par, tmb%linmat%s%isfvctr_par, &
               meth_overlap, blocksize, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_, &
               ntheta=ntheta, istheta=istheta, theta=theta)
      else
          call loewdin_charge_analysis_core(CHARGE_ANALYSIS_LOEWDIN, bigdft_mpi%iproc, bigdft_mpi%nproc, &
               tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, tmb%linmat%s%isfvctr, &
               tmb%linmat%s%nfvctr_par, tmb%linmat%s%isfvctr_par, &
               meth_overlap, blocksize, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_)
      end if
    
      call f_release_routine()
    
    
    end subroutine loewdin_charge_analysis




    subroutine loewdin_charge_analysis_core(method, iproc, nproc, norb, norbp, isorb, &
            norb_par, isorb_par, meth_overlap, blocksize, smats, smatl, atoms, kernel, ovrlp, &
               ntheta, istheta, theta)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   DENSE_FULL, SPARSE_TASKGROUP, SPARSE_FULL, &
                                   deallocate_matrices, matrices_null, allocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: uncompress_matrix2, matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups, &
                              transform_sparse_matrix
      use matrix_operations, only: overlapPowerGeneral
      use io, only: write_partial_charges
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: method, iproc, nproc, norb, norbp, isorb, meth_overlap, blocksize
      integer,dimension(0:nproc-1),intent(in) :: norb_par, isorb_par
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(atoms_data),intent(in) :: atoms
      type(matrices),intent(inout) :: kernel
      type(matrices),intent(inout) :: ovrlp
      integer,intent(in),optional :: ntheta, istheta
      real(kind=8),dimension(:,:),intent(in),optional :: theta !dimension (atoms%astruct%nat,ntheta)

      ! Local variables
      integer :: ierr, iorb, iat, ind, ist, ishift, ispin, iiorb
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:,:,:),allocatable :: proj_mat
      real(kind=8),dimension(:,:),allocatable :: weight_matrix, weight_matrixp, proj_ovrlp_half
      real(kind=8),dimension(:),allocatable :: charge_per_atom, proj_ovrlp_half_compr
      real(kind=8),dimension(:),allocatable :: weight_matrix_compr_tg, weight_matrix_compr
      real(kind=8) :: mean_error, max_error
      integer,parameter :: DENSE=101, SPARSE=102
      integer :: imode=SPARSE
      logical :: optionals_present

      call f_routine(id='loewdin_charge_analysis_core')
      
      ! Check the arguments
      select case (method)
      case(CHARGE_ANALYSIS_LOEWDIN)
          if (iproc==0) call yaml_map('Method','Loewdin')
      case(CHARGE_ANALYSIS_MULLIKEN)
          if (iproc==0) call yaml_map('Method','Mulliken')
      case default
          call f_err_throw('Wrong Method',err_name='BIGDFT_RUNTIME_ERROR')
      end select

      if (present(theta)) then
          if (.not.present(ntheta)) then
              call f_err_throw('ntheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else if (.not.present(istheta)) then
              call f_err_throw('istheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else
              optionals_present = .true.
          end if
          if (size(theta,1)/=atoms%astruct%nat) then
              call f_err_throw('wrong first dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(theta,2)/=ntheta) then
              call f_err_throw('wrong second dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      else
          optionals_present = .false.
      end if


      if (imode==DENSE) then

          call f_err_throw('Dense mode is deprecated',err_name='BIGDT_RUNTIME_ERROR')

          !!!inv_ovrlp(1) = matrices_null()
          !!!call allocate_matrices(smatl, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

          !!!ovrlp%matrix = sparsematrix_malloc_ptr(smats, iaction=DENSE_FULL, id='ovrlp%matrix')
          !!!call uncompress_matrix2(iproc, nproc, smats, &
          !!!     ovrlp%matrix_compr, ovrlp%matrix)
          !!!call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), -1, &
          !!!     imode=2, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
          !!!     ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
          !!!     max_error=max_error, mean_error=mean_error)
          !!!call f_free_ptr(ovrlp%matrix)
    
          !!!! optimize this to just change the matrix multiplication?
          !!!proj_mat = sparsematrix_malloc0(smatl,iaction=DENSE_FULL,id='proj_mat')
    
          !!!call uncompress_matrix2(iproc, nproc, smatl, kernel%matrix_compr, proj_mat)

          !!!proj_ovrlp_half=f_malloc((/norb,norbp/),id='proj_ovrlp_half')
          !!!if (norbp>0) then
          !!!   call dgemm('n', 'n', norb, norbp, &
          !!!          norb, 1.d0, &
          !!!          proj_mat(1,1,1), norb, &
          !!!          inv_ovrlp(1)%matrix(1,isorb+1,1), norb, 0.d0, &
          !!!          proj_ovrlp_half(1,1), norb)
          !!!end if
          !!!call f_free(proj_mat)
          !!!weight_matrixp=f_malloc((/norb,norbp/), id='weight_matrixp')
          !!!if (norbp>0) then
          !!!   call dgemm('n', 'n', norb, norbp, &
          !!!        norb, 1.d0, &
          !!!        inv_ovrlp(1)%matrix(1,1,1), norb, &
          !!!        proj_ovrlp_half(1,1), norb, 0.d0, &
          !!!        weight_matrixp(1,1), norb)
          !!!end if
          !!!call f_free(proj_ovrlp_half)
          !!!weight_matrix=f_malloc((/norb,norb/), id='weight_matrix')
          !!!if (nproc>1) then
          !!!   call mpi_allgatherv(weight_matrixp, norb*norbp, mpi_double_precision, weight_matrix, &
          !!!        norb*norb_par(:), norb*isorb_par, &
          !!!        mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          !!!else
          !!!   call vcopy(norb*norb,weight_matrixp(1,1),1,weight_matrix(1,1),1)
          !!!end if
          !!!call f_free(weight_matrixp)
    
          !!!charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')
    
          !!!do iorb=1,norb
          !!!    iat=smats%on_which_atom(iorb)
          !!!    charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(iorb,iorb)
          !!!end do
          !!!if (iproc==0) then
          !!!    call write_partial_charges(atoms, charge_per_atom, .true.)
          !!!    call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
          !!!    call calculate_dipole(iproc, atoms, charge_per_atom)
          !!!    call calculate_quadropole(iproc, atoms, charge_per_atom)
          !!!    call yaml_sequence_close()
          !!!end if
          !!!!!call support_function_multipoles()
    
          !!!call deallocate_matrices(inv_ovrlp(1))
          !!!call f_free(charge_per_atom)
          !!!call f_free(weight_matrix)

      else if (imode==SPARSE) then


          weight_matrix_compr_tg = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='weight_matrix_compr_tg')
          proj_ovrlp_half_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='proj_ovrlp_half_compr')

          if (method==CHARGE_ANALYSIS_LOEWDIN) then
              inv_ovrlp(1) = matrices_null()
              inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='inv_ovrlp(1)%matrix_compr')

              call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), blocksize, &
                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
                   ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
                   max_error=max_error, mean_error=mean_error)
              !call f_free_ptr(ovrlp%matrix)

              do ispin=1,smatl%nspin
                  ist = (ispin-1)*smatl%nvctrp_tg + 1
                  !not sure where exactly the problem is but this can't be called for only some mpi procs otherwise the code hangs
                  !if (norbp>0) then
                     call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                          kernel%matrix_compr(ist:), inv_ovrlp(1)%matrix_compr(ist:), proj_ovrlp_half_compr(ist:))
                  !end if

                  !if (norbp>0) then
                     call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                          inv_ovrlp(1)%matrix_compr(ist:), proj_ovrlp_half_compr(ist:), weight_matrix_compr_tg(ist:))
                  !end if
              end do
              call deallocate_matrices(inv_ovrlp(1))
          else if (method==CHARGE_ANALYSIS_MULLIKEN) then
              call transform_sparse_matrix(smats, smatl, ovrlp%matrix_compr, proj_ovrlp_half_compr, 'small_to_large')
              do ispin=1,smatl%nspin
                  ist = (ispin-1)*smatl%nvctrp_tg + 1
                  !if (norbp>0) then
                     call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                          kernel%matrix_compr(ist:), proj_ovrlp_half_compr(ist:), weight_matrix_compr_tg(ist:))
                  !end if
              end do
          end if
          call f_free(proj_ovrlp_half_compr)

          charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')

          ! Maybe this can be improved... not really necessary to gather the entire matrix
          weight_matrix_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='weight_matrix_compr')
          call gather_matrix_from_taskgroups(iproc, nproc, smatl, weight_matrix_compr_tg, weight_matrix_compr)

          if (optionals_present) then
              do ispin=1,smatl%nspin
                  ishift = (ispin-1)*smatl%nvctr
                  do iorb=1,ntheta
                      iiorb = iorb + istheta
                      iiorb = modulo(iiorb-1,smatl%nfvctr)+1
                      ind = matrixindex_in_compressed(smatl, iorb, iorb)
                      !if (iproc==0) write(*,*) 'iorb, trace charge', iorb, weight_matrix_compr(ind)
                      do iat=1,atoms%astruct%nat
                          ind = ind + ishift
                          charge_per_atom(iat) = charge_per_atom(iat) + theta(iat,iorb)*weight_matrix_compr(ind)
                       end do
                  end do
              end do
          else 
              do ispin=1,smatl%nspin
                  ishift = (ispin-1)*smatl%nvctr
                  do iorb=1,norb
                      iiorb = modulo(iorb-1,smatl%nfvctr)+1
                      iat=smats%on_which_atom(iiorb)
                      ind = matrixindex_in_compressed(smatl, iorb, iorb)
                      ind = ind + ishift
                      !if (iproc==0) write(*,*) 'iorb, trace charge', iorb, weight_matrix_compr(ind)
                      charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix_compr(ind)
                  end do
              end do
          end if

          if (iproc==0) then
              call write_partial_charges(atoms, charge_per_atom, .true.)
              call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
              call calculate_dipole(iproc, atoms, charge_per_atom)
              call calculate_quadropole(iproc, atoms, charge_per_atom)
              call yaml_sequence_close()
          end if
    
          call f_free(charge_per_atom)
          call f_free(weight_matrix_compr_tg)
          call f_free(weight_matrix_compr)

      else
          call f_err_throw('wrong value for imode',err_name='BIGDFT_RUNTIME_ERROR')
      end if

      call f_release_routine()

    end subroutine loewdin_charge_analysis_core




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
      use locreg_operations
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
          call initialize_work_arrays_sumrho(1,[tmb%lzd%Llr(ilr)],.true.,w)
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

      !call yaml_map('rxyz_center',rxyz_center)
      !call yaml_map('charge_center_elec',charge_center_elec)
      !call yaml_map('qtot',qtot)
    
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
      !call yaml_map('dipole_el',dipole_el)
      !call yaml_map('dipole_center',dipole_center)
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
      use module_interfaces, only: get_coeff, write_eigenvalues_data, write_orbital_density
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), SPARSE_FULL
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use yaml_output
      use rhopotential, only: updatePotential, sumrho_for_TMBs, corrections_for_negative_charge
      use locregs_init, only: small_to_large_locreg
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
           input%calculate_KS_residue,input%calculate_gap, energs_work, .false., input%lin%coeff_factor, &
           input%lin%pexsi_npoles, input%lin%pexsi_mumin, input%lin%pexsi_mumax, input%lin%pexsi_mu, &
           input%lin%pexsi_temperature, input%lin%pexsi_tol_charge)
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
               KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .true.)
      else if (input%write_orbitals==3) then 
          call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KS', &
               KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .false.)
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
           input%calculate_KS_residue, input%calculate_gap, energs_work, .false., input%lin%coeff_factor, &
           input%lin%pexsi_npoles, input%lin%pexsi_mumin, input%lin%pexsi_mumax, input%lin%pexsi_mu, &
           input%lin%pexsi_temperature, input%lin%pexsi_tol_charge, updatekernel=.false.)
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
      use module_interfaces, only: get_coeff, write_eigenvalues_data
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix
      use yaml_output
      use rhopotential, only: updatepotential, sumrho_for_TMBs, corrections_for_negative_charge
      use locregs_init, only: small_to_large_locreg
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
           input%calculate_KS_residue,input%calculate_gap,energs_work, .false., input%lin%coeff_factor, &
           input%lin%pexsi_npoles, input%lin%pexsi_mumin, input%lin%pexsi_mumax, input%lin%pexsi_mu, &
           input%lin%pexsi_temperature, input%lin%pexsi_tol_charge)
    
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
           input%calculate_KS_residue, input%calculate_gap, energs_work, .false., input%lin%coeff_factor, &
           input%lin%pexsi_npoles, input%lin%pexsi_mumin, input%lin%pexsi_mumax, input%lin%pexsi_mu, &
           input%lin%pexsi_temperature, input%lin%pexsi_tol_charge, updatekernel=.false.)
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



    subroutine calculate_theta(nat, rxyz, nphidim, phi, nphirdim, orbs, lzd, theta)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors
      use locreg_operations
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: nat, nphidim, nphirdim
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      real(kind=8),dimension(nphidim),intent(in) :: phi
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(nat,orbs%norbp),intent(out) :: theta

      ! Local variables
      real(kind=8),dimension(:),allocatable :: psir
      type(workarr_sumrho) :: w
      integer :: ist, istr, iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, iat, iiat, l, m
      real(kind=8),dimension(3) :: com
      real(kind=8),dimension(-1:1) :: dipole
      real(kind=8) :: weight, tt, x, y, z, r2, hxh, hyh, hzh, q, qtot, monopole, r
      real(kind=8),parameter :: sigma2=0.1d0

      call f_zero(theta)

      ! Transform the support functions to real space
      psir = f_malloc(max(nphirdim,1),id='psir')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,[lzd%Llr(ilr)],.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=nphirdim+1) then
          call f_err_throw('ERROR on process '//adjustl(trim(yaml_toa(bigdft_mpi%iproc)))//': istr/=nphirdim+1', &
               err_name='BIGDFT_RUNTIME_ERROR')
          stop
      end if

      hxh = 0.5d0*lzd%hgrids(1)
      hyh = 0.5d0*lzd%hgrids(2)
      hzh = 0.5d0*lzd%hgrids(3)


      !! METHOD 1 ####################################################
      !istr = 1
      !do iorb=1,orbs%norbp
      !    iiorb=orbs%isorb+iorb
      !    ilr=orbs%inwhichlocreg(iiorb)
      !    write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
      !    com(1:3) = 0.d0
      !    weight = 0.d0
      !    do i3=1,lzd%llr(ilr)%d%n3i
      !        ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
      !        z = ii3*hzh
      !        do i2=1,lzd%llr(ilr)%d%n2i
      !            ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
      !            y = ii2*hyh
      !            do i1=1,lzd%llr(ilr)%d%n1i
      !                ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
      !                x = ii1*hxh
      !                tt = psir(istr)**2
      !                com(1) = com(1) + x*tt
      !                com(2) = com(2) + y*tt
      !                com(3) = com(3) + z*tt
      !                weight = weight + tt
      !                istr = istr + 1
      !            end do
      !        end do
      !    end do
      !    call yaml_map('weight',weight)

      !    weight = 0.d0
      !    do iat=1,nat
      !        write(*,*) 'com, rxzy', com, rxyz(:,iat)
      !        r2 = (rxyz(1,iat)-com(1))**2 + (rxyz(2,iat)-com(2))**2 + (rxyz(3,iat)-com(3))**2
      !        tt = exp(-r2/(2*sigma2))
      !        theta(iat,iorb) = tt
      !        weight = weight + tt
      !    end do
      !    call yaml_map('weight',weight)
      !    tt = 1.d0/weight
      !    call dscal(nat, tt, theta(1,iorb), 1)

      !    call yaml_map('theta '//adjustl(trim(yaml_toa(iorb))),theta(:,iorb))

      !end do
      !! END METHOD 1 ################################################


      ! METHOD 2 ####################################################
      istr = 1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          iiat = orbs%onwhichatom(iiorb)
          !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
          com(1:3) = 0.d0
          dipole(:) = 0.d0
          weight = 0.d0
          qtot = 0.d0
          do i3=1,lzd%llr(ilr)%d%n3i
              ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
              z = ii3*hzh
              !write(*,*) 'i3, ii3, z, d', i3, ii3, z, z-rxyz(3,iiat)
              do i2=1,lzd%llr(ilr)%d%n2i
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
                  y = ii2*hyh
                  do i1=1,lzd%llr(ilr)%d%n1i
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
                      x = ii1*hxh
                      q = psir(istr)**2!*hxh*hyh*hzh
                      com(1) = com(1) + x*q
                      com(2) = com(2) + y*q
                      com(3) = com(3) + z*q
                      qtot = qtot + q
                      !r2 = (rxyz(1,iiat)-x)**2 + (rxyz(2,iiat)-y)**2 + (rxyz(3,iiat)-z)**2
                      !tt = 1.d0/r2**4
                      !!r2 = (rxyz(1,iiat)-x)**2 + (rxyz(2,iiat)-y)**2 + (rxyz(3,iiat)-z)**2
                      !!if (r2/(2.d0*0.2d0)>300) then
                      !!    tt = 0.d0
                      !!else
                      !!    tt = safe_exp(-r2/(2.d0*0.2d0))
                      !!end if
                      !!if (tt*psir(istr)**2<1.d-300) then
                      !!    q =0.d0
                      !!else
                      !!    q = tt*psir(istr)**2
                      !!end if
                      do iat=1,nat
                          !write(*,*) 'com, rxzy', com, rxyz(:,iat)
                          !r2 = (rxyz(1,iat)-com(1))**2 + (rxyz(2,iat)-com(2))**2 + (rxyz(3,iat)-com(3))**2
                          r2 = (rxyz(1,iat)-x)**2 + (rxyz(2,iat)-y)**2 + (rxyz(3,iat)-z)**2
                          if (r2/(2*sigma2)>300) then
                              tt = 0.d0
                          else
                              theta(iat,iorb) = theta(iat,iorb) + q*safe_exp(-r2/(2*sigma2))
                          end if
                      end do
                      istr = istr + 1
                  end do
              end do
          end do
          dipole(-1) = com(2)
          dipole( 0) = com(3)
          dipole( 1) = com(1)
          com(1:3) = com(1:3)/qtot
          monopole = qtot - 1.d0
          dipole(:) = dipole(:) - qtot*rxyz(1:3,iiat)
          call yaml_map('monopole',monopole)
          call yaml_map('dipole',dipole)
          call yaml_map('qtot', qtot)
          call yaml_map('rxyz(..,iiat)', rxyz(1:3,iiat))
          call yaml_map('com', com)
          !if (iiat==1 .and. iiorb==5) then
          !    theta(:,iorb) = 0.d0
          !    theta(iiat,iorb) = 1.d0
          !end if

          do iat=1,nat
              theta(iat,iorb) = 0.d0
              !x = rxyz(1,iat) - rxyz(1,iiat)
              !y = rxyz(2,iat) - rxyz(2,iiat)
              !z = rxyz(3,iat) - rxyz(3,iiat)
              x = rxyz(1,iat) - com(1)
              y = rxyz(2,iat) - com(2)
              z = rxyz(3,iat) - com(3)
              r = sqrt( x**2 + y**2 + z**2 )
              do l=0,0!1
                  do m=-l,l
                      tt = spherical_harmonic(l, m, x, y, z)*gaussian(1.d0, r)
                      if(l==1) tt = tt * dipole(m)
                      write(*,*) 'iorb, iat, m, r, rxyz, tt', iorb, iat, m, r, rxyz(:,iat), tt
                      theta(iat,iorb) = theta(iat,iorb) + tt**2
                  end do
              end do
          end do

          tt = sum(theta(1:nat,iorb))
          !write(*,*) 'tt',tt
          tt = 1.d0/tt
          call dscal(nat, tt, theta(1,iorb), 1)

          call yaml_map('theta '//adjustl(trim(yaml_toa(iorb))),theta(:,iorb))

      end do
      ! END METHOD 2 ################################################

    end subroutine calculate_theta


    function gaussian(sigma, r) result(g)
      use module_base
      implicit none
      ! Calling arguments
      real(kind=8),intent(in) :: sigma, r
      real(kind=8) :: g
      
      g = exp(-r**2/(2.d0*sigma**2))
      g = g/sqrt(2.d0*pi*sigma**2)**3
      !g = g/(sigma**3*sqrt(2.d0*pi)**3)
    end function gaussian

    !> Calculates the real spherical harmonic for given values of l, m, x, y, z.
    function spherical_harmonic(l, m, x, y, z) result(sh)
      !use module_base, only: pi => pi_param
      use module_base
      implicit none
      ! Calling arguments
      integer,intent(in) :: l, m
      real(kind=8),intent(in) :: x, y, z
      real(kind=8) :: sh

      ! Local variables
      integer,parameter :: l_max=2
      real(kind=8) :: r, r2, rnorm

      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
      if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
          err_name='BIGDFT_RUNTIME_ERROR')
      if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


      ! Normalization for a sphere of radius rmax
      select case (l)
      case (0)
          sh = 0.5d0*sqrt(1/pi)
      case (1)
          r = sqrt(x**2+y**2+z**2)
          ! fix for small r (needs proper handling later...)
          if (r==0.d0) r=1.d-20
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
          ! fix for small r2 (needs proper handling later...)
          if (r2==0.d0) r2=1.d-20
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





 !> SM: similar to support_function_multipoles. This one calculates the "gross" multipoles (i.e. without taking into account the "core" contribution)
 subroutine support_function_gross_multipoles(iproc, tmb, atoms, denspot, multipoles)
   use module_base
   use module_types
   use locreg_operations
   use yaml_output
   use multipole_base, only: lmax
   use multipole, only: multipole_analysis_core, write_multipoles_new
   use bounds, only: geocode_buffers
   
   ! Calling arguments
   integer,intent(in) :: iproc
   type(DFT_wavefunction),intent(in) :: tmb
   type(atoms_data),intent(in) :: atoms
   type(DFT_local_fields), intent(inout) :: denspot
   real(kind=8),dimension(-lmax:lmax,0:lmax,tmb%orbs%norb),optional :: multipoles
 
   integer :: ist, istr, iorb, iiorb, ilr, i, iat, iter, itype, i1, i2, i3, ii1, ii2, ii3, nl1, nl2, nl3, ii
   real(kind=8),dimension(:),allocatable :: rmax
   real(kind=8),dimension(:,:),allocatable :: centers
   real(kind=8),dimension(3) :: charge_center_elec
   real(kind=8),dimension(:),allocatable :: phir, phir_one
   real(kind=8) :: hxh, hyh, hzh, tt, x, y, z, weight
   type(workarr_sumrho) :: w
   character(len=20) :: atomname
   integer,dimension(:),allocatable :: iatype_tmp
 
   call f_routine(id='support_function_multipoles')
 
   rmax = f_malloc0(tmb%orbs%norb,id='rmax')
   phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
   phir_one = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir_one')
   phir_one = 1.d0
 
   !call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
   !call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))


   iter_loop: do iter=1,1

       multipoles = 0.d0
 
       ist=1
       istr=1
       do iorb=1,tmb%orbs%norbp
           iiorb=tmb%orbs%isorb+iorb
           ilr=tmb%orbs%inwhichlocreg(iiorb)
           iat=tmb%orbs%onwhichatom(iiorb)
           call initialize_work_arrays_sumrho(1,[tmb%lzd%Llr(ilr)],.true.,w)
           ! Transform the support function to real space
           call daub_to_isf(tmb%lzd%llr(ilr), w, tmb%psi(ist), phir(istr))
           call deallocate_work_arrays_sumrho(w)

           ! NEW: CALCULATE THE WEIGHT CENTER OF THE SUPPORT FUNCTION ############################
           hxh = 0.5d0*tmb%lzd%hgrids(1)
           hyh = 0.5d0*tmb%lzd%hgrids(2)
           hzh = 0.5d0*tmb%lzd%hgrids(3)
           ii = istr
           call geocode_buffers(tmb%lzd%Llr(ilr)%geocode, tmb%lzd%glr%geocode, nl1, nl2, nl3)
           !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
           !com(1:3,iorb) = 0.d0
           weight = 0.d0
           charge_center_elec(:) = 0.d0
           do i3=1,tmb%lzd%llr(ilr)%d%n3i
               ii3 = tmb%lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
               z = ii3*hzh
               do i2=1,tmb%lzd%llr(ilr)%d%n2i
                   ii2 = tmb%lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
                   y = ii2*hyh
                   do i1=1,tmb%lzd%llr(ilr)%d%n1i
                       ii1 = tmb%lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                       x = ii1*hxh
                       tt = phir(ii)**2
                       charge_center_elec(1) = charge_center_elec(1) + x*tt
                       charge_center_elec(2) = charge_center_elec(2) + y*tt
                       charge_center_elec(3) = charge_center_elec(3) + z*tt
                       weight = weight + tt
                       ii = ii + 1
                   end do
               end do
           end do
           !call yaml_map('weight',weight)
           charge_center_elec(1:3) = charge_center_elec(1:3)/weight
           !write(*,'(a,i4,3es13.4)') 'iorb, charge_center_elec(1:3)', iorb, charge_center_elec(1:3)
           ! ######################################################################################



           ! Calculate the charge center
           !call charge_center(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
           !     denspot%dpbox%hgrids, phir(istr), charge_center_elec)
           !write(*,*) 'ilr, tmb%lzd%llr(ilr)%locregcenter', iat, tmb%lzd%llr(ilr)%locregcenter
           atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
           ! OLD VERSION #######################################################
           if (iter==2) then
               write(*,*) 'call to calculate_multipoles_of_one_supportfunction commented'
               !!call calculate_multipoles_of_one_supportfunction(tmb%lzd%llr(ilr)%d%n1i, &
               !!     tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
               !!     denspot%dpbox%hgrids, phir(istr), tmb%lzd%llr(ilr)%locregcenter, &
               !!     multipoles(:,:,iiorb))
           ! NEW VERSION #######################################################
           else if (iter==1) then
               !write(*,*) 'SUM(phir**2)',SUM(phir**2)
               !write(*,*) 'SUM(phir)',SUM(phir)
               !write(*,*) 'SUM(phir)*sqrt(product(tmb%lzd%hgrids))',SUM(phir)*sqrt(product(tmb%lzd%hgrids))
               call multipole_analysis_core(0, 1, 1, 0, 1, 1, 1, &
                    -1,  (/1/), (/1/), (/1/), (/1/), &
                    (/tmb%lzd%Llr(ilr)%d%n1i/), (/tmb%lzd%Llr(ilr)%d%n2i/), (/tmb%lzd%Llr(ilr)%d%n3i/), &
                    (/tmb%lzd%Llr(ilr)%nsi1/), (/tmb%lzd%Llr(ilr)%nsi2/), (/tmb%lzd%Llr(ilr)%nsi3/), rmax(iiorb), &
                    tmb%lzd%hgrids, tmb%lzd%llr(ilr)%locregcenter, &
                    tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i, &
                    phir(istr), phir_one(istr), &
                    1, (/1.d0/), 0, multipoles(:,:,iiorb), rmax(iiorb), 102, matrixindex=(/1/))
           end if
           !write(*,*) 'charge_center', charge_center_elec
           ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
           istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
       end do
       if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
           write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
           stop
       end if
 
 
       if (bigdft_mpi%nproc>1) then
           call mpiallred(multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
           call mpiallred(rmax, mpi_sum, comm=bigdft_mpi%mpi_comm)
       end if
 
       if (iproc==0) then
           call yaml_sequence_open('Gross support functions moments (iter='//trim(adjustl(yaml_toa(iter)))//')')
           !do iorb=1,tmb%orbs%norb
           !    iat=tmb%orbs%onwhichatom(iorb)
           !    atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
           !    call yaml_sequence_open('number'//trim(yaml_toa(iorb))// &
           !         ' (atom number ='//trim(yaml_toa(iat))//', type = '//trim(atomname)//')')
           !    call yaml_sequence(advance='no')
           !    call yaml_map('monopole',multipoles(0,0,iorb),fmt='(es18.10)')
           !    call yaml_map('dinopole',multipoles(-1:1,1,iorb),fmt='(es18.10)')
           !    !call yaml_sequence(advance='no')
           !    !call yaml_sequence_open('net quadropole')
           !    !do i=1,3
           !    !   call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3,iorb),fmt='(es15.8)')))
           !    !end do
           !    !call yaml_sequence_close()
           !    call yaml_sequence_close()
           !end do
           !do iorb=1,tmb%orbs%norb
           !    call write_multipoles_new(1, 1, (/1/), (/'testatom'/), (/0.d0,0.d0,0.d0/), 'letssee', &
           !         multipoles(:,:,iorb), rmax, tmb%lzd%hgrids, without_normalization=.true.)
           !end do
           iatype_tmp = f_malloc(tmb%orbs%norb,id='iatype_tmp')
           centers = f_malloc((/3,tmb%orbs%norb/),id='centers')
           do iorb=1,tmb%orbs%norb
               iat = tmb%orbs%onwhichatom(iorb)
               ilr = tmb%orbs%inwhichlocreg(iorb)
               itype = atoms%astruct%iatype(iat)
               iatype_tmp(iorb) = itype
               centers(1:3,iorb) = tmb%lzd%llr(ilr)%locregcenter(1:3)
           end do
           call write_multipoles_new(tmb%orbs%norb, atoms%astruct%ntypes, iatype_tmp, &
                atoms%astruct%atomnames, centers, atoms%astruct%units, &
                multipoles)
           call f_free(centers)
           call f_free(iatype_tmp)
           call yaml_sequence_close()
       end if

   end do iter_loop
 
   call f_free(rmax)
   call f_free(phir)
   call f_free(phir_one)
   call f_release_routine()
 
 end subroutine support_function_gross_multipoles


 !! THIS IS PROBABLY WRONG
 !> SM: similar to calculate_multipoles. This one calculates the "gross" multipoles of one 
 !! support function (i.e. without taking !into account the "core" contribution
 subroutine calculate_multipoles_of_one_supportfunction(n1i, n2i, n3i, hgrids, phir, rxyz_center, &
            multipoles)
   use yaml_output
   use module_base
   use multipole_base, only: lmax
   implicit none
   ! Calling arguments
   integer,intent(in) :: n1i, n2i, n3i
   real(kind=8),dimension(3),intent(in) :: hgrids
   real(kind=8),dimension(n1i*n2i*n3i),intent(in) :: phir
   real(kind=8),dimension(3),intent(in) :: rxyz_center
   real(kind=8),dimension(-lmax:lmax,0:lmax),intent(out) :: multipoles
 
   integer :: i1, i2, i3, jj, iz, iy, ix, ii, i, j
   real(kind=8) :: q, x, y, z, qtot, ri, rj, delta_term
   real(kind=8),dimension(3) :: dipole_center, dipole_el
   real(kind=8),dimension(3,3) :: quadropole_center, quadropole_el
 
   stop 'MOST LIKELY BUGGY (DIPOLE SEEMS TO BE WRONG)'
 
   !!call yaml_map('rxyz_center',rxyz_center,fmt='(es16.6)')
   !!call yaml_map('charge_center_elec',charge_center_elec,fmt='(es16.6)')
   !!call yaml_map('sum phir',sum(phir),fmt='(es16.6)')
 
   dipole_el=0.d0
   quadropole_el=0.d0

   call f_zero(multipoles)
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
 
               q = phir(jj)!**2 !* product(hgrids)
               x = (ix*hgrids(1) - rxyz_center(1))
               y = (iy*hgrids(2) - rxyz_center(2))
               z = (iz*hgrids(3) - rxyz_center(3))

               ! Monopole part
               multipoles(0,0) = multipoles(0,0) + q
 
               ! Dipole part
               multipoles(-1,1) = multipoles(-1,1) + q*y
               multipoles( 0,1) = multipoles( 0,1) + q*z
               multipoles( 1,1) = multipoles( 1,1) + q*x

 
               !!! Quadrupole part not ready yet
               !!do i=1,3
               !!    ri=get_r(i, x, y, z)
               !!    do j=1,3
               !!        rj=get_r(j, x, y, z)
               !!        if (i==j) then
               !!            delta_term = x**2 + y**2 + z**2
               !!        else
               !!            delta_term=0.d0
               !!        end if
               !!        quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
               !!    end do
               !!end do
           end do
       end do
   end do
 
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
 
 end subroutine calculate_multipoles_of_one_supportfunction

end module postprocessing_linear
