module postprocessing_linear
  use public_enums
  implicit none

  private

  !> Public routines
  public :: loewdin_charge_analysis
  public :: loewdin_charge_analysis_core
  public :: build_ks_orbitals

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
               meth_overlap, blocksize, tmb%linmat%smmd, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_, &
               ntheta=ntheta, istheta=istheta, theta=theta)
      else
          call loewdin_charge_analysis_core(CHARGE_ANALYSIS_LOEWDIN, bigdft_mpi%iproc, bigdft_mpi%nproc, &
               tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, tmb%linmat%s%isfvctr, &
               tmb%linmat%s%nfvctr_par, tmb%linmat%s%isfvctr_par, &
               meth_overlap, blocksize, tmb%linmat%smmd, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_)
      end if
    
      call f_release_routine()
    
    
    end subroutine loewdin_charge_analysis




    subroutine loewdin_charge_analysis_core(method, iproc, nproc, norb, norbp, isorb, &
               norb_par, isorb_par, meth_overlap, blocksize, smmd, smats, smatl, atoms, kernel, ovrlp, &
               ntheta, istheta, theta)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   DENSE_FULL, SPARSE_TASKGROUP, SPARSE_FULL, &
                                   deallocate_matrices, matrices_null, allocate_matrices, &
                                   sparse_matrix_metadata
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
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smats, smatl
      type(atoms_data),intent(in) :: atoms
      type(matrices),intent(inout) :: kernel
      type(matrices),intent(inout) :: ovrlp
      integer,intent(in),optional :: ntheta, istheta
      real(kind=8),dimension(:,:),intent(in),optional :: theta !dimension (atoms%astruct%nat,ntheta)

      ! Local variables
      integer :: ierr, iorb, iat, ind, ist, ishift, ispin, iiorb
      integer, dimension(1) :: power
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

              power(1)=2
              call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
                   meth_overlap, 1, power, blocksize, &
                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
                   ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=meth_overlap<1000, &
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
              call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
                   smat_in=ovrlp%matrix_compr, lmat_out=proj_ovrlp_half_compr)
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
          call gather_matrix_from_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, &
               smatl, weight_matrix_compr_tg, weight_matrix_compr)

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
                      iat=smmd%on_which_atom(iiorb)
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
      real(kind=8),dimension(3) :: dipole
    
      dipole(1:3) = 0._gp
      do iat=1,atoms%astruct%nat
          dipole(1:3) = dipole(1:3) + &
                            (atoms%nelpsp(atoms%astruct%iatype(iat))-charge_per_atom(iat))*atoms%astruct%rxyz(1:3,iat)
      end do
    
      if (iproc==0) then
          call yaml_map('net dipole', dipole,fmt='(es12.5)')
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




    subroutine build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, order_taylor, &
               energy, energyDiff, energyold, ref_frags, frag_coeffs)
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
      use module_fragments
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
      type(system_fragment), dimension(:), pointer :: ref_frags
      logical, intent(in) :: frag_coeffs
    
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
      integer :: nstates_max, ndimcoeff
      logical :: overlap_calculated=.false. ! recalculate just to be safe
      real(kind=8), allocatable, dimension(:) :: coeff_tmp

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
          call corrections_for_negative_charge(iproc, nproc, at, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
    
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
    
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      if (.not. frag_coeffs) then
         call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
              energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
              order_taylor,input%lin%max_inversion_error,&
              input%calculate_KS_residue,input%calculate_gap, energs_work, .false., input%lin%coeff_factor, &
              input%lin%pexsi_npoles, input%lin%pexsi_mumin, input%lin%pexsi_mumax, input%lin%pexsi_mu, &
              input%lin%pexsi_temperature, input%lin%pexsi_tol_charge)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
    
         if (bigdft_mpi%iproc ==0) then
            call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
         end if
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
    
      !apply frag coeffs instead (might still be orthonormalized)
      !could do this in transfer_integrals to get naming correct and be exactly equivalent, but just want to do easiest way for now
      !assume not cdft, so can use in%frag%charge instead of modifying it
      if (frag_coeffs) then
           !don't overwrite coeffs
           ndimcoeff=size(tmb%coeff)
           coeff_tmp=f_malloc(ndimcoeff, id='coeff_tmp')
           call vcopy(ndimcoeff,tmb%coeff(1,1),1,coeff_tmp(1),1)
           call fragment_coeffs_to_kernel(iproc,input,input%frag%charge,ref_frags,tmb,KSwfn%orbs,overlap_calculated,&
                nstates_max,input%lin%constrained_dft,input%lin%kernel_restart_mode,input%lin%kernel_restart_noise)
      end if
    
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
      call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%linmat%m%nfvctr, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
                 tmb%linmat%m%nfvctr, 0.d0, phiwork_global, nvctrp)

      if (frag_coeffs) then
           call vcopy(ndimcoeff,coeff_tmp(1),1,tmb%coeff(1,1),1)
           call f_free(coeff_tmp)
      end if
      
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
          if (frag_coeffs) then
             call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KSDens', &
                  KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .true.)
          else
             call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KSDensFrag', &
                  KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .true.)
          end if
      else if (input%write_orbitals==3) then 
          if (frag_coeffs) then
              call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KSFrag', &
                   KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .false.)
          else
              call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KS', &
                   KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz, .false.)
          end if
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
          call corrections_for_negative_charge(iproc, nproc, at, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor, input%lin%max_inversion_error, &
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
    


end module postprocessing_linear
