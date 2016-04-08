module orthonormalization

  implicit none

  private

  !> Public routines
  public :: orthonormalizeLocalized
  public :: iterative_orthonormalization
  public :: orthoconstraintNonorthogonal
  public :: gramschmidt_coeff_trans
  !!public :: orthonormalize_subset
  !!public :: gramschmidt_subset
  public :: overlap_matrix

  contains

    !> extract the overlap matrix from two compressed functions
    !this is a quick wrapper to be used whan memory allocation of work arrays is not to be optimized
    subroutine overlap_matrix(phi1,nphi,lzd,orbs,collcom,smat,matrix,phi2)
      use module_base
      use module_types, only: orbitals_data,local_zone_descriptors
      use sparsematrix_base, only: sparse_matrix,matrices
      use transposed_operations, only: calculate_overlap_transposed
      use communications_base, only: TRANSPOSE_FULL,comms_linear
      use communications, only: transpose_localized
      implicit none
      integer, intent(in) :: nphi
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      type(local_zone_descriptors),intent(in) :: lzd
      real(wp),dimension(nphi),intent(in) :: phi1
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: matrix
      
      real(wp), dimension(nphi),intent(in), optional :: phi2
      !local variables
      logical :: binary
      real(wp), dimension(:), pointer :: phi1t_c, phi1t_f, sphi2t_c, sphi2t_f

      binary=present(phi2)

      ! Calculate the scalar products, i.e. the matrix <phi_ab|S_lm|phi_ab>
      phi1t_c = f_malloc_ptr(collcom%ndimind_c,id='phi1t_c')
      phi1t_f = f_malloc_ptr(7*collcom%ndimind_f,id='phi1t_f')
      call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, nphi, orbs, collcom, &
           TRANSPOSE_FULL, phi1, phi1t_c, phi1t_f, lzd)

      if (binary) then
         sphi2t_c = f_malloc_ptr(collcom%ndimind_c,id='sphit2_c')
         sphi2t_f = f_malloc_ptr(7*collcom%ndimind_f,id='sphit2_f')
         call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,&
              nphi, orbs, collcom, &
              TRANSPOSE_FULL, phi2, sphi2t_c, sphi2t_f, lzd)
      else
         sphi2t_c => phi1t_c
         sphi2t_f => phi1t_f
      end if
      call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, orbs, collcom, &
           phi1t_c, sphi2t_c, phi1t_f, sphi2t_f, smat, matrix)

      call f_free_ptr(phi1t_c)
      call f_free_ptr(phi1t_f)

      if (binary) then
         call f_free_ptr(sphi2t_c)
         call f_free_ptr(sphi2t_f)
      else
         nullify(sphi2t_c)
         nullify(sphi2t_f)
      end if

    end subroutine overlap_matrix


    !> Orthonormalized the localized orbitals
    subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, max_inversion_error, npsidim_orbs, &
               orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed)
      use module_base
      use module_types
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized, untranspose_localized
      use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                                   assignment(=), sparsematrix_malloc_ptr, SPARSE_TASKGROUP
      use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace
      use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed, &
                                       normalize_transposed
      use matrix_operations, only: overlapPowerGeneral, overlap_power_minus_one_half_parallel, check_taylor_order, &
                                   calculate_S_minus_one_half_onsite
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc,nproc,npsidim_orbs
      integer,intent(inout) :: methTransformOverlap
      real(kind=8),intent(in) :: max_inversion_error
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(in) :: ovrlp
      type(sparse_matrix),intent(in) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
      type(comms_linear),intent(in) :: collcom
      type(orthon_data),intent(in) :: orthpar
      real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
      real(kind=8),dimension(:),pointer :: psit_c, psit_f
      logical,intent(inout) :: can_use_transposed
    
      ! Local variables
      integer :: it, istat, iall
      real(kind=8), dimension(:),allocatable :: psittemp_c, psittemp_f, norm
      character(len=*), parameter :: subname='orthonormalizeLocalized'
      real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
      real(kind=8) :: mean_error, max_error
      logical :: ovrlp_associated, inv_ovrlp_associated
      type(matrices) :: ovrlp_
      type(matrices),dimension(1) :: inv_ovrlp_half_
      integer :: ii, i, ispin
      integer, dimension(1) :: power
    
    
      call f_routine(id='orthonormalizeLocalized')

    
      inv_ovrlp_half_(1) = matrices_null()
      !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_(1))
      inv_ovrlp_half_(1)%matrix_compr = sparsematrix_malloc_ptr(inv_ovrlp_half, &
          iaction=SPARSE_TASKGROUP, id='inv_ovrlp_half_(1)%matrix_compr')
    
    
    
      if(.not.can_use_transposed) then
          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
               TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.
          !!do i=1,collcom%ndimind_c
          !!    write(750+iproc,'(a,2i8,es14.5)') 'i, mod(i-1,ndimind_c/2)+1, val', i, mod(i-1,collcom%ndimind_c/2)+1, psit_c(i)
          !!end do
      end if
    
      ovrlp_ = matrices_null()
      !call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
      ovrlp_%matrix_compr = sparsematrix_malloc_ptr(ovrlp, &
          iaction=SPARSE_TASKGROUP, id='ovrlp_%matrix_compr')
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
      !!ii=0
      !!do ispin=1,ovrlp%nspin
      !!    do i=1,ovrlp%nvctr
      !!        ii=ii+1
      !!        write(930+iproc,*) 'ii, i, val', ii, i, ovrlp_%matrix_compr(ii)
      !!    end do
      !!end do
    
    
      if (methTransformOverlap==-2) then
          call calculate_S_minus_one_half_onsite(iproc, nproc, bigdft_mpi%mpi_comm, &
               orbs%norb, orbs%onwhichatom, &
               ovrlp, inv_ovrlp_half, ovrlp_, inv_ovrlp_half_(1))
      else if (methTransformOverlap==-1) then
          !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
          call overlap_power_minus_one_half_parallel(iproc, nproc, 0, ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_(1))
      else
          power(1)=-2
          call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
               methTransformOverlap, 1, power, &
               orthpar%blocksize_pdgemm, &
               imode=1, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
               ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
               verbosity=0, &
               check_accur=.true., mean_error=mean_error, max_error=max_error)!!, &
          !if (iproc==0) call yaml_map('max error',max_error)
          !if (iproc==0) call yaml_map('mean error',mean_error)
          call check_taylor_order(iproc, mean_error, max_inversion_error, methTransformOverlap)
          !!ii=0
          !!do ispin=1,inv_ovrlp_half%nspin
          !!    do i=1,inv_ovrlp_half%nvctr
          !!        ii=ii+1
          !!        write(1930+iproc,*) 'ii, i, val', ii, i, inv_ovrlp_half_%matrix_compr(ii)
          !!    end do
          !!end do
      end if
    
      call deallocate_matrices(ovrlp_)
    
      psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
      psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')
    
      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
    
      if (methTransformOverlap==-1) then
          ! this is only because overlap_power_minus_one_half_parallel still needs the entire array without taskgroups... to be improved
      call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      else
      call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      end if
      !if (iproc==0) then
      !    do i=1,size(inv_ovrlp_half_(1)%matrix_compr)
      !        write(*,*) 'i, val', i, inv_ovrlp_half_(1)%matrix_compr(i)
      !    end do
      !end if
    
    
      norm = f_malloc(orbs%norb,id='norm')
      call normalize_transposed(iproc, nproc, orbs, ovrlp%nspin, collcom, psit_c, psit_f, norm)
    
      call f_free(norm)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)
    
      call f_free(psittemp_c)
      call f_free(psittemp_f)
    
      !call f_free_ptr(inv_ovrlp_half%matrix_compr)
    
      call deallocate_matrices(inv_ovrlp_half_(1))
    
      call f_release_routine()
    
    end subroutine orthonormalizeLocalized

    !> Orthogonalization where the "additional" support functions are handled separately
    subroutine iterative_orthonormalization(iproc, nproc, verbosity, iorder, at, nspin, sf_per_type, tmb)
      use module_base
      use module_types, only: DFT_wavefunction
      use module_atoms, only: atoms_data
      use ao_inguess, only: aoig_data, aoig_data_null, aoig_set
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, verbosity, iorder
      type(atoms_data),intent(in) :: at
      integer,intent(in) :: nspin
      integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type !< number of support functions per atom type
      type(DFT_wavefunction),intent(inout) :: tmb
    
      ! Local variables
      integer,dimension(:),allocatable :: minorbs_type, maxorbs_type
      logical,dimension(:),allocatable :: type_covered
      integer :: iortho, itype, jj, iat, inl
      logical :: finished
      integer, dimension(:,:), allocatable :: nl_default
      type(aoig_data),dimension(:),allocatable :: aoig_default
    
      allocate(aoig_default(at%astruct%nat))
      do iat=1,at%astruct%nat
          aoig_default(iat)=aoig_data_null()
      end do
      nl_default=f_malloc((/0.to.3,1.to.at%astruct%nat/),id='nl_default')
      do iat=1,at%astruct%nat
         itype = at%astruct%iatype(iat)
         aoig_default(iat) = aoig_set(at%nzatom(itype), at%nelpsp(itype), &
                             at%astruct%input_polarization(iat), nspin)
         nl_default(:,iat)=aoig_default(iat)%nl(:)
      end do
      deallocate(aoig_default)
    
      if (iproc==0) call yaml_map('orthonormalization of input guess','generalized')
      maxorbs_type = f_malloc(at%astruct%ntypes,id='maxorbs_type')
      minorbs_type = f_malloc(at%astruct%ntypes,id='minorbs_type')
      type_covered = f_malloc(at%astruct%ntypes,id='type_covered')
      minorbs_type(1:at%astruct%ntypes)=0
      iortho=0
      ortho_loop: do
          finished=.true.
          type_covered=.false.
          do iat=1,at%astruct%nat
              itype=at%astruct%iatype(iat)
              if (type_covered(itype)) cycle
              type_covered(itype)=.true.
              jj=nl_default(0,iat)+3*nl_default(1,iat)+5*nl_default(2,iat)+7*nl_default(3,iat)
              maxorbs_type(itype)=jj
              !should not enter in the conditional below due to the raise of the exception above
              if (jj<sf_per_type(at%astruct%iatype(iat))) then
                  finished=.false.
                  increase_count: do inl=1,4
                     if (nl_default(inl,iat)==0) then
                        nl_default(inl,iat)=1
                        !call f_err_throw('InputguessLinear: Should not be here',&
                        !     err_name='BIGDFT_RUNTIME_ERROR')
                        exit increase_count
                     end if
                  end do increase_count
              end if
          end do
          if (iortho>0) then
              call gramschmidt_subset(iproc, nproc, verbosity, iorder, tmb%npsidim_orbs, &
                   tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%s, &
                   tmb%linmat%l, tmb%collcom, tmb%orthpar, &
                   tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
          end if
          call orthonormalize_subset(iproc, nproc, verbosity, iorder, tmb%npsidim_orbs, &                                  
               tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%s, &
               tmb%linmat%l, tmb%collcom, tmb%orthpar, &
               tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
          if (finished) exit ortho_loop
          iortho=iortho+1
          minorbs_type(1:at%astruct%ntypes)=maxorbs_type(1:at%astruct%ntypes)+1
      end do ortho_loop
    
      call f_free(maxorbs_type)
      call f_free(minorbs_type)
      call f_free(type_covered)
      call f_free(nl_default)
    
    end subroutine iterative_orthonormalization

    !> Orthonormalize a subset of orbitals
    subroutine orthonormalize_subset(iproc, nproc, verbosity, methTransformOverlap, npsidim_orbs, &
               orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
               lphi, psit_c, psit_f, can_use_transposed)
      use module_base
      use module_types
      !use module_interfaces
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized, untranspose_localized
      use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                                   sparsematrix_malloc, assignment(=), SPARSE_FULL, SPARSE_TASKGROUP, &
                                   sparsematrix_malloc_ptr
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed, &
                                       normalize_transposed
      use matrix_operations, only: overlapPowerGeneral, overlap_power_minus_one_half_parallel
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc,nproc,verbosity,methTransformOverlap,npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(atoms_data),intent(in) :: at
      integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(in) :: ovrlp
      type(sparse_matrix),intent(in) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
      type(comms_linear),intent(in) :: collcom
      type(orthon_data),intent(in) :: orthpar
      real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
      real(kind=8),dimension(:),pointer :: psit_c, psit_f
      logical,intent(inout) :: can_use_transposed
    
      ! Local variables
      integer :: it, istat, iall, iorb, jorb, iat, jat, ii, itype
      logical :: iout, jout
      integer, dimension(1) :: power
      integer,dimension(:),allocatable :: icount_norb, jcount_norb
      real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm, tmparr
      !type(sparse_matrix) :: inv_ovrlp_half
      character(len=*),parameter :: subname='orthonormalize_subset'
      real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
      real(kind=8) :: max_error, mean_error
      type(matrices) :: ovrlp_
      type(matrices),dimension(1) :: inv_ovrlp_half_
    
      call f_routine(id='orthonormalize_subset')
    
      if (iproc==0 .and. verbosity>1) then
          call yaml_sequence_open('Loewdin orthonormalization for the following orbitals')
          do itype=1,at%astruct%ntypes
              if (minorbs_type(itype)<=maxorbs_type(itype)) then
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(flow=.true.)
                  call yaml_map('atom type',adjustl(trim(at%astruct%atomnames(itype))))
                  call yaml_map('first orbital',minorbs_type(itype))
                  call yaml_map('last orbital',maxorbs_type(itype))
                  call yaml_mapping_close()
              end if
          end do
          call yaml_sequence_close()
      end if
    
      !call nullify_sparse_matrix(inv_ovrlp_half)
      !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
      !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
      !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
    
      inv_ovrlp_half_(1) = matrices_null()
      !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_(1))
      inv_ovrlp_half_(1)%matrix_compr = sparsematrix_malloc_ptr(inv_ovrlp_half, &
          iaction=SPARSE_TASKGROUP, id='inv_ovrlp_half_(1)%matrix_compr')
    
    
      if(.not.can_use_transposed) then
          !!if(associated(psit_c)) then
          !!    call f_free_ptr(psit_c)
          !!end if
          !!if(associated(psit_f)) then
          !!    call f_free_ptr(psit_f)
          !!end if
          !!psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
          !!psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')
    
          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
               TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.
    
      end if
    
      ovrlp_ = matrices_null()
      !call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
      ovrlp_%matrix_compr = sparsematrix_malloc_ptr(ovrlp, &
          iaction=SPARSE_TASKGROUP, id='ovrlp%matrix_compr')
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
      ! This can then be deleted if the transition to the new type has been completed.
    
      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      icount_norb = f_malloc(at%astruct%nat,id='icount_norb')
      jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
      icount_norb=0
      do iorb=1,orbs%norb
          iat=orbs%onwhichatom(iorb)
          icount_norb(iat)=icount_norb(iat)+1
          if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
              icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
              iout=.true.
          else
              iout=.false.
          end if
          jcount_norb=0
          do jorb=1,orbs%norb
              jat=orbs%onwhichatom(jorb)
              jcount_norb(jat)=jcount_norb(jat)+1
              if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
                  jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
                  jout=.true.
              else
                  jout=.false.
              end if
              ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
              !!if (iproc==0) write(444,'(a,2i7,2x,2i7,3x,2l4,3x,3i6)') 'iorb, jorb, iat, jat, iout, jout, icount_norb(iat), minorbs_type(at%iatype(iat)), maxorbs_type(at%iatype(iat))', &
              !!                                                         iorb, jorb, iat, jat, iout, jout, icount_norb(iat), minorbs_type(at%iatype(iat)), maxorbs_type(at%iatype(iat))
              if (ii/=0 .and. (iout .or. jout)) then
                  if (jorb==iorb) then
                      ovrlp_%matrix_compr(ii)=1.d0
                  else
                      ovrlp_%matrix_compr(ii)=0.d0
                  end if
              end if
          end do
      end do
      call f_free(icount_norb)
      call f_free(jcount_norb)
    
    
      if (methTransformOverlap==-1) then
          !ovrlp_ = matrices_null()
          !call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
          !ovrlp_%matrix_compr=ovrlp%matrix_compr
          call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, &
               ovrlp, ovrlp_, inv_ovrlp_half, inv_ovrlp_half_(1))
          !call deallocate_matrices(ovrlp_)
      else
          nullify(inv_ovrlp_null)
          ! do sparse.. check later
          !ovrlp%matrix_compr=ovrlp_%matrix_compr
          !tmparr = sparsematrix_malloc(ovrlp,iaction=SPARSE_FULL,id='tmparr')
          !call vcopy(ovrlp%nvctr*ovrlp%nspin, ovrlp_%matrix_compr(1), 1, tmparr(1), 1)
          !call extract_taskgroup_inplace(ovrlp, ovrlp_)
          power(1)=-2
          call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
               methTransformOverlap, 1, power, &
               orthpar%blocksize_pdsyev, &
               imode=1, check_accur=.true., &
               ovrlp_mat=ovrlp_, inv_ovrlp_mat=inv_ovrlp_half_, &
               ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
               max_error=max_error, mean_error=mean_error)
          !call vcopy(ovrlp%nvctr*ovrlp%nspin, tmparr(1), 1, ovrlp_%matrix_compr(1), 1)
          !call f_free(tmparr)
      end if
    
      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      icount_norb = f_malloc(at%astruct%nat,id='icount_norb')
      jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
      icount_norb=0
      do iorb=1,orbs%norb
          iat=orbs%onwhichatom(iorb)
          icount_norb(iat)=icount_norb(iat)+1
          if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
              icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
              iout=.true.
          else
              iout=.false.
          end if
          jcount_norb=0
          do jorb=1,orbs%norb
              jat=orbs%onwhichatom(jorb)
              jcount_norb(jat)=jcount_norb(jat)+1
              if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
                  jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
                  jout=.true.
              else
                  jout=.false.
              end if
              ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
              if (ii/=0 .and. (iout .or. jout)) then
                  if (jorb==iorb) then
                      ovrlp_%matrix_compr(ii)=1.d0
                  else
                      ovrlp_%matrix_compr(ii)=0.d0
                  end if
              end if
          end do
      end do
      call f_free(icount_norb)
      call f_free(jcount_norb)
    
      psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
      psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')
    
      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
    
      !inv_ovrlp_half_%matrix_compr = inv_ovrlp_half%matrix_compr
      call build_linear_combination_transposed(collcom, inv_ovrlp_half, inv_ovrlp_half_(1), &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
    
    
    
      call deallocate_matrices(ovrlp_)
    
      call deallocate_matrices(inv_ovrlp_half_(1))
    
    
      norm = f_malloc(orbs%norb,id='norm')
      call normalize_transposed(iproc, nproc, orbs, ovrlp%nspin, collcom, psit_c, psit_f, norm)
      call f_free(norm)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, & 
           TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)
    
      call f_free(psittemp_c)
      call f_free(psittemp_f)
      !!call f_free_ptr(inv_ovrlp_half%matrix_compr)
    
      call f_release_routine()
    
    end subroutine orthonormalize_subset


    !> Can still tidy this up more when tmblarge is removed
    !! use sparsity of density kernel for all inverse quantities
    subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
               correction_orthoconstraint, linmat, lphi, lhphi, lagmat, lagmat_, psit_c, psit_f, &
               hpsit_c, hpsit_f, &
               can_use_transposed, overlap_calculated, experimental_mode, calculate_inverse, norder_taylor, max_inversion_error, &
               npsidim_orbs_small, lzd_small, hpsi_noprecond, wt_philarge, wt_hphi, wt_hpsinoprecond)
      use module_base
      use module_types
      use yaml_output
      use communications_base, only: work_transpose, TRANSPOSE_POST, TRANSPOSE_FULL, TRANSPOSE_GATHER
      use communications, only: transpose_localized, untranspose_localized
      use sparsematrix_base, only: matrices_null, allocate_matrices, deallocate_matrices, sparsematrix_malloc, &
                                   sparsematrix_malloc_ptr, DENSE_FULL, DENSE_MATMUL, SPARSE_FULL, SPARSEMM_SEQ, &
                                   assignment(=), SPARSE_TASKGROUP
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: uncompress_matrix, &
                              sequential_acces_matrix_fast2, transform_sparse_matrix, &
                              gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace, &
                              uncompress_matrix_distributed2, &
                              matrix_matrix_mult_wrapper
      use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
      use matrix_operations, only: overlapPowerGeneral, check_taylor_order
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp, npsidim_orbs_small
      type(local_zone_descriptors),intent(in) :: lzd, lzd_small
      !type(orbitals_Data),intent(in) :: orbs
      type(orbitals_Data),intent(inout) :: orbs !temporary inout
      type(comms_linear),intent(in) :: collcom
      type(orthon_data),intent(in) :: orthpar
      integer,intent(in) :: correction_orthoconstraint
      real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
      real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
      type(sparse_matrix),intent(in) :: lagmat
      type(matrices),intent(out) :: lagmat_
      real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: hpsit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: hpsit_f
      real(kind=8),dimension(:),pointer :: psit_c, psit_f
      logical,intent(inout) :: can_use_transposed, overlap_calculated
      type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
      logical,intent(in) :: experimental_mode, calculate_inverse
      integer,intent(inout) :: norder_taylor
      real(kind=8),intent(in) :: max_inversion_error
      real(kind=8),dimension(npsidim_orbs_small),intent(out) :: hpsi_noprecond
      type(work_transpose),intent(inout) :: wt_philarge
      type(work_transpose),intent(out) :: wt_hphi, wt_hpsinoprecond
    
      ! Local variables
      real(kind=8) :: max_error, mean_error
      real(kind=8),dimension(:),allocatable :: tmp_mat_compr, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra
      integer,dimension(:),allocatable :: ipiv
      type(matrices),dimension(1) :: inv_ovrlp_
      integer :: ist, ispin
      integer, dimension(1) :: power
      real(8),dimension(:),allocatable :: inv_ovrlp_seq, lagmat_large, tmpmat, tmparr
      real(8),dimension(:,:),allocatable :: lagmatp, inv_lagmatp
      integer,dimension(2) :: irowcol
      real(kind=8),dimension(:),pointer :: matrix_local
      integer,parameter :: GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: data_strategy_main=SUBMATRIX!GLOBAL_MATRIX
      integer,parameter :: verbosity = 0
      !type(work_transpose) :: wt_
    
      call f_routine(id='orthoconstraintNonorthogonal')
    
    
    
      !inv_ovrlp_seq = sparsematrix_malloc(linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_seq')
      !inv_lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='inv_lagmatp')
      !lagmatp = sparsematrix_malloc(linmat%l, iaction=DENSE_MATMUL, id='lagmatp')
      !!inv_ovrlp_(1) = matrices_null()
      !!inv_ovrlp_(1)%matrix_compr = sparsematrix_malloc_ptr(linmat%l,iaction=SPARSE_FULL,id='inv_ovrlp_(1)%matrix_compr')
    
      if (calculate_inverse) then
          ! Invert the overlap matrix
          if (iproc==0 .and. verbosity>0) call yaml_map('calculation of S^-1','direct calculation')
          !!tmparr = sparsematrix_malloc(linmat%s,iaction=SPARSE_FULL,id='tmparr')
          !!call vcopy(linmat%s%nvctr*linmat%s%nspin, linmat%ovrlp_%matrix_compr(1), 1, tmparr(1), 1)
          !!call extract_taskgroup_inplace(linmat%s, linmat%ovrlp_)
          power(1)=1
          call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
               norder_taylor, 1, power, -1, &
               imode=1, ovrlp_smat=linmat%s, inv_ovrlp_smat=linmat%l, &
               ovrlp_mat=linmat%ovrlp_, inv_ovrlp_mat=linmat%ovrlppowers_(3), &
               verbosity=0, &
               check_accur=.true., max_error=max_error, mean_error=mean_error)
          !!call vcopy(linmat%s%nvctr*linmat%s%nspin, tmparr(1), 1, linmat%ovrlp_%matrix_compr(1), 1)
          !!call f_free(tmparr)
          call check_taylor_order(iproc, mean_error, max_inversion_error, norder_taylor)
      else
          if (iproc==0 .and. verbosity>0) call yaml_map('calculation of S^-1','from memory')
      end if
    
      ! Gather together the data (has been posted in getLocalizedBasis
      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_GATHER, lphi, psit_c, psit_f, lzd, wt_philarge)
      can_use_transposed=.true.
    
      ! Calculate <phi_alpha|g_beta>
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat, lagmat_)
      !call gather_matrix_from_taskgroups_inplace(iproc, nproc, lagmat, lagmat_)
    
      lagmat_large = sparsematrix_malloc(linmat%l, iaction=SPARSE_TASKGROUP, id='lagmat_large')
    
      call symmetrize_matrix()
    
    
      ! Apply S^-1
      !!call sequential_acces_matrix_fast2(linmat%l, linmat%ovrlppowers_(3)%matrix_compr, inv_ovrlp_seq)
      ! Transform the matrix to the large sparsity pattern (necessary for the following uncompress_matrix_distributed)
      if (correction_orthoconstraint==0) then
          if (data_strategy_main==GLOBAL_MATRIX) then
              stop 'deprecated'
              !!!call transform_sparse_matrix(iproc, linmat%m, linmat%l, SPARSE_TASKGROUP, 'small_to_large', &
              !!!     smat_in=lagmat_%matrix_compr, lmat_out=lagmat_large)
          end if
          if (iproc==0) call yaml_map('correction orthoconstraint',.true.)
          !!call uncompress_matrix_distributed2(iproc, linmat%l, DENSE_MATMUL, lagmat_large, lagmatp)
          !!call sparsemm(linmat%l, inv_ovrlp_seq, lagmatp, inv_lagmatp)
          !!write(*,*) 'iproc, sum(inv_lagmatp)', iproc, sum(inv_lagmatp)
          !!call compress_matrix_distributed(iproc, nproc, linmat%l, DENSE_MATMUL, &
          !!     inv_lagmatp, lagmat_large)
          do ispin=1,linmat%l%nspin
              ist=(ispin-1)*linmat%l%nvctrp_tg+1
              call matrix_matrix_mult_wrapper(iproc, nproc, linmat%l, &
                   linmat%ovrlppowers_(3)%matrix_compr(ist:), lagmat_large(ist:), lagmat_large(ist:))
          end do
      end if
      if (data_strategy_main==SUBMATRIX) then
          call transform_sparse_matrix(iproc, linmat%m, linmat%l, SPARSE_TASKGROUP, 'large_to_small', &
               lmat_in=lagmat_large, smat_out=lagmat_%matrix_compr)
      end if
      call f_free(lagmat_large)
      !call f_free(inv_ovrlp_seq)
      !call f_free(lagmatp)
      !call f_free(inv_lagmatp)
    
      !!tmparr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmparr')
      !!call vcopy(lagmat%nvctr, lagmat_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, lagmat, lagmat_)
      call build_linear_combination_transposed(collcom, lagmat, lagmat_, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
      !!call vcopy(lagmat%nvctr, tmparr(1), 1, lagmat_%matrix_compr(1), 1)
      !!call f_free(tmparr)
    
    
      ! Start the untranspose process (will be gathered together in
      ! calculate_energy_and_gradient_linear)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_POST, hpsit_c, hpsit_f, lhphi, lzd, wt_hphi)
    
      ! The symmetrized Lagrange multiplier matrix has now the wrong sign
      call dscal(lagmat%nvctrp_tg*lagmat%nspin, -1.d0, lagmat_%matrix_compr(1), 1)

      overlap_calculated=.false.
      
    
      ! @NEW apply S^-1 to the gradient
      hpsit_tmp_c = f_malloc(collcom%ndimind_c,id='psit_tmp_c')
      hpsit_tmp_f = f_malloc(7*collcom%ndimind_f,id='psit_tmp_f')
      hphi_nococontra = f_malloc(npsidim_orbs,id='hphi_nococontra')
      call build_linear_combination_transposed(collcom, linmat%l, linmat%ovrlppowers_(3), &
           hpsit_c, hpsit_f, .true., hpsit_tmp_c, hpsit_tmp_f, iproc)
    
      ! Start the untranspose process (will be gathered together in
      ! getLocalizedBasis)
      ! Pass hphi_nococontra even if it is not used
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_POST, hpsit_tmp_c, hpsit_tmp_f, hphi_nococontra, lzd, wt_hpsinoprecond)
      !!call large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs, lzd_small, lzd, &
      !!     orbs, hphi_nococontra, hpsi_noprecond)
      ! END @NEW
    
      !call deallocate_matrices(inv_ovrlp_(1))
      call f_free(hpsit_tmp_c)
      call f_free(hpsit_tmp_f)
      call f_free(hphi_nococontra)
    
      call f_release_routine()
    
    
      contains
    
    
        subroutine symmetrize_matrix()
          implicit none
          integer :: ishift, itg, iitg
          integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
          integer,parameter :: comm_strategy=GET
          integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
          integer :: iorb, jorb, ii, ii_trans, irow, jcol, info, lwork, jj, ispin, iseg, i
          integer :: isegstart, isegend, ierr
    
    
          call f_routine(id='symmetrize_matrix')
    
          ! Just to check the consistency
          if (data_strategy_main/=data_strategy) then
              stop 'data_strategy_main/=data_strategy'
          end if
    
          if (lagmat%nvctrp>0) then
              isegstart = lagmat%istsegline(lagmat%isfvctr+1)
              isegend = lagmat%istsegline(lagmat%isfvctr+lagmat%nfvctrp) + &
                        lagmat%nsegline(lagmat%isfvctr+lagmat%nfvctrp)-1
          else
              isegstart = 1
              isegend = 0
          end if
          if (data_strategy==GLOBAL_MATRIX) then
              stop 'symmetrize_matrix: option GLOBAL_MATRIX is deprecated'
              !!matrix_local = f_malloc_ptr(max(lagmat%nvctrp,1),id='matrix_local')
              !!tmp_mat_compr = sparsematrix_malloc(lagmat,iaction=SPARSE_FULL,id='tmp_mat_compr')
              !!call vcopy(lagmat%nvctr*lagmat%nspin, lagmat_%matrix_compr(1), 1, tmp_mat_compr(1), 1)
              !!do ispin=1,lagmat%nspin
              !!    ishift=(ispin-1)*lagmat%nvctr
              !!    if (isegend>=isegstart) then
              !!        !$omp parallel default(none) &
              !!        !$omp shared(isegstart,isegend,ishift,lagmat,matrix_local,tmp_mat_compr) &
              !!        !$omp private(iseg,ii,i,irowcol,ii_trans)
              !!        !$omp do
              !!        do iseg=isegstart,isegend
              !!            ii=lagmat%keyv(iseg)
              !!            ! A segment is always on one line, therefore no double loop
              !!            do i=lagmat%keyg(1,1,iseg),lagmat%keyg(2,1,iseg)
              !!               irowcol = orb_from_index(lagmat, i)
              !!               ii_trans = matrixindex_in_compressed(lagmat,lagmat%keyg(1,2,iseg),i)
              !!               matrix_local(ii-lagmat%isvctr) = -0.5d0*tmp_mat_compr(ii+ishift)-0.5d0*tmp_mat_compr(ii_trans+ishift)
              !!               ii=ii+1
              !!            end do
              !!        end do
              !!        !$omp end do
              !!        !$omp end parallel
              !!    end if
              !!    if (nproc>1) then
              !!        !!call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
              !!        !!     lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
              !!        !!     bigdft_mpi%mpi_comm, ierr)
              !!        if (comm_strategy==ALLGATHERV) then
              !!            call mpi_allgatherv(matrix_local(1), lagmat%nvctrp, mpi_double_precision, &
              !!                 lagmat_%matrix_compr(ishift+1), lagmat%nvctr_par, lagmat%isvctr_par, mpi_double_precision, &
              !!                 bigdft_mpi%mpi_comm, ierr)
              !!            call f_free_ptr(matrix_local)
              !!        else if (comm_strategy==GET) then
              !!            !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, lagmat%nvctrp, matrix_local, &
              !!            !!     lagmat%nvctr_par, lagmat%isvctr_par, lagmat%nvctr, lagmat_%matrix_compr(ishift+1:ishift+lagmat%nvctr))
              !!            call mpi_get_to_allgatherv(matrix_local(1), lagmat%nvctrp, &
              !!                 lagmat_%matrix_compr(ishift+1), &
              !!                 lagmat%nvctr_par, lagmat%isvctr_par, bigdft_mpi%mpi_comm)
              !!        else
              !!            stop 'symmetrize_matrix: wrong communication strategy'
              !!        end if
              !!    else
              !!        call vcopy(lagmat%nvctr, matrix_local(1), 1, lagmat_%matrix_compr(ishift+1), 1)
              !!    end if
              !!    if (ispin==lagmat%nspin) call f_free_ptr(matrix_local)
              !!end do
          else if (data_strategy==SUBMATRIX) then
              ! Directly use the large sparsity pattern as this one is used later
              ! for the matrix vector multiplication
              tmp_mat_compr = sparsematrix_malloc(linmat%l,iaction=SPARSE_TASKGROUP,id='tmp_mat_compr')
              call transform_sparse_matrix(iproc, linmat%m, linmat%l, SPARSE_TASKGROUP, 'small_to_large', &
                   smat_in=lagmat_%matrix_compr, lmat_out=tmp_mat_compr)
              do ispin=1,lagmat%nspin
                  ishift=(ispin-1)*linmat%l%nvctrp_tg
                  !$omp parallel default(none) &
                  !$omp shared(linmat,lagmat_large,tmp_mat_compr,ishift) &
                  !$omp private(iseg,ii,i,ii_trans)
                  !$omp do
                  !do iseg=linmat%l%iseseg_tg(1),linmat%l%iseseg_tg(2)
                  do iseg=linmat%l%istartendseg_local(1),linmat%l%istartendseg_local(2)
                      ii = linmat%l%keyv(iseg)
                      ! A segment is always on one line, therefore no double loop
                      do i=linmat%l%keyg(1,1,iseg),linmat%l%keyg(2,1,iseg) !this is too much, but for the moment ok
                          ii_trans = matrixindex_in_compressed(linmat%l,linmat%l%keyg(1,2,iseg),i)
                          lagmat_large(ii+ishift-linmat%l%isvctrp_tg) = &
                              - 0.5d0*tmp_mat_compr(ii+ishift-linmat%l%isvctrp_tg) &
                              - 0.5d0*tmp_mat_compr(ii_trans+ishift-linmat%l%isvctrp_tg)
                          ii=ii+1
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
              end do
          else
              stop 'symmetrize_matrix: wrong data strategy'
          end if
    
          call f_free(tmp_mat_compr)
          call f_release_routine()
    
        end subroutine symmetrize_matrix
    
    
    end subroutine orthoconstraintNonorthogonal


    subroutine gramschmidt_coeff_trans(iproc,nproc,norbu,norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices
      implicit none
    
      integer, intent(in) :: iproc, nproc, norbu, norb
      type(orbitals_data), intent(in) :: basis_orbs
      type(sparse_matrix),intent(in) :: basis_overlap
      type(matrices),intent(inout) :: basis_overlap_mat
      real(kind=8),dimension(basis_overlap%nfvctr,basis_orbs%norb),intent(inout) :: coeff
    
      integer :: iorb, jtmb, corb, ierr, isorb, ieorb, ispin
      real(kind=8), dimension(:,:), allocatable :: ovrlp_coeff, coeff_tmp, coeff_trans, coeff_transp
    
      real(kind=4) :: tr0, tr1
      real(kind=8) :: time0, time1, time2, time3, time4, time5
    
      call f_routine(id='gramschmidt_coeff_trans')
    
      time0=0.0d0
      time1=0.0d0
      time2=0.0d0
      time3=0.0d0
      time4=0.0d0
      time5=0.0d0
    
      call cpu_time(tr0)
    
      spin_loop: do ispin=1,basis_overlap%nspin
    
          if (ispin==1) then
              isorb=1
              ieorb=norbu
          else if (ispin==2) then
              isorb=norbu+1
              ieorb=norb
          end if
    
          coeff_transp=f_malloc((/(ieorb-isorb+1),basis_overlap%nfvctrp/),id='coeff_transp')
          !$omp parallel do default(private) shared(coeff,coeff_transp,isorb,ieorb,basis_overlap)
          do iorb=isorb,ieorb
             do jtmb=1,basis_overlap%nfvctrp
                coeff_transp(iorb-isorb+1,jtmb) = coeff(jtmb+basis_overlap%isfvctr,iorb)
             end do
          end do
          !$omp end parallel do
    
          call cpu_time(tr1)
          time0=time0+real(tr1-tr0,kind=8)
    
          ! orthonormalizing all iorb<corb wrt corb (assume original vectors were normalized)
          do corb=ieorb-isorb+1,1,-1
    
    
             call cpu_time(tr0)
    
             coeff_tmp=f_malloc((/basis_overlap%nfvctr,1/),id='coeff_tmp')
             ! calculate relevant part of cSc
             if (basis_overlap%nfvctrp>0) then
                call dgemm('n', 't', basis_overlap%nfvctr, 1, basis_overlap%nfvctrp, 1.d0, &
                     basis_overlap_mat%matrix(1,basis_overlap%isfvctr+1,1), basis_overlap%nfvctr, &
                     coeff_transp(corb,1), ieorb-isorb+1, 0.d0, &
                     coeff_tmp(1,1), basis_overlap%nfvctr)
             else
                !!call to_zero(corb,coeff_tmp(1,1)) !!!LG: is this a typo?
                call f_zero(coeff_tmp)
             end if
    
             if (nproc>1) then
                call mpiallred(coeff_tmp, mpi_sum, comm=bigdft_mpi%mpi_comm)
             end if
    
             call cpu_time(tr1)
             time1=time1+real(tr1-tr0,kind=8)
    
             ovrlp_coeff=f_malloc((/corb,1/),id='ovrlp_coeff')
             if (basis_overlap%nfvctrp>0) then
                call dgemm('n', 'n', corb, 1, basis_overlap%nfvctrp, 1.d0, &
                     coeff_transp(1,1), ieorb-isorb+1, &
                     coeff_tmp(1+basis_overlap%isfvctr,1), basis_overlap%nfvctr, 0.d0, &
                     ovrlp_coeff(1,1), corb)
             else
                call f_zero(ovrlp_coeff)
             end if
    
             if (nproc>1) then
                call mpiallred(ovrlp_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
             end if
             call f_free(coeff_tmp)
    
             call cpu_time(tr0)
             time2=time2+real(tr0-tr1,kind=8)
    
             ! (c_corb S c_iorb) * c_corb
             coeff_tmp=f_malloc((/corb,basis_overlap%nfvctrp/),id='coeff_tmp')
             if (basis_overlap%nfvctrp>0) then
                call dgemm('n', 'n', corb-1, basis_overlap%nfvctrp, 1, 1.d0, &
                     ovrlp_coeff(1,1), corb, &
                     coeff_transp(corb,1), ieorb-isorb+1, 0.d0, &
                     coeff_tmp(1,1), corb)
             end if
    
             call cpu_time(tr1)
             time3=time3+real(tr1-tr0,kind=8)
    
             ! sum and transpose coeff for allgatherv
             !$omp parallel do default(private) shared(coeff_transp,coeff_tmp,isorb,corb,basis_overlap,ovrlp_coeff)
             do iorb=isorb,corb-1
                do jtmb=1,basis_overlap%nfvctrp
                   coeff_transp(iorb,jtmb) = coeff_transp(iorb,jtmb) - coeff_tmp(iorb,jtmb)/ovrlp_coeff(corb,1)
                end do
             end do
             !$omp end parallel do
    
             do jtmb=1,basis_overlap%nfvctrp
                coeff_transp(corb,jtmb) = coeff_transp(corb,jtmb)/sqrt(ovrlp_coeff(corb,1))
             end do
    
             call cpu_time(tr0)
             time4=time4+real(tr0-tr1,kind=8)
    
             call f_free(ovrlp_coeff)
             call f_free(coeff_tmp)
          end do
    
          call cpu_time(tr0)
    
          coeff_trans=f_malloc((/(ieorb-isorb+1),basis_overlap%nfvctr/),id='coeff_tmp')
          if(nproc > 1) then
             call mpi_allgatherv(coeff_transp(1,1), basis_overlap%nfvctrp*(ieorb-isorb+1), mpi_double_precision, coeff_trans(1,1), &
                  (ieorb-isorb+1)*basis_overlap%nfvctr_par(:), (ieorb-isorb+1)*basis_overlap%isfvctr_par, &
                  mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          else
             call vcopy(basis_overlap%nfvctrp*(ieorb-isorb+1),coeff_transp(1,1),1,coeff_trans(1,1),1)
          end if
          call f_free(coeff_transp)
    
          call cpu_time(tr1)
          time5=time5+real(tr1-tr0,kind=8)
    
          ! untranspose coeff
          !$omp parallel do default(private) shared(coeff,coeff_trans,isorb,ieorb,basis_overlap)
          do jtmb=1,basis_overlap%nfvctr
             do iorb=isorb,ieorb
                coeff(jtmb,iorb) = coeff_trans(iorb-isorb+1,jtmb)
             end do
          end do
          !$omp end parallel do
    
          call cpu_time(tr0)
          time0=time0+real(tr0-tr1,kind=8)
    
          call f_free(coeff_trans)
    
          !if (iproc==0) print*,'Time in gramschmidt_coeff',time0,time1,time2,time3,time4,time5,&
          !     time0+time1+time2+time3+time4+time5
    
      end do spin_loop
    
      call f_release_routine()
    
    
    end subroutine gramschmidt_coeff_trans






    subroutine gramschmidt_subset(iproc, nproc, verbosity, methTransformOverlap, npsidim_orbs, &
               orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
               lphi, psit_c, psit_f, can_use_transposed)
      use module_base
      use module_types
      !use module_interfaces, exceptThisOne => gramschmidt_subset
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized, untranspose_localized
      use sparsematrix_base, only: sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                                   SPARSE_TASKGROUP, assignment(=), sparsematrix_malloc_ptr
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace
      use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc,nproc,verbosity,methTransformOverlap,npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(atoms_data),intent(in) :: at
      integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(in) :: ovrlp
      type(sparse_matrix),intent(in) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
      type(comms_linear),intent(in) :: collcom
      type(orthon_data),intent(in) :: orthpar
      real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
      real(kind=8),dimension(:),pointer :: psit_c, psit_f
      logical,intent(inout) :: can_use_transposed
    
      ! Local variables
      integer :: it, istat, iall, iorb, jorb, iat, jat, ii, itype
      logical :: iout, jout
      integer,dimension(:),allocatable :: icount_norb, jcount_norb
      real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
      !type(sparse_matrix) :: inv_ovrlp_half
      character(len=*),parameter :: subname='gramschmidt_subset'
      type(matrices) :: ovrlp_
    
      call f_routine('gramschmidt_subset')
    
      if (iproc==0 .and. verbosity>1) then
          call yaml_sequence_open('Gram-Schmidt orthogonalization for the following orbitals')
          do itype=1,at%astruct%ntypes
              if (minorbs_type(itype)<=maxorbs_type(itype)) then
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(flow=.true.)
                  call yaml_map('atom type',adjustl(trim(at%astruct%atomnames(itype))))
                  call yaml_map('first orbital',minorbs_type(itype))
                  call yaml_map('last orbital',maxorbs_type(itype))
                  call yaml_mapping_close()
              end if
          end do
          call yaml_sequence_close()
      end if
    
    
      !call nullify_sparse_matrix(inv_ovrlp_half)
      !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
      !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
      !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
      !!inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')
    
    
      if(.not.can_use_transposed) then
          !!if(associated(psit_c)) then
          !!    call f_free_ptr(psit_c)
          !!end if
          !!if(associated(psit_f)) then
          !!    call f_free_ptr(psit_f)
          !!end if
          !!psit_c = f_malloc_ptr(sum(collcom%nrecvcounts_c),id='psit_c')
          !!psit_f = f_malloc_ptr(7*sum(collcom%nrecvcounts_f),id='psit_f')
    
          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
               TRANSPOSE_FULL, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.
    
      end if
    
    
      ovrlp_ = matrices_null()
      !call allocate_matrices(ovrlp, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
      ovrlp_%matrix_compr = sparsematrix_malloc_ptr(ovrlp, &
          iaction=SPARSE_TASKGROUP, id='ovrlp%matrix_compr')
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp, ovrlp_)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp, ovrlp_)
      ! This can then be deleted if the transition to the new type has been completed.
      !ovrlp%matrix_compr=ovrlp_%matrix_compr
    
      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      icount_norb = f_malloc0(at%astruct%nat,id='icount_norb')
      jcount_norb = f_malloc(at%astruct%nat,id='jcount_norb')
      do iorb=1,orbs%norb
          iat=orbs%onwhichatom(iorb)
          icount_norb(iat)=icount_norb(iat)+1
          if (icount_norb(iat)<minorbs_type(at%astruct%iatype(iat)) .or. &
              icount_norb(iat)>maxorbs_type(at%astruct%iatype(iat))) then
              iout=.true.
          else
              iout=.false.
          end if
          call f_zero(jcount_norb)
          do jorb=1,orbs%norb
              jat=orbs%onwhichatom(jorb)
              jcount_norb(jat)=jcount_norb(jat)+1
              !!if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat)) .or. &
              !!    jcount_norb(jat)>maxorbs_type(at%astruct%iatype(jat))) then
              if (jcount_norb(jat)<minorbs_type(at%astruct%iatype(jat))) then
                  jout=.true.
              else
                  jout=.false.
              end if
              ii=matrixindex_in_compressed(ovrlp,jorb,iorb)
              if (ii/=0) then
                  if (iout) then
                      ovrlp_%matrix_compr(ii)=0.d0
                  else
                      if (jout) then
                          ovrlp_%matrix_compr(ii)=-ovrlp_%matrix_compr(ii)
                      else
                          ovrlp_%matrix_compr(ii)=0.d0
                      end if
                  end if
              end if
              !!if (iout .or. jout) then
              !!    if (jorb==iorb) then
              !!        ovrlp%matrix_compr(ii)=1.d0
              !!    else
              !!        ovrlp%matrix_compr(ii)=0.d0
              !!    end if
              !!end if
          end do
      end do
      call f_free(icount_norb)
      call f_free(jcount_norb)
    
    
      !!if (methTransformOverlap==-1) then
      !!    call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, orbs, ovrlp, inv_ovrlp_half)
      !!else
      !!    call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, methTransformOverlap, orthpar%blocksize_pdsyev, &
      !!        orthpar%blocksize_pdgemm, orbs%norb, ovrlp, inv_ovrlp_half)
      !!end if
    
      !!! For the "higher" TMBs: delete off-diagonal elements and
      !!! set diagonal elements to 1
      !!allocate(icount_norb(at%nat),stat=istat)
      !!call memocc(istat,icount_norb,'icount_norb',subname)
      !!allocate(jcount_norb(at%nat),stat=istat)
      !!call memocc(istat,jcount_norb,'jcount_norb',subname)
      !!do iorb=1,orbs%norb
      !!    iat=orbs%onwhichatom(iorb)
      !!    icount_norb(iat)=icount_norb(iat)+1
      !!    if (icount_norb(iat)<minorbs_type(at%iatype(iat)) .or. &
      !!        icount_norb(iat)>maxorbs_type(at%iatype(iat))) then
      !!        iout=.true.
      !!    else
      !!        iout=.false.
      !!    end if
      !!    do jorb=1,orbs%norb
      !!        jat=orbs%onwhichatom(jorb)
      !!        jcount_norb(jat)=jcount_norb(jat)+1
      !!        if (jcount_norb(jat)>maxorbs_type(at%iatype(jat))) then
      !!            jout=.true.
      !!        else
      !!            jout=.false.
      !!        end if
      !!        ii=ovrlp%matrixindex_in_compressed(jorb,iorb)
      !!        if (iout .or. jout) then
      !!            if (jorb==iorb) then
      !!                ovrlp%matrix_compr(ii)=1.d0
      !!            else
      !!                ovrlp%matrix_compr(ii)=0.d0
      !!            end if
      !!        end if
      !!    end do
      !!end do
      !!iall=-product(shape(icount_norb))*kind(icount_norb)
      !!deallocate(icount_norb, stat=istat)
      !!call memocc(istat, iall, 'icount_norb', subname)
      !!iall=-product(shape(jcount_norb))*kind(jcount_norb)
      !!deallocate(jcount_norb, stat=istat)
      !!call memocc(istat, iall, 'jcount_norb', subname)
    
      psittemp_c = f_malloc(sum(collcom%nrecvcounts_c),id='psittemp_c')
      psittemp_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psittemp_f')
    
      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
      !!call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
      !!     psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
    
    
    
      !inv_ovrlp_ = matrices_null()
      !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='ovrlp_', mat=ovrlp_)
      !@WARNING CHECK THIS
      !!ovrlp_%matrix_compr = inv_ovrlp_half%matrix_compr
      call build_linear_combination_transposed(collcom, ovrlp, ovrlp_, &
           psittemp_c, psittemp_f, .false., psit_c, psit_f, iproc)
      !call deallocate_matrices(ovrlp_)
    
    
      call deallocate_matrices(ovrlp_)
    
    
      norm = f_malloc(orbs%norb,id='norm')
      !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
    
      call f_free(norm)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
           TRANSPOSE_FULL, psit_c, psit_f, lphi, lzd)
    
      call f_free(psittemp_c)
      call f_free(psittemp_f)
      !!call f_free_ptr(inv_ovrlp_half%matrix_compr)
    
      call f_release_routine()
    
    end subroutine gramschmidt_subset



    !!subroutine gramschmidt_coeff(iproc,nproc,norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
    !!  use module_base
    !!  use module_types
    !!  use sparsematrix_base, only: sparse_matrix, matrices
    !!  implicit none
    !!
    !!  integer, intent(in) :: iproc, nproc, norb
    !!  type(orbitals_data), intent(in) :: basis_orbs
    !!  type(sparse_matrix),intent(inout) :: basis_overlap
    !!  type(matrices),intent(inout) :: basis_overlap_mat
    !!  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
    !!
    !!  integer :: iorb, jtmb, corb, ierr
    !!  real(kind=8), dimension(:,:), allocatable :: ovrlp_coeff, coeff_tmp, coeff_trans
    !!  real(kind=4) :: tr0, tr1
    !!  real(kind=8) :: time0, time1, time2, time3, time4, time5
    !!
    !!  call f_routine(id='gramschmidt_coeff')
    !!
    !!  time0=0.0d0
    !!  time1=0.0d0
    !!  time2=0.0d0
    !!  time3=0.0d0
    !!  time4=0.0d0
    !!  time5=0.0d0
    !!
    !!  ! orthonormalizing all iorb<corb wrt corb (assume original vectors were normalized)
    !!  do corb=norb,1,-1
    !!     ovrlp_coeff=f_malloc((/corb,1/),id='ovrlp_coeff')
    !!     coeff_tmp=f_malloc((/corb,basis_orbs%norbp/),id='coeff_tmp')
    !!     ! calculate relevant part of cSc
    !!
    !!     call cpu_time(tr0)
    !!     if (basis_orbs%norbp>0) then
    !!        call dgemm('t', 'n', corb, basis_orbs%norbp, basis_orbs%norb, 1.d0, &
    !!             coeff(1,1), basis_orbs%norb, &
    !!             basis_overlap_mat%matrix(1,basis_orbs%isorb+1,1), basis_orbs%norb, 0.d0, &
    !!             coeff_tmp(1,1), corb)
    !!
    !!        call dgemm('n', 'n', corb, 1, basis_orbs%norbp, 1.d0, &
    !!             coeff_tmp(1,1), corb, &
    !!             coeff(basis_orbs%isorb+1,corb), basis_orbs%norb, 0.d0, &
    !!             ovrlp_coeff(1,1), corb)
    !!     else
    !!        call f_zero(ovrlp_coeff)
    !!     end if
    !!
    !!     if (nproc>1) then
    !!        call mpiallred(ovrlp_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
    !!     end if
    !!
    !!     call cpu_time(tr1)
    !!     time1=time1+real(tr1-tr0,kind=8)
    !!
    !!     ! (c_corb S c_iorb) * c_corb
    !!     if (basis_orbs%norbp>0) then
    !!        call dgemm('n', 't', corb-1, basis_orbs%norbp, 1, 1.d0, &
    !!             ovrlp_coeff(1,1), corb, &
    !!             coeff(1+basis_orbs%isorb,corb), basis_orbs%norb, 0.d0, &
    !!             coeff_tmp(1,1), corb)
    !!     end if
    !!     call cpu_time(tr0)
    !!     time2=time2+real(tr0-tr1,kind=8)
    !!
    !!     ! sum and transpose coeff for allgatherv
    !!     !$omp parallel do default(private) shared(coeff,coeff_tmp,corb,basis_orbs,ovrlp_coeff)
    !!     do iorb=1,corb-1
    !!        do jtmb=1,basis_orbs%norbp
    !!           coeff_tmp(iorb,jtmb) = coeff(jtmb+basis_orbs%isorb,iorb) - coeff_tmp(iorb,jtmb)/ovrlp_coeff(corb,1)
    !!        end do
    !!     end do
    !!     !$omp end parallel do
    !!
    !!     do jtmb=1,basis_orbs%norbp
    !!        coeff_tmp(corb,jtmb) = coeff(jtmb+basis_orbs%isorb,corb)/sqrt(ovrlp_coeff(corb,1))
    !!     end do
    !!     call cpu_time(tr1)
    !!     time3=time3+real(tr1-tr0,kind=8)
    !!
    !!     call f_free(ovrlp_coeff)
    !!     coeff_trans=f_malloc((/corb,basis_orbs%norb/),id='coeff_tmp')
    !!     if(nproc > 1) then
    !!        call mpi_allgatherv(coeff_tmp(1,1), basis_orbs%norbp*corb, mpi_double_precision, coeff_trans(1,1), &
    !!           corb*basis_orbs%norb_par(:,0), corb*basis_orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    !!     else
    !!        call vcopy(basis_orbs%norbp*corb,coeff_tmp(1,1),1,coeff_trans(1,1),1)
    !!     end if
    !!     call f_free(coeff_tmp)
    !!
    !!     call cpu_time(tr0)
    !!     time4=time4+real(tr0-tr1,kind=8)
    !!
    !!     ! untranspose coeff
    !!     !$omp parallel do default(private) shared(coeff,coeff_trans,corb,basis_orbs)
    !!     do jtmb=1,basis_orbs%norb
    !!        do iorb=1,corb
    !!           coeff(jtmb,iorb) = coeff_trans(iorb,jtmb)
    !!        end do
    !!     end do
    !!     !$omp end parallel do
    !!
    !!     call cpu_time(tr1)
    !!     time5=time5+real(tr1-tr0,kind=8)
    !!
    !!     call f_free(coeff_trans)
    !!  end do
    !!
    !!  !if (iproc==0) print*,'Time in gramschmidt_coeff',time0,time1,time2,time3,time4,time5,&
    !!  !     time0+time1+time2+time3+time4+time5
    !!
    !!  call f_release_routine()
    !!
    !!end subroutine gramschmidt_coeff


end module orthonormalization
