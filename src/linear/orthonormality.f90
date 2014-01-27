!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(inout) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall
  !integer :: irow, ii, iorb, jcol, jorb
  real(kind=8), dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparseMatrix) :: inv_ovrlp_half
  character(len=*), parameter :: subname='orthonormalizeLocalized'
  !real(kind=8), dimension(orbs%norb,orbs%norb) :: tempmat
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparsematrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)

  do it=1,orthpar%nItOrtho

      if(.not.can_use_transposed) then
          if(associated(psit_c)) then
              iall=-product(shape(psit_c))*kind(psit_c)
              deallocate(psit_c, stat=istat)
              call memocc(istat, iall, 'psit_c', subname)
          end if
          if(associated(psit_f)) then
              iall=-product(shape(psit_f))*kind(psit_f)
              deallocate(psit_f, stat=istat)
              call memocc(istat, iall, 'psit_f', subname)
          end if
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)

          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.

      end if
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)
      !!do ii=1,ovrlp%nvctr
      !!   irow = ovrlp%orb_from_index(ii,1)
      !!   jcol = ovrlp%orb_from_index(ii,2)
      !!   Tempmat(irow,jcol)=ovrlp%matrix_compr(ii)
      !!end do
      !!if (iproc==0) then
      !!    do irow=1,orbs%norb
      !!        do jcol=1,orbs%norb
      !!            write(333,'(2i8,es12.5,2i10)') irow, jcol, tempmat(jcol, irow), &
      !!            orbs%onwhichatom(irow), orbs%onwhichatom(jcol)
      !!        end do
      !!    end do
      !!end if

      if (methTransformOverlap==-1) then
          call overlap_power_minus_one_half_parallel(iproc, nproc, 0, orbs, ovrlp, inv_ovrlp_half)
      else
          nullify(inv_ovrlp_null)
          call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, orthpar%blocksize_pdgemm, orbs%norb, &
               ovrlp%matrix, inv_ovrlp_null, error, orbs, ovrlp, inv_ovrlp_half)
      end if

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
      call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      allocate(norm(orbs%norb), stat=istat)
      call memocc(istat, norm, 'norm', subname)
      call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

      iall=-product(shape(norm))*kind(norm)
      deallocate(norm, stat=istat)
      call memocc(istat, iall, 'norm', subname)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)
  end do

  !call deallocate_sparseMatrix(inv_ovrlp_half, subname)
  iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)

end subroutine orthonormalizeLocalized


! can still tidy this up more when tmblarge is removed
! use sparsity of density kernel for all inverse quantities
subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
           correction_orthoconstraint, linmat, lphi, lhphi, lagmat, psit_c, psit_f, hpsit_c, hpsit_f, &
           can_use_transposed, overlap_calculated, experimental_mode)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintNonorthogonal
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp
  type(local_zone_descriptors),intent(in) :: lzd
  !type(orbitals_Data),intent(in) :: orbs
  type(orbitals_Data),intent(inout) :: orbs !temporary inout
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
  type(SparseMatrix),intent(inout) :: lagmat
  real(kind=8),dimension(:),pointer :: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
  logical,intent(in) :: experimental_mode

  ! Local variables
  integer :: istat, iall, iorb, jorb, ii, ii_trans, matrixindex_in_compressed, irow, jcol, info, lwork, jj
  !type(SparseMatrix) :: tmp_mat
  real(kind=8),dimension(:),allocatable :: tmp_mat_compr, lagmat_tmp_compr, work
  character(len=*),parameter :: subname='orthoconstraintNonorthogonal'
  real(kind=8),dimension(:,:),allocatable :: tmp_mat, tmp_mat2, tmp_mat3
  integer,dimension(:),allocatable :: ipiv

  ! removed option for correction orthoconstrain for now
  !if (correction_orthoconstraint==0) stop 'correction_orthoconstraint not working'

  if(.not. can_use_transposed) then
      allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psit_c, 'psit_c', subname)

      allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psit_f, 'psit_f', subname)

      call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
      can_use_transposed=.true.
  end if

  ! It is assumed that this routine is called with the transposed gradient ready if it is associated...
  if(.not.associated(hpsit_c)) then
      allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
 
      allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
 
     call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lhphi, hpsit_c, hpsit_f, lzd)
  end if

  call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)

  !call nullify_sparseMatrix(tmp_mat)
  !call sparse_copy_pattern(lagmat,tmp_mat,iproc,subname)
  !allocate(tmp_mat%matrix_compr(tmp_mat%nvctr), stat=istat)
  !call memocc(istat, tmp_mat%matrix_compr, 'tmp_mat%matrix_compr', subname)

  allocate(tmp_mat_compr(lagmat%nvctr), stat=istat) ! save cf doing sparsecopy
  call memocc(istat, tmp_mat_compr, 'tmp_mat_compr', subname)
call timing(iproc,'misc','ON')
  do ii=1,lagmat%nvctr
     iorb = lagmat%orb_from_index(1,ii)
     jorb = lagmat%orb_from_index(2,ii)
     ii_trans=matrixindex_in_compressed(lagmat,jorb, iorb)
     tmp_mat_compr(ii)=-0.5d0*lagmat%matrix_compr(ii)-0.5d0*lagmat%matrix_compr(ii_trans)
     ! SM: This is a hack, should use another variable
     !if (.false..or.correction_orthoconstraint==2) then
     if (.false..and.correction_orthoconstraint==2) then
         if (iproc==0 .and. ii==1) write(*,*) 'only normalization constraint'
         if (iorb/=jorb) then
             tmp_mat_compr(ii)=0.d0
         end if
     end if
     if (experimental_mode) then
         if (iorb==jorb) then
             if (iproc==0 .and. iorb==1) then
                 call yaml_warning('EXPERIMENTAL: modify eval')
                 call yaml_newline()
             end if
             orbs%eval(iorb)=lagmat%matrix_compr(ii)
         end if
     end if
  end do

  ! NEW: reactivate correction for non-orthogonality ##########
  if (correction_orthoconstraint==0) then
      ! WARNING: it is mandatory that the overlap matrix has been calculated before
      !!call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, linmat%ovrlp)
      if (iproc==0) write(*,*) 'correction orthoconstraint'
      allocate(linmat%ovrlp%matrix(orbs%norb,orbs%norb))
      call uncompressMatrix(iproc,linmat%ovrlp)
      allocate(tmp_mat(orbs%norb,orbs%norb))
      allocate(tmp_mat2(orbs%norb,orbs%norb))
      call to_zero(orbs%norb**2, tmp_mat(1,1))
      do ii=1,lagmat%nvctr
         irow = lagmat%orb_from_index(1,ii)
         jcol = lagmat%orb_from_index(2,ii)
         tmp_mat(irow,jcol)=tmp_mat_compr(ii)
      end do
      allocate(ipiv(orbs%norb))
      lwork=10*orbs%norb
      allocate(work(lwork))
      allocate(tmp_mat3(orbs%norb,orbs%norb))
      tmp_mat3=linmat%ovrlp%matrix

      call dgetrf(orbs%norb, orbs%norb, linmat%ovrlp%matrix, orbs%norb, ipiv, info)
      call dgetri(orbs%norb, linmat%ovrlp%matrix, orbs%norb, ipiv, work, lwork, info)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, linmat%ovrlp%matrix, orbs%norb, tmp_mat3, orbs%norb, 0.d0, tmp_mat2, orbs%norb)
      !!if (iproc==0) then
      !!  do iorb=1,orbs%norb
      !!    do jorb=1,orbs%norb
      !!      write(*,'(a,2i8,es14.5)') 'iorb, jorb, tmp_mat2(iorb,jorb)', iorb, jorb, tmp_mat2(iorb,jorb)
      !!    end do
      !!  end do
      !!end if

      ! This is the original
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, linmat%ovrlp%matrix, orbs%norb, &
           tmp_mat, orbs%norb, 0.d0, tmp_mat2, orbs%norb)

      !!! Test
      !!call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, tmp_mat, orbs%norb, linmat%ovrlp%matrix, orbs%norb, 0.d0, tmp_mat2, orbs%norb)

      do jj=1,lagmat%nvctr
         irow = lagmat%orb_from_index(1,jj)
         jcol = lagmat%orb_from_index(2,jj)
         !!if (iproc==0) write(*,'(a,3i8,2es16.6)') 'jj, irow, jcol, tmp_mat_compr(jj), tmp_mat2(irow,jcol)', &
         !!                                          jj, irow, jcol, tmp_mat_compr(jj), tmp_mat2(irow,jcol)
         tmp_mat_compr(jj)=tmp_mat2(irow,jcol)
      end do
      deallocate(linmat%ovrlp%matrix)
      deallocate(tmp_mat)
      deallocate(tmp_mat2)
      deallocate(ipiv)
      deallocate(work)
  end if
  !! ##########################################################

  allocate(lagmat_tmp_compr(lagmat%nvctr), stat=istat) ! save cf doing sparsecopy
  call memocc(istat, lagmat_tmp_compr, 'lagmat_tmp_compr', subname)

  call vcopy(lagmat%nvctr,lagmat%matrix_compr(1),1,lagmat_tmp_compr(1),1) ! need to keep a copy
  call vcopy(lagmat%nvctr,tmp_mat_compr(1),1,lagmat%matrix_compr(1),1)

  iall=-product(shape(tmp_mat_compr))*kind(tmp_mat_compr)
  deallocate(tmp_mat_compr, stat=istat)
  call memocc(istat, iall, 'tmp_mat_compr', subname)

call timing(iproc,'misc','OF')
  call build_linear_combination_transposed(collcom, lagmat, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)
  !call build_linear_combination_transposed(collcom, tmp_mat, psit_c, psit_f, .false., hpsit_c, hpsit_f, iproc)

  !!! TEST ORTHOGONALITY OF GRADIENT AND TMBs ##############################
  !!call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, hpsit_c, psit_f, hpsit_f, linmat%ovrlp)
  !!allocate(linmat%ovrlp%matrix(orbs%norb,orbs%norb))
  !!call uncompressMatrix(iproc,linmat%ovrlp)
  !!if (iproc==0) then
  !!  do iorb=1,orbs%norb
  !!    do jorb=1,orbs%norb
  !!      write(*,'(a,2i8,es16.6)') 'iorb, jorb, linmat%ovrlp%matrix(jorb,iorb)', iorb, jorb, linmat%ovrlp%matrix(jorb,iorb)
  !!    end do
  !!  end do
  !!end if
  !!deallocate(linmat%ovrlp%matrix)
  !!! END TEST #############################################################

  call vcopy(lagmat%nvctr,lagmat_tmp_compr(1),1,lagmat%matrix_compr(1),1)

  iall=-product(shape(lagmat_tmp_compr))*kind(lagmat_tmp_compr)
  deallocate(lagmat_tmp_compr, stat=istat)
  call memocc(istat, iall, 'lagmat_tmp_compr', subname)

  !call deallocate_sparseMatrix(tmp_mat, subname)

  call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, hpsit_c, hpsit_f, lhphi, lzd)

  overlap_calculated=.false.

end subroutine orthoconstraintNonorthogonal


subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: mpisource, mpidest, istsource, istdest, ncount, tag
  integer,dimension(6),intent(out) :: comarr


  ! From which MPI process shall the orbital be sent.
  comarr(1)=mpisource

  ! Starting index on the sending process.
  comarr(2)=istsource

  ! Amount of datat to be sent
  comarr(3)=ncount

  ! To which MPI process shall the orbital be sent.
  comarr(4)=mpidest

  ! Starting index on the receiving process.
  comarr(5)=istdest

  ! Tag for the communication
  comarr(6)=tag

  ! comarr(7): this entry is used as request for the mpi_isend.

  ! comarr(8): this entry is used as request for the mpi_irecv.

end subroutine setCommsParameters


!S^-1 exact only works for symmetric matrices
!BOTH sparse matrices must be present together and inv_ovrlp should be nullified pointer, NOT inv_ovrlp_smat%matrix
!when sparse matrices present, check is performed to see whether %matrix is allocated so that its allocated status remains unchanged
!contents of %matrix not guaranteed to be correct though - inv_ovrlp_smat%can_use_dense set accordingly
!power: -2 -> S^-1/2, 2 -> S^1/2, 1 -> S^-1
subroutine overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, norb, ovrlp, inv_ovrlp, error, orbs, &
     ovrlp_smat, inv_ovrlp_smat)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlapPowerGeneral
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, norb, power
  real(kind=8),dimension(:,:),pointer :: ovrlp
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp
  real(kind=8),intent(out) :: error
  type(orbitals_data), optional, intent(in) :: orbs
  type(sparseMatrix), optional, intent(inout) :: ovrlp_smat, inv_ovrlp_smat
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i
  integer :: matrixindex_in_compressed
  real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated
  logical, parameter :: check_accuracy=.false. ! move to input.perf


  call f_routine(id='overlapPowerGeneral')
  call timing(iproc,'lovrlp^-1     ','ON')

  if (power/=1.and.power/=2.and.power/=-2) stop 'Error in overlapPowerGeneral'

  if (present(inv_ovrlp_smat)) then
     ovrlp_allocated=associated(ovrlp_smat%matrix)
     inv_ovrlp_allocated=associated(inv_ovrlp_smat%matrix)
  end if

  if(iorder==0) then
     ! exact inversion
     if (present(inv_ovrlp_smat)) then
        if (.not. ovrlp_allocated) ovrlp_smat%matrix=f_malloc_ptr((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/),&
             id='ovrlp_smat%matrix')
        if (.not.(ovrlp_smat%can_use_dense.and.ovrlp_allocated)) call uncompressMatrix(iproc,ovrlp_smat)
        if (.not. inv_ovrlp_allocated) inv_ovrlp_smat%matrix=f_malloc_ptr((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%nfvctr/),&
             id='inv_ovrlp_smat%matrix')
        call vcopy(norb*norb,ovrlp_smat%matrix(1,1),1,inv_ovrlp_smat%matrix(1,1),1)
        if (.not. ovrlp_allocated) call f_free_ptr(ovrlp_smat%matrix)
        inv_ovrlp => inv_ovrlp_smat%matrix
     else
        call vcopy(norb*norb,ovrlp(1,1),1,inv_ovrlp(1,1),1)
     end if

     if (power==1) then
        if (blocksize<0) then
           call overlap_minus_one_exact_serial(norb,inv_ovrlp)
        else
           stop 'check if working - upper half may not be filled'
           call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
           call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
        end if
     else if (power==2) then
        call overlap_plus_minus_one_half_exact(norb,blocksize,.true.,inv_ovrlp,orbs)
     else if (power==-2) then
        call overlap_plus_minus_one_half_exact(norb,blocksize,.false.,inv_ovrlp,orbs)
     end if
     if (present(inv_ovrlp_smat)) then
        call compress_matrix_for_allreduce(iproc,inv_ovrlp_smat)
        nullify(inv_ovrlp)
        if (.not. inv_ovrlp_allocated) then
           call f_free_ptr(inv_ovrlp_smat%matrix)
        else
           inv_ovrlp_smat%can_use_dense=.true.
        end if
     end if
  else if(iorder>0) then
     ! Taylor expansion
     if (present(inv_ovrlp_smat)) then
        call first_order_taylor_sparse(power,ovrlp_smat,inv_ovrlp_smat)
        if (inv_ovrlp_smat%parallel_compression==1) then
           call mpiallred(inv_ovrlp_smat%matrix_compr(1), inv_ovrlp_smat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
        else if (inv_ovrlp_smat%parallel_compression==2) then
           ! could move this til after if we have higher order Taylor for sparse
           call mpi_allgatherv(inv_ovrlp_smat%matrix_comprp, inv_ovrlp_smat%nvctrp, mpi_double_precision, &
                inv_ovrlp_smat%matrix_compr, inv_ovrlp_smat%nvctr_par(:), inv_ovrlp_smat%isvctr_par, &
                mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
           call f_free_ptr(inv_ovrlp_smat%matrix_comprp)
        end if
        if (iorder>1) then
           if (.not. inv_ovrlp_allocated) inv_ovrlp_smat%matrix=f_malloc_ptr((/norb,norb/),id='inv_ovrlp_smat%matrix')
           call uncompressMatrix(iproc,inv_ovrlp_smat)
           inv_ovrlp => inv_ovrlp_smat%matrix
           if (.not. ovrlp_allocated) ovrlp_smat%matrix=f_malloc_ptr((/norb,norb/),id='ovrlp_smat%matrix')
           call uncompressMatrix(iproc,ovrlp_smat)
           ovrlp => ovrlp_smat%matrix
        end if
     end if
     if (.not.present(inv_ovrlp_smat).or.iorder>1) then 
        if (present(orbs).and.nproc>1) then
           if (norb/=orbs%norb) stop 'Error with orbs%norb in overlapPowerGeneral'
           norbp=orbs%norbp
           isorb=orbs%isorb
           inv_ovrlpp=f_malloc_ptr((/norb,norbp/),id='inv_ovrlpp')
           if (present(inv_ovrlp_smat)) call vcopy(norb*norbp,inv_ovrlp_smat%matrix(1,isorb+1),1,inv_ovrlpp(1,1),1)
        else
           norbp=norb
           isorb=0
           inv_ovrlpp => inv_ovrlp
        end if
        if (.not.present(inv_ovrlp_smat)) call first_order_taylor_dense(norb,isorb,norbp,power,ovrlp(1,isorb+1),inv_ovrlpp)
     end if

     ! add sparse here once we have sparse matrix multiply
     if (iorder>1) then
        ovrlpminone=f_malloc_ptr((/norb,norb/),id='ovrlpminone')
        if (present(orbs).and.nproc>1) then
           ovrlpminonep=f_malloc_ptr((/norb,norbp/),id='ovrlpminonep')
        else
           ovrlpminonep => ovrlpminone
        end if

        call matrix_minus_identity_dense(norb,isorb,norbp,ovrlp(1,isorb+1),ovrlpminonep)

        ovrlppoweroldp=f_malloc_ptr((/norb,norbp/),id='ovrlppoweroldp')
        call vcopy(norb*norbp,ovrlpminonep(1,1),1,ovrlppoweroldp(1,1),1)

        if(nproc > 1.and.present(orbs)) then
           call mpi_allgatherv(ovrlpminonep, norb*norbp, mpi_double_precision, ovrlpminone, &
                norb*orbs%norb_par(:,0), norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
           call f_free_ptr(ovrlpminonep)
        else
           nullify(ovrlpminonep)
        end if

        ovrlppowerp=f_malloc_ptr((/norb,norbp/),id='ovrlppowerp')

        if (power==1) then
           factor=-1.0d0
        else if (power==2) then
           factor=0.5d0
        else if (power==-2) then
           factor=-0.5d0
        end if
        do i=2,iorder
           call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlpminone(1,1), &
                norb, ovrlppoweroldp(1,1), norb, 0.d0, ovrlppowerp(1,1), norb)
           factor=newfactor(power,i,factor)
           call daxpy(norb*norbp,factor,ovrlppowerp,1,inv_ovrlpp,1)
           call vcopy(norb*norbp,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1),1)
        end do

        call f_free_ptr(ovrlppowerp)
        call f_free_ptr(ovrlppoweroldp)
        call f_free_ptr(ovrlpminone)
     end if

     if (.not. present(inv_ovrlp_smat).or.iorder>1) then
        if(nproc > 1.and.present(orbs)) then
           call mpi_allgatherv(inv_ovrlpp, norb*norbp, mpi_double_precision, inv_ovrlp, &
                norb*orbs%norb_par(:,0), norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
           call f_free_ptr(inv_ovrlpp)
        else
           nullify(inv_ovrlpp)
        end if
        if (present(inv_ovrlp_smat)) then
           call compress_matrix_for_allreduce(iproc,inv_ovrlp_smat)
           if (.not. inv_ovrlp_allocated) then
              call f_free_ptr(inv_ovrlp_smat%matrix)
           else
              inv_ovrlp_smat%can_use_dense=.true.
           end if
           if (.not. ovrlp_allocated) then
              call f_free_ptr(ovrlp_smat%matrix)
           end if
        end if
     end if
  end if

  ! check how accurate calculation was
  if (check_accuracy) then
     if (present(inv_ovrlp_smat)) then
        if (.not.inv_ovrlp_allocated) inv_ovrlp_smat%matrix=f_malloc_ptr((/norb,norb/),id='inv_ovrlp_smat%matrix')
        if (.not.(inv_ovrlp_smat%can_use_dense.and.inv_ovrlp_allocated)) call uncompressMatrix(iproc,inv_ovrlp_smat)
        if (.not.ovrlp_allocated) ovrlp_smat%matrix=f_malloc_ptr((/norb,norb/),id='ovrlp_smat%matrix')
        if (.not.(ovrlp_smat%can_use_dense.and.ovrlp_allocated)) call uncompressMatrix(iproc,ovrlp_smat)

        call check_accuracy_overlap_minus_one(iproc,norb,power,ovrlp_smat%matrix,inv_ovrlp_smat%matrix,error)
        if (.not.ovrlp_allocated) call f_free_ptr(ovrlp_smat%matrix)
        if (.not.inv_ovrlp_allocated) call f_free_ptr(inv_ovrlp_smat%matrix)
     else
        call check_accuracy_overlap_minus_one(iproc,norb,power,ovrlp,inv_ovrlp,error)
     end if
     !if (iproc==0) print*,'Accuracy of inverse overlap calculation: power, order, error',power,iorder,error
  else
     error=-1.0d0
  end if

  call timing(iproc,'lovrlp^-1     ','OF')
  call f_release_routine()

end subroutine overlapPowerGeneral


pure function newfactor(power,order,factor)
  implicit none
  integer, intent(in) :: power, order
  real(kind=8), intent(in) :: factor
  real(kind=8) :: newfactor

  if (power==1) then
     newfactor=-factor
  else if (power==2) then
     newfactor=factor*(1.5d0-real(order,kind=8))/real(order,kind=8)
  else if (power==-2) then
     newfactor=factor*(0.5d0-real(order,kind=8))/real(order,kind=8)
  end if

end function newfactor


subroutine matrix_minus_identity_dense(norb,isorb,norbp,matinp,matoutp)
  implicit none
  integer,intent(in) :: norb, isorb, norbp
  real(kind=8),dimension(norb,norbp),intent(in) :: matinp
  real(kind=8),dimension(norb,norbp),intent(out) :: matoutp

  integer :: iorb, jorb, iiorb

  !$omp parallel do default(private) shared(matinp,matoutp,norb,isorb,norbp)
  do iorb=1,norbp
     iiorb=isorb+iorb
     do jorb=1,norb
        if(iiorb==jorb) then
           matoutp(jorb,iorb) = matinp(jorb,iorb) - 1.0d0
        else
           matoutp(jorb,iorb) = matinp(jorb,iorb)
        end if
     end do
  end do
  !$omp end parallel do

end subroutine matrix_minus_identity_dense

subroutine first_order_taylor_sparse(power,ovrlp,inv_ovrlp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: power
  type(sparseMatrix),intent(in) :: ovrlp
  type(sparseMatrix),intent(out) :: inv_ovrlp

  integer :: ii,iii,iorb,jorb,ii_inv,matrixindex_in_compressed,ierr,iii_inv

  if (power/=1 .and. power/=2 .and. power/=-2) stop 'Error in first_order_taylor_sparse'

  if (inv_ovrlp%parallel_compression==0.or.bigdft_mpi%nproc==1) then
     ! inv_ovrlp can be less sparse than ovrlp, so pad with zeros first
     call to_zero(inv_ovrlp%nvctr,inv_ovrlp%matrix_compr(1))
     if (power==1) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctr
           iorb = ovrlp%orb_from_index(1,ii)
           jorb = ovrlp%orb_from_index(2,ii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 2.0d0 - ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = -ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     else if (power==2) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctr
           iorb = ovrlp%orb_from_index(1,ii)
           jorb = ovrlp%orb_from_index(2,ii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0*ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     else
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctr
           iorb = ovrlp%orb_from_index(1,ii)
           jorb = ovrlp%orb_from_index(2,ii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = -0.5d0*ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     end if
  else if (inv_ovrlp%parallel_compression==1) then
     call to_zero(inv_ovrlp%nvctr,inv_ovrlp%matrix_compr(1))
     if (power==1) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctrp
           iii=ii+ovrlp%isvctr
           iorb = ovrlp%orb_from_index(1,iii)
           jorb = ovrlp%orb_from_index(2,iii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 2.0d0 - ovrlp%matrix_compr(iii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = -ovrlp%matrix_compr(iii)
           end if
        end do
        !$omp end parallel do
     else if (power==2) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctrp
           iii=ii+ovrlp%isvctr
           iorb = ovrlp%orb_from_index(1,iii)
           jorb = ovrlp%orb_from_index(2,iii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(iii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = 0.5d0*ovrlp%matrix_compr(iii)
           end if
        end do
        !$omp end parallel do
     else
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii=1,ovrlp%nvctrp
           iii=ii+ovrlp%isvctr
           iorb = ovrlp%orb_from_index(1,iii)
           jorb = ovrlp%orb_from_index(2,iii)
           ii_inv = matrixindex_in_compressed(inv_ovrlp,iorb,jorb)
           if(iorb==jorb) then
              inv_ovrlp%matrix_compr(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(iii)
           else
              inv_ovrlp%matrix_compr(ii_inv) = -0.5d0*ovrlp%matrix_compr(iii)
           end if
        end do
        !$omp end parallel do
     end if
  else
     inv_ovrlp%matrix_comprp=f_malloc_ptr((inv_ovrlp%nvctrp),id='inv_ovrlp%matrix_comprp')
     if (power==1) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii_inv=1,inv_ovrlp%nvctrp
           iii_inv=ii_inv+inv_ovrlp%isvctr
           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
           if (ii==0) then
              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
           else if(iorb==jorb) then
              inv_ovrlp%matrix_comprp(ii_inv) = 2.0d0 - ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_comprp(ii_inv) = -ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     else if (power==2) then
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii_inv=1,inv_ovrlp%nvctrp
           iii_inv=ii_inv+inv_ovrlp%isvctr
           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
           if (ii==0) then
              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
           else if(iorb==jorb) then
              inv_ovrlp%matrix_comprp(ii_inv) = 0.5d0 + 0.5d0*ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_comprp(ii_inv) = 0.5d0*ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     else
        !$omp parallel do default(private) shared(inv_ovrlp,ovrlp)
        do ii_inv=1,inv_ovrlp%nvctrp
           iii_inv=ii_inv+inv_ovrlp%isvctr
           iorb = inv_ovrlp%orb_from_index(1,iii_inv)
           jorb = inv_ovrlp%orb_from_index(2,iii_inv)
           ii = matrixindex_in_compressed(ovrlp,iorb,jorb)
           if (ii==0) then
              inv_ovrlp%matrix_comprp(ii_inv) = 0.0d0
           else if(iorb==jorb) then
              inv_ovrlp%matrix_comprp(ii_inv) = 1.5d0 - 0.5d0*ovrlp%matrix_compr(ii)
           else
              inv_ovrlp%matrix_comprp(ii_inv) = -0.5d0*ovrlp%matrix_compr(ii)
           end if
        end do
        !$omp end parallel do
     end if
  end if

end subroutine first_order_taylor_sparse

subroutine first_order_taylor_dense(norb,isorb,norbp,power,ovrlpp,inv_ovrlpp)
use module_base
  implicit none
  integer,intent(in) :: norb, isorb, norbp, power
  real(kind=8),dimension(norb,norbp),intent(in) :: ovrlpp
  real(kind=8),dimension(norb,norbp),intent(out) :: inv_ovrlpp

  integer :: iorb, jorb, iiorb

  if (power==1) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 2.d0 - ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = -ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else if (power==2) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 0.5d0 + 0.5d0*ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = 0.5d0*ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else if (power==-2) then
     !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
     do iorb=1,norbp
        iiorb=isorb+iorb
        do jorb=1,norb
           if(iiorb==jorb) then
              inv_ovrlpp(jorb,iorb) = 1.5d0 - 0.5d0*ovrlpp(jorb,iorb)
           else
              inv_ovrlpp(jorb,iorb) = -0.5d0*ovrlpp(jorb,iorb)
           end if
        end do
     end do
     !$omp end parallel do
  else
     stop 'Error in first_order_taylor_dense'
  end if

end subroutine first_order_taylor_dense


subroutine overlap_minus_one_exact_serial(norb,inv_ovrlp)
  use module_base
  implicit none
  integer,intent(in) :: norb
  real(kind=8),dimension(norb,norb),intent(inout) :: inv_ovrlp

  integer :: info, iorb, jorb

  call dpotrf('u', norb, inv_ovrlp(1,1), norb, info)
  if(info/=0) then
     write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
     stop
  end if
  call dpotri('u', norb, inv_ovrlp(1,1), norb, info)
  if(info/=0) then
     write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
     stop
  end if

  ! fill lower half
  !$omp parallel do default(private) shared(inv_ovrlp,norb)
  do iorb=1,norb
     do jorb=1,iorb-1
        inv_ovrlp(iorb,jorb)=inv_ovrlp(jorb,iorb)
     end do
  end do
  !$omp end parallel do 

end subroutine overlap_minus_one_exact_serial


subroutine overlap_plus_minus_one_half_exact(norb,blocksize,plusminus,inv_ovrlp_half,orbs)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: norb,blocksize
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_half
  logical, intent(in) :: plusminus
  type(orbitals_data), optional, intent(in) :: orbs

  integer :: info, iorb, jorb, ierr, iiorb, isorb, norbp, lwork, jjorb
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: tempArr, orig_ovrlp
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_halfp
  real(kind=8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
  real(kind=8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
  real(dp) :: temp, error
  real(dp), allocatable, dimension(:) :: temp_vec
  logical, parameter :: symmetric=.true.
  logical, parameter :: check_lapack=.false.

  eval=f_malloc(norb,id='eval')
  if(blocksize>0) then
     call dsyev_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, min(blocksize,norb), bigdft_mpi%mpi_comm, 'v', 'l', norb, &
          inv_ovrlp_half(1,1), norb, eval(1), info)
     if(info/=0) then
        write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
     end if
  else
     if (symmetric) then
        if (check_lapack) then
           orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
           call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
        end if
        work=f_malloc(1000,id='work')
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, -1, info)
        lwork = int(work(1))
        call f_free(work)
        work=f_malloc(lwork,id='work')
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)
        if (check_lapack) then
           tempArr=f_malloc((/norbp,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           if (norbp>0) call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norbp, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0) print*,'LAPACK error',error
           call f_free_ptr(inv_ovrlp_halfp)
           call f_free(orig_ovrlp)
        end if
     else
        if (check_lapack) then
           orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
           call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
        end if
        work=f_malloc(1000,id='work')
        call dgeev( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
             norb, WORK, -1, info )
        lwork = nint(work(1))
        call f_free(work)
        call f_free(work)
        vl=f_malloc((/norb,norb/),id='vl')
        vr=f_malloc((/norb,norb/),id='vr')
        eval1=f_malloc(norb,id='eval1')
        call DGEEV( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
             norb, WORK, LWORK, info )
        call vcopy(norb*norb,vl(1,1),1,inv_ovrlp_half(1,1),1)
        call f_free(eval1)
        call f_free(vr)
        call f_free(vl)
        temp_vec=f_malloc(norb,id='temp_vec')
        do iorb=1,norb
           do jorb=iorb+1,norb
              if (eval(jorb) < eval(iorb)) then
                 temp = eval(iorb)
                 temp_vec = inv_ovrlp_half(:,iorb)
                 eval(iorb) = eval(jorb)
                 eval(jorb) = temp
                 inv_ovrlp_half(:,iorb) = inv_ovrlp_half(:,jorb)
                 inv_ovrlp_half(:,jorb) = temp_vec
              end if
           end do
        end do
        call f_free(temp_vec)
        if (check_lapack) then
           tempArr=f_malloc((/norbp,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           if (norbp>0) call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norbp, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0) print*,'LAPACK error',error
           call f_free_ptr(inv_ovrlp_halfp)
           call f_free(orig_ovrlp)
        end if
     end if

     if(info/=0) then
        write(*,'(a,i0,2x,i0)') 'ERROR in dsyev (overlap_plus_minus_one_half_exact), info, norb=', info, norb
        stop
     end if
     call f_free(work)
  end if

  ! Calculate S^{-1/2}. 
  ! First calculate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
  ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
  if (present(orbs).and.bigdft_mpi%nproc>1.and.blocksize<0) then
     if (norb/=orbs%norb) stop 'Error with orbs%norb in overlap_plus_minus_one_half_exact'
     norbp=orbs%norbp
     isorb=orbs%isorb
  else
     norbp=norb
     isorb=0
  end if
  tempArr=f_malloc((/norbp,norb/), id='tempArr')
  !$omp parallel do default(private) shared(tempArr,inv_ovrlp_half,eval,plusminus,norb,norbp,isorb)
  do iorb=1,norb
     do jorb=1,norbp
        jjorb=jorb+isorb
        if (plusminus) then
           tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*sqrt(abs(eval(iorb)))
        else
           tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
        end if
     end do
  end do
  !$omp end parallel do

  call f_free(eval)
  inv_ovrlp_halfp=f_malloc_ptr((/norb,norbp/), id='inv_ovrlp_halfp')
  ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
  ! This will give S^{+/-1/2}.
  if(blocksize<0) then
     if (norbp>0) call dgemm('n', 't', norb, norbp, norb, 1.d0, inv_ovrlp_half, &
          norb, tempArr, norbp, 0.d0, inv_ovrlp_halfp, norb)
     if (present(orbs).and.bigdft_mpi%nproc>1) then
        call mpi_allgatherv(inv_ovrlp_halfp, norb*norbp, mpi_double_precision, inv_ovrlp_half, &
                   norb*orbs%norb_par(:,0), norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
     end if
  else
     call dgemm_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, bigdft_mpi%mpi_comm, 'n', 't', norb, norb, norb, &
          1.d0, inv_ovrlp_half, norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
     call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
  end if

  call f_free_ptr(inv_ovrlp_halfp)
  call f_free(tempArr)

end subroutine overlap_plus_minus_one_half_exact



subroutine check_accuracy_overlap_minus_one(iproc,norb,power,ovrlp,inv_ovrlp,error)
  use module_base
  implicit none
  integer,intent(in) :: iproc, norb, power
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
  real(kind=8),intent(out) :: error

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2

  tmp=f_malloc((/norb,norb/),id='tmp')
  if (power==1) then
     call dgemm('n', 'n', norb, norb, norb, 1.d0, inv_ovrlp(1,1), &
          norb, ovrlp(1,1), norb, 0.d0, tmp(1,1), norb)
     call deviation_from_unity(iproc, norb, tmp, error)
  else if (power==2) then
     call dgemm('n', 'n', norb, norb, norb, 1.d0, inv_ovrlp(1,1), &
          norb, inv_ovrlp(1,1), norb, 0.d0, tmp(1,1), norb)
     call max_matrix_diff(iproc, norb, tmp, ovrlp, error)
     error=0.5d0*error
  else if (power==-2) then
     call dgemm('n', 'n', norb, norb, norb, 1.d0, inv_ovrlp(1,1), &
          norb, inv_ovrlp(1,1), norb, 0.d0, tmp(1,1), norb)
     tmp2=f_malloc((/norb,norb/),id='tmp2')
     call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlp(1,1), &
          norb, tmp(1,1), norb, 0.d0, tmp2(1,1), norb)
     call deviation_from_unity(iproc, norb, tmp2, error)
     error=0.5d0*error
     call f_free(tmp2)
  else
     stop 'Error in check_accuracy_overlap_minus_one'
  end if
  call f_free(tmp)

end subroutine check_accuracy_overlap_minus_one


subroutine max_matrix_diff(iproc, norb, mat1, mat2, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: mat1, mat2
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=1,norb
        error=abs(mat1(jorb,iorb)-mat2(jorb,iorb))
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff

subroutine deviation_from_unity(iproc, norb, ovrlp, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: ovrlp
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=1,norb
        if(iorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_unity

subroutine deviation_from_symmetry(iproc, norb, ovrlp, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: ovrlp
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norb
     do jorb=iorb+1,norb
        error=abs(ovrlp(jorb,iorb)-ovrlp(iorb,jorb))
        deviation=max(error,deviation)
     end do
  end do
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_symmetry

subroutine overlap_power_minus_one_half_parallel(iproc, nproc, meth_overlap, orbs, ovrlp, inv_ovrlp_half)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, meth_overlap
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, ii, iseg, ind, matrixindex_in_compressed
  real(kind=8) :: error
  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
  logical,dimension(:),allocatable :: in_neighborhood
  character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'

  call timing(iproc,'lovrlp^-1/2par','ON')

  allocate(in_neighborhood(orbs%norb), stat=istat)
  call memocc(istat, in_neighborhood, 'in_neighborhood', subname)
  call to_zero(inv_ovrlp_half%nvctr, inv_ovrlp_half%matrix_compr(1))

  !DEBUG
  !if (iproc==0) then
  !   jjorb=0
  !   do jorb=1,orbs%norb
  !      do korb=1,orbs%norb
  !         ind = ovrlp%matrixindex_in_compressed(korb, jorb)
  !         if (ind>0) then
  !            print*, korb,jorb,ovrlp%matrix_compr(ind)
  !         else
  !            print*, korb,jorb,0.0d0
  !         end if
  !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
  !      end do
  !   end do
  !end if
  !call mpi_barrier(bigdft_mpi%mpi_comm,istat)

  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     !ilr=orbs%inwhichlocreg(iiorb)
     ! We are at the start of a new atom
     ! Count all orbitals that are in the neighborhood

     iseg=inv_ovrlp_half%istsegline(iiorb)
     iend=iiorb*orbs%norb
     n=0
     in_neighborhood(:)=.false.
     do 
        do i=inv_ovrlp_half%keyg(1,iseg),inv_ovrlp_half%keyg(2,iseg)
           ii=i-(iiorb-1)*orbs%norb
           in_neighborhood(ii)=.true.
           n=n+1
        end do
        iseg=iseg+1
        if (iseg>inv_ovrlp_half%nseg) exit
        if (inv_ovrlp_half%keyg(1,iseg)>iend) exit
     end do

     allocate(ovrlp_tmp(n,n), stat=istat)
     call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
     call to_zero(n*n, ovrlp_tmp(1,1))

     jjorb=0
     do jorb=1,orbs%norb
        if (.not.in_neighborhood(jorb)) cycle
        jjorb=jjorb+1
        kkorb=0
        do korb=1,orbs%norb
           if (.not.in_neighborhood(korb)) cycle
           kkorb=kkorb+1
           ind = matrixindex_in_compressed(ovrlp,korb, jorb)
           if (ind>0) then
              ovrlp_tmp(kkorb,jjorb)=ovrlp%matrix_compr(ind)
           else
              ovrlp_tmp(kkorb,jjorb)=0.d0
           end if
           !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
        end do
     end do
          
     allocate(ovrlp_tmp_inv_half(n,n))
     call memocc(istat, ovrlp_tmp_inv_half, 'ovrlp_tmp_inv_half', subname)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'ovrlp_tmp',n,iiorb
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp(:,jorb)
     !end do
     !end if

     ! Calculate S^-1/2 for the small overlap matrix
     call overlapPowerGeneral(iproc, nproc, meth_overlap, -2, -8, n, ovrlp_tmp, ovrlp_tmp_inv_half, error)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'inv_ovrlp_tmp',n,iiorb,error
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp_inv_half(:,jorb)
     !end do
     !end if

     jjorb=0
     do jorb=1,orbs%norb
        if (.not.in_neighborhood(jorb)) cycle
        jjorb=jjorb+1
        kkorb=0
        if (jorb==iiorb) then
           do korb=1,orbs%norb
              if (.not.in_neighborhood(korb)) cycle
              kkorb=kkorb+1
              ind = matrixindex_in_compressed(inv_ovrlp_half,korb,jorb)
              if (ind>0) then
                 inv_ovrlp_half%matrix_compr(ind)=ovrlp_tmp_inv_half(kkorb,jjorb)
                 !if (iiorb==orbs%norb) print*,'problem here?!',iiorb,kkorb,jjorb,korb,jorb,ind,ovrlp_tmp_inv_half(kkorb,jjorb)
              end if
              !write(1300+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
           end do
           exit !no need to keep looping
        end if
     end do

     iall=-product(shape(ovrlp_tmp_inv_half))*kind(ovrlp_tmp_inv_half)
     deallocate(ovrlp_tmp_inv_half, stat=istat)
     call memocc(istat, iall, 'ovrlp_tmp_inv_half', subname)

     iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
     deallocate(ovrlp_tmp, stat=istat)
     call memocc(istat, iall, 'ovrlp_tmp', subname)
  end do

  call mpiallred(inv_ovrlp_half%matrix_compr(1), inv_ovrlp_half%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  iall=-product(shape(in_neighborhood))*kind(in_neighborhood)
  deallocate(in_neighborhood, stat=istat)
  call memocc(istat, iall, 'in_neighborhood', subname)

  !if (iproc==0) then
  !   jjorb=0
  !   do jorb=1,orbs%norb
  !      do korb=1,orbs%norb
  !         ind = inv_ovrlp_half%matrixindex_in_compressed(korb, jorb)
  !         if (ind>0) then
  !            print*, korb,jorb,inv_ovrlp_half%matrix_compr(ind)
  !         else
  !            print*, korb,jorb,0.0d0
  !         end if
  !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
  !      end do
  !   end do
  !end if
  !call mpi_finalize(ind)
  !stop

  call timing(iproc,'lovrlp^-1/2par','OF')

end subroutine overlap_power_minus_one_half_parallel

subroutine orthonormalize_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalize_subset
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(inout) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii, matrixindex_in_compressed
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparseMatrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='orthonormalize_subset'
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparsematrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)

  do it=1,orthpar%nItOrtho

      if(.not.can_use_transposed) then
          if(associated(psit_c)) then
              iall=-product(shape(psit_c))*kind(psit_c)
              deallocate(psit_c, stat=istat)
              call memocc(istat, iall, 'psit_c', subname)
          end if
          if(associated(psit_f)) then
              iall=-product(shape(psit_f))*kind(psit_f)
              deallocate(psit_f, stat=istat)
              call memocc(istat, iall, 'psit_f', subname)
          end if
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)

          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.

      end if
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)

      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      allocate(icount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,icount_norb,'icount_norb',subname)
      allocate(jcount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,jcount_norb,'jcount_norb',subname)
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
                      ovrlp%matrix_compr(ii)=1.d0
                  else
                      ovrlp%matrix_compr(ii)=0.d0
                  end if
              end if
          end do
      end do
      iall=-product(shape(icount_norb))*kind(icount_norb)
      deallocate(icount_norb, stat=istat)
      call memocc(istat, iall, 'icount_norb', subname)
      iall=-product(shape(jcount_norb))*kind(jcount_norb)
      deallocate(jcount_norb, stat=istat)
      call memocc(istat, iall, 'jcount_norb', subname)


      if (methTransformOverlap==-1) then
          call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, orbs, ovrlp, inv_ovrlp_half)
      else
          nullify(inv_ovrlp_null)
          call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, orthpar%blocksize_pdsyev, orbs%norb, &
               ovrlp%matrix, inv_ovrlp_null, error, orbs, ovrlp, inv_ovrlp_half)
      end if

      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      allocate(icount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,icount_norb,'icount_norb',subname)
      allocate(jcount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,jcount_norb,'jcount_norb',subname)
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
                      ovrlp%matrix_compr(ii)=1.d0
                  else
                      ovrlp%matrix_compr(ii)=0.d0
                  end if
              end if
          end do
      end do
      iall=-product(shape(icount_norb))*kind(icount_norb)
      deallocate(icount_norb, stat=istat)
      call memocc(istat, iall, 'icount_norb', subname)
      iall=-product(shape(jcount_norb))*kind(jcount_norb)
      deallocate(jcount_norb, stat=istat)
      call memocc(istat, iall, 'jcount_norb', subname)

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)

      call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
           psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      allocate(norm(orbs%norb), stat=istat)
      call memocc(istat, norm, 'norm', subname)
      call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

      iall=-product(shape(norm))*kind(norm)
      deallocate(norm, stat=istat)
      call memocc(istat, iall, 'norm', subname)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)
  end do

  !call deallocate_sparseMatrix(inv_ovrlp_half, subname)
  iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)

end subroutine orthonormalize_subset



subroutine gramschmidt_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gramschmidt_subset
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparseMatrix),intent(inout) :: ovrlp
  type(sparseMatrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(collective_comms),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii, matrixindex_in_compressed
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparseMatrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='gramschmidt_subset'

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparsematrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)

  do it=1,orthpar%nItOrtho

      if(.not.can_use_transposed) then
          if(associated(psit_c)) then
              iall=-product(shape(psit_c))*kind(psit_c)
              deallocate(psit_c, stat=istat)
              call memocc(istat, iall, 'psit_c', subname)
          end if
          if(associated(psit_f)) then
              iall=-product(shape(psit_f))*kind(psit_f)
              deallocate(psit_f, stat=istat)
              call memocc(istat, iall, 'psit_f', subname)
          end if
          allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, psit_c, 'psit_c', subname)
          allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, psit_f, 'psit_f', subname)

          call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, lphi, psit_c, psit_f, lzd)
          can_use_transposed=.true.

      end if
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)

      ! For the "higher" TMBs: delete off-diagonal elements and
      ! set diagonal elements to 1
      allocate(icount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,icount_norb,'icount_norb',subname)
      allocate(jcount_norb(at%astruct%nat),stat=istat)
      call memocc(istat,jcount_norb,'jcount_norb',subname)
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
                      ovrlp%matrix_compr(ii)=0.d0
                  else
                      if (jout) then
                          ovrlp%matrix_compr(ii)=-ovrlp%matrix_compr(ii)
                      else
                          ovrlp%matrix_compr(ii)=0.d0
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
      iall=-product(shape(icount_norb))*kind(icount_norb)
      deallocate(icount_norb, stat=istat)
      call memocc(istat, iall, 'icount_norb', subname)
      iall=-product(shape(jcount_norb))*kind(jcount_norb)
      deallocate(jcount_norb, stat=istat)
      call memocc(istat, iall, 'jcount_norb', subname)


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

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
      call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
      !!call build_linear_combination_transposed(collcom, inv_ovrlp_half, &
      !!     psittemp_c, psittemp_f, .true., psit_c, psit_f, iproc)
      call build_linear_combination_transposed(collcom, ovrlp, &
           psittemp_c, psittemp_f, .false., psit_c, psit_f, iproc)
      allocate(norm(orbs%norb), stat=istat)
      call memocc(istat, norm, 'norm', subname)
      !call normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)

      iall=-product(shape(norm))*kind(norm)
      deallocate(norm, stat=istat)
      call memocc(istat, iall, 'norm', subname)
      call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, lphi, lzd)

      iall=-product(shape(psittemp_c))*kind(psittemp_c)
      deallocate(psittemp_c, stat=istat)
      call memocc(istat, iall, 'psittemp_c', subname)
      iall=-product(shape(psittemp_f))*kind(psittemp_f)
      deallocate(psittemp_f, stat=istat)
      call memocc(istat, iall, 'psittemp_f', subname)
  end do

  !call deallocate_sparseMatrix(inv_ovrlp_half, subname)
  iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)

end subroutine gramschmidt_subset
