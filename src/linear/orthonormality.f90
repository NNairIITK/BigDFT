!> @file
!! Orthonormalization
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Orthonormalized the localized orbitals
subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed, foe_obj)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalizeLocalized
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: compress_matrix, uncompress_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed
  type(foe_data),intent(in) :: foe_obj

  ! Local variables
  integer :: it, istat, iall
  !integer :: irow, ii, iorb, jcol, jorb
  real(kind=8), dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*), parameter :: subname='orthonormalizeLocalized'
  !real(kind=8), dimension(orbs%norb,orbs%norb) :: tempmat
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error
  logical :: ovrlp_associated, inv_ovrlp_associated

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')

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
          !!nullify(inv_ovrlp_null)
          !!if (methTransformOverlap>1 .and.  .not.associated(foe_obj%nsegline)) then
          !!    ! this is a bit cheating, to be improved
          !!    ovrlp_associated=associated(ovrlp%matrix)
          !!    inv_ovrlp_associated=associated(inv_ovrlp_half%matrix)
          !!    if (.not.ovrlp_associated) then
          !!        ovrlp%matrix=f_malloc_ptr((/orbs%norb,orbs%norb/),id='ovrlp%matrix')
          !!    end if
          !!    if (.not.inv_ovrlp_associated) then
          !!        inv_ovrlp_half%matrix=f_malloc_ptr((/orbs%norb,orbs%norb/),id='inv_ovrlp_half%matrix')
          !!    end if
          !!    call uncompress_matrix(iproc, ovrlp)
          !!    call uncompress_matrix(iproc, inv_ovrlp_half)
          !!    call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, &
          !!         orthpar%blocksize_pdgemm, orbs%norb, orbs, &
          !!         imode=2, check_accur=.false.,ovrlp=ovrlp%matrix, &
          !!         inv_ovrlp=inv_ovrlp_half%matrix, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half)
          !!    call compress_matrix(iproc, ovrlp)
          !!    call compress_matrix(iproc, inv_ovrlp_half)
          !!    if (.not.ovrlp_associated) call f_free_ptr(ovrlp%matrix)
          !!    if (.not.inv_ovrlp_associated) call f_free_ptr(inv_ovrlp_half%matrix)
          !!else
              call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, &
                   orthpar%blocksize_pdgemm, orbs%norb, orbs, &
                   imode=1, check_accur=.false.,ovrlp=ovrlp%matrix, &
                   inv_ovrlp=inv_ovrlp_null, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half)!!, &
                   !!foe_nseg=foe_obj%nseg, foe_kernel_nsegline=foe_obj%nsegline, &
                   !!foe_istsegline=foe_obj%istsegline, foe_keyg=foe_obj%keyg)
          !!end if
      end if

      allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psittemp_c, 'psittemp_c', subname)
      allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psittemp_f, 'psittemp_f', subname)

      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
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

  !call deallocate_sparse_matrix(inv_ovrlp_half, subname)
  !!iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  !!deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  !!call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)
  call f_free_ptr(inv_ovrlp_half%matrix_compr)

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
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: uncompress_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp
  type(local_zone_descriptors),intent(in) :: lzd
  !type(orbitals_Data),intent(in) :: orbs
  type(orbitals_Data),intent(inout) :: orbs !temporary inout
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  integer,intent(in) :: correction_orthoconstraint
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
  real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
  type(sparse_matrix),intent(inout) :: lagmat
  real(kind=8),dimension(:),pointer :: psit_c, psit_f, hpsit_c, hpsit_f
  logical,intent(inout) :: can_use_transposed, overlap_calculated
  type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
  logical,intent(in) :: experimental_mode

  ! Local variables
  integer :: istat, iall, iorb, jorb, ii, ii_trans, irow, jcol, info, lwork, jj
  !type(sparse_matrix) :: tmp_mat
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

  !call nullify_sparse_matrix(tmp_mat)
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
      linmat%ovrlp%matrix=f_malloc_ptr((/orbs%norb,orbs%norb/),id='linmat%ovrlp%matrix')
      call uncompress_matrix(iproc,linmat%ovrlp)
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
      call f_free_ptr(linmat%ovrlp%matrix)
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
  !!call uncompress_matrix(iproc,linmat%ovrlp)
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

  !call deallocate_sparse_matrix(tmp_mat, subname)

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
subroutine overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, norb, orbs, imode, check_accur, ovrlp, inv_ovrlp, error, &
     ovrlp_smat, inv_ovrlp_smat)!!, &
     !!foe_nseg, foe_kernel_nsegline, foe_istsegline, foe_keyg)
  use module_base
  use module_types
  use module_interfaces, except_this_one => overlapPowerGeneral
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: compress_matrix, uncompress_matrix, transform_sparse_matrix
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, iorder, blocksize, norb, power
  type(orbitals_data),intent(in) :: orbs
  integer,intent(in) :: imode
  logical,intent(in) :: check_accur
  real(kind=8),dimension(:,:),pointer,optional :: ovrlp
  real(kind=8),dimension(:,:),pointer,optional :: inv_ovrlp
  type(sparse_matrix), optional, intent(inout) :: ovrlp_smat, inv_ovrlp_smat
  real(kind=8),intent(out),optional :: error
  !!integer,intent(in),optional :: foe_nseg
  !!integer,dimension(:),intent(in),optional :: foe_kernel_nsegline, foe_istsegline
  !!integer,dimension(:,:),intent(in),optional :: foe_keyg
  
  ! Local variables
  integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i
  integer :: matrixindex_in_compressed, nmaxvalk
  real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
  real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp, ovrlp_local, inv_ovrlp_local
  real(kind=8) :: factor, newfactor
  logical :: ovrlp_allocated, inv_ovrlp_allocated

  ! new for sparse taylor
  integer :: nout, nseq, nmaxsegk, nmaxval
  integer,dimension(:),allocatable :: ivectorindex
  integer,dimension(:,:),pointer :: onedimindices
  integer,dimension(:,:,:),allocatable :: istindexarr
  real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse, ovrlpminone_sparse_seq, ovrlp_large_compr
  real(kind=8),dimension(:),allocatable :: invovrlp_compr_seq
  real(kind=8),dimension(:,:),allocatable :: ovrlpminoneoldp, invovrlpp, ovrlp_largep
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2


  !!write(*,*) 'iorder',iorder


  call f_routine(id='overlapPowerGeneral')
  call timing(iproc,'lovrlp^-1     ','ON')


  ! Perform a check of the arguments

  if (imode/=SPARSE .and. imode/=DENSE) stop 'wrong imode'

  if (imode==DENSE) then
      if (present(ovrlp)) then
          if (.not.associated(ovrlp)) stop 'ovrlp not associated'
      else
          stop 'ovrlp not present'
      end if
      if (present(inv_ovrlp)) then
          if (.not.associated(inv_ovrlp)) stop 'inv_ovrlp not associated'
      else
          stop 'inv_ovrlp not present'
      end if
  else if (imode==SPARSE .and. (iorder>1 .or. check_accur)) then
      if (.not.present(ovrlp_smat)) stop 'ovrlp_smat not present'
      if (.not.present(inv_ovrlp_smat)) stop 'inv_ovrlp_smat not present'
      !!if (.not.present(foe_nseg)) stop 'foe_nseg not present'
      !!if (.not.present(foe_kernel_nsegline)) stop 'foe_kernel_nsegline not present'
      !!if (.not.present(foe_istsegline)) stop 'foe_istsegline not present'
      !!if (.not.present(foe_keyg)) stop 'foe_keyg not present'
  end if
  
  if (check_accur) then
      if (.not.present(error)) stop 'error not present'
  end if

  if (power/=-2 .and. power/=1 .and. power/=2) stop 'wrong value of power'

  if (nproc/=1 .and. nproc/=bigdft_mpi%nproc) stop 'wrong value of nproc'

  ! Decide whether this routine is called in parallel or in serial.
  ! If parallel, take the default values from orbs, otherwise adjust them.
  if (nproc>1) then
      norbp=orbs%norbp
      isorb=orbs%isorb
  else
      norbp=norb
      isorb=0
  end if

  sparse_dense: if (imode==DENSE) then
      if (iorder==0) then
          call vcopy(norb*norb,ovrlp(1,1),1,inv_ovrlp(1,1),1)
          if (power==1) then
             if (blocksize<0) then
                call overlap_minus_one_exact_serial(norb,inv_ovrlp)
             else
                stop 'check if working - upper half may not be filled'
                call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
                call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp(1,1), norb)
             end if
          else if (power==2) then
              if (nproc>1) then
                 call overlap_plus_minus_one_half_exact(norb,blocksize,.true.,inv_ovrlp,orbs)
             else
                 call overlap_plus_minus_one_half_exact(norb,blocksize,.true.,inv_ovrlp)
             end if
          else if (power==-2) then
              if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.false.,inv_ovrlp,orbs)
              else
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.false.,inv_ovrlp)
              end if
          end if
      else
            ovrlpminone=f_malloc_ptr((/norb,norb/),id='ovrlpminone')
            if (nproc>1) then
               ovrlpminonep=f_malloc_ptr((/norb,norbp/),id='ovrlpminonep')
            else
               ovrlpminonep => ovrlpminone
            end if

            call matrix_minus_identity_dense(norb,isorb,norbp,ovrlp(1,isorb+1),ovrlpminonep)

            ovrlppoweroldp=f_malloc_ptr((/norb,norbp/),id='ovrlppoweroldp')
            call vcopy(norb*norbp,ovrlpminonep(1,1),1,ovrlppoweroldp(1,1),1)

            if(nproc > 1) then
               call mpi_allgatherv(ovrlpminonep, norb*norbp, mpi_double_precision, ovrlpminone, &
                    norb*orbs%norb_par(:,0), norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
               !call f_free_ptr(ovrlpminonep)
            !else
            !   nullify(ovrlpminonep)
            end if

            ovrlppowerp=f_malloc_ptr((/norb,norbp/),id='ovrlppowerp')
            if (nproc>1) then
               inv_ovrlpp=f_malloc_ptr((/norb,norbp/),id='inv_ovrlpp')
           else
               inv_ovrlpp => inv_ovrlp
            end if

            if (power==1) then
               factor=-1.0d0
            else if (power==2) then
               factor=0.5d0
            else if (power==-2) then
               factor=-0.5d0
            end if
            call first_order_taylor_dense(norb,isorb,norbp,power,ovrlp(1,isorb+1),inv_ovrlpp)
            do i=2,iorder
               call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlpminone(1,1), &
                    norb, ovrlppoweroldp(1,1), norb, 0.d0, ovrlppowerp(1,1), norb)
               factor=newfactor(power,i,factor)
               call daxpy(norb*norbp,factor,ovrlppowerp,1,inv_ovrlpp,1)
               call vcopy(norb*norbp,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1),1)
            end do

            if(nproc > 1) then
               call mpi_allgatherv(inv_ovrlpp, norb*norbp, mpi_double_precision, inv_ovrlp, &
                    norb*orbs%norb_par(:,0), norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
               call f_free_ptr(inv_ovrlpp)
            else
               nullify(inv_ovrlpp)
            end if

            !!if (iorder>1) then
            !!   call f_free_ptr(inv_ovrlpp)
            !!end if

            if(nproc > 1) then
               call f_free_ptr(ovrlpminonep)
            else
               nullify(ovrlpminonep)
            end if

            call f_free_ptr(ovrlppowerp)
            call f_free_ptr(ovrlppoweroldp)
            call f_free_ptr(ovrlpminone)
        end if

        if (check_accur) then
            call check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,error)
        end if
  else if (imode==SPARSE) then
      if (iorder==0) then
          !!if (iproc==0) call yaml_warning('The compressed matrix will not be filled! You should know what you do.')
          ovrlp_local=f_malloc_ptr((/norb,norb/),id='ovrlp_local')
          inv_ovrlp_local=f_malloc_ptr((/norb,norb/),id='inv_ovrlp_local')
          call uncompress_matrix(iproc, ovrlp_smat, outmat=ovrlp_local)
          !!write(*,*) ovrlp_smat%matrix_compr
          !!write(*,*) '==============='
          !!write(*,*) ovrlp_local
          call vcopy(norb*norb,ovrlp_local(1,1),1,inv_ovrlp_local(1,1),1)
          !!do iorb=1,orbs%norb
          !!    do jorb=1,orbs%norb
          !!        inv_ovrlp_local(jorb,iorb)=0.5d0*(ovrlp_local(jorb,iorb)+ovrlp_local(iorb,jorb))
          !!    end do
          !!end do
          if (power==1) then
             if (blocksize<0) then
                call overlap_minus_one_exact_serial(norb,inv_ovrlp_local)
             else
                stop 'check if working - upper half may not be filled'
                call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp_local(1,1), norb)
                call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', norb, inv_ovrlp_local(1,1), norb)
             end if
          else if (power==2) then
              if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.true.,inv_ovrlp_local,orbs)
              else
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.true.,inv_ovrlp_local)
              end if
          else if (power==-2) then
              if (nproc>1) then
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.false.,inv_ovrlp_local,orbs)
              else
                  call overlap_plus_minus_one_half_exact(norb,blocksize,.false.,inv_ovrlp_local)
              end if
          end if
          !!! These two lines can be deleted as soon as the tests are stabilized ##########
          !!inv_ovrlp_smat%matrix=inv_ovrlp
          call compress_matrix(iproc, inv_ovrlp_smat, inmat=inv_ovrlp_local)
          call f_free_ptr(ovrlp_local)
          call f_free_ptr(inv_ovrlp_local)
          ! #############################################################################
      else
          inv_ovrlpp=f_malloc_ptr((/norb,norbp/),id='inv_ovrlpp')
          if (iorder>1) then
              ovrlp_compr_seq = f_malloc(inv_ovrlp_smat%smmm%nseq,id='ovrlp_compr_seq')
              ovrlpminone_sparse_seq = f_malloc(inv_ovrlp_smat%smmm%nseq,id='ovrlpminone_sparse_seq')
           end if
            ovrlpminone_sparse=f_malloc(inv_ovrlp_smat%nvctr,id='ovrlpminone_sparse')
            ovrlpminonep=f_malloc_ptr((/norb,norbp/),id='ovrlpminone_sparse')
            ovrlpminoneoldp=f_malloc((/norb,norbp/),id='ovrlpminone_sparse')
            invovrlpp=f_malloc((/norb,norbp/),id='ovrlpminone_sparse')
            ovrlp_large_compr=f_malloc(inv_ovrlp_smat%nvctr,id='ovrlp_large_compr')
            call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                 ovrlp_smat%matrix_compr, ovrlp_large_compr, 'small_to_large')
            call matrix_minus_identity_sparse(norb, inv_ovrlp_smat, ovrlp_large_compr, ovrlpminone_sparse)
          if (iorder>1) then
              !call sequential_acces_matrix(norb, norbp, isorb, inv_ovrlp_smat%smmm%nseg, &
              !     inv_ovrlp_smat%smmm%nsegline, inv_ovrlp_smat%smmm%istsegline, inv_ovrlp_smat%smmm%keyg, &
              !     inv_ovrlp_smat, ovrlpminone_sparse, inv_ovrlp_smat%smmm%nseq, &
              !     inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
              !     ovrlpminone_sparse_seq)
              call sequential_acces_matrix_fast(inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%nvctr, &
                   inv_ovrlp_smat%smmm%indices_extract_sequential, ovrlpminone_sparse, ovrlpminone_sparse_seq)
              call extract_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                   inv_ovrlp_smat, ovrlpminone_sparse, ovrlpminoneoldp)
           end if
            if (power==1) then
               factor=-1.0d0
            else if (power==2) then
               factor=0.5d0
            else if (power==-2) then
               factor=-0.5d0
            end if
            if (norbp>0) then
                 call extract_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                      inv_ovrlp_smat, ovrlp_large_compr, ovrlpminonep)
                call first_order_taylor_dense(norb,isorb,norbp,power,ovrlpminonep,invovrlpp)
            end if
            do i=2,iorder
               call sparsemm(inv_ovrlp_smat%smmm%nseq, ovrlpminone_sparse_seq, ovrlpminoneoldp, ovrlpminonep, &
                    norb, norbp, inv_ovrlp_smat%smmm%ivectorindex, &
                    inv_ovrlp_smat%smmm%nout, inv_ovrlp_smat%smmm%onedimindices)
               factor=newfactor(power,i,factor)
               call daxpy(norb*norbp,factor,ovrlpminonep,1,invovrlpp,1)
               call vcopy(norb*norbp,ovrlpminonep(1,1),1,ovrlpminoneoldp(1,1),1)
            end do
            !!call to_zero(inv_ovrlp_smat%nvctr, inv_ovrlp_smat%matrix_compr(1))
            call compress_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                 inv_ovrlp_smat, invovrlpp, inv_ovrlp_smat%matrix_compr)
            !!call mpiallred(inv_ovrlp_smat%matrix_compr(1), inv_ovrlp_smat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)

            if (iorder>1) then
                call f_free(ovrlpminone_sparse_seq)
            end if
            call f_free(ovrlpminone_sparse)
            call f_free_ptr(ovrlpminonep)
            call f_free(ovrlpminoneoldp)
            if (iorder>1) then
                if (.not.check_accur) call f_free(ovrlp_compr_seq)
                !!if (.not.check_accur) call f_free(istindexarr)
                !!if (.not.check_accur) call f_free(ivectorindex)
                !!if (.not.check_accur) call f_free_ptr(onedimindices)
            end if
            if (.not.check_accur) call f_free(invovrlpp)
            if (.not.check_accur) call f_free(ovrlp_large_compr)
            call f_free_ptr(inv_ovrlpp)
      end if
      if (check_accur) then
            ! HERE STARTS LINEAR CHECK ##########################
            if (iorder<=1) then
                ovrlp_compr_seq = f_malloc(inv_ovrlp_smat%smmm%nseq,id='ovrlp_compr_seq')
             end if
             if (iorder<1) then
                invovrlpp=f_malloc((/norb,norbp/),id='ovrlpminone_sparse')
                ovrlp_large_compr=f_malloc(inv_ovrlp_smat%nvctr,id='ovrlp_large_compr')
                call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                     ovrlp_smat%matrix_compr, ovrlp_large_compr, 'small_to_large')
             end if
            invovrlp_compr_seq=f_malloc(inv_ovrlp_smat%smmm%nseq,id='ovrlp_large_compr_seq')
            ovrlp_largep=f_malloc((/norb,norbp/),id='ovrlp_largep')
            call extract_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                 inv_ovrlp_smat, ovrlp_large_compr, ovrlp_largep)
            !call sequential_acces_matrix(norb, norbp, isorb, inv_ovrlp_smat%smmm%nseg, &
            !     inv_ovrlp_smat%smmm%nsegline, inv_ovrlp_smat%smmm%istsegline, inv_ovrlp_smat%smmm%keyg, &
            !     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, inv_ovrlp_smat%smmm%nseq, &
            !     inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
            !     invovrlp_compr_seq)
            call sequential_acces_matrix_fast(inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%nvctr, &
                 inv_ovrlp_smat%smmm%indices_extract_sequential, inv_ovrlp_smat%matrix_compr, invovrlp_compr_seq)
            if (power==1) then
                 !!call sequential_acces_matrix(norb, norbp, isorb, foe_nseg, &
                 !!     foe_kernel_nsegline, foe_istsegline, foe_keyg, &
                 !!     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
                 !!     invovrlp_compr_seq)
                 call check_accur_overlap_minus_one_sparse(iproc, nproc, norb, norbp, isorb, &
                      inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                      inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                      invovrlp_compr_seq, ovrlp_largep, power, &
                      error)
             else if (power==2) then
                !!call sequential_acces_matrix(norb, norbp, isorb, foe_nseg, &
                !!     foe_kernel_nsegline, foe_istsegline, foe_keyg, &
                !!     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
                !!     invovrlp_compr_seq)
                call extract_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, invovrlpp)
                call check_accur_overlap_minus_one_sparse(iproc, nproc, norb, norbp, isorb, &
                     inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                     inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                     invovrlp_compr_seq, invovrlpp, power, &
                     error,cmatp=ovrlp_largep)
             else if (power==-2) then
                !!call sequential_acces_matrix(norb, norbp, isorb, foe_nseg, &
                !!     foe_kernel_nsegline, foe_istsegline, foe_keyg, &
                !!     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
                !!     invovrlp_compr_seq)
                !call sequential_acces_matrix(norb, norbp, isorb, inv_ovrlp_smat%smmm%nseg, &
                !     inv_ovrlp_smat%smmm%nsegline, inv_ovrlp_smat%smmm%istsegline, inv_ovrlp_smat%smmm%keyg, &
                !     inv_ovrlp_smat, ovrlp_large_compr, inv_ovrlp_smat%smmm%nseq, &
                !     inv_ovrlp_smat%smmm%nmaxsegk, inv_ovrlp_smat%smmm%nmaxvalk, &
                !     ovrlp_compr_seq)
                call sequential_acces_matrix_fast(inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%nvctr, &
                     inv_ovrlp_smat%smmm%indices_extract_sequential, ovrlp_large_compr, ovrlp_compr_seq)
                call extract_matrix_distributed(iproc, nproc, norb, norbp, orbs%isorb_par, &
                     inv_ovrlp_smat, inv_ovrlp_smat%matrix_compr, invovrlpp)
                call check_accur_overlap_minus_one_sparse(iproc, nproc, norb, norbp, isorb, &
                     inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                     inv_ovrlp_smat%smmm%ivectorindex, inv_ovrlp_smat%smmm%onedimindices, &
                     invovrlp_compr_seq, invovrlpp, power, &
                     error, &
                     ovrlp_compr_seq)
             else
                 stop 'wrong power'
             end if
             call f_free(invovrlp_compr_seq)
             call f_free(ovrlp_largep)
             call f_free(ovrlp_compr_seq)
             call f_free(invovrlpp)
             call f_free(ovrlp_large_compr)
            ! HERE END LINEAR CHECK #############################
        end if
  end if sparse_dense


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
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none
  integer, intent(in) :: power
  type(sparse_matrix),intent(in) :: ovrlp
  type(sparse_matrix),intent(out) :: inv_ovrlp

  integer :: ii,iii,iorb,jorb,ii_inv,ierr,iii_inv

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
  logical, parameter :: check_lapack=.true.


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
           !!do iorb=1,norb
           !!   do jorb=1,norb
           !!      write(2000+bigdft_mpi%iproc,'(a,3i8,es16.7)') 'iproc, iorb, jorb, val', bigdft_mpi%iproc, iorb, jorb, inv_ovrlp_half(jorb,iorb)
           !!   end do
           !!end do
        call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)
        if (check_lapack) then
           tempArr=f_malloc((/norb,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0.and.abs(error)>1.0d-8) then
              print*,'LAPACK error for dsyev in overlap_plus_minus_one_half_exact',error
              open(99,file='dsyev_input.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                 end do
              end do
              close(99)
              open(99,file='dsyev_output.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                 end do
              end do
              close(99)
              call mpi_finalize(bigdft_mpi%mpi_comm)
              stop
           end if
           call f_free_ptr(inv_ovrlp_halfp)
           call f_free(orig_ovrlp)
        end if
     else
        if (check_lapack) then
           orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
           call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
        end if
        work=f_malloc(1000,id='work')
        eval1=f_malloc(norb,id='eval1')
        vl=f_malloc((/norb,norb/),id='vl')
        vr=f_malloc((/norb,norb/),id='vr')
        call dgeev( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
             norb, WORK, -1, info )
        lwork = nint(work(1))
        call f_free(work)
        work=f_malloc(lwork,id='work')
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
           tempArr=f_malloc((/norb,norb/), id='tempArr')
           do iorb=1,norb
              do jorb=1,norb
                 tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
              end do
           end do
           inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
           call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
           call f_free(tempArr)
           call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, error)
           if (bigdft_mpi%iproc==0.and.abs(error)>1.0d-8) then
              print*,'LAPACK error for dgeev in overlap_plus_minus_one_half_exact',error
              open(99,file='dgeev_input.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                 end do
              end do
              close(99)
              open(99,file='dgeev_output.txt')
              do iorb=1,norb
                 do jorb=1,norb
                   write(99,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                 end do
              end do
              close(99)
              call mpi_finalize(bigdft_mpi%mpi_comm)
              stop
           end if
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
  !write(*,*) 'iproc, present(orbs), norb, orbs%norb', bigdft_mpi%iproc, present(orbs), norb, orbs%norb
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



subroutine check_accur_overlap_minus_one_sparse(iproc, nproc, norb, norbp, isorb, nseq, nout, &
           ivectorindex, onedimindices, amat_seq, bmatp, power, &
           error, &
           dmat_seq, cmatp)
  use module_base
  implicit none
  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseq, nout, power
  integer,dimension(nseq),intent(in) :: ivectorindex
  integer,dimension(4,nout) :: onedimindices
  real(kind=8),dimension(nseq),intent(in) :: amat_seq
  real(kind=8),dimension(norb,norbp),intent(in) :: bmatp
  real(kind=8),intent(out) :: error
  real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
  real(kind=8),dimension(norb,norbp),intent(in),optional :: cmatp

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  tmpp=f_malloc0((/norb,norbp/),id='tmpp')
  if (power==1) then
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(nseq, amat_seq, bmatp, tmpp, &
          norb, norbp, ivectorindex, nout, onedimindices)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, error)
  else if (power==2) then
      if (.not.present(cmatp)) stop 'cmatp not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(nseq, amat_seq, bmatp, tmpp, &
          norb, norbp, ivectorindex, nout, onedimindices)
     call max_matrix_diff_parallel(iproc, norb, norbp, tmpp, cmatp, error)
     error=0.5d0*error
  else if (power==-2) then
     if (.not.present(dmat_seq)) stop 'dmat_seq not present'
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
     !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call sparsemm(nseq, amat_seq, bmatp, tmpp, &
          norb, norbp, ivectorindex, nout, onedimindices)
     tmp2p=f_malloc0((/norb,norbp/),id='tmp2p')
     !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
     !!     norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
     call sparsemm(nseq, dmat_seq, tmpp, tmp2p, &
          norb, norbp, ivectorindex, nout, onedimindices)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, error)
     error=0.5d0*error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one_sparse'
  end if
  call f_free(tmpp)

end subroutine check_accur_overlap_minus_one_sparse




subroutine check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,error)
  use module_base
  implicit none
  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, power
  real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
  real(kind=8),intent(out) :: error

  real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
  real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
  integer :: ierr, i,j

  tmpp=f_malloc((/norb,norbp/),id='tmpp')
  if (power==1) then
     call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
          norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, error)
  else if (power==2) then
     call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
          norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     call max_matrix_diff_parallel(iproc, norb, norbp, tmpp, ovrlp(1,isorb+1), error)
     error=0.5d0*error
  else if (power==-2) then
     call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
          norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
     tmp2p=f_malloc((/norb,norbp/),id='tmp2p')
     call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
          norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
     call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, error)
     error=0.5d0*error
     call f_free(tmp2p)
  else
     stop 'Error in check_accur_overlap_minus_one'
  end if
  call f_free(tmpp)

end subroutine check_accur_overlap_minus_one


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


subroutine max_matrix_diff_parallel(iproc, norb, norbp, mat1, mat2, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, norb, norbp
  real(8),dimension(norb,norbp),intent(in):: mat1, mat2
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, jorb, ierr
  real(8):: error

  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norbp
     !$omp parallel default(private) shared(norb, mat1, mat2, deviation)
     !$omp do reduction(max:deviation)
     do jorb=1,norb
        error=abs(mat1(jorb,iorb)-mat2(jorb,iorb))
        deviation=max(error,deviation)
     end do
     !$omp end do
     !$omp end parallel
  end do
  call mpiallred(deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
  call timing(iproc,'dev_from_unity','OF') 

end subroutine max_matrix_diff_parallel



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


subroutine deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, ovrlp, deviation)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb, norbp, isorb
  real(8),dimension(norb,norbp),intent(in):: ovrlp
  real(8),intent(out):: deviation

  ! Local variables
  integer:: iorb, iiorb, jorb, ierr
  real(8):: error


  call timing(iproc,'dev_from_unity','ON') 
  deviation=0.d0
  do iorb=1,norbp
     iiorb=iorb+isorb
     !$omp parallel default(private) shared(norb, iiorb, ovrlp, iorb, deviation)
     !$omp do reduction(max:deviation)
     do jorb=1,norb
        if(iiorb==jorb) then
           error=abs(ovrlp(jorb,iorb)-1.d0)
        else
           error=abs(ovrlp(jorb,iorb))
        end if
        deviation=max(error,deviation)
     end do
     !$omp end do
     !$omp end parallel
  end do
  if (nproc>1) then
      call mpiallred(deviation, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
  call timing(iproc,'dev_from_unity','OF') 

end subroutine deviation_from_unity_parallel


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
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, meth_overlap
  type(orbitals_data),intent(in) :: orbs
  type(sparse_matrix),intent(in) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, ii, iseg, ind
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
     !!call overlapPowerGeneral(iproc, nproc, meth_overlap, -2, -8, n, orbs, imode=2, check_accur=.true.,&
     !!     ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)
     call overlapPowerGeneral(iproc, 1, meth_overlap, -2, -8, n, orbs, imode=2, check_accur=.true.,&
          ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)


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

  if (nproc>1)then
      call mpiallred(inv_ovrlp_half%matrix_compr(1), inv_ovrlp_half%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
  end if

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

!> Orthonormalize a subset of orbitals
subroutine orthonormalize_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthonormalize_subset
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='orthonormalize_subset'
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_null
  real(kind=8) :: error

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')

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
          call overlap_power_minus_one_half_parallel(iproc, nproc, methTransformOverlap, &
               orbs, ovrlp, inv_ovrlp_half)
      else
          nullify(inv_ovrlp_null)
          ! do sparse.. check later
          call overlapPowerGeneral(iproc, nproc, methTransformOverlap, -2, &
               orthpar%blocksize_pdsyev, orbs%norb, orbs,&
               imode=1, check_accur=.true., ovrlp=ovrlp%matrix, inv_ovrlp=inv_ovrlp_null, &
               error=error, ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half)
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

      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)

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

  !call deallocate_sparse_matrix(inv_ovrlp_half, subname)
  !!iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  !!deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  !!call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)
  call f_free_ptr(inv_ovrlp_half%matrix_compr)

end subroutine orthonormalize_subset



subroutine gramschmidt_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
           orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
           lphi, psit_c, psit_f, can_use_transposed)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => gramschmidt_subset
  use communications, only: transpose_localized, untranspose_localized
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
  type(local_zone_descriptors),intent(in) :: lzd
  type(sparse_matrix),intent(inout) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
  type(comms_linear),intent(in) :: collcom
  type(orthon_data),intent(in) :: orthpar
  real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  logical,intent(inout) :: can_use_transposed

  ! Local variables
  integer :: it, istat, iall, iorb, jorb, iat, jat, ii
  logical :: iout, jout
  integer,dimension(:),allocatable :: icount_norb, jcount_norb
  real(kind=8),dimension(:),allocatable :: psittemp_c, psittemp_f, norm
  !type(sparse_matrix) :: inv_ovrlp_half
  character(len=*),parameter :: subname='gramschmidt_subset'

  if(orthpar%nItOrtho>1) write(*,*) 'WARNING: might create memory problems...'

  !call nullify_sparse_matrix(inv_ovrlp_half)
  !call sparse_copy_pattern(inv_ovrlp, inv_ovrlp_half, iproc, subname)
  !!allocate(inv_ovrlp_half%matrix_compr(inv_ovrlp_half%nvctr), stat=istat)
  !!call memocc(istat, inv_ovrlp_half%matrix_compr, 'inv_ovrlp_half%matrix_compr', subname)
  inv_ovrlp_half%matrix_compr=f_malloc_ptr(inv_ovrlp_half%nvctr,id='inv_ovrlp_half%matrix_compr')

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

      call vcopy(sum(collcom%nrecvcounts_c), psit_c(1), 1, psittemp_c(1), 1)
      call vcopy(7*sum(collcom%nrecvcounts_f), psit_f(1), 1, psittemp_f(1), 1)
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

  !call deallocate_sparse_matrix(inv_ovrlp_half, subname)
  !!iall=-product(shape(inv_ovrlp_half%matrix_compr))*kind(inv_ovrlp_half%matrix_compr)
  !!deallocate(inv_ovrlp_half%matrix_compr, stat=istat)
  !!call memocc(istat, iall, 'inv_ovrlp_half%matrix_compr', subname)
  call f_free_ptr(inv_ovrlp_half%matrix_compr)

end subroutine gramschmidt_subset





subroutine diagonalize_localized(iproc, nproc, orbs, ovrlp, inv_ovrlp_half)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(sparse_matrix),intent(in) :: ovrlp
  type(sparse_matrix),intent(inout) :: inv_ovrlp_half

  ! Local variables
  integer :: iend, i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
  integer :: iiorb, ierr, ii, iseg, ind, lwork
  real(kind=8) :: error
  real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
  real(kind=8),dimension(:),allocatable :: eval, work
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
          
     allocate(ovrlp_tmp_inv_half(n,n),stat=istat)
     call memocc(istat, ovrlp_tmp_inv_half, 'ovrlp_tmp_inv_half', subname)

     !if (iiorb==orbs%norb) then
     !print*,''
     !print*,'ovrlp_tmp',n,iiorb
     !do jorb=1,n
     !print*,jorb,ovrlp_tmp(:,jorb)
     !end do
     !end if

     call vcopy(n*n, ovrlp_tmp(1,1), 1, ovrlp_tmp_inv_half(1,1), 1)
     ! Diagonalize the small matrix
     lwork=n**2
     allocate(work(lwork),stat=istat)
     call memocc(istat,work,'work',subname)
     allocate(eval(n),stat=istat)
     call memocc(istat,eval,'eval',subname)
     call dsyev('v', 'l', n, ovrlp_tmp_inv_half, n, eval, work, lwork, istat)
     if (istat/=0) stop 'error in dsyev'
     iall=-product(shape(work))*kind(work)
     deallocate(work, stat=istat)
     call memocc(istat, iall, 'work', subname)
     iall=-product(shape(eval))*kind(eval)
     deallocate(eval, stat=istat)
     call memocc(istat, iall, 'eval', subname)

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

  call mpiallred(inv_ovrlp_half%matrix_compr(1), inv_ovrlp_half%nvctr, mpi_sum, bigdft_mpi%mpi_comm)

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

end subroutine diagonalize_localized









  subroutine extract_matrix_distributed(iproc, nproc, norb, norbp, isorb_par, smat, matrix_compr, matrixp)
    use module_base
    use module_types
    use sparsematrix_base, only: sparse_matrix
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, nproc, norb, norbp
    integer,dimension(0:nproc-1),intent(in) :: isorb_par
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctr),intent(in) :: matrix_compr
    real(kind=8),dimension(norb,norbp),intent(out) :: matrixp

    ! Local variables
    integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


       if (norbp>0) then

           call to_zero(norb*norbp,matrixp(1,1))

           isegstart=smat%istsegline(isorb_par(iproc)+1)
           if (isorb_par(iproc)+norbp<norb) then
               isegend=smat%istsegline(isorb_par(iproc+1)+1)-1
           else
               isegend=smat%nseg
           end if
           !$omp parallel do default(private) &
           !$omp shared(isegstart, isegend, smat, norb, isorb_par, matrixp, matrix_compr, iproc)
           do iseg=isegstart,isegend
               ii=smat%keyv(iseg)-1
               do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                   ii=ii+1
                   iiorb = (jorb-1)/norb + 1
                   jjorb = jorb - (iiorb-1)*norb
                   matrixp(jjorb,iiorb-isorb_par(iproc)) = matrix_compr(ii)
               end do
           end do
           !$omp end parallel do
       end if

   end subroutine extract_matrix_distributed




   subroutine matrix_minus_identity_sparse(norb, smat, ovrlp_compr, ovrlpminone_compr)
     use module_base
     use module_types
     use sparsematrix_base, only: sparse_matrix
     implicit none

     ! Calling arguments
     integer,intent(in) :: norb
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctr),intent(in) :: ovrlp_compr
     real(kind=8),dimension(smat%nvctr),intent(out) :: ovrlpminone_compr

     ! Local variables
     integer :: ii, iseg, jorb, iiorb, jjorb

     !$omp parallel do default(private) &
     !$omp shared(smat, norb, ovrlpminone_compr, ovrlp_compr)
     do iseg=1,smat%nseg
         ii=smat%keyv(iseg)-1
         do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
             iiorb = (jorb-1)/norb + 1
             jjorb = jorb - (iiorb-1)*norb
             ii=ii+1
             if (iiorb==jjorb) then
                 ovrlpminone_compr(ii)=ovrlp_compr(ii)-1.d0
             else
                 ovrlpminone_compr(ii)=ovrlp_compr(ii)
             end if
         end do
     end do
     !$omp end parallel do

   end subroutine matrix_minus_identity_sparse



   subroutine compress_matrix_distributed(iproc, nproc, norb, norbp, isorb_par, &
              smat, matrixp, matrix_compr)
     use module_base
     use module_types
     use sparsematrix_base, only: sparse_matrix
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, norb, norbp
     integer,dimension(0:nproc-1),intent(in) :: isorb_par
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(norb,norbp),intent(in) :: matrixp
     real(kind=8),dimension(smat%nvctr),intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, ierr


     call to_zero(smat%nvctr, matrix_compr(1))

     if (norbp>0) then
         isegstart=smat%istsegline(isorb_par(iproc)+1)
         if (isorb_par(iproc)+norbp<norb) then
             isegend=smat%istsegline(isorb_par(iproc+1)+1)-1
         else
             isegend=smat%nseg
         end if
         !$omp parallel default(none) &
         !$omp shared(isegstart, isegend, matrixp, smat, norb, matrix_compr, isorb_par, iproc) &
         !$omp private(iseg, ii, jorb, iiorb, jjorb)
         !$omp do
         do iseg=isegstart,isegend
             ii=smat%keyv(iseg)-1
             do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                 ii=ii+1
                 iiorb = (jorb-1)/norb + 1
                 jjorb = jorb - (iiorb-1)*norb
                 matrix_compr(ii)=matrixp(jjorb,iiorb-isorb_par(iproc))
             end do
         end do
         !$omp end do
         !$omp end parallel
     end if

     call mpiallred(matrix_compr(1), smat%nvctr, mpi_sum, bigdft_mpi%mpi_comm)

  end subroutine compress_matrix_distributed
