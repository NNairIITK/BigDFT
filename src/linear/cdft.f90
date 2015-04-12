!> @file
!>  Contains files for constrained DFT
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> CDFT: calculates the weight matrix w_ab via the expression S^1/2PS^1/2, where S is the overlap of the whole system
!! CDFT: and P is a projector matrix onto the tmbs of the desired fragment
!! CDFT: for standalone CDFT calculations, assuming one charged fragment, for transfer integrals assuming two fragments
!! CDFT: where we constrain the difference - should later generalize this
subroutine calculate_weight_matrix_lowdin_wrapper(cdft,tmb,input,ref_frags,calculate_overlap_matrix,meth_overlap)
  use module_base
  use module_types
  use module_fragments
  use constrained_dft
  use module_interfaces
  implicit none
  integer, intent(in) :: meth_overlap
  type(cdft_data), intent(inout) :: cdft
  type(input_variables),intent(in) :: input
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

  integer :: nfrag_charged
  !real(kind=gp), dimension(:,:), pointer :: ovrlp_half
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'

  call f_routine('calculate_weight_matrix_lowdin_wrapper')

  ! wrapper here so we can modify charge and avoid restructuring the code as much whilst still using the routine elsewhere
  if (.not. input%lin%calc_transfer_integrals) then
     cdft%charge=ref_frags(input%frag%frag_index(cdft%ifrag_charged(1)))%nelec-input%frag%charge(cdft%ifrag_charged(1))
  end if

  if (input%lin%calc_transfer_integrals) then
     nfrag_charged=2
  else
     nfrag_charged=1
  end if

  !ovrlp_half=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
  call calculate_weight_matrix_lowdin(cdft%weight_matrix,cdft%weight_matrix_,nfrag_charged,cdft%ifrag_charged,tmb,input%frag,&
       ref_frags,calculate_overlap_matrix,.true.,meth_overlap)
  !call f_free_ptr(ovrlp_half)

  call f_release_routine()

end subroutine calculate_weight_matrix_lowdin_wrapper


subroutine calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
     ref_frags,calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap)
  use module_base
  use module_types
  use module_fragments
  use module_interfaces, except_this_one => calculate_weight_matrix_lowdin
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: matrices, sparse_matrix, sparsematrix_malloc_ptr, &
                               DENSE_FULL, assignment(=), &
                               allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix, &
                          gather_matrix_from_taskgroups_inplace, uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed
  use matrix_operations, only: overlapPowerGeneral
  implicit none
  type(sparse_matrix), intent(inout) :: weight_matrix
  type(matrices), intent(inout) :: weight_matrix_
  type(fragmentInputParameters),intent(in) :: input_frag
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  integer, intent(in) :: nfrag_charged, meth_overlap
  integer, dimension(2), intent(in) :: ifrag_charged
  !local variables
  integer :: ifrag,iorb,ifrag_ref,isforb,ierr
  real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
  real(kind=gp) :: max_error, mean_error
  type(matrices),dimension(1) :: inv_ovrlp

  call f_routine(id='calculate_weight_matrix_lowdin')

  call allocate_matrices(tmb%linmat%s, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

  if (calculate_overlap_matrix) then
     if(.not.tmb%can_use_transposed) then
         !if(.not.associated(tmb%psit_c)) then
         !    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         !end if
         !if(.not.associated(tmb%psit_f)) then
         !    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         !end if
         call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
          tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
  end if   

  if (calculate_ovrlp_half) then
     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, &
          tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
     ! Maybe not clean here to use twice tmb%linmat%s, but it should not
     ! matter as dense is used
     call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/2/), &
          tmb%orthpar%blocksize_pdsyev, &
          imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%s, &
          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
          check_accur=.false., max_error=max_error, mean_error=mean_error)
          !check_accur=.true., max_error=max_error, mean_error=mean_error)
     call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  end if

  proj_mat=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='proj_mat')
  !call f_zero(tmb%orbs%norb**2,proj_mat(1,1))

  isforb=0
  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)
     if (ifrag==ifrag_charged(1)) then
        do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
        end do
     end if
     if (nfrag_charged==2) then
        if (ifrag==ifrag_charged(2)) then
           do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
              proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
           end do
        end if
     end if
     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do

  proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
            tmb%orbs%norb, 1.d0, &
            proj_mat(1,1), tmb%orbs%norb, &
            inv_ovrlp(1)%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, 0.d0, &
            proj_ovrlp_half(1,1), tmb%orbs%norb)
  end if
  call f_free(proj_mat)
  weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
          tmb%orbs%norb, 1.d0, &
          inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, &
          proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
          weight_matrixp(1,1), tmb%orbs%norb)
  end if
  call f_free(proj_ovrlp_half)
  weight_matrix_%matrix = sparsematrix_malloc_ptr(weight_matrix,iaction=DENSE_FULL,id='weight_matrix_%matrix')
  if (bigdft_mpi%nproc>1) then
     call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix_%matrix, &
          tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
          mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      if (weight_matrix%nspin/=1) then
          stop 'NEED TO FIX THE SPIN HERE: calculate_weight_matrix_lowdin'
      end if
     call vcopy(tmb%orbs%norb*tmb%orbs%norb*weight_matrix%nspin,weight_matrixp(1,1),1,weight_matrix_%matrix(1,1,1),1)
  end if
  call f_free(weight_matrixp)
  call compress_matrix(bigdft_mpi%iproc,weight_matrix,weight_matrix_%matrix,weight_matrix_%matrix_compr)
  call f_free_ptr(weight_matrix_%matrix)
  call deallocate_matrices(inv_ovrlp(1))
  call f_release_routine()

end subroutine calculate_weight_matrix_lowdin


!EXPERIMENTAL and temporary, currently deactivated
subroutine calculate_weight_matrix_lowdin_gradient_fd(weight_matrix,weight_matrix_,ifrag_charged,tmb,input_frag,&
     ref_frags,calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,cdft_grad)
  use module_base
  use module_types
  use module_fragments
  use module_interfaces
  use sparsematrix_base, only: matrices, sparse_matrix, sparsematrix_malloc_ptr, &
                               DENSE_FULL, assignment(=), &
                               allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix, &
                          gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace, &
                          uncompress_matrix2
   use matrix_operations, only: overlapPowerGeneral
  implicit none
  type(sparse_matrix), intent(inout) :: weight_matrix
  type(matrices), intent(inout) :: weight_matrix_
  type(fragmentInputParameters),intent(in) :: input_frag
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  integer, intent(in) :: meth_overlap
  integer, dimension(2), intent(in) :: ifrag_charged
  real(kind=8),dimension(tmb%npsidim_orbs),intent(out) :: cdft_grad
  !local variables
  integer :: ifrag,iorb,ifrag_ref,isforb,ierr
  real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp!, ovrlp_half
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
  real(kind=gp) :: max_error, mean_error
  type(matrices),dimension(1) :: inv_ovrlp
  type(matrices),dimension(1) :: ovrlp_half
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f, matrix_compr, psi_orig
  !real(kind=8),dimension(:),pointer :: psitlarge_c, psitlarge_f
  real(kind=8),dimension(:),pointer :: weight_matrix_tmp
  real(kind=8), dimension(:), allocatable :: calc_grad, fd_grad
  integer :: nfrag_charged, jorb, i, ncount, istart, iiorb, ilr
  real(kind=8) :: ddot, trkw, trkw_new, trkw_old
  real(kind=8), parameter :: h=0.00001d0
  logical, parameter :: forward=.true. !forward or centered diff
  call f_routine(id='calculate_weight_matrix_lowdin')

  !might be better to make this a parameter
  if (ifrag_charged(2)==0) then
     nfrag_charged=1
  else
     nfrag_charged=2
  end if

      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
print*,iorb,ilr,ncount
 end do

!stop

  psi_orig=f_malloc_ptr(tmb%npsidim_orbs,id='psi_orig')
  call vcopy(tmb%npsidim_orbs,tmb%psi(1),1,psi_orig(1),1)
  !call f_zero(tmb%npsidim_orbs,cdft_grad)
  call f_zero(cdft_grad)

  !calculate Tr(KW)
  !ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
  call calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
       ref_frags,calculate_overlap_matrix,.true.,meth_overlap)
 
  !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
  call extract_taskgroup_inplace(weight_matrix, weight_matrix_)
  call calculate_kernel_and_energy(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%linmat%l,weight_matrix, &
           tmb%linmat%kernel_,weight_matrix_,trkw,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
  !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%l, tmb%linmat%kernel_)
  call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, weight_matrix, weight_matrix_)

  weight_matrix_tmp = f_malloc_ptr(weight_matrix%nvctr,id='weight_matrix_tmp')
  call vcopy(weight_matrix%nvctr,weight_matrix_%matrix_compr(1),1,weight_matrix_tmp(1),1)
  do i=1,tmb%npsidim_orbs !19000,21500!
     call vcopy(tmb%npsidim_orbs,psi_orig(1),1,tmb%psi(1),1)
     if (forward) then
       tmb%psi(i)=tmb%psi(i)+h
     else
       tmb%psi(i)=tmb%psi(i)+0.5d0*h
     end if
     tmb%can_use_transposed=.false.
     call calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
          ref_frags,.true.,.true.,meth_overlap)
     !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
     call extract_taskgroup_inplace(weight_matrix, weight_matrix_)
     call calculate_kernel_and_energy(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%linmat%l,weight_matrix, &
           tmb%linmat%kernel_,weight_matrix_,trkw_new,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%l, tmb%linmat%kernel_)
     call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, weight_matrix, weight_matrix_)

     if (forward) then
        cdft_grad(i)=(trkw_new-trkw)/h
     else
        call vcopy(tmb%npsidim_orbs,psi_orig(1),1,tmb%psi(1),1)
        tmb%psi(i)=tmb%psi(i)-0.5d0*h
        tmb%can_use_transposed=.false.
        call calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
             ref_frags,.true.,.true.,meth_overlap)
        !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
        call extract_taskgroup_inplace(weight_matrix, weight_matrix_)
        call calculate_kernel_and_energy(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%linmat%l,weight_matrix, &
              tmb%linmat%kernel_,weight_matrix_,trkw_old,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
        !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%l, tmb%linmat%kernel_)
        call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, weight_matrix, weight_matrix_)
        cdft_grad(i)=(trkw_new-trkw_old)/h
     end if


     !print*,100.0d0*i/tmb%npsidim_orbs,'% (',i,tmb%npsidim_orbs,')',cdft_grad(i),tmb%psi(i),trkw_new
     print*,100.0d0*i/tmb%npsidim_orbs,'% (',cdft_grad(i),tmb%psi(i),trkw_new,')'



  end do

if (.false.) then
     ! for ease, just calc d(S^1/2)ab/dphi_c, where c=9, b=10
  call allocate_matrices(tmb%linmat%s, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))
  call allocate_matrices(tmb%linmat%s, allocate_full=.true., matname='ovrlp_half', mat=ovrlp_half(1))
     !quick test
     print*,tmb%orbs%inwhichlocreg(1),tmb%orbs%inwhichlocreg(2),&
          tmb%lzd%llr(tmb%orbs%inwhichlocreg(1))%wfd%nvctr_c+7*tmb%lzd%llr(tmb%orbs%inwhichlocreg(1))%wfd%nvctr_f,&
          tmb%lzd%llr(tmb%orbs%inwhichlocreg(2))%wfd%nvctr_c+7*tmb%lzd%llr(tmb%orbs%inwhichlocreg(2))%wfd%nvctr_f

     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, &
          tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)


     !S^-1/2 for calc grad
             call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/-2/), &
                  tmb%orthpar%blocksize_pdsyev, &
                  imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%s, &
                  ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
                  check_accur=.false., max_error=max_error, mean_error=mean_error)       

             call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/2/), &
                  tmb%orthpar%blocksize_pdsyev, &
                  imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%s, &
                  ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=ovrlp_half, &
                  check_accur=.false., max_error=max_error, mean_error=mean_error)       


istart=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          calc_grad=f_malloc(ncount,id='calc_grad')
          call vcopy(ncount,tmb%psi(istart+ncount),1,calc_grad(1),1)
          call dscal(inv_ovrlp(1)%matrix(9,10,1)*0.5d0,calc_grad(1))

call daxpy(ncount,inv_ovrlp(1)%matrix(10,9,1)*0.5d0,tmb%psi(istart+ncount),1,calc_grad(1),1)
          if (iorb==9) exit
                   istart=istart+ncount 
          call f_free(calc_grad)
     end do

     !test instead deriv of S^1/2
      !istart=1
      !do iorb=1,1 !tmb%orbs%norbp
      !    iiorb=tmb%orbs%isorb+iorb
      !    ilr=tmb%orbs%inwhichlocreg(iiorb)
      !    ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          fd_grad=f_malloc(ncount,id='fd_grad')

          do i=0,ncount-1
             call vcopy(ncount,psi_orig(istart),1,tmb%psi(istart),1)
             tmb%psi(istart+i)=tmb%psi(istart+i)+h

             ! re-calc overlap first but only for a given row
             ! ddot is sloppy but they should both be same size (same ilr), print to check
             tmb%linmat%ovrlp_%matrix(9,10,1)=ddot(ncount,tmb%psi(istart),1,tmb%psi(istart+ncount),1)
             tmb%linmat%ovrlp_%matrix(10,9,1)=ddot(ncount,tmb%psi(istart),1,tmb%psi(istart+ncount),1)

             call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/2/), &
                  tmb%orthpar%blocksize_pdsyev, &
                  imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%s, &
                  ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
                  check_accur=.false., max_error=max_error, mean_error=mean_error)             

             fd_grad(i+1) = (inv_ovrlp(1)%matrix(9,10,1)-ovrlp_half(1)%matrix(9,10,1)+ &
                             inv_ovrlp(1)%matrix(10,9,1)-ovrlp_half(1)%matrix(10,9,1))/h

             write(30,*) fd_grad(i+1),calc_grad(i+1),tmb%psi(istart+i)!abs(fd_grad(i+1)-calc_grad(i+1))!,inv_ovrlp%matrix(9,10,1),ovrlp_half%matrix(9,10,1)

          end do
       !   istart=istart+ncount
          call f_free(fd_grad)

      !end do
          call f_free(calc_grad)
     call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  call deallocate_matrices(inv_ovrlp(1))
  call deallocate_matrices(ovrlp_half(1))
end if

  call vcopy(tmb%npsidim_orbs,psi_orig(1),1,tmb%psi(1),1)
  call f_free_ptr(psi_orig)
  call vcopy(weight_matrix%nvctr,weight_matrix_tmp(1),1,weight_matrix_%matrix_compr(1),1)
  call f_free_ptr(weight_matrix_tmp)
  !call f_free(ovrlp_half)
    !cdft%charge is always constant (as is lagmult in this loop) so could in theory be ignored as in optimize_coeffs
    !cdft%lag_mult*(trkw - cdft%charge)

end subroutine calculate_weight_matrix_lowdin_gradient_fd



!EXPERIMENTAL, currently deactivated
subroutine calculate_weight_matrix_lowdin_gradient(weight_matrix,weight_matrix_,ifrag_charged,tmb,input_frag,ref_frags,&
     calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,psitlarge_c,psitlarge_f)
  use module_base
  use module_types
  use module_fragments
  use module_interfaces
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: matrices, sparse_matrix, sparsematrix_malloc_ptr, &
                               DENSE_FULL, assignment(=), &
                               allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          gather_matrix_from_taskgroups_inplace, uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
  use matrix_operations, only: overlapPowerGeneral
  implicit none
  type(sparse_matrix), intent(inout) :: weight_matrix
  type(matrices), intent(inout) :: weight_matrix_
  type(fragmentInputParameters),intent(in) :: input_frag
  type(dft_wavefunction), intent(inout) :: tmb
  logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  integer, intent(in) :: meth_overlap
  integer, dimension(2), intent(in) :: ifrag_charged
  real(kind=8),dimension(:),pointer :: psitlarge_c, psitlarge_f
  !local variables
  integer :: ifrag,iorb,ifrag_ref,isforb,ierr
  real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp, ovrlp_half
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
  real(kind=gp) :: max_error, mean_error
  type(matrices),dimension(1) :: inv_ovrlp
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f, matrix_compr
  !real(kind=8),dimension(:),pointer :: psitlarge_c, psitlarge_f
  real(kind=8),dimension(:,:),pointer :: weight_matrix_tmp
  integer :: nfrag_charged, jorb
real(kind=8) :: ddot

  call f_routine(id='calculate_weight_matrix_lowdin')

  !might be better to make this a parameter
  if (ifrag_charged(2)==0) then
     nfrag_charged=1
  else
     nfrag_charged=2
  end if


  call allocate_matrices(tmb%linmat%s, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

  if (calculate_overlap_matrix) then
     if(.not.tmb%can_use_transposed) then
         !if(.not.associated(tmb%psit_c)) then
         !    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         !end if
         !if(.not.associated(tmb%psit_f)) then
         !    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         !end if
         call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
          tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
  end if   

  !calculate s^-1/2 and generate s^1/2 from it
  if (calculate_ovrlp_half) then !always do so as not sure how this would make sense otherwise?  for that matter how does it make sense anyway?  presumably always true?
     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, &
          tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
!print*,'g',ddot(tmb%orbs%norb*tmb%orbs%norb, tmb%linmat%ovrlp_%matrix(1,1,1), 1, tmb%linmat%ovrlp_%matrix(1,1,1), 1)
     ! Maybe not clean here to use twice tmb%linmat%s, but it should not
     ! matter as dense is used
     call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/-2/), &
          tmb%orthpar%blocksize_pdsyev, &
          imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%s, &
          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
          check_accur=.false., max_error=max_error, mean_error=mean_error)
!print*,'f',ddot(tmb%orbs%norb*tmb%orbs%norb, inv_ovrlp%matrix(1,1,1), 1, inv_ovrlp%matrix(1,1,1), 1)
     ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ovrlp_half')
     !forget parallelism for now
          call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, &
            tmb%orbs%norb, 1.d0, &
            tmb%linmat%ovrlp_%matrix(1,1,1), tmb%orbs%norb, &
            inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, 0.d0, &
            ovrlp_half(1,1), tmb%orbs%norb)
!print*,'e',ddot(tmb%orbs%norb*tmb%orbs%norb, ovrlp_half(1,1), 1, ovrlp_half(1,1), 1)
     !call f_free_ptr(tmb%linmat%ovrlp_%matrix)

  end if

  proj_mat=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='proj_mat')

  !call f_zero(tmb%orbs%norb**2,proj_mat(1,1))
  isforb=0
  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)
     if (ifrag==ifrag_charged(1)) then
        do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
        end do
     end if
     if (nfrag_charged==2) then
        if (ifrag==ifrag_charged(2)) then
           do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
              proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
           end do
        end if
     end if
     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do
!print*,'h',ddot(tmb%orbs%norb*tmb%orbs%norb, proj_mat(1,1), 1, proj_mat(1,1), 1)
  !calculate 0.5*KS^1/2P
  proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
            tmb%orbs%norb, 1.d0, &
            ovrlp_half(1,1), tmb%orbs%norb, &
            proj_mat(1,tmb%orbs%isorb+1), tmb%orbs%norb, 0.d0, &
            proj_ovrlp_half(1,1), tmb%orbs%norb)
  end if
!print*,'d',ddot(tmb%orbs%norb*tmb%orbs%norbp, proj_ovrlp_half(1,1), 1, proj_ovrlp_half(1,1), 1) !?!?ARGH 
  call f_free(proj_mat)
  !call f_free(ovrlp_half)

  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%l, &
       tmb%linmat%kernel_%matrix_compr, tmb%linmat%kernel_%matrix)

  weight_matrix_tmp = f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norbp/),id='weight_matrix_tmp')
  weight_matrixp = f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='weight_matrixp')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
          tmb%orbs%norb, 1.d0, &
          tmb%linmat%kernel_%matrix(1,1,1), tmb%orbs%norb, &
          proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
          !weight_matrixp(1,1), tmb%orbs%norb)
          weight_matrix_tmp(1,1), tmb%orbs%norb)
  end if
  call f_free_ptr(tmb%linmat%kernel_%matrix)
  call f_free(proj_ovrlp_half)
!print*,'c',ddot(tmb%orbs%norb*tmb%orbs%norbp, weight_matrix_tmp(1,1), 1, weight_matrix_tmp(1,1), 1)
  ! multiply by S^-1/2

  !weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
          tmb%orbs%norb, 1.d0, &
          inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, &
          !ovrlp_half(1,1), tmb%orbs%norb, &
          weight_matrix_tmp(1,1), tmb%orbs%norb, 0.d0, &
          weight_matrixp(1,1), tmb%orbs%norb)
  end if
!print*,'b',ddot(tmb%orbs%norb*tmb%orbs%norbp, weight_matrixp(1,1), 1, weight_matrixp(1,1), 1)
  call f_free_ptr(weight_matrix_tmp)
  !call f_free(ovrlp_half)!*
  !do iorb=1,tmb%orbs%norbp
  !   jorb=tmb%orbs%isorb+iorb
  !   weight_matrixp(jorb,iorb)=0.5d0*weight_matrixp(jorb,iorb)
  !end do

  weight_matrix_%matrix = sparsematrix_malloc_ptr(weight_matrix,iaction=DENSE_FULL,id='weight_matrix_%matrix')
  if (bigdft_mpi%nproc>1) then
     call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix_%matrix, &
          tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
          mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      if (weight_matrix%nspin/=1) then
          stop 'NEED TO FIX THE SPIN HERE: calculate_weight_matrix_lowdin'
      end if
     call vcopy(tmb%orbs%norb*tmb%orbs%norb*weight_matrix%nspin,weight_matrixp(1,1),1,weight_matrix_%matrix(1,1,1),1)
  end if
  call f_free(weight_matrixp)

  !add transpose
  weight_matrix_tmp = f_malloc_ptr((/weight_matrix%nfvctr,weight_matrix%nfvctr/),id='weight_matrix_tmp')
  do iorb=1,weight_matrix%nfvctr
    do jorb=1,weight_matrix%nfvctr
       weight_matrix_tmp(jorb,iorb)=weight_matrix_%matrix(iorb,jorb,1)+weight_matrix_%matrix(jorb,iorb,1)
    end do
  end do


!pre mutliply instead?
!          call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, &
!            tmb%orbs%norb, 1.d0, &
!            ovrlp_half(1,1), tmb%orbs%norb, &
!            weight_matrix_tmp(1,1), tmb%orbs%norb, 0.d0, &
!            weight_matrix_%matrix(1,1,1), tmb%orbs%norb)

!post-multiply by S
!          call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, &
!            tmb%orbs%norb, 1.d0, &
!            weight_matrix_%matrix(1,1,1), tmb%orbs%norb, &
!            tmb%linmat%ovrlp_%matrix(1,1,1), tmb%orbs%norb, 0.d0, &
!            !ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
!            !inv_ovrlp%matrix(1,1,1), tmb%orbs%norb, 0.d0, &
!            weight_matrix_tmp(1,1), tmb%orbs%norb)


  call f_free(ovrlp_half)
call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  call vcopy(weight_matrix%nfvctr*weight_matrix%nfvctr,weight_matrix_tmp(1,1),1,weight_matrix_%matrix(1,1,1),1)
  call f_free_ptr(weight_matrix_tmp)

  !want to preserve original weight matrix not gradient version
  matrix_compr=f_malloc_ptr(weight_matrix%nvctr,id='matrix_compr')
  call vcopy(weight_matrix%nvctr,weight_matrix_%matrix_compr(1),1,matrix_compr(1),1)
  call compress_matrix(bigdft_mpi%iproc,weight_matrix,weight_matrix_%matrix,weight_matrix_%matrix_compr)
  call f_free_ptr(weight_matrix_%matrix)
  call deallocate_matrices(inv_ovrlp(1)) ! can move this


  !psitlarge_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='psitlarge_c')
  !psitlarge_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='psitlarge_f')
     !if (tmb%ham_descr%npsidim_orbs > 0)  call f_zero(tmb%ham_descr%npsidim_orbs,tmb%ham_descr%psi(1))
     !call small_to_large_locreg(bigdft_mpi%iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
     !     tmb%orbs, tmb%psi, tmb%ham_descr%psi)
!assuming we still have tmb%ham_descr%psi
!print*,size(psitlarge_c),size(psitlarge_f),tmb%ham_descr%collcom%ndimind_c,7*tmb%ham_descr%collcom%ndimind_f
          call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, &
               tmb%ham_descr%collcom, TRANSPOSE_FULL, tmb%ham_descr%psi, psitlarge_c, psitlarge_f, tmb%ham_descr%lzd)
!print*,'1',ddot(tmb%ham_descr%collcom%ndimind_c, psitlarge_c(1), 1, psitlarge_c(1), 1),ddot(7*tmb%ham_descr%collcom%ndimind_f, psitlarge_f(1), 1, psitlarge_f(1), 1)
  hpsittmp_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsittmp_c')
  hpsittmp_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsittmp_f')
!print*,'a',ddot(weight_matrix%nvctr, weight_matrix_%matrix_compr(1), 1, weight_matrix_%matrix_compr(1), 1)

      if(tmb%ham_descr%collcom%ndimind_c>0) &
          call vcopy(tmb%ham_descr%collcom%ndimind_c, psitlarge_c(1), 1, hpsittmp_c(1), 1)
      if(tmb%ham_descr%collcom%ndimind_f>0) &
          call vcopy(7*tmb%ham_descr%collcom%ndimind_f, psitlarge_f(1), 1, hpsittmp_f(1), 1)

  call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%m, weight_matrix_)
  call build_linear_combination_transposed(tmb%ham_descr%collcom, &
       tmb%linmat%m, weight_matrix_, hpsittmp_c, hpsittmp_f, .true., psitlarge_c, psitlarge_f, bigdft_mpi%iproc)
!print*,'2',ddot(tmb%ham_descr%collcom%ndimind_c, psitlarge_c(1), 1, psitlarge_c(1), 1),ddot(7*tmb%ham_descr%collcom%ndimind_f, psitlarge_f(1), 1, psitlarge_f(1), 1)
  call f_free_ptr(hpsittmp_c)
  call f_free_ptr(hpsittmp_f)

  !call f_free_ptr(psitlarge_c)
  !call f_free_ptr(psitlarge_f)

  !copy back original weight matrix
  call vcopy(weight_matrix%nvctr,matrix_compr(1),1,weight_matrix_%matrix_compr(1),1)
  call f_free_ptr(matrix_compr)

  call f_release_routine()

end subroutine calculate_weight_matrix_lowdin_gradient



!> CDFT: calculates the weight matrix w_ab given w(r)
!! for the moment putting densities in global box and ignoring parallelization
subroutine calculate_weight_matrix_using_density(iproc,cdft,tmb,at,input,GPU,denspot)
  use module_base
  use module_types
  use constrained_dft, only: cdft_data
  use module_interfaces, except_this_one => calculate_weight_matrix_using_density
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only : matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace
  use transposed_operations, only: calculate_overlap_transposed
  use potential, only: full_local_potential
  implicit none
  integer,intent(in) :: iproc
  type(cdft_data), intent(inout) :: cdft
  type(atoms_data), intent(in) :: at
  type(input_variables),intent(in) :: input
  type(dft_wavefunction), intent(inout) :: tmb
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU

  !integer :: iorb, jorb
  real(kind=gp),dimension(:),allocatable :: hpsit_c, hpsit_f
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  type(energy_terms) :: energs
  character(len=*),parameter :: subname='calculate_weight_matrix_using_density'
  type(matrices) :: weight_

  energs=energy_terms_null()
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(bigdft_mpi%iproc,bigdft_mpi%nproc,&
       denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),max(denspot%dpbox%nscatterarr(:,2),1),cdft%weight_function, &
       tmb%ham_descr%comgp%nrecvbuf,tmb%ham_descr%comgp%recvbuf,tmb%ham_descr%comgp,tmb%ham_descr%lzd)

  allocate(confdatarrtmp(tmb%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

  call small_to_large_locreg(bigdft_mpi%iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%orbs, tmb%psi, tmb%ham_descr%psi)

  if (tmb%ham_descr%npsidim_orbs > 0) call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

  call full_local_potential(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%lzd,2, &
       denspot%dpbox,denspot%xc,cdft%weight_function,denspot%pot_work,tmb%ham_descr%comgp)
  call LocalHamiltonianApplication(bigdft_mpi%iproc,bigdft_mpi%nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
       energs,input%SIC,GPU,2,denspot%xc,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
       potential=cdft%weight_function,comgp=tmb%ham_descr%comgp)

  deallocate(confdatarrtmp)

  print *,'CDFT Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc

  call f_free_ptr(denspot%pot_work)


  ! calculate w_ab
  ! Calculate the matrix elements <phi|H|phi>.
  if(.not.tmb%ham_descr%can_use_transposed) then
      !!if(associated(tmb%ham_descr%psit_c)) then
      !!    call f_free_ptr(tmb%ham_descr%psit_c)
      !!end if
      !!if(associated(tmb%ham_descr%psit_f)) then
      !!    call f_free_ptr(tmb%ham_descr%psit_f)
      !!end if

      !!tmb%ham_descr%psit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='tmb%ham_descr%psit_c')
      !!tmb%ham_descr%psit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='tmb%ham_descr%psit_f')
      call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs, &
           tmb%ham_descr%collcom,TRANSPOSE_FULL,tmb%ham_descr%psi,tmb%ham_descr%psit_c,tmb%ham_descr%psit_f,tmb%ham_descr%lzd)
      tmb%ham_descr%can_use_transposed=.true.
  end if

  hpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
  hpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
  call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,  &
       tmb%ham_descr%collcom,TRANSPOSE_FULL,tmb%hpsi,hpsit_c,hpsit_f,tmb%ham_descr%lzd)

  weight_ = matrices_null()
  call allocate_matrices(tmb%linmat%m, allocate_full=.false., &
       matname='weight_', mat=weight_)

  call calculate_overlap_transposed(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%collcom, &
       tmb%ham_descr%psit_c,hpsit_c,tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, weight_)
  call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc,bigdft_mpi%nproc, tmb%linmat%m, weight_)
  ! This can then be deleted if the transition to the new type has been completed.
  cdft%weight_matrix_%matrix_compr=weight_%matrix_compr
  call deallocate_matrices(weight_)

  call f_free(hpsit_c)
  call f_free(hpsit_f)

  ! debug
  !allocate(cdft%weight_matrix%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  !call memocc(istat, cdft%weight_matrix%matrix, 'cdft%weight_matrix%matrix', subname)
  !call uncompress_matrix(bigdft_mpi%iproc,cdft%weight_matrix)
  !do iorb=1,tmb%orbs%norb
  !   do jorb=1,tmb%orbs%norb
  !      write(87,*) iorb,jorb,cdft%weight_matrix%matrix(iorb,jorb)
  !   end do
  !end do
  !deallocate(cdft%weight_matrix%matrix, stat=istat)
  !call memocc(istat, iall, 'cdft%weight_matrix%matrix', subname)
  ! end debug

end subroutine calculate_weight_matrix_using_density


!> CDFT: calculates the weight function w(r)
!! for the moment putting densities in global box and ignoring parallelization
subroutine calculate_weight_function(in,ref_frags,cdft,ndimrho_all_fragments,rho_all_fragments,tmb,atoms,rxyz,denspot)
  use module_base
  use module_types
  use module_fragments
  use constrained_dft
  implicit none
  type(input_variables), intent(in) :: in
  type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags
  type(cdft_data), intent(inout) :: cdft
  integer, intent(in) :: ndimrho_all_fragments
  real(kind=wp),dimension(ndimrho_all_fragments),intent(in) :: rho_all_fragments
  type(dft_wavefunction), intent(in) :: tmb
  type(atoms_data), intent(in) :: atoms ! just for plotting
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz ! just for plotting
  type(DFT_local_fields), intent(in) :: denspot ! just for plotting

  integer :: ifrag, ifrag_charged, ref_frag_charged, iorb_start, ipt

  ! unecessary recalculation here, but needs generalizing anyway
  iorb_start=1
  do ifrag=1,in%frag%nfrag
     if (in%frag%charge(ifrag)/=0) then
        ifrag_charged=ifrag
        exit
     end if
     iorb_start=iorb_start+ref_frags(in%frag%frag_index(ifrag))%fbasis%forbs%norb
  end do
  ref_frag_charged=in%frag%frag_index(ifrag_charged)

  ! check both densities are the same size (for now calculating fragment density in global box)
  if (ndimrho_all_fragments /= cdft%ndim_dens) stop 'Fragment density dimension does not match global density dimension'

  ! only need to calculate the fragment density for the fragment with a charge (stored in frag%fbasis%density)
  ! ideally would use already calculated kernel here, but charges could vary so would need an actual fragment
  ! rather than a reference fragment - could point to original fragment but modify nks...
  ! also converting input charge to integer which might not be the case
  call calculate_fragment_density(ref_frags(ref_frag_charged),ndimrho_all_fragments,tmb,iorb_start,&
       int(in%frag%charge(ifrag_charged)),atoms,rxyz,denspot)

  ! the weight function becomes the fragment density of ifrag_charged divided by the sum of all fragment densities
  ! using a cutoff of 1.e-6 for fragment density and 1.e-12 for denominator
  do ipt=1,ndimrho_all_fragments
      if (denspot%rhov(ipt) >= 1.0e-12) then
         cdft%weight_function(ipt)=ref_frags(ref_frag_charged)%fbasis%density(ipt)/denspot%rhov(ipt)
      else
         cdft%weight_function(ipt)=0.0_gp
      end if
  end do

  ! COME BACK TO THIS
  cdft%charge=ref_frags(ref_frag_charged)%nelec-in%frag%charge(ref_frag_charged)

  ! plot the weight function, fragment density and initial total density (the denominator) to check

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'fragment_density.cube', &
       atoms,rxyz,denspot%dpbox,1,ref_frags(ref_frag_charged)%fbasis%density)

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'weight_function.cube', &
       atoms,rxyz,denspot%dpbox,1,cdft%weight_function)

  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'initial_density.cube', &
       atoms,rxyz,denspot%dpbox,1,denspot%rhov)

  ! deallocate fragment density here as we don't need it any more
  if (associated(ref_frags(ref_frag_charged)%fbasis%density)) call f_free_ptr(ref_frags(ref_frag_charged)%fbasis%density)

end subroutine calculate_weight_function
