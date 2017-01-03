!> @file
!! Optimize the coefficients
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Coefficients are defined for Ntmb KS orbitals so as to maximize the number
!!  of orthonormality constraints. This should speedup the convergence by
!!  reducing the effective number of degrees of freedom.
subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, sd_fit_curve, &
    factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, num_extra)
  use module_base
  use module_types
  !use module_interfaces, fake_name => optimize_coeffs
  use diis_sd_optimization
  use yaml_output
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, itmax, itout, it_scc, it_cdft
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(orbitals_data),intent(inout):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  type(DIIS_obj), intent(inout) :: ldiis_coeff
  real(kind=gp),intent(in):: fnrm_crit
  real(kind=gp),intent(out):: fnrm
  real(kind=gp), intent(inout) :: energy
  logical, intent(in) :: sd_fit_curve
  real(kind=gp), intent(in) :: factor
  integer, optional, intent(in) :: num_extra
  logical, optional, intent(in) :: reorder

  ! Local variables
  integer:: iorb, jorb, iiorb, ierr, it, itlast, ispin
  integer, dimension(1) :: iproc_arr, ncomp
  real(kind=gp),dimension(:,:),allocatable:: grad, grad_cov_or_coeffp !coeffp, grad_cov
  real(kind=gp),dimension(:),allocatable:: mat_coeff_diag, occup_tmp
  real(kind=gp) :: tt, ddot, energy0, pred_e


  !!if (present(num_extra) .and. tmb%linmat%m%nspin==2) then
  !!    stop 'ERROR: optimize_coeffs not yet implemented for nspin=2 and num_extra present!'
  !!end if

  call f_routine(id='optimize_coeffs')


  if (present(num_extra)) then
     if (num_extra > 0) then
        occup_tmp=f_malloc(orbs%norb,id='occup_tmp')
        ! make a copy of 'correct' occup
        call vcopy(orbs%norb, orbs%occup(1), 1, occup_tmp(1), 1)
     end if
  end if


  if (ldiis_coeff%idsx == 0 .and. sd_fit_curve) then
     ! calculate initial energy for SD line fitting and printing (maybe don't need to (re)calculate kernel here?)
     !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
     !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
          tmb%linmat%kernel_, tmb%linmat%ham_, energy0,&
          tmb%coeff,orbs,tmb%orbs,.true.)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
     !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
  else
     energy0=energy
  end if

  grad=f_malloc((/tmb%linmat%m%nfvctr,orbs%norbp/), id='grad')
  grad_cov_or_coeffp=f_malloc((/tmb%linmat%m%nfvctr,orbs%norbp/), id='grad_cov_or_coeffp')

  if (iproc==0) then
      call yaml_newline()
      call yaml_sequence_open('expansion coefficients optimization',label=&
           'it_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
           trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//&
           '_'//trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))))
  end if


  do it=1,itmax

      if (iproc==0) then
          call yaml_newline()
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_comment('it coeff:'//yaml_toa(it,fmt='(i6)'),hfill='-')
      end if

     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(5300+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
     !!    end do
     !!end do

     if (present(num_extra)) then
        if (num_extra > 0) then
           if (tmb%linmat%m%nspin==1) then
               tt = 2.0_gp
           else if (tmb%linmat%m%nspin==2) then
               tt = 1.0_gp
           end if
           do iorb=1,orbs%norb
              orbs%occup(iorb)=tt
           end do
        end if
     end if

     call calculate_coeff_gradient(iproc,nproc,tmb,order_taylor,max_inversion_error,orbs,grad_cov_or_coeffp,grad)

     if (present(num_extra)) then
        if (num_extra > 0) then
           call vcopy(orbs%norb, occup_tmp(1), 1, orbs%occup(1), 1)
        end if
     end if

     ! scale the gradient by a factor of 2 in case of spin polarized calculation
     ! in order to make the gradient analogous to the case of no polarization
     ! (in that case the gradient is included in the kernel).
     if (tmb%linmat%l%nspin==2) then
         call dscal(orbs%norbp*tmb%linmat%m%nfvctr, 2.d0, grad, 1)
     end if

     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(5400+10*iproc+ispin,'(a,2i8,2es18.7)') 'iiorb, jorb, coeff(jorb,iiorb),grad(jorb,iorb)', iiorb, jorb, tmb%coeff(jorb,iiorb), grad(jorb,iorb)
     !!    end do
     !!end do

     ! Precondition the gradient (only making things worse...)
     !call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, grad)

     call timing(iproc,'dirmin_sddiis','ON')

     !For fnrm, we only sum on the occupied KS orbitals
     tt=0.d0
     do iorb=1,orbs%norbp
         tt=tt+ddot(tmb%linmat%m%nfvctr, grad_cov_or_coeffp(1,iorb), 1, grad(1,iorb), 1)
     end do

     if (nproc > 1) then
        call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     end if
     fnrm=2.0_gp*tt

     !scale the gradient - useful for fragments/constrained
     call dscal(tmb%linmat%m%nfvctr*orbs%norbp,factor,grad,1)

     !write(*,*) 'present(num_extra), ldiis_coeff%idsx', present(num_extra), ldiis_coeff%idsx
     if (ldiis_coeff%idsx > 0) then !do DIIS
        !TO DO: make sure DIIS works
        ldiis_coeff%mids=mod(ldiis_coeff%ids,ldiis_coeff%idsx)+1
        ldiis_coeff%ids=ldiis_coeff%ids+1

        if (orbs%norbp>0) call vcopy(tmb%linmat%m%nfvctr*orbs%norbp,tmb%coeff(1,orbs%isorb+1), &
            1,grad_cov_or_coeffp(1,1),1)

        iproc_arr(1)=iproc
        ncomp(1)=tmb%linmat%m%nfvctr*orbs%norbp
        call diis_opt(iproc,nproc,1,0,1,iproc_arr,ncomp,tmb%linmat%m%nfvctr*orbs%norbp, &
             grad_cov_or_coeffp,grad,ldiis_coeff) 
     else  !steepest descent with curve fitting for line minimization
        call timing(iproc,'dirmin_sddiis','OF')
        if (sd_fit_curve) call find_alpha_sd(iproc,nproc,ldiis_coeff%alpha_coeff,tmb,orbs,&
             grad_cov_or_coeffp,grad,energy0,fnrm,pred_e,order_taylor)
        call timing(iproc,'dirmin_sddiis','ON')
        do iorb=1,orbs%norbp
           iiorb = orbs%isorb + iorb
           do jorb=1,tmb%linmat%m%nfvctr
              grad_cov_or_coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*grad(jorb,iorb)
           end do
        end do
     end if

     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(4500+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, grad_cov_or_coeffp(jorb,iorb)', iiorb, jorb, grad_cov_or_coeffp(jorb,iorb)
     !!    end do
     !!end do

     call timing(iproc,'dirmin_sddiis','OF')
     

     call timing(iproc,'dirmin_allgat','ON')
     if(nproc > 1) then 
        call mpi_allgatherv(grad_cov_or_coeffp, tmb%linmat%m%nfvctr*orbs%norbp, mpi_double_precision, tmb%coeff, &
           tmb%linmat%m%nfvctr*orbs%norb_par(:,0), tmb%linmat%m%nfvctr*orbs%isorb_par, mpi_double_precision, &
           bigdft_mpi%mpi_comm, ierr)
     else
        call vcopy(tmb%linmat%m%nfvctr*orbs%norb,grad_cov_or_coeffp(1,1),1,tmb%coeff(1,1),1)
     end if

     call timing(iproc,'dirmin_allgat','OF')

     fnrm=sqrt(fnrm/dble(orbs%norb))
     ! Multiply the gradient by sqrt(2) to make it homologous to the case of
     ! no spin polarization (since with polarization norb is doubled).
     if (tmb%linmat%l%nspin==2) then
         fnrm=fnrm*sqrt(2.d0)
     end if


     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(4600+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
     !!    end do
     !!end do

     !! experimenting with calculating cHc and diagonalizing
     !if (present(num_extra).and.present(reorder)) then
     !   call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
     !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
     !   call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !   if (reorder) call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, &
     !        tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !else if (present(reorder)) then
     !   call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !   if (reorder) call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !else
     !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, .false.)
     !end if

     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(5200+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
     !!    end do
     !!end do


     ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
     ! instead of twice could add some criterion to check accuracy?
     ! first need to check if we're at least close to normalized, otherwise reorthonormalize_coeff may explode (should maybe incorporate this into reorthonormalize_coeff)

     ! should really check this first, for now just normalize anyway
     !mat_coeff_diag=f_malloc(orbs%norb,id='mat_coeff_diag')
     !call calculate_coeffMatcoeff_diag(tmb%linmat%ovrlp_%matrix,tmb%orbs,orbs,tmb%coeff,mat_coeff_diag)
     !do iorb=1,orbs%norb
     !   call dscal(tmb%orbs%norb,1.0d0/sqrt(mat_coeff_diag(iorb)),tmb%coeff(1,iorb),1)
     !end do
     !call f_free(mat_coeff_diag)

     call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, order_taylor, &
          tmb%orbs, tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, orbs)
     !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     !!!!!!!!!!!!!!!!!!!!!!!!
     !can't put coeffs directly in ksorbs%eval as intent in, so change after - problem with orthonormality of coeffs so adding extra
     !call find_eval_from_coeffs(iproc, nproc, tmb%orthpar%methTransformOverlap, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, &
     !     tmb%coeff, tmb%orbs%eval, .true., .true.)
     !ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
     !call order_coeffs_by_energy(orbs%norb,tmb%orbs%norb,tmb%coeff,tmb%orbs%eval,ipiv)
     !call f_free(ipiv)
     !!!!!!!!!!!!!!!!!!!!!!!!

     !!do iorb=1,orbs%norbp
     !!    iiorb=orbs%isorb+iorb
     !!    if (orbs%spinsgn(iiorb)>0.d0) then
     !!        ispin=1
     !!    else
     !!        ispin=2
     !!    end if
     !!    do jorb=1,tmb%linmat%m%nfvctr
     !!        write(5100+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
     !!    end do
     !!end do

     !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
     !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m,&
          tmb%linmat%kernel_, tmb%linmat%ham_, energy,&
          tmb%coeff,orbs,tmb%orbs,.true.)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)

     !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
     !write(127,*) ldiis_coeff%alpha_coeff,energy
     !close(127)

     ! can't check ebs only, need to check Etot in linearScaling, but if we do it>1 without curve fitting will still need to update here
     !if (ldiis_coeff%idsx == 0 .and. (.not. sd_fit_curve) .and. energy0/=0.0_gp) then ! only update alpha after first iteration
     !   ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
     !   if ((energy-energy0)<0.d0 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
     !      ldiis_coeff%alpha_coeff=1.1d0*ldiis_coeff%alpha_coeff
     !   else if (ldiis_coeff%alpha_coeff > 1.7d-3) then
     !      ldiis_coeff%alpha_coeff=0.5d0*ldiis_coeff%alpha_coeff
     !   end if
     !   if (iproc==0) print*,'EBSdiff,alpha',energy-energy0,ldiis_coeff%alpha_coeff,energy,energy0
     !end if

     !if (iproc==0) write(*,*) ''
     if (sd_fit_curve .and. ldiis_coeff%idsx == 0) then
        !!if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha, pred E, diff',&
        !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff,pred_e,pred_e-energy
        if (iproc==0) then
            call yaml_map('method','DminSD')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
            call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
            call yaml_map('predicted energy',pred_e,fmt='(es24.17)')
            call yaml_map('D',pred_e-energy,fmt='(es10.3)')
        end if
     else if (ldiis_coeff%idsx == 0) then
        !!if (iproc==0) write(*,'(a,I4,2x,4(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha',&
        !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff
        if (iproc==0) then
            call yaml_map('method','DminSD')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
            call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
        end if
     else
        !!if (iproc==0) write(*,'(a,I4,2x,3(ES16.6e3,2x))')'DminDIIS: it, fnrm, ebs, ebsdiff',&
        !!     it,fnrm,energy0,energy-energy0
        if (iproc==0) then
            call yaml_map('method','DminDIIS')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
        end if
     end if

     energy0=energy

     if (iproc==0) then
         call yaml_mapping_close()
         call yaml_flush_document()
         !call bigdft_utils_flush(unit=6)
     end if


     itlast=it
     if (fnrm<fnrm_crit) exit

  end do

  !!if (iproc==0) then
  !!    call yaml_newline()
  !!    call yaml_sequence(label='final_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
  !!         trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//'_'//&
  !!         trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))),advance='no')
  !!    call yaml_mapping_open(flow=.true.)
  !!    call yaml_comment('iter:'//yaml_toa(itlast,fmt='(i6)'),hfill='-')
  !!    call yaml_map('iter',itlast,fmt='(i6)')
  !!    call yaml_map('fnrm',fnrm,fmt='(es9.2)')
  !!    call yaml_map('eBS',energy0,fmt='(es24.17)')
  !!    call yaml_map('D',energy-energy0,fmt='(es10.3)')
  !!    call yaml_mapping_close()
  !!    call bigdft_utils_flush(unit=6)
  !!end if


  if (iproc==0) then
      call yaml_sequence_close()
      call yaml_newline()
  end if


  call f_free(grad_cov_or_coeffp)
  call f_free(grad)

  if (present(num_extra)) then
     if (num_extra > 0) then
        call f_free(occup_tmp)
     end if
  end if


  call f_release_routine()

end subroutine optimize_coeffs


!> subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
!! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
subroutine coeff_weight_analysis(iproc, nproc, input, ksorbs, tmb, ref_frags)
  use module_base
  use module_types
  use module_fragments
  use constrained_dft
  use yaml_output
  use sparsematrix_base, only: sparse_matrix, matrices, sparse_matrix_null, deallocate_sparse_matrix, &
                               sparsematrix_malloc_ptr, DENSE_FULL, SPARSE_FULL, assignment(=), &
                               matrices_null, allocate_matrices, deallocate_matrices,copy_sparse_matrix
  use sparsematrix, only: uncompress_matrix, uncompress_matrix2
  use matrix_operations, only: overlapPowerGeneral
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc
  type(orbitals_data), intent(in) :: ksorbs
  type(dft_wavefunction), intent(inout) :: tmb
  type(input_variables),intent(in) :: input
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

  integer :: iorb, istat, iall, ifrag
  integer, dimension(2) :: ifrag_charged
  integer, dimension(1) :: power
  !real(kind=8), dimension(:,:,:), allocatable :: weight_coeff
  real(kind=8), dimension(:,:), allocatable :: weight_coeff_diag
  real(kind=8), dimension(:,:), pointer :: ovrlp_half
  real(kind=8), dimension(:), allocatable :: weight_sum
  real(kind=8) :: max_error, mean_error, occsum
  type(sparse_matrix) :: weight_matrix
  type(matrices) :: weight_matrix_
  type(matrices),dimension(1) :: inv_ovrlp
  character(len=256) :: subname='coeff_weight_analysis'

  call timing(iproc,'weightanalysis','ON')

  call f_routine('coeff_weight_analysis')

  !call nullify_sparse_matrix(weight_matrix)
  weight_matrix=sparse_matrix_null()
  weight_matrix_=matrices_null()
  !call sparse_copy_pattern(tmb%linmat%m, weight_matrix, iproc, subname)
  call copy_sparse_matrix(tmb%linmat%m, weight_matrix)
  !weight_matrix_%matrix_compr=f_malloc_ptr(weight_matrix%nvctr,id='weight_matrix%matrix_compr')
  weight_matrix_%matrix_compr=sparsematrix_malloc_ptr(weight_matrix,iaction=SPARSE_FULL,id='weight_matrix%matrix_compr')

  inv_ovrlp(1) = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

  !weight_coeff=f_malloc((/ksorbs%norb,ksorbs%norb,input%frag%nfrag/), id='weight_coeff')
  weight_coeff_diag=f_malloc((/ksorbs%norb,input%frag%nfrag/), id='weight_coeff')
  !ovrlp_half=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
  tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
  power(1)=2
  call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
       input%lin%order_taylor, 1, power, &
       tmb%orthpar%blocksize_pdsyev, imode=2, &
       ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
       ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
       max_error=max_error, mean_error=mean_error)
  call f_free_ptr(tmb%linmat%ovrlp_%matrix)

  do ifrag=1,input%frag%nfrag
     ifrag_charged(1)=ifrag
     call calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,1,ifrag_charged,tmb,input%frag,ref_frags,&
          .false.,.true.,input%lin%order_taylor)
     weight_matrix_%matrix = sparsematrix_malloc_ptr(weight_matrix,iaction=DENSE_FULL,id='weight_matrix%matrix')
     call uncompress_matrix(iproc,nproc,weight_matrix,weight_matrix_%matrix_compr,weight_matrix_%matrix)
     !call calculate_coeffMatcoeff(nproc,weight_matrix%matrix,tmb%orbs,ksorbs,tmb%coeff,weight_coeff(1,1,ifrag))
     call calculate_coeffMatcoeff_diag(weight_matrix_%matrix,tmb%orbs,ksorbs,tmb%coeff,weight_coeff_diag(1,ifrag))
     call f_free_ptr(weight_matrix_%matrix)
  end do
  !call f_free_ptr(ovrlp_half)

  if (iproc==0) call yaml_sequence_open('Weight analysis',flow=.true.)
  if (iproc==0) call yaml_newline()
  if (iproc==0) call yaml_comment ('coeff, occ, eval, frac for each frag')
  ! only care about diagonal elements
  weight_sum=f_malloc(input%frag%nfrag,id='weight_sum')
  weight_sum=0.0d0
  do iorb=1,ksorbs%norb
     if (iproc==0) then
         call yaml_mapping_open(flow=.true.)
         call yaml_map('iorb',iorb,fmt='(i4)')
         call yaml_map('occ',KSorbs%occup(iorb),fmt='(f6.4)')
         call yaml_map('eval',tmb%orbs%eval(iorb),fmt='(f10.6)')
     end if
     do ifrag=1,input%frag%nfrag
        if (iproc==0) call yaml_map('frac',weight_coeff_diag(iorb,ifrag),fmt='(f6.4)')
        weight_sum(ifrag)=weight_sum(ifrag)+KSorbs%occup(iorb)*weight_coeff_diag(iorb,ifrag)
     end do
     if (iproc==0) call yaml_mapping_close()
     if (iproc==0) call yaml_newline()
  end do
     
  !sum for each fragment
  occsum=sum(KSorbs%occup(1:ksorbs%norb))
  if (iproc==0) then
     call yaml_mapping_open(flow=.true.)
     call yaml_map('sum',occsum,fmt='(f8.4)')
  end if
  do ifrag=1,input%frag%nfrag
     if (iproc==0) call yaml_map('frac',weight_sum(ifrag),fmt='(f6.4)')
  end do
  if (iproc==0) call yaml_mapping_close()
  if (iproc==0) call yaml_newline()
  call f_free(weight_sum)

  if (iproc==0) call yaml_sequence_close()

  call deallocate_matrices(inv_ovrlp(1))

  call deallocate_sparse_matrix(weight_matrix)
  call deallocate_matrices(weight_matrix_)
  call f_free(weight_coeff_diag)
  !call f_free(weight_coeff)

  call f_release_routine()

  call timing(iproc,'weightanalysis','OF')


end subroutine coeff_weight_analysis


!!! subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
!!! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
!!subroutine find_eval_from_coeffs(iproc, nproc, meth_overlap, ksorbs, basis_orbs, ham, ovrlp, coeff, eval, calc_overlap, diag)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use sparsematrix_base, only: sparse_matrix
!!  implicit none
!!
!!  ! Calling arguments
!!  integer, intent(in) :: iproc, nproc, meth_overlap
!!  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
!!  type(sparse_matrix),intent(in) :: ham
!!  type(sparse_matrix),intent(inout) :: ovrlp
!!  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
!!  real(kind=8),dimension(ksorbs%norb),intent(inout) :: eval
!!  logical, intent(in) :: diag, calc_overlap
!!
!!  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
!!  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff,  ovrlp_coeff
!!  real(kind=8) :: offdiagsum, offdiagsum2, coeff_orthog_threshold
!!  character(len=256) :: subname='reordering_coeffs'
!!
!!  coeff_orthog_threshold=1.0d-3
!!
!!  allocate(ham_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
!!  call memocc(istat, ham_coeff, 'ham_coeff', subname)
!!
!!  call calculate_coeffMatcoeff(ham%matrix,basis_orbs,ksorbs,coeff,ham_coeff)
!!
!!  if (calc_overlap) then
!!     allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
!!     call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
!!     call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
!!  end if
!!
!!  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
!!  offdiagsum=0.0d0
!!  offdiagsum2=0.0d0
!!  do iorb=1,ksorbs%norb
!!     do jorb=1,ksorbs%norb
!!        if (iorb==jorb) then
!!           eval(iorb)=ham_coeff(iorb,iorb)
!!        else
!!           offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
!!           if (calc_overlap) offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
!!        end if
!!     end do
!!  end do
!!  offdiagsum=offdiagsum/(ksorbs%norb**2-ksorbs%norb)
!!  if (calc_overlap) offdiagsum2=offdiagsum2/(ksorbs%norb**2-ksorbs%norb)
!!  if (calc_overlap.and.iproc==0) print*,''
!!  if (calc_overlap) then
!!     if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2
!!  else
!!     if (iproc==0) print*,'offdiagsum (ham):',offdiagsum
!!  end if
!!
!!  ! if coeffs are too far from orthogonality
!!  if (calc_overlap .and. offdiagsum2>coeff_orthog_threshold) then
!!     call reorthonormalize_coeff(iproc, nproc, ksorbs%norb, -8, -8, meth_overlap, basis_orbs, ovrlp, coeff, ksorbs)
!!  end if
!!
!!  if (diag.or.offdiagsum>1.0d-2) then
!!     ! diagonalize within the space of occ+extra
!!     if (.not.calc_overlap) then
!!        ! assume ovrlp_coeff is orthgonal for now
!!        allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
!!        call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
!!        call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
!!        do iorb=1,ksorbs%norb
!!           do jorb=1,ksorbs%norb
!!              if (iorb==jorb) then
!!                 ovrlp_coeff(iorb,jorb)=1.0d0
!!              else
!!                 ovrlp_coeff(iorb,jorb)=0.0d0
!!              end if
!!           end do
!!        end do
!!     end if
!!     ! diagonalize within the space of occ+extra
!!     call diagonalizeHamiltonian2(iproc, ksorbs%norb, ham_coeff, ovrlp_coeff, eval)
!!     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb/),id='coeff_tmp')
!!
!!     ! multiply new eigenvectors by coeffs
!!     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb, ksorbs%norb, 1.d0, coeff(1,1), &
!!          basis_orbs%norb, ham_coeff, ksorbs%norb, 0.d0, coeff_tmp, basis_orbs%norb)
!! 
!!     call vcopy(basis_orbs%norb*(ksorbs%norb),coeff_tmp(1,1),1,coeff(1,1),1)
!!     call f_free(coeff_tmp)
!!  end if
!!
!!  if (calc_overlap.or.diag) then
!!     iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
!!     deallocate(ovrlp_coeff,stat=istat)
!!     call memocc(istat,iall,'ovrlp_coeff',subname)
!!  end if
!!
!!  iall=-product(shape(ham_coeff))*kind(ham_coeff)
!!  deallocate(ham_coeff,stat=istat)
!!  call memocc(istat,iall,'ham_coeff',subname)
!!
!!end subroutine find_eval_from_coeffs


subroutine calculate_coeffMatcoeff(nproc,matrix,basis_orbs,ksorbs,coeff,mat_coeff)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: nproc
  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(in) :: matrix
  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
  real(kind=8),dimension(ksorbs%norb,ksorbs%norb),intent(inout) :: mat_coeff

  integer :: iall, istat, ierr
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp
  character(len=256) :: subname='calculate_coeffMatcoeff'

  coeff_tmp = f_malloc((/ basis_orbs%norbp, ksorbs%norb /),id='coeff_tmp')

  if (basis_orbs%norbp>0) then
     call dgemm('n', 'n', basis_orbs%norbp, ksorbs%norb, basis_orbs%norb, 1.d0, matrix(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
     call dgemm('t', 'n', ksorbs%norb, ksorbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, mat_coeff, ksorbs%norb)
  else
     call f_zero(mat_coeff)
  end if

  call f_free(coeff_tmp)

  if (nproc>1) then
      call mpiallred(mat_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

end subroutine calculate_coeffMatcoeff


!> Same as above but only calculating diagonal elements
subroutine calculate_coeffMatcoeff_diag(matrix,basis_orbs,ksorbs,coeff,mat_coeff_diag)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(in) :: matrix
  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
  real(kind=8),dimension(ksorbs%norb),intent(inout) :: mat_coeff_diag

  integer :: iall, istat, ierr, iorb
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp
  real(kind=8), dimension(:), allocatable :: mat_coeff_diagp
  logical, parameter :: allgather=.true.


  if (allgather) then
     mat_coeff_diagp=f_malloc((/ksorbs%norbp/), id='mat_coeff_diagp')
  else
     call f_zero(mat_coeff_diag)
  end if
  if (ksorbs%norbp>0) then
     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norbp/), id='coeff_tmp')
     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norbp, basis_orbs%norb, 1.d0, matrix(1,1), &
          basis_orbs%norb, coeff(1,ksorbs%isorb+1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
     if (allgather) then
        do iorb=1,ksorbs%norbp
           call dgemm('t', 'n', 1, 1, basis_orbs%norb, 1.d0, coeff(1,ksorbs%isorb+iorb), basis_orbs%norb, &
                coeff_tmp(1,iorb), basis_orbs%norb, 0.d0, mat_coeff_diagp(iorb), ksorbs%norb)
        end do
     else
        do iorb=1,ksorbs%norbp
           call dgemm('t', 'n', 1, 1, basis_orbs%norb, 1.d0, coeff(1,ksorbs%isorb+iorb), basis_orbs%norb, &
                coeff_tmp(1,iorb), basis_orbs%norb, 0.d0, mat_coeff_diag(ksorbs%isorb+iorb), ksorbs%norb)
        end do
     end if
     call f_free(coeff_tmp)
  end if

  if (bigdft_mpi%nproc>1) then
     if (allgather) then
        call mpi_allgatherv(mat_coeff_diagp, ksorbs%norbp, mpi_double_precision, mat_coeff_diag, &
             ksorbs%norb_par(:,0), ksorbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call mpiallred(mat_coeff_diag, mpi_sum, comm=bigdft_mpi%mpi_comm)
     end if
  else
     if (allgather) then
        call vcopy(ksorbs%norb, mat_coeff_diagp(1), 1, mat_coeff_diag(1), 1)
     end if
  end if

  if (allgather) then
     call f_free(mat_coeff_diagp)
 end if

end subroutine calculate_coeffMatcoeff_diag
 
!> not really fragment related so prob should be moved - reorders coeffs by eval
  ! output ipiv in case want to use it for something else
  subroutine order_coeffs_by_energy(nstate,ntmb,coeff,eval,ipiv)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nstate, ntmb
  real(kind=gp), dimension(ntmb,nstate), intent(inout) :: coeff
  real(kind=gp), dimension(nstate), intent(inout) :: eval
    integer, dimension(nstate), intent(out) :: ipiv

  integer :: itmb, jorb
    !integer, allocatable, dimension(:) :: ipiv
  real(gp), dimension(:), allocatable :: tmp_array
  real(gp), dimension(:,:), allocatable :: tmp_array2

    !ipiv=f_malloc(nstate,id='coeff_final')
  tmp_array=f_malloc(nstate,id='tmp_array')

  do itmb=1,nstate
     tmp_array(itmb)=-eval(itmb)
  end do

  call sort_positions(nstate,tmp_array,ipiv)

  do itmb=1,nstate
     eval(itmb)=-tmp_array(ipiv(itmb))
  end do

  call f_free(tmp_array)

  tmp_array2=f_malloc((/ntmb,nstate/),id='tmp_array2')

  do jorb=1,nstate
     do itmb=1,ntmb
        tmp_array2(itmb,jorb)=coeff(itmb,jorb)
     end do
  end do

  do jorb=1,nstate
     do itmb=1,ntmb
        coeff(itmb,jorb)=tmp_array2(itmb,ipiv(jorb))
     end do
  end do

  call f_free(tmp_array2)
    !call f_free(ipiv)

end subroutine order_coeffs_by_energy

!!! experimental - for num_extra see if extra states are actually lower in energy and reorder (either by sorting or diagonalization)
!!! could improve efficiency by only calculating cSc or cHc up to norb (i.e. ksorbs%norb+extra)
!!subroutine reordering_coeffs(iproc, nproc, num_extra, ksorbs, basis_orbs, ham, ovrlp, coeff, reorder)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use sparsematrix_base, only: sparse_matrix
!!  implicit none
!!
!!  ! Calling arguments
!!  integer, intent(in) :: iproc, nproc, num_extra
!!  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
!!  type(sparse_matrix),intent(in) :: ham, ovrlp
!!  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
!!  logical, intent(in) :: reorder
!!
!!  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
!!  integer, allocatable, dimension(:) :: ipiv
!!  real(gp), dimension(:), allocatable :: tmp_array, eval
!!  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff, tmp_array2, ham_coeff_small, ovrlp_coeff_small, ovrlp_coeff
!!  real(kind=8) :: offdiagsum, offdiagsum2
!!  character(len=256) :: subname='reordering_coeffs'
!!
!!  allocate(ham_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
!!  call memocc(istat, ham_coeff, 'ham_coeff', subname)
!!
!!  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
!!  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
!!
!!  if (basis_orbs%norbp>0) then
!!     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ham%matrix(basis_orbs%isorb+1,1), &
!!          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
!!     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
!!          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ham_coeff, basis_orbs%norb)
!!  else
!!     call to_zero(basis_orbs%norb**2, ham_coeff(1,1))
!!  end if
!!
!!  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
!!  deallocate(coeff_tmp,stat=istat)
!!  call memocc(istat,iall,'coeff_tmp',subname)
!!
!!  if (nproc>1) then
!!      call mpiallred(ham_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!  end if
!!
!!  allocate(ovrlp_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
!!  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
!!
!!  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
!!  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
!!
!!  if (basis_orbs%norbp>0) then
!!     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ovrlp%matrix(basis_orbs%isorb+1,1), &
!!          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
!!     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
!!          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, basis_orbs%norb)
!!  else
!!     call to_zero(basis_orbs%norb**2, ovrlp_coeff(1,1))
!!  end if
!!
!!  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
!!  deallocate(coeff_tmp,stat=istat)
!!  call memocc(istat,iall,'coeff_tmp',subname)
!!
!!  if (nproc>1) then
!!      call mpiallred(ovrlp_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!  end if
!!
!!  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
!!  offdiagsum=0.0d0
!!  offdiagsum2=0.0d0
!!  do iorb=1,ksorbs%norb+num_extra!basis_orbs%norb
!!     do jorb=1,ksorbs%norb+num_extra!basis_orbs%norb
!!        if (iorb==jorb) cycle
!!        offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
!!        offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
!!     end do
!!  end do
!!  offdiagsum=offdiagsum/(basis_orbs%norb**2-basis_orbs%norb)
!!  offdiagsum2=offdiagsum2/(basis_orbs%norb**2-basis_orbs%norb)
!!  if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2
!!
!!  ! sort the states - really need just ks+extra not all, otherwise sloshing!
!!  ipiv=f_malloc(ksorbs%norb+num_extra,id='coeff_final')
!!  tmp_array=f_malloc(ksorbs%norb+num_extra,id='tmp_array')
!!  do itmb=1,ksorbs%norb+num_extra
!!     tmp_array(itmb)=-ham_coeff(itmb,itmb)
!!  end do
!!
!!  do itmb=1,ksorbs%norb+num_extra
!!     ipiv(itmb)=itmb
!!  end do
!!  !call sort_positions(ksorbs%norb+num_extra,tmp_array,ipiv)
!!  do jtmb=1,ksorbs%norb+num_extra
!!     ham_coeff(jtmb,jtmb)=-tmp_array(ipiv(jtmb))
!!  end do
!!
!!  tmp_array2=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='tmp_array2')
!!  do jtmb=1,ksorbs%norb+num_extra
!!     do itmb=1,basis_orbs%norb
!!        tmp_array2(itmb,jtmb)=coeff(itmb,jtmb)
!!     end do
!!  end do
!!
!!  do jtmb=1,ksorbs%norb+num_extra
!!     do itmb=1,basis_orbs%norb
!!        coeff(itmb,jtmb)=tmp_array2(itmb,ipiv(jtmb))
!!     end do
!!  end do
!!  call f_free(tmp_array2)
!!
!!  if (reorder) then
!!     ! diagonalize within the space of occ+extra
!!     eval=f_malloc((/ksorbs%norb+num_extra/),id='eval')
!!     ham_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ham_coeff_small')
!!     ovrlp_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ovrlp_coeff_small')
!!     ! assume ovrlp_coeff is orthgonal for now
!!     do iorb=1,ksorbs%norb+num_extra
!!        do jorb=1,ksorbs%norb+num_extra
!!           ham_coeff_small(iorb,jorb)=ham_coeff(iorb,jorb)
!!           if (iorb==jorb) then
!!              ovrlp_coeff_small(iorb,jorb)=1.0d0
!!           else
!!              ovrlp_coeff_small(iorb,jorb)=0.0d0
!!           end if
!!        end do
!!     end do
!!     ! diagonalize within the space of occ+extra
!!     call diagonalizeHamiltonian2(iproc, ksorbs%norb+num_extra, ham_coeff_small, ovrlp_coeff_small, eval)
!!     call f_free(ovrlp_coeff_small)
!!     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='coeff_tmp')
!!
!!     ! multiply new eigenvectors by coeffs
!!     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb+num_extra, ksorbs%norb+num_extra, 1.d0, coeff(1,1), &
!!          basis_orbs%norb, ham_coeff_small, ksorbs%norb+num_extra, 0.d0, coeff_tmp, basis_orbs%norb)
!! 
!!     call f_free(ham_coeff_small)
!!     call vcopy(basis_orbs%norb*(ksorbs%norb+num_extra),coeff_tmp(1,1),1,coeff(1,1),1)
!!     call f_free(coeff_tmp)
!!
!!     do iorb=1,basis_orbs%norb
!!        if (iorb<=ksorbs%norb) then
!!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
!!                ksorbs%occup(iorb),ham_coeff(iorb,iorb),eval(iorb)
!!        else if (iorb<=ksorbs%norb+num_extra) then
!!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
!!                0.0d0,ham_coeff(iorb,iorb),eval(iorb)
!!        !else
!!        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
!!        !        0.0d0,ham_coeff(iorb,iorb)
!!        end if
!!     end do
!!     call f_free(eval)
!!   else
!!      do iorb=1,basis_orbs%norb
!!        if (iorb<=ksorbs%norb) then
!!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
!!                ksorbs%occup(iorb),ham_coeff(iorb,iorb)
!!        else if (iorb<=ksorbs%norb+num_extra) then
!!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
!!                0.0d0,ham_coeff(iorb,iorb)
!!        !else
!!        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
!!        !        0.0d0,ham_coeff(iorb,iorb)
!!        end if
!!     end do
!!  end if
!!
!!  call f_free(ipiv)
!!  call f_free(tmp_array)
!!  iall=-product(shape(ham_coeff))*kind(ham_coeff)
!!  deallocate(ham_coeff,stat=istat)
!!  call memocc(istat,iall,'ham_coeff',subname)
!!
!!end subroutine reordering_coeffs


subroutine find_alpha_sd(iproc,nproc,alpha,tmb,orbs,coeffp,grad,energy0,fnrm,pred_e,taylor_order)
  use module_base
  use module_types
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
  implicit none
  integer, intent(in) :: iproc, nproc, taylor_order
  real(kind=gp), intent(inout) :: alpha
  type(DFT_wavefunction) :: tmb
  type(orbitals_data), intent(in) :: orbs
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(inout) :: coeffp
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(in) :: grad
  real(kind=gp), intent(in) :: energy0, fnrm
  real(kind=gp), intent(out) :: pred_e
  integer :: iorb, iiorb, jorb, ierr
  real(kind=gp) :: energy1, a, b, c, alpha_old
  real(kind=gp),dimension(:,:),allocatable :: coeff_tmp

  call timing(iproc,'dirmin_sdfit','ON')

  ! take an initial step to get 2nd point
  coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
  do iorb=1,orbs%norbp
     iiorb = orbs%isorb + iorb
     do jorb=1,tmb%linmat%m%nfvctr
        coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-alpha*grad(jorb,iorb)
     end do
  end do

  if(nproc > 1) then 
     call mpi_allgatherv(coeffp, tmb%linmat%m%nfvctr*orbs%norbp, mpi_double_precision, coeff_tmp, &
          tmb%linmat%m%nfvctr*orbs%norb_par(:,0), tmb%linmat%m%nfvctr*orbs%isorb_par, &
          mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(tmb%linmat%m%nfvctr*orbs%norb,coeffp(1,1),1,coeff_tmp(1,1),1)
  end if

  ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
  ! instead of twice could add some criterion to check accuracy?
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, taylor_order, tmb%orbs, &
       tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, coeff_tmp, orbs)
  !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, &
  !     tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, coeff_tmp, orbs)
  !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
  !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m,&
       tmb%linmat%kernel_, tmb%linmat%ham_, energy1,&
       coeff_tmp,orbs,tmb%orbs,.true.)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
  !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
  call f_free(coeff_tmp)

  ! find ideal alpha using both points
  alpha_old=alpha
  a=fnrm/alpha_old+(energy1-energy0)/alpha_old**2
  b=-fnrm
  c=energy0
  alpha=-0.5_gp*b/a
  ! don't update if we have found a maximum, or negative alpha is predicted
  ! do something better here - don't just want to do the same thing twice, so at least check if energy has decreased
  if (alpha<0.0_gp .or. a<0.0_gp) alpha=alpha_old
  pred_e=a*alpha**2+b*alpha+c

  !open(127)
  !write(127,*) '#',a,b,c,(energy1-energy0)/alpha_old,b-(energy1-energy0)/alpha_old
  !write(127,*) 0.0_gp,energy0
  !write(127,*) alpha_old,energy1

  call timing(iproc,'dirmin_sdfit','OF')

end subroutine find_alpha_sd


subroutine calculate_kernel_and_energy(iproc,nproc,denskern,ham,denskern_mat,ham_mat, &
           energy,coeff,orbs,tmb_orbs,calculate_kernel)
  use module_base
  use module_types
  use module_interfaces, only: calculate_density_kernel
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: extract_taskgroup_inplace, gather_matrix_from_taskgroups_inplace
  implicit none
  integer, intent(in) :: iproc, nproc
  type(sparse_matrix), intent(in) :: ham
  type(sparse_matrix), intent(in) :: denskern
  type(matrices),intent(in) :: ham_mat
  type(matrices),intent(out) :: denskern_mat
  logical, intent(in) :: calculate_kernel
  real(kind=gp), intent(out) :: energy
  type(orbitals_data), intent(in) :: orbs, tmb_orbs
  real(kind=gp), dimension(denskern%nfvctr,tmb_orbs%norb), intent(in) :: coeff

  integer :: iorb, jorb, ind_ham, ind_denskern, ierr, iorbp, is, ie, ispin

  if (calculate_kernel) then 
     !!call extract_taskgroup_inplace(denskern, denskern_mat)
     call calculate_density_kernel(iproc, nproc, .true., orbs, tmb_orbs, coeff, denskern, denskern_mat)
     !call gather_matrix_from_taskgroups_inplace(iproc, nproc, denskern, denskern_mat)
     !denskern%matrix_compr = denskern_mat%matrix_compr
  end if

  call timing(iproc,'calc_energy','ON')
  energy=0.0_gp
  do iorbp=1,tmb_orbs%norbp
     iorb=iorbp+tmb_orbs%isorb
     if (tmb_orbs%spinsgn(iorb)>0.d0) then
         ! spin up support function or non-polarized case
         is=1
         ie=tmb_orbs%norbu
         ispin=1
     else
         ! spin down support function
         is=tmb_orbs%norbu+1
         ie=tmb_orbs%norb
         ispin=2
     end if
     !$omp parallel default(private) shared(is,ie,iorb,denskern,ham,denskern_mat,ham_mat,tmb_orbs,energy,ispin)
     !$omp do reduction(+:energy)
     !do jorb=1,tmb_orbs%norb
     do jorb=is,ie
        ind_ham = matrixindex_in_compressed(ham,iorb,jorb)
        ind_denskern = matrixindex_in_compressed(denskern,jorb,iorb)
        if (ind_ham==0.or.ind_denskern==0) cycle
        energy = energy + &
            denskern_mat%matrix_compr(ind_denskern-denskern%isvctrp_tg)*ham_mat%matrix_compr(ind_ham-ham%isvctrp_tg)
            !!write(*,'(a,5i8,2es16.7)') 'iorb, jorb, ispin, ind_denskern, ind_ham, val_denskern, val_ham', &
            !!    iorb, jorb, ispin, mod(ind_denskern-denskern%isvctrp_tg-1,denskern%nvctr)+1, mod(ind_ham-ham%isvctrp_tg-1,ham%nvctr)+1, denskern_mat%matrix_compr(ind_denskern-denskern%isvctrp_tg), ham_mat%matrix_compr(ind_ham-ham%isvctrp_tg)
     end do
     !$omp end do
     !$omp end parallel
  end do
  if (nproc>1) then
     call mpiallred(energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if
  call timing(iproc,'calc_energy','OF')

end subroutine calculate_kernel_and_energy


!> calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
!! then grad=S^-1grad_cov
subroutine calculate_coeff_gradient(iproc,nproc,tmb,order_taylor,max_inversion_error,KSorbs,grad_cov,grad)
  use module_base
  use module_types
  use module_interfaces, only: calculate_density_kernel
  use sparsematrix_base, only: matrices, sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: extract_taskgroup_inplace, gather_matrix_from_taskgroups_inplace
  use matrix_operations, only: overlapPowerGeneral, check_taylor_order
  use parallel_linalg, only: dgesv_parallel
  implicit none

  integer, intent(in) :: iproc, nproc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(DFT_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: KSorbs
  real(gp), dimension(tmb%linmat%m%nfvctr,KSorbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp

  integer :: iorb, iiorb, info, ierr, ispin
  integer, dimension(1) :: power
  real(gp),dimension(:,:,:),allocatable :: sk, skh, skhp
  real(gp),dimension(:,:),pointer :: inv_ovrlp
  real(kind=gp), dimension(:,:), allocatable:: grad_full
  character(len=*),parameter:: subname='calculate_coeff_gradient'
  type(matrices),dimension(1) :: inv_ovrlp_

  integer :: itmp, itrials, jorb
  integer :: ncount1, ncount_rate, ncount_max, ncount2
  real(kind=4) :: tr0, tr1
  real(kind=8) :: deviation, time, max_error, mean_error, maxerror


  call f_routine(id='calculate_coeff_gradient')
  call timing(iproc,'dirmin_lagmat1','ON')

  inv_ovrlp_(1) = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_(1))


  ! we have the kernel already, but need it to not contain occupations so recalculate here
  ! don't want to lose information in the compress/uncompress process - ideally need to change sparsity pattern of kernel
  !call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')

  !!do iorb=1,KSorbs%norbp
  !!    iiorb=KSorbs%isorb+iorb
  !!    if (KSorbs%spinsgn(iiorb)>0.d0) then
  !!        ispin=1
  !!    else
  !!        ispin=2
  !!    end if
  !!    do jorb=1,tmb%linmat%m%nfvctr
  !!        write(4900+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
  !!    end do
  !!end do

  !!do iiorb=1,KSorbs%norb
  !!    if (KSorbs%spinsgn(iiorb)>0.d0) then
  !!        ispin=1
  !!    else
  !!        ispin=2
  !!    end if
  !!    do jorb=1,tmb%linmat%m%nfvctr
  !!        write(4800+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
  !!    end do
  !!end do

 !call uncompress_matrix(iproc,tmb%linmat%denskern)
  !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
  call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, &
       tmb%coeff, tmb%linmat%l, tmb%linmat%kernel_, .true.)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)

  sk=f_malloc0((/tmb%linmat%m%nfvctrp,tmb%linmat%m%nfvctr,tmb%linmat%m%nspin/), id='sk')

  ! calculate I-S*K - first set sk to identity
  do ispin=1,tmb%linmat%m%nspin
      do iorb=1,tmb%linmat%m%nfvctrp
          iiorb=mod(tmb%linmat%m%isfvctr+iorb-1,tmb%linmat%m%nfvctr)+1
          sk(iorb,iiorb,ispin) = 1.d0
      end do 
  end do

  if (tmb%linmat%m%nfvctrp>0) then
      do ispin=1,tmb%linmat%m%nspin
          call dgemm('t', 'n', tmb%linmat%m%nfvctrp, tmb%linmat%m%nfvctr, tmb%linmat%m%nfvctr, -1.d0, &
               tmb%linmat%ovrlp_%matrix(1,tmb%linmat%m%isfvctr+1,ispin), tmb%linmat%m%nfvctr, &
               tmb%linmat%kernel_%matrix(1,1,ispin), tmb%linmat%m%nfvctr, 1.d0, sk(1,1,ispin), tmb%linmat%m%nfvctrp)
      end do
  end if
  !!do ispin=1,tmb%linmat%m%nspin
  !!    do iorb=1,tmb%linmat%m%nfvctr
  !!        do jorb=1,tmb%linmat%m%nfvctr
  !!            write(6300+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, tmb%linmat%ovrlp_%matrix(jorb,iorb,ispin)', iorb, jorb, tmb%linmat%ovrlp_%matrix(jorb,iorb,ispin)
  !!            write(6400+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, tmb%linmat%kernel_%matrix(jorb,iorb,ispin)', iorb, jorb, tmb%linmat%kernel_%matrix(jorb,iorb,ispin)
  !!        end do
  !!    end do
  !!end do
  !!do ispin=1,tmb%linmat%m%nspin
  !!    do iorb=1,tmb%linmat%m%nfvctr
  !!        do jorb=1,tmb%linmat%m%nfvctrp
  !!            write(6100+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, sk(jorb,iorb,ispin)', iorb, jorb, sk(jorb,iorb,ispin)
  !!        end do
  !!    end do
  !!end do

  ! coeffs and therefore kernel will change, so no need to keep it
  call f_free_ptr(tmb%linmat%kernel_%matrix)

  skhp=f_malloc((/tmb%linmat%m%nfvctr,max(tmb%linmat%m%nfvctrp,1),tmb%linmat%m%nspin/), id='skhp')

  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
  if (tmb%linmat%m%nfvctrp>0) then
      do ispin=1,tmb%linmat%m%nspin
          call dgemm('t', 't', tmb%linmat%m%nfvctr, tmb%linmat%m%nfvctrp, tmb%linmat%m%nfvctr, &
               1.d0, tmb%linmat%ham_%matrix(1,1,ispin), tmb%linmat%m%nfvctr, &
               sk(1,1,ispin), tmb%linmat%m%nfvctrp, 0.d0, skhp(1,1,ispin), tmb%linmat%m%nfvctr)
      end do
  end if
  !!do ispin=1,tmb%linmat%m%nspin
  !!    do iorb=1,tmb%linmat%m%nfvctrp
  !!        do jorb=1,tmb%linmat%m%nfvctr
  !!            write(6200+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, skhp(jorb,iorb,ispin)', iorb, jorb, skhp(jorb,iorb,ispin)
  !!        end do
  !!    end do
  !!end do

  call f_free(sk)

  skh=f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr,tmb%linmat%m%nspin/), id='skh')

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_lagmat2','ON')

  ! gather together
  if(nproc > 1) then
      do ispin=1,tmb%linmat%m%nspin
         call mpi_allgatherv(skhp(1,1,ispin), tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctrp, &
            mpi_double_precision, skh(1,1,ispin), tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctr_par, &
            tmb%linmat%m%nfvctr*tmb%linmat%m%isfvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      end do
  else
     call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctrp*tmb%linmat%m%nspin,skhp(1,1,1),1,skh(1,1,1),1)
  end if

  !!do ispin=1,tmb%linmat%m%nspin
  !!    do iorb=1,tmb%linmat%m%nfvctr
  !!        do jorb=1,tmb%linmat%m%nfvctr
  !!            write(6000+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, skh(jorb,iorb,ispin)', iorb, jorb, skh(jorb,iorb,ispin)
  !!        end do
  !!    end do
  !!end do

  call timing(iproc,'dirmin_lagmat2','OF')
  call timing(iproc,'dirmin_lagmat1','ON')

  call f_free(skhp)

  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
  if (KSorbs%norbp>0) then
     !call dgemm('t', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
     !     tmb%orbs%norb, tmb%coeff(1,KSorbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,1), tmb%orbs%norb)
     do iorb=1,KSorbs%norbp
         iiorb=KSorbs%isorb+iorb
         if (KSorbs%spinsgn(iiorb)>0.d0) then
             ispin=1
         else
             ispin=2
         end if
         call dgemm('t', 'n', tmb%linmat%m%nfvctr, 1, tmb%linmat%m%nfvctr, 1.d0, skh(1,1,ispin), &
              tmb%linmat%m%nfvctr, tmb%coeff(1,KSorbs%isorb+iorb), tmb%linmat%m%nfvctr, &
              0.d0, grad_cov(1,iorb), tmb%linmat%m%nfvctr)
     end do
  end if

  call f_free(skh)

  ! multiply by f_i to get grad_i^a
  do iorb=1,KSorbs%norbp
     iiorb=KSorbs%isorb+iorb
     grad_cov(:,iorb)=grad_cov(:,iorb)*KSorbs%occup(iiorb)
  end do

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_dgesv','ON')

  ! Solve the linear system ovrlp*grad=grad_cov
  if(tmb%orthpar%blocksize_pdsyev<0) then
     call timing(iproc,'dirmin_dgesv','OF')
     inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     power(1)=1
     call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
          order_taylor, 1, power, -8, &
          imode=2, &
          ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, check_accur=.true., &
          max_error=max_error, mean_error=mean_error)
     call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)

     !!!DEBUG checking S^-1 etc.
     !!!test dense version of S^-1
     !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     !!if (iproc==0)write(*,*) ''
     !!do itmp=0,20
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, itmp, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
     !!           inv_ovrlp, error, tmb%orbs)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do
     !!call f_free_ptr(inv_ovrlp)

     !!! test S^+/-1/2
     !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     !!if (iproc==0)write(*,*) ''
     !!do itmp=0,20
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, itmp, 2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
     !!           inv_ovrlp, error, tmb%orbs)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do

     !!if (iproc==0)write(*,*) ''
     !!do itmp=0,20
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, itmp, -2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
     !!           inv_ovrlp, error, tmb%orbs)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do

     !!!test sparse version of S^-1
     !!deallocate(tmb%linmat%ovrlp%matrix)

     !!if (iproc==0)write(*,*) ''
     !!call f_free_ptr(inv_ovrlp)
     !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
     !!do itmp=0,2
     !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, 1, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
     !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do
     !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)

     !!!test sparse version of S^1/2
     !!if (iproc==0)write(*,*) ''
     !!nullify(inv_ovrlp)
     !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
     !!do itmp=0,2
     !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, 1, 2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
     !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do
     !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)

     !!!test sparse version of S^-1/2
     !!if (iproc==0)write(*,*) ''
     !!nullify(inv_ovrlp)
     !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
     !!do itmp=0,2
     !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr0)
     !!   call system_clock(ncount1,ncount_rate,ncount_max)
     !!   maxerror=0.0d0
     !!   do itrials=1,50
     !!      call overlapPowerGeneral(iproc, nproc, 1, -2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
     !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
     !!      maxerror=max(maxerror,error)
     !!   end do
     !!   call mpi_barrier(mpi_comm_world,ierr)
     !!   call cpu_time(tr1)
     !!   call system_clock(ncount2,ncount_rate,ncount_max)
     !!   if (iproc==0) then
     !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
     !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
     !!   end if
     !!end do
     !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)

     !!call mpi_finalize(bigdft_mpi%mpi_comm)
     !!stop
     !!END DEBUG

     call timing(iproc,'dirmin_dgesv','ON')
     if (KSorbs%norbp>0) then
         !call dgemm('n', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, inv_ovrlp_%matrix(1,1,1), &
         !     tmb%orbs%norb, grad_cov(1,1), tmb%orbs%norb, 0.d0, grad(1,1), tmb%orbs%norb)
         do iorb=1,KSorbs%norbp
             iiorb=KSorbs%isorb+iorb
             if (KSorbs%spinsgn(iiorb)>0.d0) then
                 ispin=1
             else
                 ispin=2
             end if
             call dgemm('n', 'n', tmb%linmat%m%nfvctr, 1, tmb%linmat%m%nfvctr, 1.d0, inv_ovrlp_(1)%matrix(1,1,ispin), &
                  tmb%linmat%m%nfvctr, grad_cov(1,iorb), tmb%linmat%m%nfvctr, 0.d0, grad(1,iorb), tmb%linmat%m%nfvctr)
         end do
     end if
     call f_free_ptr(inv_ovrlp)
  else
     info = 0 ! needed for when some processors have orbs%norbp=0
     grad_full=f_malloc((/tmb%linmat%m%nfvctr,KSorbs%norb/),id='grad_full')
     ! do allgather instead of allred so we can keep grad as per proc
     if(nproc > 1) then 
        call mpi_allgatherv(grad_cov, tmb%linmat%m%nfvctr*KSorbs%norbp, mpi_double_precision, grad_full, &
             tmb%linmat%m%nfvctr*KSorbs%norb_par(:,0), tmb%linmat%m%nfvctr*KSorbs%isorb_par, &
             mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call vcopy(tmb%linmat%m%nfvctr*KSorbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
     end if

     call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
          tmb%linmat%m%nfvctr, KSorbs%norb, tmb%linmat%ovrlp_%matrix, tmb%linmat%m%nfvctr, &
          grad_full, tmb%linmat%m%nfvctr, info)
     call vcopy(tmb%linmat%m%nfvctr*KSorbs%norbp,grad_full(1,KSorbs%isorb+1),1,grad(1,1),1)

     call f_free(grad_full)
     if(info/=0) then
        write(*,'(a,i0)') 'ERROR in dgesv: info=',info
        stop
     end if
  end if

  call deallocate_matrices(inv_ovrlp_(1))

  call timing(iproc,'dirmin_dgesv','OF')
  call f_release_routine()

end subroutine calculate_coeff_gradient


!!!> calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
!!!! then grad=S^-1grad_cov
!!subroutine calculate_coeff_gradient_extra(iproc,nproc,num_extra,tmb,order_taylor,max_inversion_error,KSorbs,grad_cov,grad)
!!  use module_base
!!  use module_types
!!  use module_interfaces, only: calculate_density_kernel
!!  use sparsematrix_base, only: matrices, sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
!!                               matrices_null, allocate_matrices, deallocate_matrices
!!  use sparsematrix, only: extract_taskgroup_inplace, gather_matrix_from_taskgroups_inplace
!!  use matrix_operations, only: overlapPowerGeneral, check_taylor_order
!!  use parallel_linalg, only: dgesv_parallel
!!  implicit none
!!
!!  integer, intent(in) :: iproc, nproc, num_extra
!!  integer,intent(inout) :: order_taylor
!!  real(kind=8),intent(in) :: max_inversion_error
!!  type(DFT_wavefunction), intent(inout) :: tmb
!!  type(orbitals_data), intent(in) :: KSorbs
!!  real(gp), dimension(tmb%linmat%m%nfvctr,tmb%orbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp
!!
!!  integer :: iorb, iiorb, info, ierr, ispin
!!  integer, dimension(1) :: power
!!  real(gp),dimension(:,:,:),allocatable :: sk, skhp, skh
!!  real(gp),dimension(:,:),pointer ::  inv_ovrlp
!!  real(kind=gp), dimension(:), allocatable:: occup_tmp
!!  real(kind=gp), dimension(:,:), allocatable:: grad_full
!!  real(kind=gp) :: max_error, mean_error
!!  type(matrices),dimension(1) :: inv_ovrlp_
!!
!!  !!if (tmb%linmat%m%nspin==2) stop 'ERROR: calculate_coeff_gradient_extra not yet implemented for npsin==2'
!!
!!  call f_routine(id='calculate_coeff_gradient_extra')
!!  call timing(iproc,'dirmin_lagmat1','ON')
!!
!!  inv_ovrlp_(1) = matrices_null()
!!  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_(1))
!!
!!
!!  occup_tmp=f_malloc(tmb%orbs%norb,id='occup_tmp')
!!  call vcopy(tmb%orbs%norb,tmb%orbs%occup(1),1,occup_tmp(1),1)
!!
!!  call f_zero(tmb%orbs%norb,tmb%orbs%occup(1))
!!  do iorb=1,KSorbs%norb+num_extra
!!     tmb%orbs%occup(iorb)=1.0d0
!!  end do
!!
!!  ! we have the kernel already, but need it to not contain occupations so recalculate here
!!  !call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
!!  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
!!  !call uncompress_matrix(iproc,tmb%linmat%denskern)
!!  !!call calculate_density_kernel_uncompressed (iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%kernel_%matrix)
!!  !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
!!  call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, &
!!       tmb%coeff, tmb%linmat%l, tmb%linmat%kernel_, .true.)
!!  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
!!
!!
!!  call vcopy(tmb%orbs%norb,occup_tmp(1),1,tmb%orbs%occup(1),1)
!!  call f_free(occup_tmp)
!!
!!  sk=f_malloc0((/tmb%linmat%l%nfvctrp,tmb%linmat%l%nfvctr,tmb%linmat%l%nspin/), id='sk')
!!
!!  ! calculate I-S*K - first set sk to identity
!!  do ispin=1,tmb%linmat%l%nspin
!!     do iorb=1,tmb%linmat%l%nfvctrp
!!        iiorb=mod(tmb%linmat%l%isfvctr+iorb-1,tmb%linmat%l%nfvctr)+1
!!        sk(iorb,iiorb,ispin) = 1.d0
!!     end do 
!!  end do
!!
!!
!!  if (tmb%linmat%l%nfvctrp>0) then
!!     do ispin=1,tmb%linmat%l%nspin
!!        call dgemm('t', 'n', tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctr, -1.d0, &
!!             tmb%linmat%ovrlp_%matrix(1,tmb%linmat%l%isfvctr+1,ispin), tmb%linmat%l%nfvctr, &
!!             tmb%linmat%kernel_%matrix(1,1,ispin), tmb%linmat%l%nfvctr, 1.d0, sk(1,1,ispin), tmb%linmat%l%nfvctrp)
!!      end do
!!  end if
!!
!!
!!  ! coeffs and therefore kernel will change, so no need to keep it
!!  call f_free_ptr(tmb%linmat%kernel_%matrix)
!!
!!  skhp=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctrp,tmb%linmat%l%nspin/), id='skhp')
!!
!!  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
!!  if (tmb%linmat%l%nfvctrp>0) then
!!     do ispin=1,tmb%linmat%l%nspin
!!        call dgemm('t', 't', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
!!             1.d0, tmb%linmat%ham_%matrix(1,1,ispin), &
!!             tmb%linmat%l%nfvctr, sk(1,1,ispin), tmb%linmat%l%nfvctrp, 0.d0, skhp(1,1,ispin), tmb%linmat%l%nfvctr)
!!      end do
!!  end if
!!
!!
!!  call f_free(sk)
!!
!!  skh=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr,tmb%linmat%l%nspin/), id='skh')
!!
!!  call timing(iproc,'dirmin_lagmat1','OF')
!!  call timing(iproc,'dirmin_lagmat2','ON')
!!
!!  ! gather together
!!  if(nproc > 1) then
!!     do ispin=1,tmb%linmat%l%nspin
!!        call mpi_allgatherv(skhp(1,1,ispin), tmb%linmat%l%nfvctr*tmb%linmat%l%nfvctrp, &
!!            mpi_double_precision, skh(1,1,ispin), &
!!            tmb%linmat%l%nfvctr*tmb%linmat%l%nfvctr_par(:), &
!!            tmb%linmat%l%nfvctr*tmb%linmat%l%isfvctr_par, &
!!            mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!     end do
!!  else
!!     call vcopy(tmb%linmat%l%nfvctrp*tmb%linmat%l%nfvctr*tmb%linmat%l%nspin,skhp(1,1,1),1,skh(1,1,1),1)
!!  end if
!!
!!
!!  call timing(iproc,'dirmin_lagmat2','OF')
!!  call timing(iproc,'dirmin_lagmat1','ON')
!!
!!  call f_free(skhp)
!!
!!  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
!!  if (tmb%linmat%l%nfvctrp>0) then
!!     do iorb=1,tmb%linmat%l%nfvctrp
!!        iiorb=tmb%linmat%l%isfvctr+iorb
!!        if (tmb%orbs%spinsgn(iiorb)>0.d0) then
!!            ispin=1
!!        else
!!            ispin=2
!!        end if
!!        call dgemm('t', 'n', tmb%linmat%l%nfvctr, 1, tmb%linmat%l%nfvctr, 1.d0, skh(1,1,ispin), &
!!             tmb%linmat%l%nfvctr, tmb%coeff(1,tmb%linmat%l%isfvctr+iorb), tmb%linmat%l%nfvctr, &
!!             0.d0, grad_cov(1,iorb), tmb%linmat%l%nfvctr)
!!      end do
!!  end if
!!
!!  call f_free(skh)
!!
!!  ! multiply by f_i to get grad_i^a
!!  do iorb=1,tmb%orbs%norbp
!!     iiorb=tmb%orbs%isorb+iorb
!!     grad_cov(:,iorb)=grad_cov(:,iorb)*tmb%orbs%occup(iiorb)
!!  end do
!!
!!  call timing(iproc,'dirmin_lagmat1','OF')
!!  call timing(iproc,'dirmin_dgesv','ON')
!!
!!
!!  info = 0 ! needed for when some processors have orbs%norbp=0
!!  ! Solve the linear system ovrlp*grad=grad_cov
!!  if(tmb%orthpar%blocksize_pdsyev<0) then
!!     !! keep the covariant gradient to calculate fnrm correctly
!!     !call vcopy(tmb%orbs%norb*tmb%orbs%norbp,grad_cov,1,grad,1)
!!     !if (tmb%orbs%norbp>0) then
!!     !   ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
!!     !   call dgesv(tmb%orbs%norb, tmb%orbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
!!     !        grad(1,1), tmb%orbs%norb, info)
!!     !   call f_free(ipiv)
!!     !end if
!!     !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
!!     power(1)=1
!!     call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
!!          order_taylor, 1, power, -8, &
!!          imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
!!          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, check_accur=.true., &
!!          max_error=max_error, mean_error=mean_error)
!!     call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
!!
!!     if (tmb%orbs%norbp>0) then
!!        call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%orbs%norbp, tmb%linmat%l%nfvctr, 1.d0, inv_ovrlp_(1)%matrix(1,1,1), &
!!             tmb%linmat%l%nfvctr, grad_cov(1,1), tmb%linmat%l%nfvctr, 0.d0, grad(1,1), tmb%linmat%l%nfvctr)
!!     end if
!!     !!call f_free_ptr(inv_ovrlp)
!!  else
!!      grad_full=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_full')
!!      ! do allgather instead of allred so we can keep grad as per proc
!!      if(nproc > 1) then 
!!         call mpi_allgatherv(grad_cov, tmb%linmat%l%nfvctr*tmb%orbs%norbp, &
!!              mpi_double_precision, grad_full, &
!!              tmb%linmat%l%nfvctr*tmb%orbs%norb_par(:,0), &
!!              tmb%linmat%l%nfvctr*tmb%orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!      else
!!         call vcopy(tmb%linmat%l%nfvctr*tmb%orbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
!!      end if
!!      !call mpiallred(grad(1,1), tmb%orbs%norb*tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!
!!      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
!!           tmb%orbs%norb, tmb%orbs%norb, tmb%linmat%ovrlp_%matrix, tmb%orbs%norb, grad_full, tmb%orbs%norb, info)
!!
!!      call vcopy(tmb%linmat%l%nfvctr*tmb%orbs%norbp,grad_full(1,tmb%orbs%isorb+1),1,grad(1,1),1)
!!
!!      call f_free(grad_full)
!!
!!
!!  end if
!!
!!  if(info/=0) then
!!      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
!!      stop
!!  end if
!!
!!  call deallocate_matrices(inv_ovrlp_(1))
!!
!!  call timing(iproc,'dirmin_dgesv','OF')
!!  call f_release_routine()
!!
!!end subroutine calculate_coeff_gradient_extra


subroutine precondition_gradient_coeff(ntmb, norb, ham, ovrlp, grad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: ntmb, norb
  real(8),dimension(ntmb,ntmb),intent(in):: ham, ovrlp
  real(8),dimension(ntmb,norb),intent(inout):: grad
  
  ! Local variables
  integer:: iorb, itmb, jtmb, info, istat, iall
  complex(8),dimension(:,:),allocatable:: mat
  complex(8),dimension(:,:),allocatable:: rhs
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='precondition_gradient_coeff'
  
  mat = f_malloc((/ ntmb, ntmb /),id='mat')
  rhs = f_malloc((/ ntmb, norb /),id='rhs')
  
  ! Build the matrix to be inverted
  do itmb=1,ntmb
      do jtmb=1,ntmb
          mat(jtmb,itmb) = cmplx(ham(jtmb,itmb)+.5d0*ovrlp(jtmb,itmb),0.d0,kind=8)
      end do
      mat(itmb,itmb)=mat(itmb,itmb)+cmplx(0.d0,-1.d-1,kind=8)
      !mat(itmb,itmb)=mat(itmb,itmb)-cprec
  end do
  do iorb=1,norb
      do itmb=1,ntmb
          rhs(itmb,iorb)=cmplx(grad(itmb,iorb),0.d0,kind=8)
      end do
  end do
  
  ipiv = f_malloc(ntmb,id='ipiv')
  
  call zgesv(ntmb, norb, mat(1,1), ntmb, ipiv, rhs(1,1), ntmb, info)
  if(info/=0) then
      stop 'ERROR in zgesv'
  end if
  !call vcopy(nel, rhs(1), 1, grad(1), 1)
  do iorb=1,norb
      do itmb=1,ntmb
          grad(itmb,iorb)=real(rhs(itmb,iorb))
      end do
  end do
  
  call f_free(ipiv)
  call f_free(mat)
  call f_free(rhs)

end subroutine precondition_gradient_coeff



subroutine initialize_DIIS_coeff(isx, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: isx
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  character(len=*),parameter:: subname='initialize_DIIS_coeff'
    
  ldiis%isx=isx
  ldiis%is=0
  ldiis%switchSD=.false.
  ldiis%trmin=1.d100
  ldiis%trold=1.d100
  ldiis%alpha_coeff=0.1d0

end subroutine initialize_DIIS_coeff


subroutine allocate_DIIS_coeff(tmb, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(in):: tmb
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  integer:: ii, istat
  character(len=*),parameter:: subname='allocate_DIIS_coeff'

  ldiis%mat = f_malloc_ptr((/ldiis%isx,ldiis%isx,tmb%orbs%norbp/),id='ldiis%mat')

  ii=ldiis%isx*tmb%orbs%norb*tmb%orbs%norbp
  ldiis%phiHist = f_malloc_ptr(ii,id='ldiis%phiHist')
  ldiis%hphiHist = f_malloc_ptr(ii,id='ldiis%hphiHist')

end subroutine allocate_DIIS_coeff
