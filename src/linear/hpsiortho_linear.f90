!> @file
!! H|psi> and orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, kernel_compr, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
           energy_increased, tmb, lhphiold, tmblarge, overlap_calculated, &
           energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
           hpsi_noprecond)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(DFT_wavefunction),target,intent(inout):: tmblarge, tmb
  real(8),dimension(tmblarge%mad%nvctr),target,intent(in) :: kernel_compr
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(8),dimension(tmb%orbs%norb),intent(inout) :: fnrmOldArr
  real(8),dimension(tmb%orbs%norbp),intent(inout) :: alpha
  real(8),intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(8),intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphiold
  logical,intent(inout):: overlap_calculated
  type(energy_terms),intent(in) :: energs
  real(8),dimension(:),pointer:: hpsit_c, hpsit_f
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
  real(kind=8),dimension(tmb%orbs%npsidim_orbs),optional,intent(out) :: hpsi_noprecond

  ! Local variables
  integer :: iorb, jorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, iseg, jjorb, i
  real(kind=8) :: ddot, tt, eval_zero
  character(len=*),parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8),dimension(:),pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8),dimension(:,:),allocatable :: fnrmOvrlpArr, fnrmArr
  ! real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
  ! real(kind=8),dimension(:,:),allocatable :: gnrmArr
  real(kind=8),dimension(:),allocatable :: lagmat_compr, hpsi_conf, hpsi_tmp
  real(kind=8),dimension(:),pointer :: kernel_compr_tmp


  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      allocate(hpsi_conf(tmb%orbs%npsidim_orbs), stat=istat)
      call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
      call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmblarge%hpsi, hpsi_conf)
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
          call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
          ist=ist+ncount
      end do
  end if

  allocate(fnrmOvrlpArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! by default no quick exit
  energy_increased=.false.


  allocate(hpsittmp_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
  allocate(hpsittmp_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)

  if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     target_function==TARGET_FUNCTION_IS_HYBRID) then
      if(sum(tmblarge%collcom%nrecvcounts_c)>0) &
          call dcopy(sum(tmblarge%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmblarge%collcom%nrecvcounts_f)>0) &
          call dcopy(7*sum(tmblarge%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          allocate(kernel_compr_tmp(tmblarge%mad%nvctr), stat=istat)
          call memocc(istat, kernel_compr_tmp, 'kernel_compr_tmp', subname)
          call vcopy(tmblarge%mad%nvctr, kernel_compr(1), 1, kernel_compr_tmp(1), 1)
          ii=0
          do iseg=1,tmblarge%mad%nseg
              do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
                  ii=ii+1
                  iiorb = (jorb-1)/tmblarge%orbs%norb + 1
                  jjorb = jorb - (iiorb-1)*tmblarge%orbs%norb
                  if(iiorb==jjorb) then
                      kernel_compr_tmp(ii)=0.d0
                  else
                      kernel_compr_tmp(ii)=kernel_compr(ii)
                  end if
              end do
          end do
      else
          kernel_compr_tmp => kernel_compr
      end if

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then

          ist=1
          do iorb=tmblarge%orbs%isorb+1,tmblarge%orbs%isorb+tmblarge%orbs%norbp
              ilr=tmblarge%orbs%inwhichlocreg(iorb)
              ii=0
              do iseg=1,tmblarge%mad%nseg
                  do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/tmblarge%orbs%norb + 1
                      jjorb = jorb - (iiorb-1)*tmblarge%orbs%norb
                      if(iiorb==jjorb .and. iiorb==iorb) then
                          ncount=tmblarge%lzd%llr(ilr)%wfd%nvctr_c+7*tmblarge%lzd%llr(ilr)%wfd%nvctr_f
                          call dscal(ncount, kernel_compr(ii), tmblarge%hpsi(ist), 1)
                          ist=ist+ncount
                      end if
                  end do
              end do
          end do
          call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, tmblarge%hpsi, hpsit_c, hpsit_f, tmblarge%lzd)
          call build_linear_combination_transposed(tmblarge%orbs%norb, kernel_compr_tmp, tmblarge%collcom, &
               tmblarge%mad, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
          iall=-product(shape(kernel_compr_tmp))*kind(kernel_compr_tmp)
          deallocate(kernel_compr_tmp, stat=istat)
          call memocc(istat, iall, 'kernel_compr_tmp', subname)
      else

          call build_linear_combination_transposed(tmblarge%orbs%norb, kernel_compr_tmp, tmblarge%collcom, &
               tmblarge%mad, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)

      end if

  end if



  allocate(lagmat_compr(tmblarge%mad%nvctr), stat=istat)
  call memocc(istat, lagmat_compr, 'lagmat_compr', subname)

  call orthoconstraintNonorthogonal(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%mad, &
       tmblarge%collcom, tmblarge%orthpar, correction_orthoconstraint, tmblarge%psi, tmblarge%hpsi, &
       lagmat_compr, tmblarge%psit_c, tmblarge%psit_f, hpsit_c, hpsit_f, tmblarge%can_use_transposed, overlap_calculated)

  call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmblarge%hpsi, tmb%hpsi)


  if (present(hpsi_noprecond)) call dcopy(tmb%orbs%npsidim_orbs, tmb%hpsi, 1, hpsi_noprecond, 1)

  ! Calculate trace (or band structure energy, resp.)
  trH=0.d0
  ii=0
  do iseg=1,tmblarge%mad%nseg
      do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
          iiorb = (jorb-1)/tmb%orbs%norb + 1
          jjorb = jorb - (iiorb-1)*tmblarge%orbs%norb
          ii=ii+1
          if (iiorb==jjorb) then
              trH = trH + lagmat_compr(ii)
          end if
      end do
  end do


  iall=-product(shape(lagmat_compr))*kind(lagmat_compr)
  deallocate(lagmat_compr, stat=istat)
  call memocc(istat, iall, 'lagmat_compr', subname)


  ! trH is now the total energy (name is misleading, correct this)
  ! Multiply by 2 because when minimizing trace we don't have kernel
  ! if(iproc==0)print *,'trH,energs',trH,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
  if(tmb%orbs%nspin==1 .and. target_function/= TARGET_FUNCTION_IS_ENERGY) trH=2.d0*trH
  trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp


  ! Cycle if the trace increased (steepest descent only)
  if(.not. ldiis%switchSD .and. ldiis%isx==0) then
      if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
          if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
              'WARNING: the target function is larger than its minimal value reached so far:',trH,' > ', ldiis%trmin
          if(iproc==0) write(*,'(1x,a)') 'Decrease step size and restart with previous TMBs'
          energy_increased=.true.
      end if
  end if

  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, tmb%hpsi(ist), 1, lhphiold(ist), 1)
      fnrmArr(iorb,1)=ddot(ncount, tmb%hpsi(ist), 1, tmb%hpsi(ist), 1)
      fnrmOldArr(iorb)=ddot(ncount, lhphiold(ist), 1, lhphiold(ist), 1)
      ist=ist+ncount
  end do


  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.
  fnrm=0.d0
  fnrmMax=0.d0
  do iorb=1,tmb%orbs%norbp
      fnrm=fnrm+fnrmArr(iorb,1)
      if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
      if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      ! Adapt step size for the steepest descent minimization.
          tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
          if(tt>.6d0 .and. trH<trHold) then
              alpha(iorb)=alpha(iorb)*1.1d0
          else
              alpha(iorb)=alpha(iorb)*.6d0
          end if
          !!alpha(iorb)=min(alpha(iorb),1.5d0)
      end if
  end do
  call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)
  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmb%orbs%npsidim_orbs, tmb%hpsi, 1, lhphiold, 1)

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)',advance='no') 'Preconditioning... '
  end if
 

  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      allocate(hpsi_tmp(tmb%orbs%npsidim_orbs), stat=istat)
      call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
  end if
  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr = tmb%orbs%inWhichLocreg(iiorb)
      ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(target_function==TARGET_FUNCTION_IS_HYBRID) then
          tt=ddot(ncnt, hpsi_conf(ist), 1, tmb%hpsi(ist), 1)
          tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
          do i=ist,ist+ncnt-1
              hpsi_tmp(i)=tt*hpsi_conf(i)
          end do
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               nit_precond, hpsi_tmp(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
          call daxpy(ncnt, -tt, hpsi_conf(ist), 1, tmb%hpsi(ist), 1)
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               nit_precond, tmb%hpsi(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               0.d0, iorb, eval_zero)
          call daxpy(ncnt, 1.d0, hpsi_tmp(ist), 1, tmb%hpsi(ist), 1)
      else
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               nit_precond, tmb%hpsi(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
      end if
      ist=ist+ncnt
  end do


  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      iall=-product(shape(hpsi_conf))*kind(hpsi_conf)
      deallocate(hpsi_conf, stat=istat)
      call memocc(istat, iall, 'hpsi_conf', subname)
      iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
      deallocate(hpsi_tmp, stat=istat)
      call memocc(istat, iall, 'hpsi_tmp', subname)
  end if


  if(iproc==0) then
      write(*,'(a)') 'done.'
  end if


  call timing(iproc,'eglincomms','ON')
  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  alpha_mean=tt/dble(tmb%orbs%norb)
  alpha_max=maxval(alpha)
  call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  call timing(iproc,'eglincomms','OF')

  iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
  deallocate(fnrmOvrlpArr, stat=istat)
  call memocc(istat, iall, 'fnrmOvrlpArr', subname)

  iall=-product(shape(fnrmArr))*kind(fnrmArr)
  deallocate(fnrmArr, stat=istat)
  call memocc(istat, iall, 'fnrmArr', subname)

  iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
  deallocate(hpsittmp_c, stat=istat)
  call memocc(istat, iall, 'hpsittmp_c', subname)
  iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
  deallocate(hpsittmp_f, stat=istat)
  call memocc(istat, iall, 'hpsittmp_f', subname)


end subroutine calculate_energy_and_gradient_linear



subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, tmblarge, &
           lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, psidiff)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb, tmblarge
  real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: lphiold
  real(kind=8),intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8),dimension(tmb%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  real(kind=8),dimension(tmb%orbs%npsidim_orbs),optional,intent(out) :: psidiff
  
  ! Local variables
  integer :: istat, iall, i
  character(len=*),parameter :: subname='hpsitopsi_linear'


  call DIISorSD(iproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold)
  if(iproc==0) then
      if(ldiis%isx>0) then
          write(*,'(1x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
              ldiis%icountDIISFailureCons, ', total failures=', ldiis%icountDIISFailureTot
      else
          write(*,'(1x,2(a,es9.3),a,i0)') 'SD informations: mean alpha=', alpha_mean, ', max alpha=', alpha_max,&
          ', consecutive successes=', ldiis%icountSDSatur
      end if
  end if

  ! Improve the orbitals, depending on the choice made above.
  if (present(psidiff)) call dcopy(tmb%orbs%npsidim_orbs, tmb%psi, 1, psidiff, 1)
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, tmb, ldiis, alpha)
  else
      if(iproc==0) write(*,'(1x,a)') 'no improvement of the orbitals, recalculate gradient'
  end if
  if (present(psidiff)) then
      do i=1,tmb%orbs%npsidim_orbs
          psidiff(i)=tmb%psi(i)-psidiff(i)
      end do 
  end if

  ! The transposed quantities can now not be used any more...
  if(tmb%can_use_transposed) then
      iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
      deallocate(tmb%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%psit_c', subname)
      iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
      deallocate(tmb%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%psit_f', subname)
      tmb%can_use_transposed=.false.
  end if

  if(.not.ldiis%switchSD) then
      if(iproc==0) then
           write(*,'(1x,a)',advance='no') 'Orthonormalization... '
      end if
      ! Give tmblarge%mad since this is the correct matrix description
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orbs, tmb%lzd, &
           tmblarge%mad, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear

