!> @file
!! H|psi> and orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, kernel, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, &
           fnrmMax, alpha_mean, alpha_max, energy_increased, tmb, lhphi, lhphiold, &
           tmblarge, lhphilarge, overlap_calculated, ovrlp, energs, hpsit_c, hpsit_f)
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(DFT_wavefunction),target,intent(inout):: tmblarge, tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in) :: kernel
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(8),dimension(tmb%orbs%norb),intent(inout) :: fnrmOldArr
  real(8),dimension(tmb%orbs%norbp),intent(inout) :: alpha
  real(8),intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(8),intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(8),dimension(tmblarge%orbs%npsidim_orbs),intent(inout):: lhphilarge
  real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphi, lhphiold
  logical,intent(inout):: overlap_calculated
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: ovrlp
  type(energy_terms),intent(in) :: energs
  real(8),dimension(:),pointer:: hpsit_c, hpsit_f

  ! Local variables
  integer :: iorb, jorb, iiorb, ilr, ncount, korb, ierr, ist, ncnt, istat, iall
  real(kind=8) :: ddot, tt, eval_zero
  character(len=*),parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8),dimension(:),pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8),dimension(:,:),allocatable :: fnrmOvrlpArr, fnrmArr, lagmat


  allocate(lagmat(tmblarge%orbs%norb,tmblarge%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)
  allocate(fnrmOvrlpArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! by default no quick exit
  energy_increased=.false.

  !!call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)

  if(tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      if(.not. tmblarge%can_use_transposed) then
          allocate(tmblarge%psit_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, tmblarge%psit_c, 'tmblarge%psit_c', subname)
          allocate(tmblarge%psit_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, tmblarge%psit_f, 'tmblarge%psit_f', subname)
          call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, &
               tmblarge%psi, tmblarge%psit_c, tmblarge%psit_f, tmblarge%lzd)
          tmblarge%can_use_transposed=.true.


      end if
      allocate(hpsittmp_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
      allocate(hpsittmp_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)
      if(sum(tmblarge%collcom%nrecvcounts_c)>0) &
          call dcopy(sum(tmblarge%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmblarge%collcom%nrecvcounts_f)>0) &
          call dcopy(7*sum(tmblarge%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)
      call build_linear_combination_transposed(tmblarge%orbs%norb, kernel, tmblarge%collcom, &
           hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
      deallocate(hpsittmp_c, stat=istat)
      call memocc(istat, iall, 'hpsittmp_c', subname)
      iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
      deallocate(hpsittmp_f, stat=istat)
      call memocc(istat, iall, 'hpsittmp_f', subname)
  end if


  call orthoconstraintNonorthogonal(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, tmblarge%mad, &
       tmblarge%collcom, tmblarge%orthpar, tmblarge%wfnmd%bpo, tmblarge%wfnmd%bs, tmblarge%psi, lhphilarge, lagmat, ovrlp, &
       tmblarge%psit_c, tmblarge%psit_f, hpsit_c, hpsit_f, tmblarge%can_use_transposed, overlap_calculated)


  call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, lhphilarge, lhphi)


  ! Calculate trace (or band structure energy, resp.)
  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      trH=0.d0
      do jorb=1,tmb%orbs%norb
          do korb=1,tmb%orbs%norb
              tt = kernel(korb,jorb)*lagmat(korb,jorb)
              trH = trH + tt
          end do
      end do
  else
      trH=0.d0
      do jorb=1,tmb%orbs%norb
          trH = trH + lagmat(jorb,jorb)
      end do
  end if

  ! trH is now the total energy (name is misleading, correct this)
  if(tmb%orbs%nspin==1) trH=2.d0*trH
  trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp



  ! Cycle if the trace increased (steepest descent only)
  if(.not. ldiis%switchSD .and. ldiis%isx==0) then
      if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
          if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
              'WARNING: the target function is larger than it minimal value reached so far:',trH,' > ', ldiis%trmin
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
      if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphi(ist), 1, lhphiold(ist), 1)
      fnrmArr(iorb,1)=ddot(ncount, lhphi(ist), 1, lhphi(ist), 1)
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
  call dcopy(tmb%orbs%npsidim_orbs, lhphi, 1, lhphiold, 1)
  trHold=trH

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)',advance='no') 'Preconditioning... '
  end if
 
  !!call project_gradient(iproc, nproc, tmb, tmb%psi, lhphi)

  !!call get_both_gradients(iproc, nproc, tmb%lzd, tmb%orbs, lhphi, gnrm_in, gnrm_out)

  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr = tmb%orbs%inWhichLocreg(iiorb)
      ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmb%wfnmd%bs%nit_precond, lhphi(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
      ist=ist+ncnt
  end do

  if(iproc==0) then
      write(*,'(a)') 'done.'
  end if

      call timing(iproc,'eglincomms','ON') ! lr408t
  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  alpha_mean=tt/dble(tmb%orbs%norb)
  alpha_max=maxval(alpha)
  call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  call timing(iproc,'eglincomms','OF') ! lr408t

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)

  iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
  deallocate(fnrmOvrlpArr, stat=istat)
  call memocc(istat, iall, 'fnrmOvrlpArr', subname)

  iall=-product(shape(fnrmArr))*kind(fnrmArr)
  deallocate(fnrmArr, stat=istat)
  call memocc(istat, iall, 'fnrmArr', subname)


end subroutine calculate_energy_and_gradient_linear



subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lhphi, lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb
  real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: lhphi, lphiold
  real(kind=8),intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8),dimension(tmb%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  
  ! Local variables
  integer :: istat, iall
  real(kind=8),dimension(:,:),allocatable :: ovrlp
  character(len=*),parameter :: subname='hpsitopsi_linear'


  allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)


  call DIISorSD(iproc, nproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold)
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
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, nproc, it, tmb, ldiis, lhphi, alpha)
  else
      if(iproc==0) write(*,'(1x,a)') 'no improvement of the orbitals, recalculate gradient'
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
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
           tmb%orbs, tmb%op, tmb%comon, tmb%lzd, &
           tmb%mad, tmb%collcom, tmb%orthpar, tmb%wfnmd%bpo, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)

  end if

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear
