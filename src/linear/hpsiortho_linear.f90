!> @file
!! H|psi> and orthonormalization
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, kernel_compr, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, &
           fnrmMax, alpha_mean, alpha_max, energy_increased, tmb, lhphi, lhphiold, &
           tmblarge, lhphilarge, overlap_calculated, energs, hpsit_c, hpsit_f)
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
  real(8),dimension(tmblarge%mad%nvctr),target,intent(in) :: kernel_compr
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(8),dimension(tmb%orbs%norb),intent(inout) :: fnrmOldArr
  real(8),dimension(tmb%orbs%norbp),intent(inout) :: alpha
  real(8),intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(8),intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(8),dimension(tmblarge%orbs%npsidim_orbs),intent(inout):: lhphilarge
  real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphi, lhphiold
  logical,intent(inout):: overlap_calculated
  type(energy_terms),intent(in) :: energs
  real(8),dimension(:),pointer:: hpsit_c, hpsit_f

  ! Local variables
  integer :: iorb, jorb, iiorb, ilr, ncount, korb, ierr, ist, ncnt, istat, iall, ii, iseg, jjorb
  real(kind=8) :: ddot, tt, eval_zero, fnrm_old
  character(len=*),parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8),dimension(:),pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8),dimension(:,:),allocatable :: fnrmOvrlpArr, fnrmArr
  real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
  real(kind=8),dimension(:,:),allocatable :: gnrmArr
  real(kind=8),dimension(:),allocatable :: lagmat_compr, hpsi_tmp
  real(kind=8),dimension(:),pointer :: kernel_compr_tmp

  allocate(fnrmOvrlpArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! by default no quick exit
  energy_increased=.false.

  !!call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)

  allocate(hpsittmp_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
  allocate(hpsittmp_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)

  if(tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_HYBRID) then
      if(sum(tmblarge%collcom%nrecvcounts_c)>0) &
          call dcopy(sum(tmblarge%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmblarge%collcom%nrecvcounts_f)>0) &
          call dcopy(7*sum(tmblarge%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_HYBRID) then
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

      !!call build_linear_combination_transposed(tmblarge%orbs%norb, kernel_compr, tmblarge%collcom, &
      !!     tmblarge%mad, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      call build_linear_combination_transposed(tmblarge%orbs%norb, kernel_compr_tmp, tmblarge%collcom, &
           tmblarge%mad, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)

      if (tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_HYBRID) then
          iall=-product(shape(kernel_compr_tmp))*kind(kernel_compr_tmp)
          deallocate(kernel_compr_tmp, stat=istat)
          call memocc(istat, iall, 'kernel_compr_tmp', subname)

          allocate(hpsi_tmp(tmblarge%orbs%npsidim_orbs), stat=istat)
          call memocc(istat, hpsi_tmp, 'hpsi_tmp', subname)
          call untranspose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, hpsit_c, hpsit_f, hpsi_tmp, tmblarge%lzd)

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
                          call daxpy(ncount, kernel_compr(ii), lhphilarge(ist), 1, hpsi_tmp(ist), 1)
                          ist=ist+ncount
                      end if
                  end do
              end do
          end do
          call dcopy(tmblarge%orbs%npsidim_orbs, hpsi_tmp(1), 1, lhphilarge(1), 1)
          call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)

          iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
          deallocate(hpsi_tmp, stat=istat)
          call memocc(istat, iall, 'hpsi_tmp', subname)
      end if
  end if


  !!do istat=1,size(hpsit_c)
  !!    write(1000+iproc,*) istat, hpsit_c(istat)
  !!end do


  allocate(lagmat_compr(tmblarge%mad%nvctr), stat=istat)
  call memocc(istat, lagmat_compr, 'lagmat_compr', subname)

  call orthoconstraintNonorthogonal(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, tmblarge%mad, &
       tmblarge%collcom, tmblarge%orthpar, tmblarge%wfnmd%bpo, tmblarge%wfnmd%bs, tmblarge%psi, lhphilarge, lagmat_compr, &
       tmblarge%psit_c, tmblarge%psit_f, hpsit_c, hpsit_f, tmblarge%can_use_transposed, overlap_calculated)


  call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, lhphilarge, lhphi)



  !!! Calculate trace (or band structure energy, resp.)
  !!if (tmb%wfnmd%bs%target_function == TARGET_FUNCTION_IS_TRACE) then
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
  !!else
  !!    trH=0.d0
  !!    ii=0
  !!    do iseg=1,tmblarge%mad%nseg
  !!        do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
  !!            ii=ii+1
  !!            trH = trH + kernel_compr(ii)*lagmat_compr(ii)
  !!        end do
  !!    end do
  !!end if



  iall=-product(shape(lagmat_compr))*kind(lagmat_compr)
  deallocate(lagmat_compr, stat=istat)
  call memocc(istat, iall, 'lagmat_compr', subname)


  ! trH is now the total energy (name is misleading, correct this)
  ! Multiply by 2 because when minimizing trace we don't have kernel
  if(tmb%orbs%nspin==1 .and. tmb%wfnmd%bs%target_function/= TARGET_FUNCTION_IS_ENERGY) trH=2.d0*trH
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
  !trHold=trH
  !if (iproc==0) write(*,'(a,2es16.6)') 'BEFORE: fnrm, fnrmmax',fnrm,fnrmmax

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)',advance='no') 'Preconditioning... '
  end if
 
  !!call project_gradient(iproc, nproc, tmb, tmb%psi, lhphi)

  !!call get_both_gradients(iproc, nproc, tmb%lzd, tmb%orbs, lhphi, gnrm_in, gnrm_out)

  !!tt=ddot(tmb%orbs%npsidim_orbs, lhphi,1 , lhphi, 1)
  !!call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!if (iproc==0) write(*,'(a,es20.12)') 'before precond: tt',tt

  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr = tmb%orbs%inWhichLocreg(iiorb)
      ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!if (tmb%confdatarr(iorb)%prefac>=1.d-2) then
      !!    tt=tmb%confdatarr(iorb)%prefac
      !!else
      !!    tt=0.d0
      !!end if
      if(tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_HYBRID) then
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmb%wfnmd%bs%nit_precond, lhphi(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               0.d0, iorb, eval_zero)
      else
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmb%wfnmd%bs%nit_precond, lhphi(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
      end if
      !call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,tmb%lzd%hgrids(1),tmb%lzd%hgrids(2),&
      !      tmb%lzd%hgrids(3),tmb%wfnmd%bs%nit_precond,lhphi,tmb%confdatarr,&
      !      gnrm,gnrm_zero)
      ist=ist+ncnt
  end do

  !!tt=ddot(tmb%orbs%npsidim_orbs, lhphi,1 , lhphi, 1)
  !!call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!if (iproc==0) write(*,'(a,es20.12)') 'after precond: tt',tt

  if(iproc==0) then
      write(*,'(a)') 'done.'
  end if

  !allocate(gnrmArr(tmb%orbs%norb,2), stat=istat)
  !call memocc(istat, gnrmArr, 'gnrmArr', subname)
  !ist=1
  !do iorb=1,tmb%orbs%norbp
  !    iiorb=tmb%orbs%isorb+iorb
  !    ilr=tmb%orbs%inwhichlocreg(iiorb)
  !    ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !    gnrmArr(iorb,1)=ddot(ncount, lhphi(ist), 1, lhphi(ist), 1)
  !    ist=ist+ncount
  !end do
  !gnrm_old=gnrm
  !gnrm=0.d0
  !gnrmMax=0.d0
  !do iorb=1,tmb%orbs%norbp
  !   gnrm=gnrm+gnrmArr(iorb,1)
  !    if(gnrmArr(iorb,1)>gnrmMax) gnrmMax=gnrmArr(iorb,1)
  !end do
  !call mpiallred(gnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !call mpiallred(gnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  !gnrm=sqrt(gnrm/dble(tmb%orbs%norb))
  !gnrmMax=sqrt(gnrmMax)
  !if (iproc==0) write(*,'(a,3es16.6)') 'AFTER: gnrm, gnrmmax, gnrm/gnrm_old',gnrm,gnrmmax,gnrm/gnrm_old
  !iall=-product(shape(gnrmArr))*kind(gnrmArr)
  !deallocate(gnrmArr, stat=istat)
  !call memocc(istat, iall, 'gnrmArr', subname)



      call timing(iproc,'eglincomms','ON') ! lr408t
  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  alpha_mean=tt/dble(tmb%orbs%norb)
  alpha_max=maxval(alpha)
  call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  call timing(iproc,'eglincomms','OF') ! lr408t


  !!if (iproc==0 .and. tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
  !!    do istat=1,tmb%orbs%norb
  !!        do iall=1,tmb%orbs%norb
  !!            write(333,*) istat,iall,lagmat(istat,iall)
  !!        end do
  !!    end do 
  !!end if


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
           lhphi, lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb, tmblarge
  real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: lhphi, lphiold
  real(kind=8),intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8),dimension(tmb%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  
  ! Local variables
  integer :: istat, iall
  character(len=*),parameter :: subname='hpsitopsi_linear'



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
      ! Give tmblarge%mad since this is the correct matrix description
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
           tmb%orbs, tmb%op, tmb%comon, tmb%lzd, &
           tmblarge%mad, tmb%collcom, tmb%orthpar, tmb%wfnmd%bpo, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)

  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear


