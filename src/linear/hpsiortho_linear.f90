subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           variable_locregs, tmbopt, kernel, &
           confdatarr, ldiis, lhphiopt, lphioldopt, lhphioldopt, consecutive_rejections, fnrmArr, &
           fnrmOvrlpArr, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, meanAlpha, ovrlp)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  logical,intent(in):: variable_locregs
  type(DFT_wavefunction),intent(inout):: tmbopt
  real(8),dimension(tmbopt%orbs%norb,tmbopt%orbs%norb),intent(in):: kernel
  type(confpot_data), dimension(tmbopt%orbs%norbp),intent(in):: confdatarr
  type(localizedDIISParameters),intent(inout):: ldiis
  real(8),dimension(tmbopt%wfnmd%nphi),intent(inout):: lhphiopt
  real(8),dimension(tmbopt%wfnmd%nphi),intent(inout):: lphioldopt, lhphioldopt
  integer,intent(inout):: consecutive_rejections
  real(8),dimension(tmbopt%orbs%norb,2),intent(inout):: fnrmArr, fnrmOvrlpArr
  real(8),dimension(tmbopt%orbs%norb),intent(inout):: fnrmOldArr
  real(8),dimension(tmbopt%orbs%norbp),intent(inout):: alpha
  real(8),intent(out):: trH, trHold, fnrm, fnrmMax, meanAlpha
  real(8),dimension(tmbopt%orbs%norb,tmbopt%orbs%norb),intent(in):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, istart, ncount, korb, nvctr_c, nvctr_f, ierr, ind2, ncnt, istat, iall
  real(8):: tt1, tt2, tt3, tt4, tt5,  timecommunp2p, timecommuncoll, timecompress, ddot, tt, eval_zero
  character(len=*),parameter:: subname='calculate_energy_and_gradient_linear'
  real(8),dimension(:,:),allocatable:: lagmat

  allocate(lagmat(tmbopt%orbs%norb,tmbopt%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)



  if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      call allocateSendBufferOrtho(tmbopt%comon, subname)
      call allocateRecvBufferOrtho(tmbopt%comon, subname)
      ! Extract the overlap region from the orbitals phi and store them in tmbopt%comon%sendBuf.
      call extractOrbital3(iproc, nproc, tmbopt%orbs, max(tmbopt%orbs%npsidim_orbs,tmbopt%orbs%npsidim_comp), &
           tmbopt%orbs%inWhichLocreg, tmbopt%lzd, tmbopt%op, &
           lhphiopt, tmbopt%comon%nsendBuf, tmbopt%comon%sendBuf)
      call postCommsOverlapNew(iproc, nproc, tmbopt%orbs, tmbopt%op, tmbopt%lzd, lhphiopt, tmbopt%comon, tt1, tt2)
      call collectnew(iproc, nproc, tmbopt%comon, tmbopt%mad, tmbopt%op, tmbopt%orbs, tmbopt%lzd, tmbopt%comon%nsendbuf, &
           tmbopt%comon%sendbuf, tmbopt%comon%nrecvbuf, tmbopt%comon%recvbuf, tt3, tt4, tt5)
      call build_new_linear_combinations(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon%nrecvbuf, &
           tmbopt%comon%recvbuf, kernel, .true., lhphiopt)
      call deallocateRecvBufferOrtho(tmbopt%comon, subname)
      call deallocateSendBufferOrtho(tmbopt%comon, subname)
  end if
  call orthoconstraintNonorthogonal(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon, tmbopt%mad, ovrlp, &
       tmbopt%orthpar%methTransformOverlap, tmbopt%orthpar%blocksize_pdgemm, tmbopt%psi, lhphiopt, lagmat)


  ! Calculate trace (or band structure energy, resp.)
  if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      trH=0.d0
      do jorb=1,tmbopt%orbs%norb
          do korb=1,tmbopt%orbs%norb
              trH = trH + kernel(korb,jorb)*lagmat(korb,jorb)
          end do
      end do
  else
      trH=0.d0
      do jorb=1,tmbopt%orbs%norb
          trH = trH + lagmat(jorb,jorb)
      end do
  end if


  ! Cycle if the trace increased (steepest descent only)
  if(.not. ldiis%switchSD .and. ldiis%isx==0) then
       if(trH > trHold + 1.d-8*abs(trHold)) then
           consecutive_rejections=consecutive_rejections+1
           if(iproc==0) write(*,'(1x,a,es9.2,a)') 'WARNING: the trace increased by ', 100.d0*(trH-trHold)/abs(trHold), '%.'
           if(consecutive_rejections<=3) then
               ! If the trace increased three times consecutively, do not decrease the step size any more and go on.
               alpha=alpha*.6d0
               if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
                   if(iproc==0) write(*,'(1x,a)') 'Reject orbitals, reuse the old ones and decrease step size.'
                   call dcopy(size(tmbopt%psi), lphioldopt, 1, tmbopt%psi, 1)
               else
                   ! It is not possible to use the old orbitals since the locregs might have changed.
                   if(iproc==0) write(*,'(1x,a)') 'Decrease step size, but accept new orbitals'
               end if
           else
               consecutive_rejections=0
           end if
       else
           consecutive_rejections=0
       end if
  end if





  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  istart=1
  do iorb=1,tmbopt%orbs%norbp
      if(.not.variable_locregs .or. tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
          iiorb=tmbopt%orbs%isorb+iorb
          ilr=tmbopt%orbs%inWhichLocreg(iiorb)
          ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
          if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphioldopt(istart), 1)
          fnrmArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphiopt(istart), 1)
      else
          ! Here the angle between the current and the old gradient cannot be determined since
          ! the locregs might have changed, so we assign to fnrmOvrlpArr a fake value of 1.d0
          iiorb=tmbopt%orbs%isorb+iorb
          ilr=tmbopt%orbs%inWhichLocreg(iiorb)
          ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
          if(it>1) fnrmOvrlpArr(iorb,1)=1.d0
          fnrmArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphiopt(istart), 1)
      end if
      istart=istart+ncount
  end do

  ! Keep the gradient for the next iteration.
  if(it>1) then
      call dcopy(tmbopt%orbs%norbp, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
  end if

  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.
  do iorb=1,tmbopt%orbs%norbp
      fnrm=fnrm+fnrmArr(iorb,1)
      if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
      if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      ! Adapt step size for the steepest descent minimization.
          tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
          if(tt>.9d0 .and. trH<trHold) then
              alpha(iorb)=alpha(iorb)*1.1d0
          else
              alpha(iorb)=alpha(iorb)*.6d0
          end if
      end if
  end do
  call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
  fnrm=sqrt(fnrm/dble(tmbopt%orbs%norb))
  fnrmMax=sqrt(fnrmMax)
  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmbopt%orbs%npsidim_orbs, lhphiopt, 1, lhphioldopt, 1)
  if(variable_locregs .and. tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
      call dcopy(max(tmbopt%orbs%npsidim_orbs,tmbopt%orbs%npsidim_comp), lhphiopt, 1, lhphioldopt, 1)
  trHold=trH

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)') 'Preconditioning.'
  end if


  ind2=1
  do iorb=1,tmbopt%orbs%norbp
      iiorb=tmbopt%orbs%isorb+iorb
      ilr = tmbopt%orbs%inWhichLocreg(iiorb)
      ncnt=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
      call choosePreconditioner2(iproc, nproc, tmbopt%orbs, tmbopt%lzd%llr(ilr), &
           tmbopt%lzd%hgrids(1), tmbopt%lzd%hgrids(2), tmbopt%lzd%hgrids(3), &
           tmbopt%wfnmd%bs%nit_precond, lhphiopt(ind2:ind2+ncnt-1), confdatarr(iorb)%potorder, &
           confdatarr(iorb)%prefac, it, iorb, eval_zero)
      ind2=ind2+ncnt
  end do


  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  meanAlpha=tt/dble(tmbopt%orbs%norb)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)



end subroutine calculate_energy_and_gradient_linear
