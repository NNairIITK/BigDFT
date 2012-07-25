subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, kernel, &
           ldiis, consecutive_rejections, fnrmOldArr, alpha, trH, trHold, fnrm, &
           fnrmMax, meanAlpha, emergency_exit, tmb, lhphi, lhphiold, &
           tmblarge, lhphilarge2, overlap_calculated, ovrlp, energs, hpsit_c, hpsit_f)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  type(DFT_wavefunction),target,intent(inout):: tmblarge, tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: kernel
  type(localizedDIISParameters),intent(inout):: ldiis
  integer,intent(inout):: consecutive_rejections
  real(8),dimension(tmb%orbs%norb),intent(inout):: fnrmOldArr
  real(8),dimension(tmb%orbs%norbp),intent(inout):: alpha
  real(8),intent(out):: trH, trHold, fnrm, fnrmMax, meanAlpha
  logical,intent(out):: emergency_exit
  real(8),dimension(:),target,intent(inout):: lhphilarge2
  real(8),dimension(:),target,intent(inout):: lhphi, lhphiold
  logical,intent(inout):: overlap_calculated
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: ovrlp
  type(energy_terms),intent(in) :: energs
  real(8),dimension(:),intent(out),pointer,optional:: hpsit_c, hpsit_f

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, istart, ncount, korb, nvctr_c, nvctr_f, ierr, ind2, ncnt, istat, iall, jlr, lorb
  real(8):: tt1, tt2, tt3, tt4, tt5,  timecommunp2p, timecommuncoll, timecompress, ddot, tt, eval_zero
  character(len=*),parameter:: subname='calculate_energy_and_gradient_linear'
  real(8),dimension(:),pointer:: hpsittmp_c, hpsittmp_f
  real(8),dimension(:,:),allocatable:: epsmat, fnrmOvrlpArr, fnrmArr, lagmat
  real(8):: closesteval, gnrm_temple
  integer:: owa, owanext

  if(tmblarge%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE .and. &
      (.not.present(hpsit_c) .or. .not.present(hpsit_f))) stop 'ERROR: transposed quantities must be present!'
  if(tmblarge%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
      nullify(hpsittmp_c)
      nullify(hpsittmp_f)
  end if



  allocate(lagmat(tmblarge%orbs%norb,tmblarge%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)
  allocate(fnrmOvrlpArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! by default no quick exit
  emergency_exit=.false.

  !!call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, lhphilarge2, hpsit_c, hpsit_f, tmblarge%lzd)

  if(tmblarge%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      if(tmblarge%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
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
      else if (tmblarge%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call allocateSendBufferOrtho(tmblarge%comon, subname)
          call allocateRecvBufferOrtho(tmblarge%comon, subname)
          ! Extract the overlap region from the orbitals phi and store them in tmblarge%comon%sendBuf.
          call extractOrbital3(iproc, nproc, tmblarge%orbs, tmblarge%orbs, &
               max(tmblarge%orbs%npsidim_orbs,tmblarge%orbs%npsidim_comp), &
               tmblarge%lzd, tmblarge%lzd, tmblarge%op, tmblarge%op, &
               lhphilarge2, tmblarge%comon%nsendBuf, tmblarge%comon%sendBuf)
          !!call postCommsOverlapNew(iproc, nproc, tmblarge%orbs, tmblarge%op, tmblarge%lzd, lhphilarge2, tmblarge%comon, tt1, tt2)
      call timing(iproc,'eglincomms','ON') ! lr408t
          call post_p2p_communication(iproc, nproc, tmblarge%comon%nsendbuf, tmblarge%comon%sendbuf, &
               tmblarge%comon%nrecvbuf, tmblarge%comon%recvbuf, tmblarge%comon)
          !!call collectnew(iproc, nproc, tmblarge%comon, tmblarge%mad, tmblarge%op, tmblarge%orbs, tmblarge%lzd, tmblarge%comon%nsendbuf, &
          !!     tmblarge%comon%sendbuf, tmblarge%comon%nrecvbuf, tmblarge%comon%recvbuf, tt3, tt4, tt5)
          call wait_p2p_communication(iproc, nproc, tmblarge%comon)
      call timing(iproc,'eglincomms','OF') ! lr408t
          call build_new_linear_combinations(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon%nrecvbuf, &
               tmblarge%comon%recvbuf, kernel, .true., lhphilarge2)
          call deallocateRecvBufferOrtho(tmblarge%comon, subname)
          call deallocateSendBufferOrtho(tmblarge%comon, subname)
      end if
  end if


  call orthoconstraintNonorthogonal(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, tmblarge%mad, &
       tmblarge%collcom, tmblarge%orthpar, tmblarge%wfnmd%bpo, tmblarge%wfnmd%bs, tmblarge%psi, lhphilarge2, lagmat, ovrlp, &
       tmblarge%psit_c, tmblarge%psit_f, hpsit_c, hpsit_f, tmblarge%can_use_transposed, overlap_calculated)

  call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, lhphilarge2, lhphi)


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
      if(iproc==0) write(*,*) 'trH, trHold,ldiis%trmin',trH, trHold,ldiis%trmin
      if(trH > ldiis%trmin) then
          consecutive_rejections=consecutive_rejections+1
          if(iproc==0) write(*,'(1x,a,es9.2,a)') 'WARNING: the trace increased by ', 100.d0*(trH-trHold)/abs(trHold), '%.'
              consecutive_rejections=0
              if(iproc==0) write(*,'(1x,a)') 'Energy grows in spite of decreased step size, will exit...'
              emergency_exit=.true.
              !call large_to_small_locreg(iproc,nproc,tmb%lzd,tmblarge%lzd,tmb%orbs,tmblarge%orbs,tmblarge%psi,tmb%psi)
      else
          consecutive_rejections=0
      end if
  end if





  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  istart=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphiold(istart), 1)
      fnrmArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphi(istart), 1)
      fnrmOldArr(iorb)=ddot(ncount, lhphiold(istart), 1, lhphiold(istart), 1)
      istart=istart+ncount
  end do



  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.
  fnrm=0.d0
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
  call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)
  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmb%orbs%npsidim_orbs, lhphi, 1, lhphiold, 1)
  trHold=trH

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)') 'Preconditioning.'
  end if


  !!call get_both_gradients(iproc, nproc, tmb%lzd, tmb%orbs, lhphi, gnrm_in, gnrm_out)

  ind2=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr = tmb%orbs%inWhichLocreg(iiorb)
      ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
               tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmb%wfnmd%bs%nit_precond, lhphi(ind2:ind2+ncnt-1), tmb%confdatarr(iorb)%potorder, &
               tmb%confdatarr(iorb)%prefac, iorb, eval_zero, tmb, kernel)
      ind2=ind2+ncnt
  end do


      call timing(iproc,'eglincomms','ON') ! lr408t
  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, mpi_comm_world, ierr)
  meanAlpha=tt/dble(tmb%orbs%norb)
  tt=maxval(alpha)
  call mpiallred(tt, 1, mpi_max, mpi_comm_world, ierr)
  if(iproc==0) write(*,'(a,es12.4)') 'max alpha',tt
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
           lhphi, lphiold, alpha, trH, meanAlpha, alphaDIIS)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  type(localizedDIISParameters),intent(inout):: ldiis
  type(DFT_wavefunction),target,intent(inout):: tmb
  real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphi, lphiold
  real(8),intent(in):: trH, meanAlpha
  real(8),dimension(tmb%orbs%norbp),intent(out):: alpha, alphaDIIS
  
  ! Local variables
  integer:: ist, iorb, iiorb, ilrlarge, ncnt, istat, iall, ilr
  real(8):: tt
  real(8),dimension(:,:),allocatable:: ovrlp
  character(len=*),parameter:: subname='hpsitopsi_linear'


  allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)


  call DIISorSD(iproc, nproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold)
  if(iproc==0) then
      if(ldiis%isx>0) then
          write(*,'(1x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
              ldiis%icountDIISFailureCons, ', total failures=', ldiis%icountDIISFailureTot
      else
          write(*,'(1x,a,es9.3,a,i0,a)') 'steepest descent informations: mean alpha=', meanAlpha, &
          ', consecutive successes=', ldiis%icountSDSatur, ', DIIS=y'
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



  do_ortho_if2: if(.not.ldiis%switchSD) then

      !tmbopt => tmb
      !lhphiopt => lhphi
      !tmbopt%confdatarr => tmb%confdatarr
      !if (tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
           tmb%orbs, tmb%op, tmb%comon, tmb%lzd, &
           tmb%mad, tmb%collcom, tmb%orthpar, tmb%wfnmd%bpo, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)

  end if do_ortho_if2

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear
