subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           tmbopt, kernel, &
           ldiis, lhphiopt, lphioldopt, lhphioldopt, consecutive_rejections, fnrmArr, &
           fnrmOvrlpArr, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, gnrm_in, gnrm_out, meanAlpha, emergency_exit, &
           tmb, lhphi, lphiold, lhphiold, &
           tmblarge2, lhphilarge2, lphilargeold2, lhphilargeold2, orbs, overlap_calculated, ovrlp)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  type(DFT_wavefunction),pointer,intent(inout):: tmbopt
  real(8),dimension(tmbopt%orbs%norb,tmbopt%orbs%norb),intent(inout):: kernel
  type(localizedDIISParameters),intent(inout):: ldiis
  real(8),dimension(:),pointer,intent(inout):: lhphiopt
  real(8),dimension(:),pointer,intent(inout):: lphioldopt, lhphioldopt
  integer,intent(inout):: consecutive_rejections
  real(8),dimension(tmbopt%orbs%norb,2),intent(inout):: fnrmArr, fnrmOvrlpArr
  real(8),dimension(tmbopt%orbs%norb),intent(inout):: fnrmOldArr
  real(8),dimension(tmbopt%orbs%norbp),intent(inout):: alpha
  real(8),intent(out):: trH, trHold, fnrm, fnrmMax, meanAlpha, gnrm_in, gnrm_out
  logical,intent(out):: emergency_exit
  type(DFT_wavefunction),target,intent(inout):: tmblarge2, tmb
  real(8),dimension(:),target,intent(inout):: lhphilarge2
  real(8),dimension(:),target,intent(inout):: lphilargeold2, lhphilargeold2, lhphi, lphiold, lhphiold
  type(orbitals_data),intent(in):: orbs
  logical,intent(inout):: overlap_calculated
  real(8),dimension(tmbopt%orbs%norb,tmbopt%orbs%norb),intent(inout):: ovrlp

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, istart, ncount, korb, nvctr_c, nvctr_f, ierr, ind2, ncnt, istat, iall, jlr, lorb
  real(8):: tt1, tt2, tt3, tt4, tt5,  timecommunp2p, timecommuncoll, timecompress, ddot, tt, eval_zero
  character(len=*),parameter:: subname='calculate_energy_and_gradient_linear'
  real(8),dimension(:),pointer:: hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f
  real(8),dimension(:,:),allocatable:: lagmat, epsmat
  real(8):: closesteval, gnrm_temple
  integer:: owa, owanext

  nullify(hpsit_c)
  nullify(hpsit_f)
  nullify(hpsittmp_c)
  nullify(hpsittmp_f)



  allocate(lagmat(tmbopt%orbs%norb,tmbopt%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)

  ! by default no quick exit
  emergency_exit=.false.




 

  if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      if(tmbopt%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
          if(.not. tmbopt%can_use_transposed) then
              allocate(tmbopt%psit_c(sum(tmbopt%collcom%nrecvcounts_c)), stat=istat)
              call memocc(istat, tmbopt%psit_c, 'tmbopt%psit_c', subname)
              allocate(tmbopt%psit_f(7*sum(tmbopt%collcom%nrecvcounts_f)), stat=istat)
              call memocc(istat, tmbopt%psit_f, 'tmbopt%psit_f', subname)
              call transpose_localized(iproc, nproc, tmbopt%orbs, tmbopt%collcom, &
                   tmbopt%psi, tmbopt%psit_c, tmbopt%psit_f, tmbopt%lzd)
              tmbopt%can_use_transposed=.true.
          end if
          allocate(hpsit_c(sum(tmbopt%collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, hpsit_c, 'hpsit_c', subname)
          allocate(hpsit_f(7*sum(tmbopt%collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, hpsit_f, 'hpsit_f', subname)
          allocate(hpsittmp_c(sum(tmbopt%collcom%nrecvcounts_c)), stat=istat)
          call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
          allocate(hpsittmp_f(7*sum(tmbopt%collcom%nrecvcounts_f)), stat=istat)
          call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)
          call transpose_localized(iproc, nproc, tmbopt%orbs, tmbopt%collcom, lhphiopt, hpsit_c, hpsit_f, tmbopt%lzd)
          call dcopy(sum(tmbopt%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
          call dcopy(7*sum(tmbopt%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)
          call build_linear_combination_transposed(tmbopt%orbs%norb, kernel, tmbopt%collcom, &
               hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f)
          iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
          deallocate(hpsittmp_c, stat=istat)
          call memocc(istat, iall, 'hpsittmp_c', subname)
          iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
          deallocate(hpsittmp_f, stat=istat)
          call memocc(istat, iall, 'hpsittmp_f', subname)
      else if (tmbopt%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call allocateSendBufferOrtho(tmbopt%comon, subname)
          call allocateRecvBufferOrtho(tmbopt%comon, subname)
          ! Extract the overlap region from the orbitals phi and store them in tmbopt%comon%sendBuf.
          call extractOrbital3(iproc, nproc, tmbopt%orbs, tmbopt%orbs, max(tmbopt%orbs%npsidim_orbs,tmbopt%orbs%npsidim_comp), &
               tmbopt%lzd, tmbopt%lzd, tmbopt%op, tmbopt%op, &
               lhphiopt, tmbopt%comon%nsendBuf, tmbopt%comon%sendBuf)
          !!call postCommsOverlapNew(iproc, nproc, tmbopt%orbs, tmbopt%op, tmbopt%lzd, lhphiopt, tmbopt%comon, tt1, tt2)
      call timing(iproc,'eglincomms','ON') ! lr408t
          call post_p2p_communication(iproc, nproc, tmbopt%comon%nsendbuf, tmbopt%comon%sendbuf, &
               tmbopt%comon%nrecvbuf, tmbopt%comon%recvbuf, tmbopt%comon)
          !!call collectnew(iproc, nproc, tmbopt%comon, tmbopt%mad, tmbopt%op, tmbopt%orbs, tmbopt%lzd, tmbopt%comon%nsendbuf, &
          !!     tmbopt%comon%sendbuf, tmbopt%comon%nrecvbuf, tmbopt%comon%recvbuf, tt3, tt4, tt5)
          call wait_p2p_communication(iproc, nproc, tmbopt%comon)
      call timing(iproc,'eglincomms','OF') ! lr408t
          call build_new_linear_combinations(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon%nrecvbuf, &
               tmbopt%comon%recvbuf, kernel, .true., lhphiopt)
          call deallocateRecvBufferOrtho(tmbopt%comon, subname)
          call deallocateSendBufferOrtho(tmbopt%comon, subname)
      end if
  end if


  call orthoconstraintNonorthogonal(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon, tmbopt%mad, &
       tmbopt%collcom, tmbopt%orthpar, tmbopt%wfnmd%bpo, tmbopt%wfnmd%bs, tmbopt%psi, lhphiopt, lagmat, ovrlp, &
       tmbopt%psit_c, tmbopt%psit_f, hpsit_c, hpsit_f, tmbopt%can_use_transposed, overlap_calculated)

  if(associated(hpsit_c)) then
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
  end if
  if(associated(hpsit_f)) then
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)
  end if

  tmbopt => tmb
  lhphiopt => lhphi
  lphioldopt => lphiold
  lhphioldopt => lhphiold
  call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, lhphilarge2, lhphi)


  ! Calculate trace (or band structure energy, resp.)
  if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      trH=0.d0
      do jorb=1,tmbopt%orbs%norb
          do korb=1,tmbopt%orbs%norb
              tt = kernel(korb,jorb)*lagmat(korb,jorb)
              trH = trH + tt
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
  !if(ldiis%isx==0) then
       if(iproc==0) write(*,*) 'trH, trHold',trH, trHold
       !if(trH > trHold + 1.d-8*abs(trHold)) then
       if(trH > ldiis%trmin) then
           consecutive_rejections=consecutive_rejections+1
           if(iproc==0) write(*,'(1x,a,es9.2,a)') 'WARNING: the trace increased by ', 100.d0*(trH-trHold)/abs(trHold), '%.'
           !!if(consecutive_rejections<=3) then
           !!    ! If the trace increased three times consecutively, do not decrease the step size any more and go on.
           !!    alpha=alpha*.6d0
           !!    if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
           !!        if(iproc==0) write(*,'(1x,a)') 'Reject orbitals, reuse the old ones and decrease step size.'
           !!        call dcopy(size(tmbopt%psi), lphioldopt, 1, tmbopt%psi, 1)
           !!        if(.not.variable_locregs) then
           !!            call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
           !!        end if
           !!    else if(.not.variable_locregs) then
           !!        if(iproc==0) write(*,'(1x,a)') 'Reject orbitals, reuse the old ones and decrease step size.'
           !!        call dcopy(size(tmbopt%psi), lphioldopt, 1, tmbopt%psi, 1)
           !!        if(.not.variable_locregs) then
           !!            call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
           !!        end if
           !!    else 
           !!        ! It is not possible to use the old orbitals since the locregs might have changed.
           !!        !if(iproc==0) write(*,'(1x,a)') 'Decrease step size, but accept new orbitals'
           !!        if(iproc==0) write(*,'(1x,a)') 'Energy grows, will exit...'
           !!        emergency_exit=.true.
           !!        if(.not.variable_locregs) then
           !!            call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
           !!        end if
           !!    end if
           !!else
               consecutive_rejections=0
               if(iproc==0) write(*,'(1x,a)') 'Energy grows in spite of decreased step size, will exit...'
               emergency_exit=.true.
               !!if(.not.variable_locregs) then
                   call large_to_small_locreg(iproc,nproc,tmb%lzd,tmblarge2%lzd,tmb%orbs,tmblarge2%orbs,tmblarge2%psi,tmb%psi)
               !!end if
           !end if
       else
           consecutive_rejections=0
       end if
  end if





  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  istart=1
  do iorb=1,tmbopt%orbs%norbp
      !!if(.not.variable_locregs .or. tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
          iiorb=tmbopt%orbs%isorb+iorb
          ilr=tmbopt%orbs%inwhichlocreg(iiorb)
          !!owa=tmbopt%orbs%onwhichatom(iiorb)
          !!if(iiorb<tmbopt%orbs%norb) then
          !!    owanext=tmbopt%orbs%onwhichatom(iiorb+1)
          !!else
          !!    owanext=tmbopt%lzd%nlr+1
          !!end if
          ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
          if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphioldopt(istart), 1)
          fnrmArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphiopt(istart), 1)
          !!write(1000+iiorb,'(2es14.5)') lagmat(iiorb,iiorb), fnrmArr(iorb,1)
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
  fnrm=0.d0
  do iorb=1,tmbopt%orbs%norbp
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
      end if
  end do
  call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
  fnrm=sqrt(fnrm/dble(tmbopt%orbs%norb))
  fnrmMax=sqrt(fnrmMax)
  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmbopt%orbs%npsidim_orbs, lhphiopt, 1, lhphioldopt, 1)
  trHold=trH

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)') 'Preconditioning.'
  end if


  call get_both_gradients(iproc, nproc, tmbopt%lzd, tmbopt%orbs, lhphiopt, gnrm_in, gnrm_out)

  ind2=1
  do iorb=1,tmbopt%orbs%norbp
      iiorb=tmbopt%orbs%isorb+iorb
      ilr = tmbopt%orbs%inWhichLocreg(iiorb)
      ncnt=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
          call choosePreconditioner2(iproc, nproc, tmbopt%orbs, tmbopt%lzd%llr(ilr), &
               tmbopt%lzd%hgrids(1), tmbopt%lzd%hgrids(2), tmbopt%lzd%hgrids(3), &
               tmbopt%wfnmd%bs%nit_precond, lhphiopt(ind2:ind2+ncnt-1), tmbopt%confdatarr(iorb)%potorder, &
               tmbopt%confdatarr(iorb)%prefac, it, iorb, eval_zero)
      ind2=ind2+ncnt
  end do



  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, mpi_comm_world, ierr)
  meanAlpha=tt/dble(tmbopt%orbs%norb)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)


end subroutine calculate_energy_and_gradient_linear



subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmblarge, tmb, tmbopt, at, rxyz, kernel, &
           lhphilarge, lphilargeold, lhphilargeold, lhphi, lphiold, lhphiold, lhphiopt, lphioldopt, &
           alpha, locregCenter, locregCenterTemp, &
           denspot, locrad, inwhichlocreg_reference, factor, trH, meanAlpha, alphaDIIS)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  type(localizedDIISParameters),intent(inout):: ldiis
  type(DFT_wavefunction),target,intent(inout):: tmblarge, tmb
  type(DFT_wavefunction),pointer,intent(inout):: tmbopt
  type(atoms_data),intent(in):: at
  real(8),dimension(3,at%nat),intent(in):: rxyz
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: kernel
  real(8),dimension(:),pointer,intent(inout):: lhphilarge, lphilargeold, lhphilargeold
  real(8),dimension(:),pointer,intent(inout):: lhphi, lphiold, lhphiold
  real(8),dimension(:),pointer,intent(inout):: lhphiopt
  real(8),dimension(:),pointer,intent(out):: lphioldopt
  real(8),dimension(3,tmb%lzd%nlr),intent(inout):: locregCenter
  real(8),dimension(3,tmb%lzd%nlr),intent(inout):: locregCenterTemp
  type(DFT_local_fields),intent(inout):: denspot
  real(8),dimension(tmb%lzd%nlr),intent(in):: locrad
  integer,dimension(tmb%orbs%norb),intent(in):: inwhichlocreg_reference
  real(8),intent(in):: factor, trH, meanAlpha
  real(8),dimension(tmb%orbs%norbp),intent(out):: alpha, alphaDIIS
  
  ! Local variables
  integer:: ist, iorb, iiorb, ilrlarge, ncnt, istat, iall, ilr
  real(8):: tt, dnrm2
  real(8),dimension(:,:),allocatable:: Umat, ovrlp
  integer,dimension(:),allocatable:: onwhichatom_reference
  real(8),dimension(:),allocatable:: locrad_tmp
  character(len=*),parameter:: subname='hpsitopsi_linear'


  allocate(Umat(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, Umat, 'Umat', subname)

  allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  allocate(onwhichatom_reference(tmb%orbs%norb), stat=istat)
  call memocc(istat, onwhichatom_reference, 'onwhichatom_reference', subname)

  allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
  call memocc(istat, locrad_tmp, 'locrad_tmp', subname)


  call DIISorSD(iproc, nproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt)
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
      call improveOrbitals(iproc, nproc, it, tmbopt, ldiis, lhphiopt, alpha)
  else
      if(iproc==0) write(*,'(1x,a)') 'no improvement of the orbitals, recalculate gradient'
  end if

  ! The transposed quantities can now not be used any more...
  if(tmbopt%can_use_transposed) then
      iall=-product(shape(tmbopt%psit_c))*kind(tmbopt%psit_c)
      deallocate(tmbopt%psit_c, stat=istat)
      call memocc(istat, iall, 'tmbopt%psit_c', subname)
      iall=-product(shape(tmbopt%psit_f))*kind(tmbopt%psit_f)
      deallocate(tmbopt%psit_f, stat=istat)
      call memocc(istat, iall, 'tmbopt%psit_f', subname)
      tmbopt%can_use_transposed=.false.
  end if



  do_ortho_if2: if(.not.ldiis%switchSD) then

      tmbopt => tmb
      lhphiopt => lhphi
      lphioldopt => lphiold
      tmbopt%confdatarr => tmb%confdatarr
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
      tmbopt%orbs, tmbopt%op, tmbopt%comon, tmbopt%lzd, &
      tmbopt%mad, tmbopt%collcom, tmbopt%orthpar, tmbopt%wfnmd%bpo, tmbopt%psi, tmbopt%psit_c, tmbopt%psit_f, &
      tmbopt%can_use_transposed)

  end if do_ortho_if2

  iall=-product(shape(Umat))*kind(Umat)
  deallocate(Umat, stat=istat)
  call memocc(istat, iall, 'Umat', subname)

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  iall=-product(shape(onwhichatom_reference))*kind(onwhichatom_reference)
  deallocate(onwhichatom_reference, stat=istat)
  call memocc(istat, iall, 'onwhichatom_reference', subname)

  iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
  deallocate(locrad_tmp, stat=istat)
  call memocc(istat, iall, 'locrad_tmp', subname)

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear
