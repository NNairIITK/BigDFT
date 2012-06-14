subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           variable_locregs, tmbopt, kernel, &
           ldiis, lhphiopt, lphioldopt, lhphioldopt, consecutive_rejections, fnrmArr, &
           fnrmOvrlpArr, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, meanAlpha, emergency_exit, &
           tmb, lhphi, lphiold, lhphiold, &
           tmblarge2, lhphilarge2, lphilargeold2, lhphilargeold2, orbs)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  logical,intent(in):: variable_locregs
  type(DFT_wavefunction),pointer,intent(inout):: tmbopt
  real(8),dimension(tmbopt%orbs%norb,tmbopt%orbs%norb),intent(inout):: kernel
  type(localizedDIISParameters),intent(inout):: ldiis
  real(8),dimension(:),pointer,intent(inout):: lhphiopt
  real(8),dimension(:),pointer,intent(inout):: lphioldopt, lhphioldopt
  integer,intent(inout):: consecutive_rejections
  real(8),dimension(tmbopt%orbs%norb,2),intent(inout):: fnrmArr, fnrmOvrlpArr
  real(8),dimension(tmbopt%orbs%norb),intent(inout):: fnrmOldArr
  real(8),dimension(tmbopt%orbs%norbp),intent(inout):: alpha
  real(8),intent(out):: trH, trHold, fnrm, fnrmMax, meanAlpha
  logical,intent(out):: emergency_exit
  type(DFT_wavefunction),target,intent(inout):: tmblarge2, tmb
  real(8),dimension(:),target,intent(inout):: lhphilarge2
  real(8),dimension(:),target,intent(inout):: lphilargeold2, lhphilargeold2, lhphi, lphiold, lhphiold
  type(orbitals_data),intent(in):: orbs

  ! Local variables
  integer:: iorb, jorb, iiorb, ilr, istart, ncount, korb, nvctr_c, nvctr_f, ierr, ind2, ncnt, istat, iall, jlr, lorb
  real(8):: tt1, tt2, tt3, tt4, tt5,  timecommunp2p, timecommuncoll, timecompress, ddot, tt, eval_zero
  character(len=*),parameter:: subname='calculate_energy_and_gradient_linear'
  real(8),dimension(:,:),allocatable:: lagmat, epsmat, ovrlp
  real(8):: closesteval, gnrm_temple



  do iorb=1,tmbopt%orbs%norb
      ilr=tmbopt%orbs%inwhichlocreg(iorb)
      do jorb=1,tmbopt%orbs%norb
          jlr=tmbopt%orbs%inwhichlocreg(jorb)
          tt = (tmbopt%lzd%llr(ilr)%locregCenter(1)-tmbopt%lzd%llr(jlr)%locregCenter(1))**2 &
              +(tmbopt%lzd%llr(ilr)%locregCenter(2)-tmbopt%lzd%llr(jlr)%locregCenter(2))**2 &
              +(tmbopt%lzd%llr(ilr)%locregCenter(3)-tmbopt%lzd%llr(jlr)%locregCenter(3))**2
          tt=sqrt(tt)
          !if(tt>10.d0) kernel(jorb,iorb)=0.d0
          if (iproc==0) write(999,'(2i8,es12.4,es15.6)') jorb,iorb,tt,kernel(jorb,iorb)
      end do
  end do
  if (iproc==0) write(999,*) '-------------------------------------------------------------------'


  allocate(lagmat(tmbopt%orbs%norb,tmbopt%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)

  ! by default no quick exit
  emergency_exit=.false.



  istart=1
  fnrm=0.d0
  do iorb=1,tmbopt%orbs%norbp
      iiorb=tmbopt%orbs%isorb+iorb
      ilr=tmbopt%orbs%inWhichLocreg(iiorb)
      ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
      fnrmArr(iorb,1)=ddot(ncount, lhphiopt(istart), 1, lhphiopt(istart), 1)
      istart=istart+ncount
  end do
  do iorb=1,tmbopt%orbs%norbp
      fnrm=fnrm+fnrmArr(iorb,1)
  end do
  call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
  fnrm=sqrt(fnrm/dble(tmbopt%orbs%norb))
  if(iproc==0) write(*,*) 'FISRT FNRM',fnrm




 

  if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      call allocateSendBufferOrtho(tmbopt%comon, subname)
      call allocateRecvBufferOrtho(tmbopt%comon, subname)
      ! Extract the overlap region from the orbitals phi and store them in tmbopt%comon%sendBuf.
      call extractOrbital3(iproc, nproc, tmbopt%orbs, tmbopt%orbs, max(tmbopt%orbs%npsidim_orbs,tmbopt%orbs%npsidim_comp), &
           tmbopt%lzd, tmbopt%lzd, tmbopt%op, tmbopt%op, &
           lhphiopt, tmbopt%comon%nsendBuf, tmbopt%comon%sendBuf)
      !!call postCommsOverlapNew(iproc, nproc, tmbopt%orbs, tmbopt%op, tmbopt%lzd, lhphiopt, tmbopt%comon, tt1, tt2)
      call post_p2p_communication(iproc, nproc, tmbopt%comon%nsendbuf, tmbopt%comon%sendbuf, &
           tmbopt%comon%nrecvbuf, tmbopt%comon%recvbuf, tmbopt%comon)
      !!call collectnew(iproc, nproc, tmbopt%comon, tmbopt%mad, tmbopt%op, tmbopt%orbs, tmbopt%lzd, tmbopt%comon%nsendbuf, &
      !!     tmbopt%comon%sendbuf, tmbopt%comon%nrecvbuf, tmbopt%comon%recvbuf, tt3, tt4, tt5)
      call wait_p2p_communication(iproc, nproc, tmbopt%comon)
      call build_new_linear_combinations(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon%nrecvbuf, &
           tmbopt%comon%recvbuf, kernel, .true., lhphiopt)
      call deallocateRecvBufferOrtho(tmbopt%comon, subname)
      call deallocateSendBufferOrtho(tmbopt%comon, subname)
  end if

  !!tmbopt => tmb
  !!lhphiopt => lhphi
  !!lphioldopt => lphiold
  !!lhphioldopt => lhphiold
  !!if(.not.variable_locregs) then
  !!    call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, lhphilarge2, lhphi)
  !!end if

  call orthoconstraintNonorthogonal(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%op, tmbopt%comon, tmbopt%mad, &
       tmbopt%collcom, tmbopt%orthpar, tmbopt%wfnmd%bpo, tmbopt%psi, lhphiopt, lagmat, &
       tmbopt%psit_c, tmbopt%psit_f, tmbopt%can_use_transposed)


  !! EXPERIMENTAL ######################################################################
  ! calculate Lagrange multiplier matrix for psi
  allocate(epsmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, epsmat, 'epsmat', subname)
  allocate(ovrlp(tmbopt%orbs%norb,tmbopt%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)
  do iorb=1,orbs%norb
     do jorb=1,orbs%norb
         tt=0.d0
         do korb=1,tmbopt%orbs%norb
             do lorb=1,tmbopt%orbs%norb
                 tt=tt+tmb%wfnmd%coeff(korb,iorb)*tmb%wfnmd%coeff(lorb,jorb)*lagmat(korb,lorb)
             end do
         end do
         epsmat(jorb,iorb)=tt
     end do
  end do
  call getOverlapMatrix2(iproc, nproc, tmbopt%lzd, tmbopt%orbs, tmbopt%comon, tmbopt%op, tmbopt%psi, tmbopt%mad, ovrlp)

  tt=0.d0
  do iorb=1,tmbopt%orbs%norb
      do jorb=1,orbs%norb
          do korb=1,tmbopt%orbs%norb
              tt=tt+tmb%wfnmd%coeff(korb,jorb)*lagmat(korb,iorb)
          end do
          do korb=1,orbs%norb
              do lorb=1,tmbopt%orbs%norb
                  tt=tt-epsmat(korb,jorb)*tmb%wfnmd%coeff(lorb,korb)*ovrlp(lorb,iorb)
              end do
          end do
      end do
  end do
  if(iproc==0) write(*,*) 'tt',tt

  iall=-product(shape(epsmat))*kind(epsmat)
  deallocate(epsmat, stat=istat)
  call memocc(istat, iall, 'epsmat', subname)
  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)
  !! ###################################################################################

      if(iproc==0) write(*,*) 'EXPERIMENTAL: precond here with confinement'
  !!if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
  !!    ind2=1
  !!    do iorb=1,tmbopt%orbs%norbp
  !!        iiorb=tmbopt%orbs%isorb+iorb
  !!        ilr = tmbopt%orbs%inWhichLocreg(iiorb)
  !!        ncnt=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
  !!        call choosePreconditioner2(iproc, nproc, tmbopt%orbs, tmbopt%lzd%llr(ilr), &
  !!             tmbopt%lzd%hgrids(1), tmbopt%lzd%hgrids(2), tmbopt%lzd%hgrids(3), &
  !!             5, lhphiopt(ind2:ind2+ncnt-1), tmbopt%confdatarr(iorb)%potorder, &
  !!             1.d-3, it, iorb, eval_zero)
  !!    end do
  !!end if


  call get_weighted_gradient(iproc, nproc, tmbopt%lzd, tmbopt%orbs, lhphiopt)
  call plot_gradient(iproc, nproc, 1000, tmbopt%lzd, tmbopt%orbs, lhphiopt)
  call plot_gradient(iproc, nproc, 2000, tmbopt%lzd, tmbopt%orbs, tmbopt%psi)


  !!tmbopt => tmb
  !!lhphiopt => lhphi
  !!lphioldopt => lphiold
  !!lhphioldopt => lhphiold
  if(.not.variable_locregs) then
      call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, lhphilarge2, lhphi)
      !call flatten_at_boundaries(tmb%lzd, tmb%orbs, lhphi)
      !!do istat=1,tmb%orbs%npsidim_orbs
      !!    write(200+iproc,*) istat, lhphi(istat)
      !!end do
      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, lhphi, lhphilarge2)
  end if

  call plot_gradient(iproc, nproc, 3000, tmbopt%lzd, tmbopt%orbs, lhphiopt)
  call plot_gradient(iproc, nproc, 4000, tmbopt%lzd, tmbopt%orbs, tmbopt%psi)




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
  !if(.not. ldiis%switchSD .and. ldiis%isx==0) then
  if(ldiis%isx==0) then
       if(trH > trHold + 1.d-8*abs(trHold)) then
           consecutive_rejections=consecutive_rejections+1
           if(iproc==0) write(*,'(1x,a,es9.2,a)') 'WARNING: the trace increased by ', 100.d0*(trH-trHold)/abs(trHold), '%.'
           if(consecutive_rejections<=3) then
               ! If the trace increased three times consecutively, do not decrease the step size any more and go on.
               alpha=alpha*.6d0
               if(tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
                   if(iproc==0) write(*,'(1x,a)') 'Reject orbitals, reuse the old ones and decrease step size.'
                   call dcopy(size(tmbopt%psi), lphioldopt, 1, tmbopt%psi, 1)
                   if(.not.variable_locregs) then
                       call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
                   end if
               else if(.not.variable_locregs) then
                   if(iproc==0) write(*,'(1x,a)') 'Reject orbitals, reuse the old ones and decrease step size.'
                   call dcopy(size(tmbopt%psi), lphioldopt, 1, tmbopt%psi, 1)
                   if(.not.variable_locregs) then
                       call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
                   end if
               else 
                   ! It is not possible to use the old orbitals since the locregs might have changed.
                   !if(iproc==0) write(*,'(1x,a)') 'Decrease step size, but accept new orbitals'
                   if(iproc==0) write(*,'(1x,a)') 'Energy grows, will exit...'
                   emergency_exit=.true.
               end if
           else
               consecutive_rejections=0
               if(iproc==0) write(*,'(1x,a)') 'Energy grows in spite of decreased step size, will exit...'
               emergency_exit=.true.
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



  ! TEST: Temple Criterium ###########################
  gnrm_temple=0.d0
  do iorb=1,tmbopt%orbs%norbp
      iiorb=tmbopt%orbs%isorb+iorb
      closesteval=1.d100
      do jorb=1,tmbopt%orbs%norb
          if(jorb/=iiorb) then
              tt=abs(lagmat(iiorb,iiorb)-lagmat(jorb,jorb))
              if(tt<=closesteval) then
                  closesteval=tt
              end if
          end if
      end do
      !write(1000+iproc,'(a,4es14.5)') 'iiorb, fnrmArr(iorb,i), closesteval, Temple', iiorb, fnrmArr(iorb,1), closesteval, fnrmArr(iorb,1)/closesteval
      !write(*,'(a,4es14.5)') 'iiorb, fnrmArr(iorb,i), closesteval, Temple', iiorb, fnrmArr(iorb,1), closesteval, fnrmArr(iorb,1)/closesteval
      !!write(*,*) 'iiorb, fnrmArr(iorb,i), closesteval, Temple', iiorb, fnrmArr(iorb,1), closesteval, fnrmArr(iorb,1)/closesteval
      gnrm_temple = gnrm_temple + fnrmArr(iorb,1)/closesteval
  end do
  call mpiallred(gnrm_temple, 1, mpi_sum, mpi_comm_world, ierr)
  gnrm_temple = gnrm_temple/dble(tmbopt%orbs%norb)
  !if(iproc==0) write(*,*) 'GNRM TEMPLE', gnrm_temple
  ! ##################################################



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
      !if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      if(it>1 .and. ldiis%isx==0) then
      ! Adapt step size for the steepest descent minimization.
          tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
          !if(tt>.9d0 .and. trH<trHold) then
          !if(tt>.7d0 .and. trH<trHold) then
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
  if(variable_locregs .and. tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
      call dcopy(max(tmbopt%orbs%npsidim_orbs,tmbopt%orbs%npsidim_comp), lhphiopt, 1, lhphioldopt, 1)
  if(trH <= trHold + 1.d-8*abs(trHold)) then
      trHold=trH
  end if

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)') 'Preconditioning.'
  end if


  !!if(iproc==0) write(*,*) 'WARNING: set prefactor manually for preconditioning!'
  ind2=1
  do iorb=1,tmbopt%orbs%norbp
      iiorb=tmbopt%orbs%isorb+iorb
      ilr = tmbopt%orbs%inWhichLocreg(iiorb)
      ncnt=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
      !if (tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
          call choosePreconditioner2(iproc, nproc, tmbopt%orbs, tmbopt%lzd%llr(ilr), &
               tmbopt%lzd%hgrids(1), tmbopt%lzd%hgrids(2), tmbopt%lzd%hgrids(3), &
               tmbopt%wfnmd%bs%nit_precond, lhphiopt(ind2:ind2+ncnt-1), tmbopt%confdatarr(iorb)%potorder, &
               tmbopt%confdatarr(iorb)%prefac, it, iorb, eval_zero)
      !!else if (tmbopt%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      !!    call choosePreconditioner2(iproc, nproc, tmbopt%orbs, tmbopt%lzd%llr(ilr), &
      !!         tmbopt%lzd%hgrids(1), tmbopt%lzd%hgrids(2), tmbopt%lzd%hgrids(3), &
      !!         tmbopt%wfnmd%bs%nit_precond, lhphiopt(ind2:ind2+ncnt-1), tmbopt%confdatarr(iorb)%potorder, &
      !!         2.d-3, it, iorb, eval_zero)
      !!end if
      ind2=ind2+ncnt
  end do

  call plot_gradient(iproc, nproc, 5000, tmbopt%lzd, tmbopt%orbs, lhphiopt)
  call plot_gradient(iproc, nproc, 6000, tmbopt%lzd, tmbopt%orbs, tmbopt%psi)

  tmbopt => tmb
  lhphiopt => lhphi
  lphioldopt => lphiold
  lhphioldopt => lhphiold
  if(.not.variable_locregs) then
      call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, lhphilarge2, lhphi)
  end if


  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, mpi_comm_world, ierr)
  meanAlpha=tt/dble(tmbopt%orbs%norb)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)


end subroutine calculate_energy_and_gradient_linear



subroutine hpsitopsi_linear(iproc, nproc, it, variable_locregs, ldiis, tmblarge, tmb, tmbopt, at, rxyz, kernel, &
           lhphilarge, lphilargeold, lhphilargeold, lhphi, lphiold, lhphiold, lhphiopt, lphioldopt, &
           alpha, locregCenter, locregCenterTemp, &
           denspot, locrad, inwhichlocreg_reference, factor, trH, meanAlpha, alphaDIIS)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  logical,intent(in):: variable_locregs
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
  !!if(iproc==0) write(*,*) 'ldiis%switchSD',ldiis%switchSD

  ! Improve the orbitals, depending on the choice made above.
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, nproc, it, variable_locregs, tmbopt, ldiis, lhphiopt, alpha)
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

  
  newgradient_if_2: if(variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      call update_confdatarr(tmblarge%lzd, tmblarge%orbs, locregCenterTemp, tmb%confdatarr)
      ! Normalize tmblarge%psi
      if(variable_locregs) then
          ist=1
          do iorb=1,tmblarge%orbs%norbp
              iiorb=tmblarge%orbs%isorb+iorb
              ilrlarge=tmblarge%orbs%inwhichlocreg(iiorb)
              ncnt=tmblarge%lzd%llr(ilrlarge)%wfd%nvctr_c+7*tmblarge%lzd%llr(ilrlarge)%wfd%nvctr_f
              tt=dnrm2(ncnt, tmblarge%psi(ist), 1)
              call dscal(ncnt, 1/tt, tmblarge%psi(ist), 1)
              ist=ist+ncnt
          end do
      end if



      ! Update tmb%confdatarr...
      call update_confdatarr(tmblarge%lzd, tmblarge%orbs, locregCenterTemp, tmb%confdatarr)
      call MLWFnew(iproc, nproc, tmblarge%lzd, tmblarge%orbs, at, tmblarge%op, &
          tmblarge%comon, tmblarge%mad, rxyz, tmb%wfnmd%bs%nit_unitary_loop, kernel, &
          tmb%confdatarr, tmb%lzd%hgrids(1), locregCenterTemp, 3.d0, tmblarge%psi, Umat, locregCenter)


      call check_locregCenters(iproc, tmb%lzd, locregCenter, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3))
      if(tmb%wfnmd%bs%nit_unitary_loop>0) then
          call update_kernel(tmb%orbs%norb, Umat, kernel)
      end if


      if(variable_locregs) then
          call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
          call destroy_new_locregs(iproc, nproc, tmb)
          !!call deallocate_auxiliary_basis_function(subname, tmb%psi, lhphi, lhphiold, lphiold)
          call update_locreg(iproc, nproc, tmblarge%lzd%nlr, locrad, inwhichlocreg_reference, locregCenter, tmblarge%lzd%glr, &
               .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmblarge%orbs, tmb%lzd, tmb%orbs, tmb%op, tmb%comon, &
               tmb%comgp, tmb%comsr, tmb%mad, tmb%collcom)
          call update_ldiis_arrays(tmb, subname, ldiis)
          !!call allocate_auxiliary_basis_function(tmb%orbs%npsidim_orbs, subname, tmb%psi, lhphi, lhphiold, lphiold)     
          call update_auxiliary_basis_function(subname, tmb%orbs%npsidim_orbs, tmb%psi, lhphi, lhphiold, lphiold)     
          call copy_basis_performance_options(tmblarge%wfnmd%bpo, tmb%wfnmd%bpo, subname)
          call copy_orthon_data(tmblarge%orthpar, tmb%orthpar, subname)
          call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmb%orbs%onwhichatom(1), 1)
          tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
      end if

      call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmblarge%psi, tmb%psi)

      call update_confdatarr(tmblarge%lzd, tmblarge%orbs, locregCenter, tmb%confdatarr)

      if(variable_locregs) then
          call vcopy(tmb%orbs%norb, tmblarge%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
          ! this communication is useless, but otherwise the wait in destroy_new_locregs makes problems... to be solved               
          call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &                            
               tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)    
          call destroy_new_locregs(iproc, nproc, tmblarge)
          !!call deallocate_auxiliary_basis_function(subname, tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
          locrad_tmp=factor*locrad
          call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, inwhichlocreg_reference, locregCenter, tmb%lzd%glr, &
               .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmb%orbs, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, &
               tmblarge%comgp, tmblarge%comsr, tmblarge%mad, tmblarge%collcom)
          call update_ldiis_arrays(tmblarge, subname, ldiis)
          !!call allocate_auxiliary_basis_function(tmblarge%orbs%npsidim_orbs, subname, tmblarge%psi, &
          !!     lhphilarge, lhphilargeold, lphilargeold)
          call update_auxiliary_basis_function(subname, tmblarge%orbs%npsidim_orbs, tmblarge%psi, &
               lhphilarge, lhphilargeold, lphilargeold)
          call copy_basis_performance_options(tmb%wfnmd%bpo, tmblarge%wfnmd%bpo, subname)
          call copy_orthon_data(tmb%orthpar, tmblarge%orthpar, subname)
          tmblarge%wfnmd%nphi=tmblarge%orbs%npsidim_orbs
          call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmblarge%orbs%onwhichatom(1), 1)
          locregCenterTemp=locregCenter
      end if


  end if newgradient_if_2


  do_ortho_if2: if(.true.) then
  !do_ortho_if2: if(.not.ldiis%switchSD) then
  !do_ortho_if2: if(ldiis%nwaitortho<=0) then

      if(.not.variable_locregs .or. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
      else
          call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)
      end if
      if(.not.variable_locregs .or. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
          tmbopt => tmb
          lhphiopt => lhphi
          lphioldopt => lphiold
      else
          tmbopt => tmblarge
          lhphiopt => lhphilarge
          lphioldopt => lphilargeold
      end if
      tmbopt%confdatarr => tmb%confdatarr
      if(iproc==0) write(*,*) 'calling orthonormalizeLocalized...'
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
           tmbopt%orbs, tmbopt%op, tmbopt%comon, tmbopt%lzd, &
           tmbopt%mad, tmbopt%collcom, tmbopt%orthpar, tmbopt%wfnmd%bpo, tmbopt%psi, tmbopt%psit_c, tmbopt%psit_f, &
           tmbopt%can_use_transposed)
      !!if(tmbopt%can_use_transposed) then
      !!!if(associated(tmbopt%psit_c)) then
      !!    iall=-product(shape(tmbopt%psit_c))*kind(tmbopt%psit_c)
      !!    deallocate(tmbopt%psit_c, stat=istat)
      !!    call memocc(istat, iall, 'tmbopt%psit_c', subname)
      !!!end if
      !!!if(associated(tmbopt%psit_f)) then
      !!    iall=-product(shape(tmbopt%psit_f))*kind(tmbopt%psit_f)
      !!    deallocate(tmbopt%psit_f, stat=istat)
      !!    call memocc(istat, iall, 'tmbopt%psit_f', subname)
      !!end if
      !!tmbopt%can_use_transposed=.false.

      if(variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
          ! This is not the ideal place for this...
          if(tmbopt%can_use_transposed) then
              tmbopt%can_use_transposed=.false.
              iall = -product(shape(tmbopt%psit_c))*kind(tmbopt%psit_c)
              deallocate(tmbopt%psit_c,stat=istat)
              call memocc(istat,iall,'tmbopt%psit_c',subname)
              iall = -product(shape(tmbopt%psit_f))*kind(tmbopt%psit_f)
              deallocate(tmbopt%psit_f,stat=istat)
              call memocc(istat,iall,'tmbopt%psit_f',subname)
          end if

          ! Optimize the locreg centers and potentially the shape of the basis functions.
          call update_confdatarr(tmblarge%lzd, tmblarge%orbs, locregCenterTemp, tmb%confdatarr)
          call MLWFnew(iproc, nproc, tmblarge%lzd, tmblarge%orbs, at, tmblarge%op, &
               tmblarge%comon, tmblarge%mad, rxyz, tmb%wfnmd%bs%nit_unitary_loop, kernel, &
               tmb%confdatarr, tmb%lzd%hgrids(1), locregCenterTemp, 3.d0, tmblarge%psi, Umat, locregCenter)

          ! Check whether the new locreg centers are ok.
          call check_locregCenters(iproc, tmb%lzd, locregCenter, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3))

          ! Update the kernel if required.
          if(tmb%wfnmd%bs%nit_unitary_loop>0) then                          
              call update_kernel(tmb%orbs%norb, Umat, kernel)
          end if

          if(variable_locregs) then
              call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
              ! this communication is useless, but otherwise the wait in destroy_new_locregs makes problems... to be solved               
              call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &                            
                   tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)    
              call destroy_new_locregs(iproc, nproc, tmb)
              !!call deallocate_auxiliary_basis_function(subname, tmb%psi, lhphi, lhphiold, lphiold)
              call update_locreg(iproc, nproc, tmblarge%lzd%nlr, locrad, &
                   inwhichlocreg_reference, locregCenter, tmblarge%lzd%glr, &
                   .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
                   tmblarge%orbs, tmb%lzd, tmb%orbs, tmb%op, tmb%comon, &
                   tmb%comgp, tmb%comsr, tmb%mad, tmb%collcom)
              call update_ldiis_arrays(tmb, subname, ldiis)
              !!call allocate_auxiliary_basis_function(tmb%orbs%npsidim_orbs, subname, tmb%psi, lhphi, lhphiold, lphiold)
              call update_auxiliary_basis_function(subname, tmb%orbs%npsidim_orbs, tmb%psi, lhphi, lhphiold, lphiold)
              call copy_basis_performance_options(tmblarge%wfnmd%bpo, tmb%wfnmd%bpo, subname)
              call copy_orthon_data(tmblarge%orthpar, tmb%orthpar, subname)
              call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmb%orbs%onwhichatom(1), 1)
              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
          end if


          !!call postCommunicationsPotential(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, tmb%comgp)
          call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
               tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)

          ! Transform back to small locreg
          call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmblarge%psi, tmb%psi)

          ! Update tmb%confdatarr...
          call update_confdatarr(tmblarge%lzd, tmblarge%orbs, locregCenter, tmb%confdatarr)
 
          ! Update the localization region if required.
          if(variable_locregs) then
              call vcopy(tmb%orbs%norb, tmblarge%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
              ! this communication is useless, but otherwise the wait in destroy_new_locregs makes problems... to be solved               
              call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &                            
                   tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)    
              call destroy_new_locregs(iproc, nproc, tmblarge)
              !!call deallocate_auxiliary_basis_function(subname, tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
              locrad_tmp=factor*locrad
              call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, inwhichlocreg_reference, locregCenter, tmb%lzd%glr, &
                   .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
                   tmb%orbs, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, &
                   tmblarge%comgp, tmblarge%comsr, tmblarge%mad, tmblarge%collcom)
              call update_ldiis_arrays(tmblarge, subname, ldiis)
              !!call allocate_auxiliary_basis_function(tmblarge%orbs%npsidim_orbs, subname, tmblarge%psi, &
              !!     lhphilarge, lhphilargeold, lphilargeold)     
              call update_auxiliary_basis_function(subname, tmblarge%orbs%npsidim_orbs, tmblarge%psi, &
                   lhphilarge, lhphilargeold, lphilargeold)     
              call copy_basis_performance_options(tmb%wfnmd%bpo, tmblarge%wfnmd%bpo, subname)
              call copy_orthon_data(tmb%orthpar, tmblarge%orthpar, subname)
              tmblarge%wfnmd%nphi=tmblarge%orbs%npsidim_orbs
              !!allocate(tmblarge%orbs%onwhichatom(tmb%orbs%norb), stat=istat)
              !!call memocc(istat, tmblarge%orbs%onwhichatom, 'tmblarge%orbs%onwhichatom', subname)
              call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmblarge%orbs%onwhichatom(1), 1)
              locregCenterTemp=locregCenter
          end if

          ! Normalize by hand
          ist=1
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              ilr=tmb%orbs%inwhichlocreg(iiorb)
              ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              tt=dnrm2(ncnt, tmb%psi(ist), 1)
              !!write(*,*) 'here iorb,tt',iorb,tt
              call dscal(ncnt, 1/tt, tmb%psi(ist), 1)
              ist=ist+ncnt
          end do

      end if

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
