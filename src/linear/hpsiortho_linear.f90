!> @file
!! H|psi> and orthonormalization (linear scaling version)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
           energy_increased, tmb, lhphiold, overlap_calculated, &
           energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
           energy_only, hpsi_small, hpsi_noprecond)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, it
  type(DFT_wavefunction), target, intent(inout):: tmb
  type(localizedDIISParameters), intent(inout) :: ldiis
  real(kind=8), dimension(tmb%orbs%norb), intent(inout) :: fnrmOldArr
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
  real(kind=8), intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(kind=8), intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
  logical,intent(inout):: overlap_calculated
  type(energy_terms), intent(in) :: energs
  real(kind=8), dimension(:), pointer:: hpsit_c, hpsit_f
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
  logical, intent(in) :: energy_only
  real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: hpsi_noprecond

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, jjorb, i
  integer :: matrixindex_in_compressed
  real(kind=8) :: ddot, tt, gnrmArr
  !real(kind=8) :: eval_zero
  character(len=*), parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8), dimension(:,:), allocatable :: fnrmOvrlpArr, fnrmArr
  real(kind=8), dimension(:), allocatable :: hpsi_conf, hpsi_tmp
  real(kind=8), dimension(:), pointer :: kernel_compr_tmp
  !type(sparseMatrix) :: lagmat
  real(kind=8), dimension(:), allocatable :: prefac
  real(wp), dimension(2) :: garray
  real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now

  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      allocate(hpsi_conf(tmb%npsidim_orbs), stat=istat)
      call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%hpsi, hpsi_conf)
      call timing(iproc,'eglincomms','ON')
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
          call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
          ist=ist+ncount
      end do
      call timing(iproc,'eglincomms','OF')
  end if

  ! by default no quick exit
  energy_increased=.false.


  allocate(hpsittmp_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
  allocate(hpsittmp_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)

  if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     target_function==TARGET_FUNCTION_IS_HYBRID) then

      if(sum(tmb%ham_descr%collcom%nrecvcounts_c)>0) &
          call dcopy(sum(tmb%ham_descr%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmb%ham_descr%collcom%nrecvcounts_f)>0) &
          call dcopy(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call timing(iproc,'eglincomms','ON')
          allocate(kernel_compr_tmp(tmb%linmat%denskern%nvctr), stat=istat)
          call memocc(istat, kernel_compr_tmp, 'kernel_compr_tmp', subname)
          call vcopy(tmb%linmat%denskern%nvctr, tmb%linmat%denskern%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !ii=0
          !do iseg=1,tmb%linmat%denskern%nseg
          !    do jorb=tmb%linmat%denskern%keyg(1,iseg), tmb%linmat%denskern%keyg(2,iseg)
          !        ii=ii+1
          !        iiorb = (jorb-1)/tmb%orbs%norb + 1
          !        jjorb = jorb - (iiorb-1)*tmb%orbs%norb
              do ii=1,tmb%linmat%denskern%nvctr
                      iiorb = tmb%linmat%denskern%orb_from_index(1,ii)
                      jjorb = tmb%linmat%denskern%orb_from_index(2,ii)
                  if(iiorb==jjorb) then
                      tmb%linmat%denskern%matrix_compr(ii)=0.d0
                  else
                      tmb%linmat%denskern%matrix_compr(ii)=kernel_compr_tmp(ii)
                  end if
              end do
          !end do

          ist=1
          do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
              ilr=tmb%orbs%inwhichlocreg(iorb)
              !ii=0
              !do iseg=1,tmb%linmat%denskern%nseg
              !    do jorb=tmb%linmat%denskern%keyg(1,iseg), tmb%linmat%denskern%keyg(2,iseg)
              !        ii=ii+1
              !        iiorb = (jorb-1)/tmb%orbs%norb + 1
              !        jjorb = jorb - (iiorb-1)*tmb%orbs%norb
              do ii=1,tmb%linmat%denskern%nvctr
                      iiorb = tmb%linmat%denskern%orb_from_index(1,ii)
                      jjorb = tmb%linmat%denskern%orb_from_index(2,ii)
                      if(iiorb==jjorb .and. iiorb==iorb) then
                          ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
                          call dscal(ncount, kernel_compr_tmp(ii), tmb%hpsi(ist), 1)
                          ist=ist+ncount
                      end if
                  !end do
              end do
          end do
          call timing(iproc,'eglincomms','OF')
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
          ! copy correct kernel back
          call vcopy(tmb%linmat%denskern%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%denskern%matrix_compr(1), 1)
          iall=-product(shape(kernel_compr_tmp))*kind(kernel_compr_tmp)
          deallocate(kernel_compr_tmp, stat=istat)
          call memocc(istat, iall, 'kernel_compr_tmp', subname)
      else
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      end if
  end if

  iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
  deallocate(hpsittmp_c, stat=istat)
  call memocc(istat, iall, 'hpsittmp_c', subname)
  iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
  deallocate(hpsittmp_f, stat=istat)
  call memocc(istat, iall, 'hpsittmp_f', subname)


  ! make lagmat a structure with same sparsity as h
  !call nullify_sparsematrix(lagmat)
  !call sparse_copy_pattern(tmb%linmat%ham, lagmat, iproc, subname)
  !allocate(lagmat%matrix_compr(lagmat%nvctr), stat=istat)
  !call memocc(istat, lagmat%matrix_compr, 'lagmat%matrix_compr', subname)

  call orthoconstraintNonorthogonal(iproc, nproc, tmb%ham_descr%lzd, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%orbs, tmb%ham_descr%collcom, tmb%orthpar, correction_orthoconstraint, tmb%linmat, tmb%ham_descr%psi, tmb%hpsi, &
       tmb%linmat%ham, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, hpsit_c, hpsit_f, tmb%ham_descr%can_use_transposed, &
       overlap_calculated)

  call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%orbs, tmb%hpsi, hpsi_small)

  if (present(hpsi_noprecond)) call dcopy(tmb%npsidim_orbs, hpsi_small, 1, hpsi_noprecond, 1)

  ! Calculate trace (or band structure energy, resp.)
  trH=0.d0
  call timing(iproc,'eglincomms','ON')
  do iorb=1,tmb%orbs%norb
     ii=matrixindex_in_compressed(tmb%linmat%ham,iorb,iorb)
     trH = trH + tmb%linmat%ham%matrix_compr(ii)
  end do
  call timing(iproc,'eglincomms','OF')
  !call deallocate_sparseMatrix(lagmat,subname)

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

  allocate(fnrmOvrlpArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb,2), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  call timing(iproc,'eglincomms','ON')
  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, hpsi_small(ist), 1, lhphiold(ist), 1)
      fnrmArr(iorb,1)=ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
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
          ! apply a threshold so that alpha never goes below around 1.d-2
          else if (alpha(iorb)>1.7d-3) then
              alpha(iorb)=alpha(iorb)*.6d0
          end if
          !!alpha(iorb)=min(alpha(iorb),1.5d0)
      end if
  end do
  call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)

  iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
  deallocate(fnrmOvrlpArr, stat=istat)
  call memocc(istat, iall, 'fnrmOvrlpArr', subname)

  iall=-product(shape(fnrmArr))*kind(fnrmArr)
  deallocate(fnrmArr, stat=istat)
  call memocc(istat, iall, 'fnrmArr', subname)

  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  alpha_mean=tt/dble(tmb%orbs%norb)
  alpha_max=maxval(alpha)
  call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)

  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmb%npsidim_orbs, hpsi_small, 1, lhphiold, 1)
  call timing(iproc,'eglincomms','OF')

  ! if energy has increased or we only wanted to calculate the energy, not gradient, we can return here
  ! rather than calculating the preconditioning for nothing
  if ((energy_increased .or. energy_only) .and. target_function/=TARGET_FUNCTION_IS_HYBRID) return

  ! Precondition the gradient.
  if(iproc==0) then
      write(*,'(a)',advance='no') 'Preconditioning... '
  end if
 

  !if (target_function==TARGET_FUNCTION_IS_HYBRID) then
  !    allocate(hpsi_tmp(tmb%npsidim_orbs), stat=istat)
  !    call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
  !end if
  !ist=1
  !do iorb=1,tmb%orbs%norbp
  !    iiorb=tmb%orbs%isorb+iorb
  !    ilr = tmb%orbs%inWhichLocreg(iiorb)
  !    ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !    if(target_function==TARGET_FUNCTION_IS_HYBRID) then
  !        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
  !        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
  !        do i=ist,ist+ncnt-1
  !            hpsi_tmp(i)=tt*hpsi_conf(i)
  !        end do
  !        call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !             tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !             nit_precond, hpsi_tmp(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !             tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
  !        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
  !        call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !             tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !             nit_precond, hpsi_small(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !             0.d0, iorb, eval_zero)
  !        call daxpy(ncnt, 1.d0, hpsi_tmp(ist), 1, hpsi_small(ist), 1)
  !    else
  !    !    call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !    !         tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !    !         nit_precond, hpsi_small(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !    !         tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
  !    end if
  !    ist=ist+ncnt
  !end do

  if(target_function==TARGET_FUNCTION_IS_HYBRID) then
     allocate(hpsi_tmp(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, hpsi_tmp, 'hpsi_tmp', subname)
     ist=1
     do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        ilr = tmb%orbs%inWhichLocreg(iiorb)
        ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
        do i=ist,ist+ncnt-1
           hpsi_tmp(i)=tt*hpsi_conf(i)
        end do
        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)

        ist=ist+ncnt
     end do

     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_tmp,tmb%confdatarr,gnrm,gnrm_zero)

     ! temporarily turn confining potential off...
     allocate(prefac(tmb%orbs%norbp),stat=istat)
     call memocc(istat, prefac, 'prefac', subname)
     prefac(:)=tmb%confdatarr(:)%prefac
     tmb%confdatarr(:)%prefac=0.0d0
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero) ! prefac should be zero
     call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_tmp(1), 1, hpsi_small(1), 1)
     ! ...revert back to correct value
     tmb%confdatarr(:)%prefac=prefac

     iall=-product(shape(prefac))*kind(prefac)
     deallocate(prefac, stat=istat)
     call memocc(istat, iall, 'prefac', subname)
     iall=-product(shape(hpsi_conf))*kind(hpsi_conf)
     deallocate(hpsi_conf, stat=istat)
     call memocc(istat, iall, 'hpsi_conf', subname)
     iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
     deallocate(hpsi_tmp, stat=istat)
     call memocc(istat, iall, 'hpsi_tmp', subname)
  else
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero)
   end if

  !sum over all the partial residues
  if (nproc > 1) then
      garray(1)=gnrm
      garray(2)=gnrm_zero
      call mpiallred(garray(1),2,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
      gnrm     =garray(1)
      gnrm_zero=garray(2)
  end if

  !if (target_function==TARGET_FUNCTION_IS_HYBRID) then
  !    iall=-product(shape(hpsi_conf))*kind(hpsi_conf)
  !    deallocate(hpsi_conf, stat=istat)
  !    call memocc(istat, iall, 'hpsi_conf', subname)
  !    iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
  !    deallocate(hpsi_tmp, stat=istat)
  !    call memocc(istat, iall, 'hpsi_tmp', subname)
  !end if

  if(iproc==0) then
      write(*,'(a)') 'done.'
  end if

  call timing(iproc,'eglincomms','ON')
  ist=1
  gnrm_old=gnrm
  gnrm=0.d0
  gnrmMax=0.d0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      gnrmArr=ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
      gnrm=gnrm+gnrmArr
      if(gnrmArr>gnrmMax) gnrmMax=gnrmArr
      ist=ist+ncount
  end do

  call mpiallred(gnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(gnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  gnrm=sqrt(gnrm/dble(tmb%orbs%norb))
  gnrmMax=sqrt(gnrmMax)
  if (iproc==0) write(*,'(a,3es16.6)') 'AFTER: gnrm, gnrmmax, gnrm/gnrm_old',gnrm,gnrmmax,gnrm/gnrm_old
  call timing(iproc,'eglincomms','OF')


end subroutine calculate_energy_and_gradient_linear



subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb,  &
           lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff)
  use module_base
  use module_types
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(localizedDIISParameters), intent(inout) :: ldiis
  type(DFT_wavefunction), target,intent(inout) :: tmb
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: lphiold
  real(kind=8), intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha, alphaDIIS
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: psidiff
  logical, intent(in) :: ortho
  
  ! Local variables
  integer :: istat, iall, i
  character(len=*), parameter :: subname='hpsitopsi_linear'

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
  if (present(psidiff)) call dcopy(tmb%npsidim_orbs, tmb%psi, 1, psidiff, 1)
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, tmb, ldiis, alpha, hpsi_small)
  else
      if(iproc==0) write(*,'(1x,a)') 'no improvement of the orbitals, recalculate gradient'
  end if
  if (present(psidiff)) then
      do i=1,tmb%npsidim_orbs
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

  if(.not.ldiis%switchSD.and.ortho) then
      if(iproc==0) then
           write(*,'(1x,a)',advance='no') 'Orthonormalization... '
      end if
      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)
  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear

