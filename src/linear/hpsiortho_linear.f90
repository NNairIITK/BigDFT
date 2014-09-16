!> @file
!! H|psi> and orthonormalization (linear scaling version)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           ldiis, fnrmOldArr, fnrm_old, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
           energy_increased, tmb, lhphiold, overlap_calculated, &
           energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
           hpsi_small, experimental_mode, correction_co_contra, hpsi_noprecond, &
           norder_taylor, max_inversion_error, method_updatekernel, precond_convol_workarrays, precond_workarrays, cdft)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  use communications, only: transpose_localized
  use sparsematrix_base, only: matrices, matrices_null, deallocate_matrices, &
                               sparsematrix_malloc_ptr, assignment(=), SPARSE_FULL
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: transform_sparse_matrix, orb_from_index
  use constrained_dft, only: cdft_data
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, it, method_updatekernel
  integer,intent(inout) :: norder_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(DFT_wavefunction), target, intent(inout):: tmb
  type(localizedDIISParameters), intent(inout) :: ldiis
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: fnrmOldArr
  real(kind=8),intent(inout) :: fnrm_old
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
  real(kind=8), intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(kind=8), intent(in):: trHold
  logical,intent(out) :: energy_increased
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
  logical, intent(inout):: overlap_calculated
  type(energy_terms), intent(in) :: energs
  real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c) :: hpsit_c
  real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f) :: hpsit_f
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
  logical, intent(in) :: experimental_mode, correction_co_contra
  real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs),intent(out) :: hpsi_noprecond
  type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
  type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
  type(cdft_data),intent(in),optional :: cdft

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, jjorb, i
  integer :: lwork, info, ishift,ispin, iseg
  real(kind=8) :: ddot, tt, fnrmOvrlp_tot, fnrm_tot, fnrmold_tot, tt2, trkw
  character(len=*), parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8), dimension(:), allocatable :: hpsi_conf
  real(kind=8), dimension(:), pointer :: kernel_compr_tmp
  real(kind=8), dimension(:), allocatable :: prefac
  real(kind=8),dimension(2) :: reducearr
  integer,dimension(2) :: irowcol
  real(wp), dimension(2) :: garray
  real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
  type(matrices) :: matrixm

  call f_routine(id='calculate_energy_and_gradient_linear')

  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      hpsi_conf = f_malloc(tmb%npsidim_orbs,id='hpsi_conf')
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%hpsi, hpsi_conf)
      call timing(iproc,'buildgrad_mcpy','ON')
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
          call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
          ist=ist+ncount
      end do
      call timing(iproc,'buildgrad_mcpy','OF')
  end if

  ! by default no quick exit
  energy_increased=.false.

  if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
      !!do i=1,size(hpsit_c)
      !!    write(4231,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_c(i), hpsit_c(i)', i, tmb%ham_descr%psit_c(i), hpsit_c(i)
      !!end do
      !!do i=1,size(hpsit_f)
      !!    write(4232,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_f(i), hpsit_f(i)', i, tmb%ham_descr%psit_f(i), hpsit_f(i)
      !!end do
  end if

  hpsittmp_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsittmp_c')
  hpsittmp_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsittmp_f')

  if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     target_function==TARGET_FUNCTION_IS_HYBRID) then
     call build_gradient(iproc, nproc, tmb, target_function, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f)
      
     !!! if(tmb%ham_descr%collcom%ndimind_c>0) &
     !!!     call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
     !!! if(tmb%ham_descr%collcom%ndimind_f>0) &
     !!!     call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)

     !!! if (target_function==TARGET_FUNCTION_IS_HYBRID) then
     !!!     call timing(iproc,'buildgrad_mcpy','ON')
     !!!     kernel_compr_tmp = sparsematrix_malloc_ptr(tmb%linmat%l,iaction=SPARSE_FULL,id='kernel_compr_tmp')
     !!!     call vcopy(tmb%linmat%l%nvctr*tmb%linmat%l%nspin, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
     !!!     do ispin=1,tmb%linmat%l%nspin
     !!!         ishift=(ispin-1)*tmb%linmat%l%nvctr
     !!!         do iseg=1,tmb%linmat%l%nseg
     !!!             ii=tmb%linmat%l%keyv(iseg)
     !!!             do i=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
     !!!                 irowcol = orb_from_index(tmb%linmat%l, i)
     !!!                 if(irowcol(1)==irowcol(2)) then
     !!!                     tmb%linmat%kernel_%matrix_compr(ii+ishift)=0.d0
     !!!                 else
     !!!                     tmb%linmat%kernel_%matrix_compr(ii+ishift)=kernel_compr_tmp(ii+ishift)
     !!!                 end if
     !!!                 ii=ii+1
     !!!             end do
     !!!         end do
     !!!     end do

     !!!     ist=1
     !!!     do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
     !!!         ilr=tmb%orbs%inwhichlocreg(iorb)
     !!!         if (tmb%orbs%spinsgn(iorb)>0.d0) then
     !!!             ispin=1
     !!!         else
     !!!             ispin=2
     !!!         end if
     !!!         do iseg=1,tmb%linmat%l%nseg
     !!!             ii=tmb%linmat%l%keyv(iseg)
     !!!             do i=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
     !!!                 irowcol = orb_from_index(tmb%linmat%l, i)
     !!!                 ishift=(ispin-1)*tmb%linmat%l%nvctr
     !!!                 if(irowcol(1)==irowcol(2) .and. irowcol(1)==iorb) then
     !!!                     ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
     !!!                     write(*,*) 'iorb, ii, ishift, ist', iorb, ii, ishift, ist
     !!!                     call dscal(ncount, kernel_compr_tmp(ii+ishift), tmb%hpsi(ist), 1)
     !!!                     ist=ist+ncount
     !!!                 end if
     !!!                 ii=ii+1
     !!!             end do
     !!!         end do
     !!!     end do
     !!!     call timing(iproc,'buildgrad_mcpy','OF')
     !!!     call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
     !!!          tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
     !!!     call build_linear_combination_transposed(tmb%ham_descr%collcom, &
     !!!          tmb%linmat%l, tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
     !!!     ! copy correct kernel back
     !!!     call vcopy(tmb%linmat%l%nvctr*tmb%linmat%l%nspin, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
     !!!     call f_free_ptr(kernel_compr_tmp)
     !!! else
     !!!     call build_linear_combination_transposed(tmb%ham_descr%collcom, &
     !!!          tmb%linmat%l, tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
     !!! end if
  end if

  ! WARNING: TO BE CHECKED!!!!
  ! For the non polarized case, a factor of two is already included in the
  ! kernel. Therefore explicitely add this factor for the polarized case 
  ! in order to make the two cases analogous.
  if (tmb%linmat%l%nspin==2 .and. target_function==TARGET_FUNCTION_IS_ENERGY) then
      if (iproc==0) call yaml_warning('multiply the gradient by 2.0, check this!')
      hpsit_c=2.d0*hpsit_c
      hpsit_f=2.d0*hpsit_f
  end if

  !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
  !!    ist=0
  !!    do ispin=1,tmb%linmat%l%nspin
  !!        do i=1,tmb%linmat%l%nvctr
  !!            ist=ist+1
  !!            write(4241,'(a,3i8,es14.6)') 'ispin, i, ist, tmb%linmat%kernel_%matrix_compr(ist)', ispin, i, ist, tmb%linmat%kernel_%matrix_compr(ist)
  !!        end do
  !!    end do
  !!end if


  if (correction_co_contra) then
      !@NEW correction for contra / covariant gradient

      if ((method_updatekernel/=UPDATE_BY_FOE .and. method_updatekernel/=UPDATE_BY_RENORMALIZATION) &
           .or. target_function/=TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.

          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
               tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
           !write(*,*) 'corr cocontra: sums(ovrlp)', &
           !    sum(tmb%linmat%ovrlp_%matrix_compr(1:tmb%linmat%s%nvctr)), sum(tmb%linmat%ovrlp_%matrix_compr(tmb%linmat%s%nvctr+1:2*tmb%linmat%s%nvctr))
      end if
      call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
      call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)

      ! Transform to the larger sparse region in order to be compatible with tmb%ham_descr%collcom.
      ! To this end use ham_.
      call transform_sparse_matrix(tmb%linmat%s, tmb%linmat%m, &
           tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ham_%matrix_compr, 'small_to_large')
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%m, tmb%linmat%ham_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)

      !@END NEW correction for contra / covariant gradient
  end if


  call f_free_ptr(hpsittmp_c)
  call f_free_ptr(hpsittmp_f)


  ! Calculate the overlap matrix if necessary
  if (.not.correction_co_contra) then
      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
           tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      tmb%can_use_transposed=.true.
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
  end if


  !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
  !!    ist=0
  !!    do iorb=1,tmb%orbs%norbp
  !!        iiorb=tmb%orbs%isorb+iorb
  !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!        ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
  !!        do i=1,ncount
  !!            ist=ist+1
  !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
  !!                write(4211,'(a,2i10,f8.1,2es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%ham_descr%psi(ist), tmb%hpsi(ist)
  !!            else
  !!                write(4212,'(a,2i10,f8.1,2es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%ham_descr%psi(ist), tmb%hpsi(ist)
  !!            end if
  !!        end do
  !!    end do
  !!    do i=1,size(hpsit_c)
  !!        write(4221,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_c(i), hpsit_c(i)', i, tmb%ham_descr%psit_c(i), hpsit_c(i)
  !!    end do
  !!    do i=1,size(hpsit_f)
  !!        write(4222,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_f(i), hpsit_f(i)', i, tmb%ham_descr%psit_f(i), hpsit_f(i)
  !!    end do
  !!end if

  call orthoconstraintNonorthogonal(iproc, nproc, tmb%ham_descr%lzd, &
       tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%orbs, tmb%ham_descr%collcom, tmb%orthpar, correction_orthoconstraint, &
       tmb%linmat, tmb%ham_descr%psi, tmb%hpsi, &
       tmb%linmat%m, tmb%linmat%ham_, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, &
       hpsit_c, hpsit_f, tmb%ham_descr%can_use_transposed, &
       overlap_calculated, experimental_mode, norder_taylor, max_inversion_error, &
       tmb%npsidim_orbs, tmb%lzd, hpsi_noprecond)

  !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
  !!    ist=0
  !!    do iorb=1,tmb%orbs%norbp
  !!        iiorb=tmb%orbs%isorb+iorb
  !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!        ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
  !!        do i=1,ncount
  !!            ist=ist+1
  !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
  !!                write(4201,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
  !!            else
  !!                write(4202,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
  !!            end if
  !!        end do
  !!    end do
  !!end if


  call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%orbs, tmb%hpsi, hpsi_small)

  !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
  !!    ist=0
  !!    do iorb=1,tmb%orbs%norbp
  !!        iiorb=tmb%orbs%isorb+iorb
  !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!        ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !!        do i=1,ncount
  !!            ist=ist+1
  !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
  !!                write(4301,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
  !!            else
  !!                write(4302,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
  !!            end if
  !!        end do
  !!    end do
  !!end if


  ! Calculate trace (or band structure energy, resp.)
  call calculate_trace()

  if (present(cdft).and..false.) then
    !only correct energy not gradient for now
    !can give tmb%orbs twice as ksorbs is only used for recalculating the kernel
    call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,cdft%weight_matrix, &
           tmb%linmat%kernel_,cdft%weight_matrix_,trkw,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
    !cdft%charge is always constant so could in theory be ignored as in optimize_coeffs
    trH = trH + cdft%lag_mult*(trkw - cdft%charge)
    !if (iproc==0) print*,'trH,trH+V(trkw-N)',trH-cdft%lag_mult*(trkw - cdft%charge),trH
  end if


  ! WARNING: TO BE CHECKED!!!!
  ! For the polarized case, a factor of two was already multiplied to the
  ! gradient (see above). Therefore now undo this again for the band structure
  ! energy in order to get the analogous result to the non-polarized case.
  if (tmb%linmat%l%nspin==2 .and. target_function==TARGET_FUNCTION_IS_ENERGY) then
      if (iproc==0) call yaml_warning('divide the band stucture energy by 2.0, check this!')
      trH=0.5d0*trH
  end if

  ! trH is now the total energy (name is misleading, correct this)
  ! Multiply by 2 because when minimizing trace we don't have kernel
  if(tmb%orbs%nspin==1 .and. target_function==TARGET_FUNCTION_IS_TRACE) trH=2.d0*trH
  trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp


  ! Determine whether the target function is increasing
  !if(.not. ldiis%switchSD .and. ldiis%isx==0) then
  if(.not. ldiis%switchSD) then
      if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
          !!if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
          !!    'WARNING: the target function is larger than its minimal value reached so far:',trH,' > ', ldiis%trmin
          if (iproc==0) then
              call yaml_newline()
              call yaml_warning('target function larger than its minimal value reached so far, &
                  &D='//trim(yaml_toa(trH-ldiis%trmin,fmt='(1es10.3)')))!//'. &
                  !&Decrease step size and restart with previous TMBs')
          end if
          !if(iproc==0) write(*,'(1x,a)') 'Decrease step size and restart with previous TMBs'
          energy_increased=.true.
      end if
  end if



  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.
  fnrm=0.d0
  fnrmMax=0.d0
  fnrmOvrlp_tot=0.d0
  fnrmOld_tot=0.d0
  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      tt = ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
      fnrm = fnrm + tt
      if(tt>fnrmMax) fnrmMax=tt
      if(it>1) then
          tt2=ddot(ncount, hpsi_small(ist), 1, lhphiold(ist), 1)
          fnrmOvrlp_tot = fnrmOvrlp_tot + tt2
          if(ldiis%isx==0 .and. .not.ldiis%switchSD) then
              ! Adapt step size for the steepest descent minimization.
              if (.not.experimental_mode) then
                  tt2=tt2/sqrt(tt*fnrmOldArr(iorb))
                  ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
                  if(tt2>.6d0 .and. trH<trHold .and. alpha(iorb)<1.8d0) then
                      alpha(iorb)=alpha(iorb)*1.1d0
                  else if (alpha(iorb)>1.7d-3) then
                      alpha(iorb)=alpha(iorb)*.6d0
                  end if
              end if
          end if
      end if
      fnrmOldArr(iorb)=tt
      ist=ist+ncount
  end do
  
  if (nproc > 1) then
     call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm)
     call mpiallred(fnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if


  if (experimental_mode .and. it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      !do iorb=1,tmb%orbs%norbp
      !    fnrmOld_tot=fnrmOld_tot+fnrmOldArr(iorb)
      !end do
      if (nproc>1) then
          call mpiallred(fnrmOvrlp_tot, 1, mpi_sum, bigdft_mpi%mpi_comm)
      end if
      !tt2=fnrmOvrlp_tot/sqrt(fnrm*fnrmOld_tot)
      tt2=fnrmOvrlp_tot/sqrt(fnrm*fnrm_old)
      ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
      if(tt2>.6d0 .and. trH<trHold .and. alpha(1)<1.8d0) then ! take alpha(1) since the value is the same for all
          alpha(:)=alpha(:)*1.1d0
      else if (alpha(1)>1.7d-3) then
          alpha(:)=alpha(:)*.6d0
      end if
  end if

  fnrm_old=fnrm ! This value will be used in th next call to this routine

  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)



  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  alpha_max=maxval(alpha)
  if (nproc > 1) then
     call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm)
     call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
  alpha_mean=tt/dble(tmb%orbs%norb)

  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call vcopy(tmb%npsidim_orbs, hpsi_small(1), 1, lhphiold(1), 1)
  call timing(iproc,'buildgrad_mcpy','OF')

  ! if energy has increased or we only wanted to calculate the energy, not gradient, we can return here
  ! rather than calculating the preconditioning for nothing
  ! Do this only for steepest descent
  if ((energy_increased) .and. target_function/=TARGET_FUNCTION_IS_HYBRID .and. ldiis%isx==0) return



  if(target_function==TARGET_FUNCTION_IS_HYBRID) then
     ist=1
     do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        ilr = tmb%orbs%inWhichLocreg(iiorb)
        ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
        call dscal(ncnt, tt, hpsi_conf(ist), 1)

        ist=ist+ncnt
     end do

     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_conf,tmb%confdatarr,gnrm,gnrm_zero, &
          precond_convol_workarrays, precond_workarrays)

     ! temporarily turn confining potential off...
     prefac = f_malloc(tmb%orbs%norbp,id='prefac')
     prefac(:)=tmb%confdatarr(:)%prefac
     tmb%confdatarr(:)%prefac=0.0d0
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero, & ! prefac should be zero
          precond_convol_workarrays, precond_workarrays)
     call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_conf(1), 1, hpsi_small(1), 1)
     ! ...revert back to correct value
     tmb%confdatarr(:)%prefac=prefac

     call f_free(prefac)
     call f_free(hpsi_conf)
  else
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero,&
          precond_convol_workarrays, precond_workarrays)
  end if

  !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
  !!    ist=0
  !!    do iorb=1,tmb%orbs%norbp
  !!        iiorb=tmb%orbs%isorb+iorb
  !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!        ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !!        do i=1,ncount
  !!            ist=ist+1
  !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
  !!                write(4401,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
  !!            else
  !!                write(4402,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
  !!            end if
  !!        end do
  !!    end do
  !!end if

  if (iproc==0) then
      call yaml_map('Preconditioning',.true.)
  end if

  call f_release_routine()


  contains

    subroutine calculate_trace()
      implicit none

      ! Local variables
      integer :: iorb, iiorb

      call f_routine(id='calculate_trace')
      call timing(iproc,'calctrace_comp','ON')

      trH=0.d0
      do iorb=1,tmb%orbs%norbp
         iiorb=tmb%orbs%isorb+iorb
         ii=matrixindex_in_compressed(tmb%linmat%m,iiorb,iiorb)
         trH = trH + tmb%linmat%ham_%matrix_compr(ii)
      end do
      call timing(iproc,'calctrace_comp','OF')
      call timing(iproc,'calctrace_comm','ON')
      if (nproc>1) then
          call mpiallred(trH, 1, mpi_sum, bigdft_mpi%mpi_comm)
      end if
      call timing(iproc,'calctrace_comm','OF')

      call f_release_routine()

    end subroutine calculate_trace

end subroutine calculate_energy_and_gradient_linear


subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_residue_ks
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: uncompress_matrix
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, num_extra
  type(dft_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: ksorbs
  real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f

  integer :: iorb, istat, iall, ierr
  real(kind=8) :: ksres_sum
  real(kind=8), dimension(:), allocatable :: ksres
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, grad_coeff
  type(sparse_matrix) :: grad_ovrlp
  type(matrices) :: grad_ovrlp_
  character(len=256) :: subname='calculate_residue_ks'


  call f_routine(id='calculate_residue_ks')

  ! want to calculate the residue of the KS states here, not just the tmbs
  ! for now just occupied, eventually would want occupied+num_extra
  ! probably better to calculate c_i^a|g_a> first but copying and pasting for now (INEFFICIENT but just testing)
  !!if(associated(hpsit_c)) then
  !!    iall=-product(shape(hpsit_c))*kind(hpsit_c)
  !!    deallocate(hpsit_c, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_c', subname)
  !!end if
  !!if(associated(hpsit_f)) then
  !!    iall=-product(shape(hpsit_f))*kind(hpsit_f)
  !!    deallocate(hpsit_f, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_f', subname)
  !!end if
  !!allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
  !!call memocc(istat, hpsit_c, 'hpsit_c', subname)
  !!allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  !!call memocc(istat, hpsit_f, 'hpsit_f', subname)

  ! should already be done in orthoconstraintnonorthogonal so can use directly
  !call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!can_use_transposed=.true.

  !call nullify_sparse_matrix(grad_ovrlp)
  grad_ovrlp=sparse_matrix_null()
  call sparse_copy_pattern(tmb%linmat%m, grad_ovrlp, iproc, subname)
  !grad_ovrlp%matrix_compr=f_malloc_ptr(grad_ovrlp%nvctr,id='grad_ovrlp%matrix_compr')
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%m, allocate_full=.false., &
       matname='grad_ovrlp_', mat=grad_ovrlp_)

  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, hpsit_c, hpsit_c, &
       hpsit_f, hpsit_f, tmb%linmat%m, grad_ovrlp_)
  !! This can then be deleted if the transition to the new type has been completed.
  !grad_ovrlp%matrix_compr=grad_ovrlp_%matrix_compr


  grad_coeff = f_malloc((/ tmb%orbs%norb, tmb%orbs%norb /),id='grad_coeff')
  coeff_tmp = f_malloc((/ tmb%orbs%norbp, max(tmb%orbs%norb, 1) /),id='coeff_tmp')

  !grad_ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_ovrlp%matrix')
  call uncompress_matrix(iproc,grad_ovrlp,grad_ovrlp_%matrix_compr,grad_ovrlp_%matrix)

  ! can change this so only go up to ksorbs%norb...
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, 1.d0, grad_ovrlp_%matrix(tmb%orbs%isorb+1,1,1), &
          tmb%orbs%norb, tmb%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp, tmb%orbs%norbp)
     call dgemm('t', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norbp, 1.d0, tmb%coeff(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, coeff_tmp, tmb%orbs%norbp, 0.d0, grad_coeff, tmb%orbs%norb)
  else
     call to_zero(tmb%orbs%norb**2, grad_coeff(1,1))
  end if

  call f_free(coeff_tmp)
  call deallocate_matrices(grad_ovrlp_)

  if (nproc>1) then
      call mpiallred(grad_coeff(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  ! now calculate sqrt(<g_i|g_i>) and mean value
  ksres = f_malloc(ksorbs%norb+num_extra,id='ksres')
  
  ksres_sum=0.0d0
  do iorb=1,ksorbs%norb+num_extra
    ksres(iorb)=dsqrt(grad_coeff(iorb,iorb))
    !ksres_sum=ksres_sum+ksres(iorb)
    ksres_sum=ksres_sum+grad_coeff(iorb,iorb)
    if (iproc==0) write(*,*) 'KS residue',iorb,ksres(iorb)!,tmb%orbs%occup(iorb)
  end do
  if (iproc==0) write(*,*) 'Average KS residue',sqrt(ksres_sum/real(ksorbs%norb+num_extra,gp))


  !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
  !     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
  ! calculate Tr[Kg]  (not recalculating kernel as don't have the correct occs here)
  !call calculate_kernel_and_energy(iproc,nproc,denskern,grad_coeff,ksres_sum,tmb%coeff,orbs,tmb%orbs,.true.)
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='grad_ovrlp_', mat=grad_ovrlp_)
  !grad_ovrlp_%matrix_compr=grad_ovrlp%matrix_compr
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,grad_ovrlp,&
       tmb%linmat%kernel_, grad_ovrlp_, &
       ksres_sum,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
  call deallocate_matrices(grad_ovrlp_)
  if (iproc==0) write(*,*) 'KS residue from trace',dsqrt(ksres_sum)/real(tmb%orbs%norb,gp) ! should update normalization as would only be occ here not extra?

  call deallocate_sparse_matrix(grad_ovrlp)

  call f_free(grad_coeff)
  call f_free(ksres)

  call f_release_routine()

end subroutine calculate_residue_ks


subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb,  &
           lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
           experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => hpsitopsi_linear
  use communications, only: transpose_localized, untranspose_localized
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(localizedDIISParameters), intent(inout) :: ldiis
  type(DFT_wavefunction), target,intent(inout) :: tmb
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: lphiold
  real(kind=8), intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha, alphaDIIS
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: psidiff
  logical, intent(in) :: ortho, experimental_mode
  real(kind=8),intent(out) :: trH_ref
  real(kind=8),dimension(tmb%linmat%l%nvctr),intent(out) :: kernel_best
  logical,intent(out) :: complete_reset

  ! Local variables
  integer :: i, iorb, ilr, ist, iiorb, ncount
  character(len=*), parameter :: subname='hpsitopsi_linear'
  real(kind=8), dimension(:), allocatable :: norm
  real(kind=8) :: ddot, dnrm2, tt

  call f_routine(id='hpsitopsi_linear')

  call DIISorSD(iproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold, trH_ref, kernel_best, complete_reset)
  if(iproc==0) then
      call yaml_newline()
      call yaml_mapping_open('Optimization',flow=.true.)
      if(ldiis%isx>0) then
          call yaml_map('algorithm','DIIS')
          call yaml_map('history length',ldiis%isx)
          call yaml_map('consecutive failures',ldiis%icountDIISFailureCons)
          call yaml_map('total failures',ldiis%icountDIISFailureTot)
      else
          call yaml_map('algorithm','SD')
          call yaml_map('mean alpha',alpha_mean,fmt='(es9.3)')
          call yaml_map('max alpha',alpha_max,fmt='(es9.3)')
          call yaml_map('consecutive successes',ldiis%icountSDSatur)
      end if
      call yaml_mapping_close()
      call yaml_newline()
  end if

  ! Improve the orbitals, depending on the choice made above.
  if (present(psidiff)) call vcopy(tmb%npsidim_orbs, tmb%psi(1), 1, psidiff(1), 1)
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, nproc, tmb, tmb%linmat%s%nspin, ldiis, alpha, hpsi_small, experimental_mode)
  else
      if (iproc==0) then
          call yaml_warning('no improvement of the orbitals, recalculate gradient')
          call yaml_newline()
      end if
  end if

  ! The transposed quantities can now not be used any more...
  if(tmb%can_use_transposed) then
      tmb%can_use_transposed=.false.
  end if


  if (.not.ortho .and. iproc==0) then
      call yaml_map('Orthogonalization',.false.)
  end if

  if(.not.ldiis%switchSD.and.ortho) then
      if (present(psidiff)) then
          do i=1,tmb%npsidim_orbs
              psidiff(i)=tmb%psi(i)-psidiff(i)
          end do 
      end if

      call orthonormalizeLocalized(iproc, nproc, order_taylor, max_inversion_error, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed, tmb%foe_obj)
      if (iproc == 0) then
          call yaml_map('Orthogonalization',.true.)
      end if
  else if (experimental_mode .or. .not.ldiis%switchSD) then
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=dnrm2(ncount, tmb%psi(ist), 1)
          tt=1.d0/tt
          call dscal(ncount, tt, tmb%psi(ist), 1)
          ist=ist+ncount
      end do
      if (iproc == 0) then
          call yaml_map('Normalization',.true.)
          !call yaml_map('Normalization',.false.)
      end if
      if (present(psidiff)) then
          do i=1,tmb%npsidim_orbs
              psidiff(i)=tmb%psi(i)-psidiff(i)
          end do 
      end if
  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

  call f_release_routine

end subroutine hpsitopsi_linear



subroutine build_gradient(iproc, nproc, tmb, target_function, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f)
  use module_base
  use module_types
  use sparsematrix_base, only: sparsematrix_malloc_ptr, SPARSE_FULL, assignment(=)
  use sparsematrix, only: orb_from_index
  use communications, only: transpose_localized
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, target_function
  type(DFT_wavefunction),intent(inout) :: tmb
  real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(inout) :: hpsit_c
  real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(inout) :: hpsit_f
  real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(out) :: hpsittmp_c !<workarray
  real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(out) :: hpsittmp_f !<workarray

  ! Local variables
  integer :: ispin, ishift, iseg, ii, i, ist, ncount, iorb, ilr, isegstart, isegend, ierr, iiorb
  integer,dimension(2) :: irowcol
  real(kind=8),dimension(:),pointer :: kernel_compr_tmp
  real(kind=8),dimension(:),pointer :: matrix_local
  integer,parameter :: ALLGATHERV=51, GET=52
  integer,parameter :: comm_strategy=GET

      call f_routine(id='build_gradient')
      call timing(iproc,'buildgrad_mcpy','ON')

      if(tmb%ham_descr%collcom%ndimind_c>0) &
          call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(tmb%ham_descr%collcom%ndimind_f>0) &
          call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          kernel_compr_tmp = sparsematrix_malloc_ptr(tmb%linmat%l,iaction=SPARSE_FULL,id='kernel_compr_tmp')
          call vcopy(tmb%linmat%l%nvctr*tmb%linmat%l%nspin, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          isegstart = tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr+1)
          isegend = tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp) + &
                    tmb%linmat%l%nsegline(tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp)-1
          matrix_local = f_malloc_ptr(tmb%linmat%l%nvctrp,id='matrix_local')
          do ispin=1,tmb%linmat%l%nspin
              ishift=(ispin-1)*tmb%linmat%l%nvctr
              !$omp parallel default(none) &
              !$omp shared(isegstart,isegend,tmb,matrix_local,kernel_compr_tmp,ishift) &
              !$omp private(iseg,ii,i,irowcol)
              !$omp do
              do iseg=isegstart,isegend
                  ii=tmb%linmat%l%keyv(iseg)
                  do i=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
                      irowcol = orb_from_index(tmb%linmat%l, i)
                      if(irowcol(1)==irowcol(2)) then
                          matrix_local(ii-tmb%linmat%l%isvctr)=0.d0
                      else
                          matrix_local(ii-tmb%linmat%l%isvctr)=kernel_compr_tmp(ii+ishift)
                      end if
                      ii=ii+1
                  end do
              end do
              !$omp end do
              !$omp end parallel
              if (nproc>1) then
                   call timing(iproc,'buildgrad_mcpy','OF')
                   call timing(iproc,'buildgrad_comm','ON')
                   !!call mpi_allgatherv(matrix_local(1), tmb%linmat%l%nvctrp, mpi_double_precision, &
                   !!     tmb%linmat%kernel_%matrix_compr(ishift+1), tmb%linmat%l%nvctr_par, &
                   !!     tmb%linmat%l%isvctr_par, mpi_double_precision, &
                   !!     bigdft_mpi%mpi_comm, ierr)
                   if (comm_strategy==ALLGATHERV) then
                       call mpi_allgatherv(matrix_local(1), tmb%linmat%l%nvctrp, mpi_double_precision, &
                            tmb%linmat%kernel_%matrix_compr(ishift+1), tmb%linmat%l%nvctr_par, &
                            tmb%linmat%l%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                       call f_free_ptr(matrix_local)
                   else if (comm_strategy==GET) then
                       !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, tmb%linmat%l%nvctrp, matrix_local, &
                       !!     tmb%linmat%l%nvctr_par, tmb%linmat%l%isvctr_par, &
                       !!     tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(ishift+1:ishift+tmb%linmat%l%nvctr))
                       call mpi_get_to_allgatherv(matrix_local(1), tmb%linmat%l%nvctrp, &
                            tmb%linmat%kernel_%matrix_compr(ishift+1), &
                            tmb%linmat%l%nvctr_par, tmb%linmat%l%isvctr_par, bigdft_mpi%mpi_comm)
                   else
                       stop 'build_gradient: wrong communication strategy'
                   end if
                   call timing(iproc,'buildgrad_comm','OF')
                   call timing(iproc,'buildgrad_mcpy','ON')
                   if (ispin==tmb%linmat%l%nspin) call f_free_ptr(matrix_local)
               else
                   call vcopy(tmb%linmat%l%nvctr, matrix_local(1), 1, &
                        tmb%linmat%kernel_%matrix_compr(ishift+1), 1)
               end if
          end do


          ist=1
          do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
              ilr=tmb%orbs%inwhichlocreg(iorb)
              if (tmb%orbs%spinsgn(iorb)>0.d0) then
                  ispin=1
              else
                  ispin=2
              end if
              iiorb = mod(iorb-1,tmb%linmat%l%nfvctr)+1 ! spin-independent index
              ishift=(ispin-1)*tmb%linmat%l%nvctr
              isegstart = tmb%linmat%l%istsegline(iiorb)
              isegend = tmb%linmat%l%istsegline(iiorb) + tmb%linmat%l%nsegline(iiorb) - 1
              do iseg=isegstart,isegend
              !do iseg=1, tmb%linmat%l%nseg
                  ii=tmb%linmat%l%keyv(iseg)
                  do i=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
                      irowcol = orb_from_index(tmb%linmat%l, i)
                      if (irowcol(2)/=iiorb) stop'irowcol(2)/=iiorb'
                      !if(irowcol(1)==irowcol(2) .and. irowcol(1)==iiorb) then
                      if(irowcol(1)==irowcol(2)) then
                          ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
                          !write(*,*) 'iorb, ii, ishift, ist', iorb, ii, ishift, ist
                          call dscal(ncount, kernel_compr_tmp(ii+ishift), tmb%hpsi(ist), 1)
                          ist=ist+ncount
                      end if
                      ii=ii+1
                  end do
              end do
          end do
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%l, tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
          ! copy correct kernel back
          call vcopy(tmb%linmat%l%nvctr*tmb%linmat%l%nspin, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
          call f_free_ptr(kernel_compr_tmp)
      else
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%l, tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      end if

      call timing(iproc,'buildgrad_mcpy','ON')
      call f_release_routine()

end subroutine build_gradient
