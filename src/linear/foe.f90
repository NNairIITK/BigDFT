!> @file
!! Fermi Operator Expansion Medthod
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Could still do more tidying - assuming all sparse matrices except for Fermi have the same pattern
subroutine foe(iproc, nproc, orbs, foe_obj, tmprtr, &
           ebs, itout, it_scc, order_taylor, &
           tmb)
  use module_base
  use module_types
  use module_interfaces, except_this_one => foe
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc,itout,it_scc, order_taylor
  type(orbitals_data),intent(in) :: orbs
  type(foe_data),intent(inout) :: foe_obj
  real(kind=8),intent(inout) :: tmprtr
  real(kind=8),intent(out) :: ebs
  type(DFT_wavefunction),intent(inout) :: tmb

  ! Local variables
  integer :: npl, istat, iall, jorb, info, ipl, i, it, ierr, ii, iiorb, jjorb, iseg, it_solver, iorb
  integer :: isegstart, isegend, iismall, iseglarge, isegsmall, is, ie, iilarge, nsize_polynomial
  integer :: iismall_ovrlp, iismall_ham
  integer,parameter :: nplx=50000
  real(kind=8),dimension(:,:),allocatable :: cc, fermip, chebyshev_polynomials
  real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value, sumn, sumnder, charge_diff, ef_interpol, ddot
  real(kind=8) :: evlow_old, evhigh_old, m, b, det, determinant, sumn_old, ef_old, bound_low, bound_up, tt
  real(kind=8) :: fscale, error, tt_ovrlp, tt_ham
  logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible, emergency_stop
  character(len=*),parameter :: subname='foe'
  real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr, SHS
  real(kind=8),dimension(4,4) :: interpol_matrix, tmp_matrix
  real(kind=8),dimension(4) :: interpol_vector, interpol_solution
  integer,dimension(4) :: ipiv
  real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion
  integer :: jproc, iorder, npl_boundaries
  logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
  real(kind=8),dimension(:,:),allocatable :: workmat
  real(kind=8) :: trace_sparse

  type(sparseMatrix),pointer :: ovrlp, ham !to be able to deallocate
  !!type(sparseMatrix),pointer :: fermi, inv_ovrlp
  integer :: irow, icol, itemp, matrixindex_in_compressed, iflag
  logical :: overlap_calculated, cycle_FOE


  allocate(tmb%linmat%inv_ovrlp_large%matrix_compr(tmb%linmat%inv_ovrlp_large%nvctr), stat=istat)
  call memocc(istat, tmb%linmat%inv_ovrlp_large%matrix_compr, 'tmb%linmat%inv_ovrlp_large%matrix_compr', subname)


  call timing(iproc, 'FOE_auxiliary ', 'ON')

  ! initialization
  interpol_solution = 0.d0


  allocate(penalty_ev(orbs%norb,orbs%norbp,2), stat=istat)
  call memocc(istat, penalty_ev, 'penalty_ev', subname)

  allocate(fermip(orbs%norb,orbs%norbp), stat=istat)
  call memocc(istat, fermip, 'fermip', subname)

  allocate(SHS(tmb%linmat%denskern_large%nvctr), stat=istat)
  call memocc(istat, SHS, 'SHS', subname)

  if (order_taylor==1) then
      ii=0
      do iseg=1,tmb%linmat%denskern_large%nseg
          do jorb=tmb%linmat%denskern_large%keyg(1,iseg),tmb%linmat%denskern_large%keyg(2,iseg)
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              ii=ii+1
              iismall = matrixindex_in_compressed(tmb%linmat%ovrlp, iiorb, jjorb)
              if (iismall>0) then
                  tt=tmb%linmat%ovrlp%matrix_compr(iismall)
              else
                  tt=0.d0
              end if
              if (iiorb==jjorb) then
                  tmb%linmat%inv_ovrlp_large%matrix_compr(ii)=1.5d0-.5d0*tt
              else
                  tmb%linmat%inv_ovrlp_large%matrix_compr(ii)=-.5d0*tt
              end if
          end do  
      end do
  else

        call timing(iproc, 'FOE_auxiliary ', 'OF')
        call overlap_minus_onehalf() ! has internal timer
        call timing(iproc, 'FOE_auxiliary ', 'ON')

  end if


  allocate(hamscal_compr(tmb%linmat%denskern_large%nvctr), stat=istat)
  call memocc(istat, hamscal_compr, 'hamscal_compr', subname)

    
  ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
  nsize_polynomial=0
  if (orbs%norbp>0) then
      isegstart=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=tmb%linmat%denskern_large%nseg
      end if
      !$omp parallel default(private) shared(isegstart, isegend, tmb, nsize_polynomial)
      !$omp do reduction(+:nsize_polynomial)
      do iseg=isegstart,isegend
          do jorb=tmb%linmat%denskern_large%keyg(1,iseg),tmb%linmat%denskern_large%keyg(2,iseg)
              nsize_polynomial=nsize_polynomial+1
          end do
      end do
      !$omp end do
      !$omp end parallel
  end if
  
  
  ! Fake allocation, will be modified later
  allocate(chebyshev_polynomials(nsize_polynomial,1),stat=istat)
  call memocc(istat,chebyshev_polynomials,'chebyshev_polynomials',subname)

  fscale=foe_obj%fscale/0.5d0 ! this will be undone in the first iteration of the following loop


  temp_loop: do itemp=1,2

      
      fscale=fscale*0.5d0 ! make the error function sharper, i.e. more "step function-like"

      evlow_old=1.d100
      evhigh_old=-1.d100
      
      if (iproc==0) call yaml_map('decay length of error function',fscale,fmt='(es10.3)')
    
          ! Don't let this value become too small.
          foe_obj%bisection_shift = max(foe_obj%bisection_shift,1.d-4)
    
          efarr(1)=foe_obj%ef-foe_obj%bisection_shift
          efarr(2)=foe_obj%ef+foe_obj%bisection_shift
          sumnarr(1)=0.d0
          sumnarr(2)=1.d100
    
          adjust_lower_bound=.true.
          adjust_upper_bound=.true.
    
          calculate_SHS=.true.
    
          if (orbs%norbp>0) then
              call to_zero(orbs%norb*orbs%norbp, fermip(1,1))
          end if
    
          if (iproc==0) then
              !call yaml_sequence(advance='no')
              call yaml_open_sequence('FOE to determine density kernel',label=&
                   'it_foe'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//&
                   trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)')))//'-'//&
                   trim(adjustl(yaml_toa(itemp,fmt='(i2.2)'))))
          end if
    
    
    
          it=0
          it_solver=0
          eval_bounds_ok=.false.
          bisection_bounds_ok=.false.
          main_loop: do 
              
              it=it+1
    
              if (iproc==0) then
                  call yaml_newline()
                  call yaml_sequence(advance='no')
                  call yaml_open_map(flow=.true.)
                  call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
              end if
              
              if (adjust_lower_bound) then
                  foe_obj%ef=efarr(1)
              else if (adjust_upper_bound) then
                  foe_obj%ef=efarr(2)
              end if
          
    
              ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
              if (foe_obj%evlow/=evlow_old .or. foe_obj%evhigh/=evhigh_old) then
                  scale_factor=2.d0/(foe_obj%evhigh-foe_obj%evlow)
                  shift_value=.5d0*(foe_obj%evhigh+foe_obj%evlow)
                  !$omp parallel default(none) private(ii,irow,icol,iismall_ovrlp,iismall_ham,tt_ovrlp,tt_ham) &
                  !$omp shared(tmb,hamscal_compr,scale_factor,shift_value)
                  !$omp do
                  do ii=1,tmb%linmat%denskern_large%nvctr
                      irow = tmb%linmat%denskern_large%orb_from_index(1,ii)
                      icol = tmb%linmat%denskern_large%orb_from_index(2,ii)
                      iismall_ovrlp = matrixindex_in_compressed(tmb%linmat%ovrlp, irow, icol)
                      iismall_ham = matrixindex_in_compressed(tmb%linmat%ham, irow, icol)
                      if (iismall_ovrlp>0) then
                          tt_ovrlp=tmb%linmat%ovrlp%matrix_compr(iismall_ovrlp)
                      else
                          tt_ovrlp=0.d0
                      end if
                      if (iismall_ham>0) then
                          tt_ham=tmb%linmat%ham%matrix_compr(iismall_ham)
                      else
                          tt_ham=0.d0
                      end if
                      hamscal_compr(ii)=scale_factor*(tt_ham-shift_value*tt_ovrlp)
                      !hamscal_compr(ii)=scale_factor*(ham%matrix_compr(ii)-shift_value*tt)
                      !hamscal_compr(ii)=scale_factor*(ham%matrix_compr(ii)-shift_value*tmb%linmat%ovrlp%matrix_compr(iismall))
                  end do
                  !$omp end do
                  !$omp end parallel
                  calculate_SHS=.true.
              else
                  calculate_SHS=.false.
              end if
              evlow_old=foe_obj%evlow
              evhigh_old=foe_obj%evhigh
    
              !!foe_obj%ef = foe_obj%evlow+1.d-4*it
    
              ! Determine the degree of the polynomial
              !npl=nint(3.0d0*(foe_obj%evhigh-foe_obj%evlow)/foe_obj%fscale)
              npl=nint(3.0d0*(foe_obj%evhigh-foe_obj%evlow)/fscale)
              npl_boundaries=nint(3.0d0*(foe_obj%evhigh-foe_obj%evlow)/1.d-3) ! max polynomial degree for given eigenvalue boundaries
              if (npl>npl_boundaries) then
                  npl=npl_boundaries
                  if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
              end if
              if (npl>nplx) stop 'npl>nplx'
    
              ! Array the holds the Chebyshev polynomials. Needs to be recalculated
              ! every time the Hamiltonian has been modified.
              if (calculate_SHS) then
                  iall=-product(shape(chebyshev_polynomials))*kind(chebyshev_polynomials)
                  deallocate(chebyshev_polynomials,stat=istat)
                  call memocc(istat,iall,'chebyshev_polynomials',subname)
                  allocate(chebyshev_polynomials(nsize_polynomial,npl),stat=istat)
                  call memocc(istat,chebyshev_polynomials,'chebyshev_polynomials',subname)
              end if
    
              if (iproc==0) then
                  call yaml_map('bisec/eval bounds',&
                       (/efarr(1),efarr(2),foe_obj%evlow,foe_obj%evhigh/),fmt='(f5.2)')
                  !call yaml_map('lower bisection bound',efarr(1),fmt='(es10.3)')
                  !call yaml_map('upper bisection bound',efarr(1),fmt='(es10.3)')
                  !call yaml_map('lower eigenvalue bound',foe_obj%evlow,fmt='(es10.3)')
                  !call yaml_map('upper eigenvalue bound',foe_obj%evhigh,fmt='(es10.3)')
                  call yaml_map('pol deg',npl,fmt='(i3)')
                  call yaml_map('eF',foe_obj%ef,fmt='(es16.9)')
                  !call yaml_map('conv crit',charge_tolerance,fmt='(es10.3)')
              end if
    
    
              allocate(cc(npl,3), stat=istat)
              call memocc(istat, cc, 'cc', subname)
    
              if (foe_obj%evlow>=0.d0) then
                  stop 'ERROR: lowest eigenvalue must be negative'
              end if
              if (foe_obj%evhigh<=0.d0) then
                  stop 'ERROR: highest eigenvalue must be positive'
              end if
    
              call timing(iproc, 'FOE_auxiliary ', 'OF')
              call timing(iproc, 'chebyshev_coef', 'ON')
    
              !call chebft(foe_obj%evlow, foe_obj%evhigh, npl, cc(1,1), foe_obj%ef, foe_obj%fscale, temperature)
              call chebft(foe_obj%evlow, foe_obj%evhigh, npl, cc(1,1), foe_obj%ef, fscale, tmprtr)
              call chder(foe_obj%evlow, foe_obj%evhigh, cc(1,1), cc(1,2), npl)
              call chebft2(foe_obj%evlow, foe_obj%evhigh, npl, cc(1,3))
              call evnoise(npl, cc(1,3), foe_obj%evlow, foe_obj%evhigh, anoise)
    
              call timing(iproc, 'chebyshev_coef', 'OF')
              call timing(iproc, 'FOE_auxiliary ', 'ON')
            
              !!if (iproc==0) then
              !!    call pltwght(npl,cc(1,1),cc(1,2),foe_obj%evlow,foe_obj%evhigh,foe_obj%ef,foe_obj%fscale,temperature)
              !!    call pltexp(anoise,npl,cc(1,3),foe_obj%evlow,foe_obj%evhigh)
              !!end if
            
            
              if (orbs%nspin==1) then
                  do ipl=1,npl
                      cc(ipl,1)=2.d0*cc(ipl,1)
                      cc(ipl,2)=2.d0*cc(ipl,2)
                      cc(ipl,3)=2.d0*cc(ipl,3)
                  end do
              end if
            
            
              call timing(iproc, 'FOE_auxiliary ', 'OF')
    
              if (calculate_SHS) then
                  ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                  if (iproc==0) call yaml_map('polynomials','recalculated')
                  call chebyshev_clean(iproc, nproc, npl, cc, orbs, foe_obj, &
                       tmb%linmat%denskern_large, hamscal_compr, &
                       tmb%linmat%inv_ovrlp_large%matrix_compr, calculate_SHS, &
                       nsize_polynomial, SHS, fermip, penalty_ev, chebyshev_polynomials, &
                       emergency_stop)
              else
                  ! The Chebyshev polynomials are already available
                  if (iproc==0) call yaml_map('polynomials','from memory')
                  call chebyshev_fast(iproc, nsize_polynomial, npl, orbs, &
                      tmb%linmat%denskern_large, chebyshev_polynomials, cc, fermip)
              end if 

             ! Check for an emergency stop, which happens if the kernel explodes, presumably due
             ! to the eigenvalue bounds being too small.
             ! mpi_lor seems not to work on certain systems...
             !call mpiallred(emergency_stop, 1, mpi_lor, bigdft_mpi%mpi_comm, ierr)
             if (emergency_stop) then
                 iflag=1
             else
                 iflag=0
             end if
             call mpiallred(iflag, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
             if (iflag>0) then
                 emergency_stop=.true.
             else
                 emergency_stop=.false.
             end if
             if (emergency_stop) then
                  eval_bounds_ok(1)=.false.
                  foe_obj%evlow=foe_obj%evlow*1.2d0
                  eval_bounds_ok(2)=.false.
                  foe_obj%evhigh=foe_obj%evhigh*1.2d0
                  if (iproc==0) then
                      call yaml_map('eval/bisection bounds ok',&
                           (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      call yaml_close_map()
                      !call bigdft_utils_flush(unit=6)
                  end if
                  iall=-product(shape(cc))*kind(cc)
                  deallocate(cc, stat=istat)
                  call memocc(istat, iall, 'cc', subname)
                  cycle main_loop
             end if
    
    
              call timing(iproc, 'FOE_auxiliary ', 'ON')
    
    
              restart=.false.
    
              ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
              ! (otherwise this has already been checked in the previous iteration).
              if (calculate_SHS) then
                  ! The penalty function must be smaller than the noise.
                  bound_low=0.d0
                  bound_up=0.d0
                  if (orbs%norbp>0) then
                      isegstart=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc)+1)
                      if (orbs%isorb+orbs%norbp<orbs%norb) then
                          isegend=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc+1)+1)-1
                      else
                          isegend=tmb%linmat%denskern_large%nseg
                      end if
                      !$omp parallel default(private) &
                      !$omp shared(isegstart, isegend, orbs, penalty_ev, tmb, bound_low, bound_up)
                      !$omp do reduction(+:bound_low,bound_up)
                      do iseg=isegstart,isegend
                          ii=tmb%linmat%denskern_large%keyv(iseg)-1
                          do jorb=tmb%linmat%denskern_large%keyg(1,iseg),tmb%linmat%denskern_large%keyg(2,iseg)
                              ii=ii+1
                              iiorb = (jorb-1)/orbs%norb + 1
                              jjorb = jorb - (iiorb-1)*orbs%norb
                              iismall = matrixindex_in_compressed(tmb%linmat%ovrlp, iiorb, jjorb)
                              if (iismall>0) then
                                  tt=tmb%linmat%ovrlp%matrix_compr(iismall)
                              else
                                  tt=0.d0
                              end if
                              bound_low = bound_low + penalty_ev(jjorb,iiorb-orbs%isorb,2)*tt
                              bound_up = bound_up +penalty_ev(jjorb,iiorb-orbs%isorb,1)*tt
                          end do  
                      end do
                      !$omp end do
                      !$omp end parallel
                  end if
    
                  allredarr(1)=bound_low
                  allredarr(2)=bound_up
                  call mpiallred(allredarr, 2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                  allredarr=abs(allredarr) !for some crazy situations this may be negative
                  anoise=10.d0*anoise
                  if (allredarr(1)>anoise) then
                      eval_bounds_ok(1)=.false.
                      foe_obj%evlow=foe_obj%evlow*1.2d0
                      restart=.true.
                  else
                      eval_bounds_ok(1)=.true.
                  end if
                  if (allredarr(2)>anoise) then
                      eval_bounds_ok(2)=.false.
                      foe_obj%evhigh=foe_obj%evhigh*1.2d0
                      restart=.true.
                  else
                      eval_bounds_ok(2)=.true.
                  end if
              end if
    
              iall=-product(shape(cc))*kind(cc)
              deallocate(cc, stat=istat)
              call memocc(istat, iall, 'cc', subname)
    
              if (restart) then
                  if (iproc==0) then
                      call yaml_map('eval/bisection bounds ok',&
                           (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      call yaml_close_map()
                      !call bigdft_utils_flush(unit=6)
                  end if
                  cycle
              end if
                  
            
    
              sumn=0.d0
              sumnder=0.d0
              if (orbs%norbp>0) then
                  !do jproc=0,nproc-1
                  !    if (iproc==jproc) then
                          isegstart=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc)+1)
                          if (orbs%isorb+orbs%norbp<orbs%norb) then
                              isegend=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc+1)+1)-1
                          else
                              isegend=tmb%linmat%denskern_large%nseg
                          end if
                          !$omp parallel default(private) shared(isegstart, isegend, orbs, fermip, tmb, sumn) 
                          !$omp do reduction(+:sumn)
                          do iseg=isegstart,isegend
                              ii=tmb%linmat%denskern_large%keyv(iseg)-1
                              do jorb=tmb%linmat%denskern_large%keyg(1,iseg),tmb%linmat%denskern_large%keyg(2,iseg)
                                  ii=ii+1
                                  iiorb = (jorb-1)/orbs%norb + 1
                                  jjorb = jorb - (iiorb-1)*orbs%norb
                                  if (jjorb==iiorb) sumn = sumn + fermip(jjorb,iiorb-orbs%isorb)
                              end do  
                          end do
                          !$omp end do
                          !$omp end parallel
                  !    end if
                  !    call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
                  !end do
              end if
    
              if (nproc>1) then
                  call mpiallred(sumn, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
              end if
    
    
              ! Make sure that the bounds for the bisection are negative and positive
              restart=.false.
              charge_diff = sumn-foe_obj%charge
              if (adjust_lower_bound) then
                  !if (iproc==0) call yaml_map('checking lower bisection bound, charge diff',charge_diff,fmt='(es9.2)')
                  if (charge_diff<=0.d0) then
                      ! Lower bound okay
                      adjust_lower_bound=.false.
                      foe_obj%bisection_shift=foe_obj%bisection_shift*9.d-1
                      sumnarr(1)=sumn
                      if (iproc==0) then
                      end if
                      restart=.true.
                      bisection_bounds_ok(1)=.true.
                  else
                      efarr(1)=efarr(1)-foe_obj%bisection_shift
                      foe_obj%bisection_shift=foe_obj%bisection_shift*1.1d0
                      restart=.true.
                      bisection_bounds_ok(1)=.false.
                  end if
              end if
              if (restart) then
                  if (iproc==0) then
                      call yaml_map('eval/bisection bounds ok',&
                           (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      call yaml_close_map()
                  end if
                  cycle
              end if
              if (adjust_upper_bound) then
                  if (charge_diff>=0.d0) then
                      ! Upper bound okay
                      adjust_upper_bound=.false.
                      foe_obj%bisection_shift=foe_obj%bisection_shift*9.d-1
                      sumnarr(2)=sumn
                      restart=.false.
                      bisection_bounds_ok(2)=.true.
                  else
                      efarr(2)=efarr(2)+foe_obj%bisection_shift
                      foe_obj%bisection_shift=foe_obj%bisection_shift*1.1d0
                      restart=.true.
                      bisection_bounds_ok(2)=.false.
                  end if
              end if
    
              if (iproc==0) then
                  call yaml_map('eval/bisection bounds ok',&
                       (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
              end if
              if (restart) then
                  if (iproc==0) then
                      call yaml_close_map()
                  end if
                  cycle
              end if
    
    
              it_solver=it_solver+1
    
              ! Check whether the system behaves reasonably.
              interpolation_possible=.true.
              if (it_solver>1) then
                  if (iproc==0) call yaml_newline()
                  if (iproc==0) call yaml_open_map('interpol check',flow=.true.)
                  if (iproc==0) call yaml_map('D eF',foe_obj%ef-ef_old,fmt='(es13.6)')
                  if (iproc==0) call yaml_map('D Tr',sumn-sumn_old,fmt='(es13.6)')
                  if (foe_obj%ef>ef_old .and. sumn<sumn_old) then
                      interpolation_possible=.false.
                  end if
                  if (foe_obj%ef<ef_old .and. sumn>sumn_old) then
                      interpolation_possible=.false.
                  end if
                  if (foe_obj%ef>ef_old .and. sumn<sumn_old .or. foe_obj%ef<ef_old .and. sumn>sumn_old) then
                      if (iproc==0) call yaml_map('interpol possible',.false.)
                  else
                      if (iproc==0) call yaml_map('interpol possible',.true.)
                  end if
                  if (iproc==0) call yaml_close_map()
                  !!call bigdft_utils_flush(unit=6)
                  if (iproc==0) call yaml_newline()
              end if
              if (.not.interpolation_possible) then
                  ! Set the history for the interpolation to zero.
                  it_solver=0
              end if
    
              ef_old=foe_obj%ef
              sumn_old=sumn
    
              ! Shift up the old results.
              if (it_solver>4) then
                  do i=1,4
                      interpol_matrix(1,i)=interpol_matrix(2,i)
                      interpol_matrix(2,i)=interpol_matrix(3,i)
                      interpol_matrix(3,i)=interpol_matrix(4,i)
                  end do
                  interpol_vector(1)=interpol_vector(2)
                  interpol_vector(2)=interpol_vector(3)
                  interpol_vector(3)=interpol_vector(4)
              end if
              !LG: if it_solver==0 this index comes out of bounds!
              ii=max(min(it_solver,4),1)
              interpol_matrix(ii,1)=foe_obj%ef**3
              interpol_matrix(ii,2)=foe_obj%ef**2
              interpol_matrix(ii,3)=foe_obj%ef
              interpol_matrix(ii,4)=1
              interpol_vector(ii)=sumn-foe_obj%charge
    
              ! Solve the linear system interpol_matrix*interpol_solution=interpol_vector
              if (it_solver>=4) then
                  do i=1,ii
                      interpol_solution(i)=interpol_vector(i)
                      tmp_matrix(i,1)=interpol_matrix(i,1)
                      tmp_matrix(i,2)=interpol_matrix(i,2)
                      tmp_matrix(i,3)=interpol_matrix(i,3)
                      tmp_matrix(i,4)=interpol_matrix(i,4)
                  end do
    
                  call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
                  if (info/=0) then
                     if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
                  end if
    
    
                  call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
                       interpol_solution(3), interpol_solution(4), foe_obj%ef, ef_interpol)
              end if
    
    
              ! Adjust the bounds for the bisection.
              if (charge_diff<0.d0) then
                  efarr(1)=foe_obj%ef
                  sumnarr(1)=sumn
              else if (charge_diff>=0.d0) then
                  efarr(2)=foe_obj%ef
                  sumnarr(2)=sumn
              end if
    
    
              ! Calculate the new Fermi energy.
              if (iproc==0) call yaml_newline()
              if (iproc==0) call yaml_open_map('Search new eF',flow=.true.)
              if (it_solver>=4 .and.  abs(sumn-foe_obj%charge)<foe_obj%ef_interpol_chargediff) then
                  det=determinant(iproc,4,interpol_matrix)
                  if (iproc==0) call yaml_map('det',det,fmt='(es10.3)')
                  if (iproc==0) call yaml_map('limit',foe_obj%ef_interpol_det,fmt='(es10.3)')
                  if(abs(det)>foe_obj%ef_interpol_det) then
                      foe_obj%ef=ef_interpol
                      if (iproc==0) call yaml_map('method','cubic interpolation')
                  else
                      ! linear interpolation
                      if (iproc==0) call yaml_map('method','linear interpolation')
                      m = (interpol_vector(4)-interpol_vector(3))/(interpol_matrix(4,3)-interpol_matrix(3,3))
                      b = interpol_vector(4)-m*interpol_matrix(4,3)
                      foe_obj%ef = -b/m
                  end if
              else
                  ! Use mean value of bisection and secant method
                  ! Secant method solution
                  foe_obj%ef = efarr(2)-(sumnarr(2)-foe_obj%charge)*(efarr(2)-efarr(1))/(sumnarr(2)-sumnarr(1))
                  ! Add bisection solution
                  foe_obj%ef = foe_obj%ef + .5d0*(efarr(1)+efarr(2))
                  ! Take the mean value
                  foe_obj%ef=.5d0*foe_obj%ef
                  if (iproc==0) call yaml_map('method','bisection / secant method')
              end if
              if (iproc==0) then
                  call yaml_close_map()
                  !!call bigdft_utils_flush(unit=6)
                  !call yaml_newline()
              end if
    
    
              if (iproc==0) then
                  call yaml_newline()
                  call yaml_map('iter',it)
                  call yaml_map('Tr(K)',sumn,fmt='(es16.9)')
                  call yaml_map('charge diff',sumn-foe_obj%charge,fmt='(es16.9)')
              end if
    
              if (iproc==0) then
                  call yaml_close_map()
                  !call bigdft_utils_flush(unit=6)
              end if
    
              if (abs(charge_diff)<charge_tolerance) then
                  if (iproc==0) call yaml_close_sequence()
                  exit
              end if
    
    
          end do main_loop
    
    
    
    
      call to_zero(tmb%linmat%denskern_large%nvctr, tmb%linmat%denskern_large%matrix_compr(1))
    
      if (orbs%norbp>0) then
          isegstart=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc)+1)
          if (orbs%isorb+orbs%norbp<orbs%norb) then
              isegend=tmb%linmat%denskern_large%istsegline(orbs%isorb_par(iproc+1)+1)-1
          else
              isegend=tmb%linmat%denskern_large%nseg
          end if
          !$omp parallel default(private) shared(isegstart, isegend, orbs, fermip, tmb)
          !$omp do
          do iseg=isegstart,isegend
              ii=tmb%linmat%denskern_large%keyv(iseg)-1
              do jorb=tmb%linmat%denskern_large%keyg(1,iseg),tmb%linmat%denskern_large%keyg(2,iseg)
                  ii=ii+1
                  iiorb = (jorb-1)/orbs%norb + 1
                  jjorb = jorb - (iiorb-1)*orbs%norb
                  tmb%linmat%denskern_large%matrix_compr(ii)=fermip(jjorb,iiorb-orbs%isorb)
              end do
          end do
          !$omp end do
          !$omp end parallel
      end if

      call timing(iproc, 'FOE_auxiliary ', 'OF')
      call timing(iproc, 'chebyshev_comm', 'ON')
    
      call mpiallred(tmb%linmat%denskern_large%matrix_compr(1), tmb%linmat%denskern_large%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    
    
      call timing(iproc, 'chebyshev_comm', 'OF')
    
    
      call overlap_minus_onehalf() !has internal timer
      call timing(iproc, 'FOE_auxiliary ', 'ON')
    
    
    
      allocate(tmb%linmat%inv_ovrlp_large%matrix(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, tmb%linmat%inv_ovrlp_large%matrix, 'tmb%linmat%inv_ovrlp_large%matrix', subname)
    
      allocate(workmat(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, workmat, 'workmat', subname)
    
      call uncompressMatrix(iproc,tmb%linmat%inv_ovrlp_large)
    
      allocate(tmb%linmat%denskern_large%matrix(orbs%norb,orbs%norb))
      call memocc(istat, tmb%linmat%denskern_large%matrix, 'tmb%linmat%denskern_large%matrix', subname)
      call uncompressMatrix(iproc,tmb%linmat%denskern_large)
    
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tmb%linmat%inv_ovrlp_large%matrix, orbs%norb, &
                 tmb%linmat%denskern_large%matrix, orbs%norb, 0.d0, workmat(1,1), orbs%norb)
      call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, workmat(1,1), orbs%norb, &
                 tmb%linmat%inv_ovrlp_large%matrix, orbs%norb, 0.d0, tmb%linmat%denskern_large%matrix(1,1), orbs%norb)
    
    
    
      call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern_large)

     iall=-product(shape(tmb%linmat%denskern_large%matrix))*kind(tmb%linmat%denskern_large%matrix)
     deallocate(tmb%linmat%denskern_large%matrix,stat=istat)
     call memocc(istat,iall,'tmb%linmat%denskern_large%matrix',subname)

     iall=-product(shape(workmat))*kind(workmat)
     deallocate(workmat,stat=istat)
     call memocc(istat,iall,'workmat',subname)

     iall=-product(shape(tmb%linmat%inv_ovrlp_large%matrix))*kind(tmb%linmat%inv_ovrlp_large%matrix)
     deallocate(tmb%linmat%inv_ovrlp_large%matrix,stat=istat)
     call memocc(istat,iall,'tmb%linmat%inv_ovrlp_large%matrix',subname)
    
  
      ! Purify the kernel
      if (iproc==0) then
          call yaml_open_sequence('Final kernel purification')
          call yaml_newline()
      end if
      overlap_calculated=.false.
      tmb%can_use_transposed=.false.
      call purify_kernel(iproc, nproc, tmb, overlap_calculated)
      if (iproc==0) then
          call yaml_close_sequence()
      end if
    
    
      ! Calculate trace(KS).
      sumn=trace_sparse(iproc, nproc, orbs, tmb%linmat%ovrlp, tmb%linmat%denskern_large)



      ! Check whether this agrees with the number of electrons. If not,
      ! calculate a new kernel with a sharper decay of the error function
      ! (correponds to a lower temperature)
      if (abs(sumn-foe_obj%charge)>1.d-6) then
          cycle_FOE=.true.
      else
          cycle_FOE=.false.
      end if
      if (iproc==0) then
          call yaml_map('trace(KS)',sumn)
          call yaml_map('need to repeat with sharper decay',cycle_FOE)
      end if
      if (.not.cycle_FOE) exit temp_loop

    

  end do temp_loop


  iall=-product(shape(tmb%linmat%inv_ovrlp_large%matrix_compr))*kind(tmb%linmat%inv_ovrlp_large%matrix_compr)
  deallocate(tmb%linmat%inv_ovrlp_large%matrix_compr,stat=istat)
  call memocc(istat,iall,'tmb%linmat%inv_ovrlp_large%matrix_compr',subname)
  


  scale_factor=1.d0/scale_factor
  shift_value=-shift_value


  ! Calculate trace(KH). Since they have the same sparsity pattern and K is
  ! symmetric, this is a simple ddot.
  ebs=ddot(tmb%linmat%denskern_large%nvctr, tmb%linmat%denskern_large%matrix_compr,1 , hamscal_compr, 1)
  ebs=ebs*scale_factor-shift_value*sumn





  if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='-')


  iall=-product(shape(chebyshev_polynomials))*kind(chebyshev_polynomials)
  deallocate(chebyshev_polynomials,stat=istat)
  call memocc(istat,iall,'chebyshev_polynomials',subname)

  iall=-product(shape(penalty_ev))*kind(penalty_ev)
  deallocate(penalty_ev, stat=istat)
  call memocc(istat, iall, 'penalty_ev', subname)

  iall=-product(shape(hamscal_compr))*kind(hamscal_compr)
  deallocate(hamscal_compr, stat=istat)
  call memocc(istat, iall, 'hamscal_compr', subname)

  iall=-product(shape(fermip))*kind(fermip)
  deallocate(fermip, stat=istat)
  call memocc(istat, iall, 'fermip', subname)

  iall=-product(shape(SHS))*kind(SHS)
  deallocate(SHS, stat=istat)
  call memocc(istat, iall, 'SHS', subname)


  call timing(iproc, 'FOE_auxiliary ', 'OF')




      contains

        subroutine overlap_minus_onehalf()
          ! Taylor approximation of S^-1/2 up to higher order
          allocate(tmb%linmat%ovrlp%matrix(orbs%norb,orbs%norb), stat=istat)
          call memocc(istat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)
          call uncompressMatrix(iproc,tmb%linmat%ovrlp)

          allocate(tmb%linmat%inv_ovrlp_large%matrix(orbs%norb,orbs%norb),stat=istat)
          call memocc(istat, tmb%linmat%inv_ovrlp_large%matrix,'tmb%linmat%inv_ovrlp_large%matrix',subname)

          call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, orbs%norb, &
               tmb%linmat%ovrlp%matrix, tmb%linmat%inv_ovrlp_large%matrix, error, orbs, check_accur=.true.)
          if (iproc==0) then
              call yaml_map('error of S^-1/2',error,fmt='(es9.2)')
          end if
          call compress_matrix_for_allreduce(iproc,tmb%linmat%inv_ovrlp_large)

          iall=-product(shape(tmb%linmat%inv_ovrlp_large%matrix))*kind(tmb%linmat%inv_ovrlp_large%matrix)
          deallocate(tmb%linmat%inv_ovrlp_large%matrix,stat=istat)
          call memocc(istat,iall,'tmb%linmat%inv_ovrlp_large%matrix',subname)

          iall=-product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
          deallocate(tmb%linmat%ovrlp%matrix,stat=istat)
          call memocc(istat,iall,'tmb%linmat%ovrlp%matrix',subname)
      end subroutine overlap_minus_onehalf


end subroutine foe




! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft(A,B,N,cc,ef,fscale,tmprtr)
  use module_base, pi => pi_param
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
  integer,intent(in) :: n
  real(8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8) :: bma, bpa, y, arg, fac, tt, erfcc
  real(kind=8),dimension(50000) :: cf
  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  if (n>50000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      if (tmprtr.eq.0.d0) then
          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
      else
          cf(k)=1.d0/(1.d0+exp( (arg-ef)*(1.d0/tmprtr) ) )
      end if
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do

end subroutine chebft



! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft2(a,b,n,cc)
  use module_base, pi => pi_param
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(kind=8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
  real(kind=8) :: tt, y, arg, fac, bma, bpa
  real(kind=8),dimension(50000) :: cf

  if (n>50000) stop 'chebft2'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  ! 3 gives broder safety zone than 4
  !tt=3.0d0*n/(b-a)
  tt=4.d0*n/(b-a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      cf(k)=exp((arg-b)*tt)
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
end subroutine chebft2


! Calculates chebychev expansion of the derivative of Fermi distribution.
subroutine chder(a,b,c,cder,n)
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(8),dimension(n),intent(in) :: c
  real(8),dimension(n),intent(out) :: cder

  ! Local variables
  integer :: j
  real(kind=8) :: con

  cder(n)=0.d0
  cder(n-1)=2*(n-1)*c(n)
  if(n>=3)then
      do j=n-2,1,-1
        cder(j)=cder(j+2)+2*j*c(j+1)
      end do
  end if
  con=2.d0/(b-a)
  do j=1,n
      cder(j)=cder(j)*con
  end do

end subroutine chder


!> Determine noise level
subroutine evnoise(npl,cc,evlow,evhigh,anoise)
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: npl
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),intent(in) :: evlow, evhigh
  real(kind=8),intent(out) :: anoise
  
  ! Local variables
  real(kind=8) :: fact, dist, ddx, cent, tt, x, chebev
  
  fact=1.d0
  dist=(fact*evhigh-fact*evlow)
  ddx=dist/(10*npl)
  cent=.5d0*(fact*evhigh+fact*evlow)
  !!tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  !!do x=ddx,.25d0*dist,ddx
  !!    tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
  !!       & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
  !!end do
  ! Rewritten version ob the above loop
  tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  x=ddx
  do 
      tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
         & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
      x=x+ddx
      if (x>=.25d0*dist) exit
  end do
  !anoise=1.d0*tt
  anoise=20.d0*tt

end subroutine evnoise



!> Calculates the error function complement with an error of less than 1.2E-7
function erfcc(x)
  implicit none

  ! Calling arguments
  real(8),intent(in) :: x
  real(8) :: erfcc

  ! Local variables
  real(8) :: z, t

  z=abs(x)
  t=1.d0/(1.+0.5d0*z)
  erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
        & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
        & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  if (x.lt.0.) erfcc=2.D0-erfcc

end function erfcc


!> Evaluates chebychev expansion
function chebev(a,b,m,x,cc)
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: a, b, x
  integer,intent(in) :: m
  real(kind=8),dimension(m),intent(in) :: cc
  real(kind=8) :: chebev
  
  ! Local variables
  integer :: j
  real(kind=8) :: d, dd, y, sv
  
  d=0.d0
  dd=0.d0
  y=2.d0*(2.d0*x-a-b)/(b-a)
  do j=m,2,-1
      sv=d
      d=y*d-dd+cc(j)
      dd=sv
  end do
  chebev= -dd + 0.5d0*(y*d+cc(1))

end function chebev
 



! plots the approximate fermi distribution
        subroutine pltwght(npl,cc,cder,evlow,evhigh,ef,fscale,tmprtr)
          implicit none

          ! Calling arguments
          integer,intent(in) :: npl
          real(kind=8),dimension(npl),intent(in) :: cc, cder
          real(kind=8),intent(in) :: evlow, evhigh, ef, fscale, tmprtr

          ! Local variables
          integer :: ic
          real(kind=8) :: ddx, x, tt, err, chebev

        open (unit=66,file='fermi',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 3'
        write(66,*) ' #LINETYPE{132}'
65        format(a,f5.2,a,i3,a)
        write(66,65) ' #TITLE{WEIGHT DISTR. for fscale=', fscale,' npl=',npl,'}'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' #XYAPT{WEIGHT DISTR.}'
        write(66,*) ' #2YAXIS{2}'
        write(66,*) ' #YLOGSCALE2'
        write(66,*) ' #YCAPT2{ERROR}'
        write(66,*) ' $'
!
! plot chebechev expansion of weight distribution function
!
        ddx=(evhigh-evlow)/(10*npl)
! number of plot p[oints
        ic=0
        !!do x=evlow,evhigh,ddx
        !!    ic=ic+1
        !!end do
        x=evlow
        do
            ic=ic+1
            x=x+ddx
            if (x>=evhigh) exit
        end do
! weight distribution
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
        !!end do
        x=evlow
        do
            write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
            x=x+ddx
            if (x>=evhigh) exit
        end do
! derivative
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
        !!end do
        x=evlow
        do
            write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
            x=x+ddx
            if (x>=evhigh) exit
        end do
! error
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    tt=tmprtr
        !!    if (tmprtr.eq.0.d0) tt=1.d-16
        !!    err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
        !!    write(66,*) x,err
        !!end do
        x=evlow
        do
            tt=tmprtr
            if (tmprtr.eq.0.d0) tt=1.d-16
            err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
            write(66,*) x,err
            x=x+ddx
            if (x>=evhigh) exit
        end do

        close(unit=66)
end subroutine pltwght




! plots the approximate fermi distribution
subroutine pltexp(anoise,npl,cc,evlow,evhigh)
        implicit none

        ! Calling arguments
        integer,intent(in) :: npl
        real(kind=8),dimension(npl),intent(in) :: cc
        real(kind=8),intent(in) :: anoise, evlow, evhigh

        ! Local variables
        integer :: ic
        real(kind=8) :: fact, ddx, tt, chebev, x

        open (unit=66,file='exp',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 2'
        write(66,*) ' #LINETYPE{12}'
        write(66,*) ' #TITLE{exponential}'
        write(66,*) ' #YLOGSCALE'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' $'
!
        fact=1.25d0
! plot chebechev expansion of weight distribution function
!
        ddx=(fact*evhigh-fact*evlow)/(10*npl)
! number of plot p[oints
        ic=0
        !!do x=fact*evlow,fact*evhigh,ddx
        !!    ic=ic+1
        !!end do
        x=fact*evlow
        do
            ic=ic+1
            x=x+ddx
            if (x>=fact*evhigh) exit
        end do
! first curve
        write(66,*) ic
        !!do x=fact*evlow,fact*evhigh,ddx
        !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
        !!    if (abs(tt).lt.anoise) tt=anoise
        !!    write(66,*) x,tt
        !!end do
        x=fact*evlow
        do
            tt=CHEBEV(evlow,evhigh,npl,x,cc)
            if (abs(tt).lt.anoise) tt=anoise
            write(66,*) x,tt
            x=x+ddx
            if (x>=fact*evhigh) exit
        end do
! second curve
        write(66,*) ic
        !!do x=fact*evhigh,fact*evlow,-ddx
        !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
        !!    if (abs(tt).lt.anoise) tt=anoise
        !!    write(66,*) fact*evhigh-(x-fact*evlow),tt
        !!end do
        x=fact*evhigh
        do
            tt=CHEBEV(evlow,evhigh,npl,x,cc)
            if (abs(tt).lt.anoise) tt=anoise
            write(66,*) fact*evhigh-(x-fact*evlow),tt
            x=x-ddx
            if (x<=fact*evlow) exit
        end do

        close(unit=66)
end subroutine pltexp



! Finds the real root of the equation ax**3 + bx**2 + cx + d which is closest to target_solution
subroutine get_roots_of_cubic_polynomial(a, b, c, d, target_solution, solution)
  use module_base
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b, c, d
  real(kind=8),intent(in) :: target_solution
  real(kind=8),intent(out) :: solution

  ! Local variables
  complex(kind=8) :: a_c, b_c, c_c, d_c, Q_c, S_c, ttp_c, ttm_c
  complex(kind=8),dimension(3) :: sol_c
  real(kind=8) :: ttmin, tt
  integer :: i

  a_c=cmplx(a,0.d0,kind=8)
  b_c=cmplx(b,0.d0,kind=8)
  c_c=cmplx(c,0.d0,kind=8)
  d_c=cmplx(d,0.d0,kind=8)

  Q_c = sqrt( (2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c)**2 - 4*(b_c**2-3*a_c*c_c)**3 )
  S_c = ( .5d0*(Q_c+2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c) )**(1.d0/3.d0)
  ttp_c = cmplx(1.d0,sqrt(3.d0),kind=8)
  ttm_c = cmplx(1.d0,-sqrt(3.d0),kind=8)
  sol_c(1) = -b_c/(3*a_c) - S_c/(3*a_c) - (b_c**2-3*a_c*c_c)/(3*a_c*S_c)
  sol_c(2) = -b_c/(3*a_c) + (S_c*ttp_c)/(6*a_c) + ttm_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
  sol_c(3) = -b_c/(3*a_c) + (S_c*ttm_c)/(6*a_c) + ttp_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
  !!if (iproc==0) then
  !!    write(*,*) 'sol 1', sol_c(1)
  !!    write(*,*) 'sol 2', sol_c(2)
  !!    write(*,*) 'sol 3', sol_c(3)
  !!end if

  ! Select the real solution that is closest to target_solution
  ttmin=1.d100
  do i=1,3
      if (abs(aimag(sol_c(i)))>1.d-14) cycle !complex solution
      tt=abs(real(sol_c(i))-target_solution)
      if (tt<ttmin) then
          ttmin=tt
          solution=real(sol_c(i))
      end if
  end do

end subroutine get_roots_of_cubic_polynomial



real(kind=8) function determinant(iproc, n, mat)
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, n
    real(kind=8),dimension(n,n),intent(in) :: mat

    ! Local variables
    integer :: i, info
    integer,dimension(n) :: ipiv
    real(kind=8),dimension(n,n) :: mat_tmp
    real(kind=8) :: sgn

    call dcopy(n**2, mat, 1, mat_tmp, 1)

    call dgetrf(n, n, mat_tmp, n, ipiv, info)
    if (info/=0) then
        if (iproc==0) write(*,'(a,i0,a)') 'ERROR in dgetrf, info=',info,'. Set determinant to zero.'
        determinant=0
        return
    end if

    determinant=1.d0
    do i=1,n
        determinant=determinant*mat_tmp(i,i)
    end do

    sgn=1.d0
    do i=1,n
        if(ipiv(i)/=i) then
            sgn=-sgn
        end if
    end do

    determinant=sgn*determinant   

end function determinant


subroutine compress_polynomial_vector(iproc, nsize_polynomial, orbs, fermi, vector, vector_compressed)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nsize_polynomial
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: fermi
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(in) :: vector
  real(kind=8),dimension(nsize_polynomial),intent(out) :: vector_compressed

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (orbs%norbp>0) then
      isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=fermi%nseg
      end if
      ii=0
      !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermi, vector, vector_compressed)
      !!$omp do
      do iseg=isegstart,isegend
          !ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              vector_compressed(ii)=vector(jjorb,iiorb-orbs%isorb)
              !write(300,*) 'ii, jjorb, iiorb-orbs%isorb', ii, jjorb, iiorb-orbs%isorb
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if
end subroutine compress_polynomial_vector



subroutine uncompress_polynomial_vector(iproc, nsize_polynomial, orbs, fermi, vector_compressed, vector)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nsize_polynomial
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial),intent(in) :: vector_compressed
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(out) :: vector

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (orbs%norbp>0) then
      call to_zero(orbs%norb*orbs%norbp, vector(1,1))
      isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=fermi%nseg
      end if
      ii=0
      !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermi, vector, vector_compressed)
      !!$omp do
      do iseg=isegstart,isegend
          !ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              vector(jjorb,iiorb-orbs%isorb)=vector_compressed(ii)
              !write(*,*) 'ii, iiorb-orbs%isorb, jjorb', ii, iiorb-orbs%isorb, jjorb
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if
end subroutine uncompress_polynomial_vector


!< Calculates the trace of the matrix product amat*bmat.
!< WARNING: It is mandatory that the sparsity pattern of amat is contained
!< within the sparsity pattern of bmat!
function trace_sparse(iproc, nproc, orbs, amat, bmat)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,  nproc
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: amat, bmat

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iilarge
  integer :: matrixindex_in_compressed, ierr
  real(kind=8) :: sumn, trace_sparse

      sumn=0.d0
      if (orbs%norbp>0) then
          isegstart=amat%istsegline(orbs%isorb_par(iproc)+1)
          if (orbs%isorb+orbs%norbp<orbs%norb) then
              isegend=amat%istsegline(orbs%isorb_par(iproc+1)+1)-1
          else
              isegend=amat%nseg
          end if
          !$omp parallel default(private) shared(isegstart, isegend, orbs, bmat, amat, sumn)
          !$omp do reduction(+:sumn)
          do iseg=isegstart,isegend
              ii=amat%keyv(iseg)-1
              do jorb=amat%keyg(1,iseg),amat%keyg(2,iseg)
                  ii=ii+1
                  iiorb = (jorb-1)/orbs%norb + 1
                  jjorb = jorb - (iiorb-1)*orbs%norb
                  iilarge = matrixindex_in_compressed(bmat, iiorb, jjorb)
                  sumn = sumn + amat%matrix_compr(ii)*bmat%matrix_compr(iilarge)
              end do  
          end do
          !$omp end do
          !$omp end parallel
      end if
      call mpiallred(sumn, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      trace_sparse = sumn

end function trace_sparse
