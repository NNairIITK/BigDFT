!> @file
!! Fermi Operator Expansion Method
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Could still do more tidying - assuming all sparse matrices except for Fermi have the same pattern
subroutine foe(iproc, nproc, tmprtr, &
           ebs, itout, it_scc, order_taylor, max_inversion_error, purification_quickreturn, foe_verbosity, &
           accuracy_level, tmb, foe_obj)
  use module_base
  use module_types
  use module_interfaces, except_this_one => foe
  use yaml_output
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, &
                               matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed, &
                          uncompress_matrix_distributed, orb_from_index
  use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                      foe_data_get_logical
  use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                         fermilevel_get_real, fermilevel_get_logical
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc,itout,it_scc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  real(kind=8),intent(in) :: tmprtr
  real(kind=8),intent(out) :: ebs
  logical,intent(in) :: purification_quickreturn
  integer,intent(in) :: foe_verbosity
  integer,intent(in) :: accuracy_level
  type(DFT_wavefunction),intent(inout) :: tmb
  type(foe_data),intent(inout) :: foe_obj

  ! Local variables
  integer :: npl, jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
  integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
  integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
  integer,parameter :: nplx=50000
  real(kind=8),dimension(:,:),allocatable :: cc, chebyshev_polynomials, cc_check, fermip_check
  real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value, sumn, sumn_check, charge_diff, ef_interpol, ddot
  real(kind=8) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
  real(kind=8) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
  logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible, emergency_stop
  real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr, SHS, fermi_check_compr
  real(kind=8),dimension(4,4) :: interpol_matrix
  real(kind=8),dimension(4) :: interpol_vector
  real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion
  logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
  real(kind=8) :: trace_sparse, temp_multiplicator, ebs_check, ef, ebsp
  integer :: irow, icol, itemp, iflag,info, ispin, isshift, imshift, ilshift, i, j
  logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
  real(kind=8),parameter :: FSCALE_LOWER_LIMIT=5.d-3
  real(kind=8),parameter :: FSCALE_UPPER_LIMIT=5.d-2
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
  real(kind=8),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
  real(kind=8),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
  real(kind=8),parameter :: CHECK_RATIO=1.25d0
  integer,parameter :: NPL_MIN=100
  type(matrices) :: inv_ovrlp
  integer,parameter :: NTEMP_ACCURATE=4
  integer,parameter :: NTEMP_FAST=1
  real(kind=8) :: degree_multiplicator
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2
  integer,parameter :: imode=SPARSE
  type(fermi_aux) :: f
  


  call f_routine(id='foe')

  if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')

  if (accuracy_level/=FOE_ACCURATE .and. accuracy_level/=FOE_FAST) then
      stop 'wrong value of accuracy_level'
  end if

  inv_ovrlp%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, &
                           iaction=SPARSE_FULL, id='inv_ovrlp%matrix_compr')


  call timing(iproc, 'FOE_auxiliary ', 'ON')


  evbounds_shrinked=.false.


  penalty_ev = f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp,2/),id='penalty_ev')
  fermip_check = f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp/),id='fermip_check')
  SHS = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='SHS')
  fermi_check_compr = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='fermi_check_compr')


  call timing(iproc, 'FOE_auxiliary ', 'OF')
  call overlap_minus_onehalf() ! has internal timer
  call timing(iproc, 'FOE_auxiliary ', 'ON')


  hamscal_compr = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='hamscal_compr')

    
  ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
  nsize_polynomial=0
  if (tmb%linmat%l%smmm%nfvctrp>0) then
      isegstart = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+1)
      isegend = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp) + &
                tmb%linmat%l%nsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp)-1
      !!isegstart=tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr+1)
      !!if (tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp<tmb%linmat%l%nfvctr) then
      !!    isegend=tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr_par(iproc+1)+1)-1
      !!else
      !!    isegend=tmb%linmat%l%nseg
      !!end if
      !$omp parallel default(private) shared(isegstart, isegend, tmb, nsize_polynomial)
      !$omp do reduction(+:nsize_polynomial)
      do iseg=isegstart,isegend
          do jorb=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
              nsize_polynomial=nsize_polynomial+1
          end do
      end do
      !$omp end do
      !$omp end parallel
  end if
  
  
  ! Fake allocation, will be modified later
  chebyshev_polynomials = f_malloc((/nsize_polynomial,1/),id='chebyshev_polynomials')


  ! try to decrease the eigenvalue spectrum a bit
  if (foe_data_get_int(foe_obj,"evbounds_isatur")>foe_data_get_int(foe_obj,"evbounds_nsatur") .and. &
      foe_data_get_int(foe_obj,"evboundsshrink_isatur")<=foe_data_get_int(foe_obj,"evboundsshrink_nsatur")) then
      do ispin=1,tmb%linmat%l%nspin
          call foe_data_set_real(foe_obj,"evlow",0.9d0*foe_data_get_real(foe_obj,"evlow",ispin),ispin)
          call foe_data_set_real(foe_obj,"evhigh",0.9d0*foe_data_get_real(foe_obj,"evhigh",ispin),ispin)
      end do
      evbounds_shrinked=.true.
  else
      evbounds_shrinked=.false.
  end if

  ! This is to distinguish whether the routine is called from get_coeff of
  ! getLocBasis, to be improved.
  if (accuracy_level==FOE_ACCURATE) then
      ntemp = NTEMP_ACCURATE
      degree_multiplicator = DEGREE_MULTIPLICATOR_ACCURATE
      temp_multiplicator = TEMP_MULTIPLICATOR_ACCURATE
  else if (accuracy_level==FOE_FAST) then
      ntemp = NTEMP_FAST
      degree_multiplicator = DEGREE_MULTIPLICATOR_FAST
      temp_multiplicator = TEMP_MULTIPLICATOR_FAST
  else
      stop 'wrong value of accuracy_level'
  end if

  fscale_new=1.d100

  ebs=0.d0

  spin_loop: do ispin=1,tmb%linmat%l%nspin

      isshift=(ispin-1)*tmb%linmat%s%nvctr
      imshift=(ispin-1)*tmb%linmat%m%nvctr
      ilshift=(ispin-1)*tmb%linmat%l%nvctr

      degree_sufficient=.true.

      fscale_new = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")

      temp_loop: do itemp=1,ntemp

          fscale = fscale_new
          fscale = max(fscale,FSCALE_LOWER_LIMIT)
          fscale = min(fscale,FSCALE_UPPER_LIMIT)
          fscale_check = CHECK_RATIO*fscale
    
          evlow_old=1.d100
          evhigh_old=-1.d100
          
          if (iproc==0) then
              call yaml_map('decay length of error function',fscale,fmt='(es10.3)')
              call yaml_map('decay length multiplicator',temp_multiplicator,fmt='(es10.3)')
              call yaml_map('polynomial degree multiplicator',degree_multiplicator,fmt='(es10.3)')
          end if
    
        
              ! Don't let this value become too small.
              call foe_data_set_real(foe_obj,"bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",ispin),1.d-4),ispin)
        
              efarr(1)=foe_data_get_real(foe_obj,"ef",ispin)-foe_data_get_real(foe_obj,"bisection_shift",ispin)
              efarr(2)=foe_data_get_real(foe_obj,"ef",ispin)+foe_data_get_real(foe_obj,"bisection_shift",ispin)

              sumnarr(1)=0.d0
              sumnarr(2)=1.d100
              call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
                   foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                   foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
              call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)
        
              adjust_lower_bound=.true.
              adjust_upper_bound=.true.
        
              calculate_SHS=.true.
        
              if (tmb%linmat%l%smmm%nfvctrp>0) then
                  call to_zero(tmb%linmat%l%nfvctr*tmb%linmat%l%smmm%nfvctrp*tmb%linmat%l%nspin,tmb%linmat%kernel_%matrixp(1,1,1))
              end if
        
              if (iproc==0) then
                  !call yaml_sequence(advance='no')
                  if (foe_verbosity>=1) then
                      call yaml_sequence_open('FOE to determine density kernel',label=&
                           'it_foe'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//&
                           trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)')))//'-'//&
                           trim(adjustl(yaml_toa(itemp,fmt='(i2.2)')))//'-'//&
                           trim(adjustl(yaml_toa(ispin,fmt='(i2.2)'))))
                  else
                      call yaml_sequence_open('FOE to determine density kernel')
                      if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='-')
                  end if
              end if
        
        
        
              it=0
              eval_bounds_ok=.false.
              bisection_bounds_ok=.false.
              main_loop: do 
                  
                  it=it+1
        
                  if (iproc==0) then
                      call yaml_newline()
                      call yaml_sequence(advance='no')
                      call yaml_mapping_open(flow=.true.)
                      if (foe_verbosity>=1) call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
                  end if
                  
    !!              if (adjust_lower_bound) then
    !!                  call foe_data_set_real(foe_obj,"ef",efarr(1))
    !!              else if (adjust_upper_bound) then
    !!                  call foe_data_set_real(foe_obj,"ef",efarr(2))
    !!              end if
              
        
                  ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
                  if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                      call scale_and_shift_hamiltonian()
                      calculate_SHS=.true.
                  else
                      calculate_SHS=.false.
                  end if
                  evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
                  evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)
        
        
                  ! Determine the degree of the polynomial
                  npl=nint(degree_multiplicator* &
                      (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))/fscale)
                  npl=max(npl,NPL_MIN)
                  !npl_check = nint(degree_multiplicator*(foe_data_get_real(foe_obj,"evhigh")-foe_data_get_real(foe_obj,"evlow"))/fscale_check)
                  !npl_check = max(npl_check,nint(real(npl,kind=8)/CHECK_RATIO)) ! this is necessary if npl was set to the minimal value
                  npl_check = nint(real(npl,kind=8)/CHECK_RATIO)
                  npl_boundaries = nint(degree_multiplicator* &
                      (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)) &
                          /foe_data_get_real(foe_obj,"fscale_lowerbound")) ! max polynomial degree for given eigenvalue boundaries
                  if (npl>npl_boundaries) then
                      npl=npl_boundaries
                      if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
                      if (iproc==0) write(*,*) 'STOP SINCE THIS WILL CREATE PROBLEMS WITH NPL_CHECK'
                      stop
                  end if
                  if (npl>nplx) stop 'npl>nplx'
        
                  ! Array the holds the Chebyshev polynomials. Needs to be recalculated
                  ! every time the Hamiltonian has been modified.
                  if (calculate_SHS) then
                      call f_free(chebyshev_polynomials)
                      chebyshev_polynomials = f_malloc((/nsize_polynomial,npl/),id='chebyshev_polynomials')
                  end if
        
                  !if (foe_verbosity>=1 .and. iproc==0) then
                  if (iproc==0) then
                      if (foe_verbosity>=1) then
                          call yaml_map('bisec/eval bounds',&
                               (/fermilevel_get_real(f,"efarr(1)"),fermilevel_get_real(f,"efarr(2)"),&
                               foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                      else
                          call yaml_map('eval bounds',&
                               (/foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                      end if
                      call yaml_map('pol deg',npl,fmt='(i3)')
                      if (foe_verbosity>=1) call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es16.9)')
                  end if
        
        
                  cc = f_malloc((/npl,3/),id='cc')
                  cc_check = f_malloc((/npl,3/),id='cc_check')
        
                  if (foe_data_get_real(foe_obj,"evlow",ispin)>=0.d0) then
                      stop 'ERROR: lowest eigenvalue must be negative'
                  end if
                  if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      stop 'ERROR: highest eigenvalue must be positive'
                  end if
        
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
                  call timing(iproc, 'chebyshev_coef', 'ON')
        
                  call chebft(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,1), &
                       foe_data_get_real(foe_obj,"ef",ispin), fscale, tmprtr)
                  call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), cc(1,1), cc(1,2), npl)
                  call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,3))
                  call evnoise(npl, cc(1,3), foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
    
                  call chebft(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,1), &
                       foe_data_get_real(foe_obj,"ef",ispin), fscale_check, tmprtr)
                  call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), &
                       cc_check(1,1), cc_check(1,2), npl_check)
                  call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,3))
        
                  call timing(iproc, 'chebyshev_coef', 'OF')
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
                
                  !!if (iproc==0) then
                  !!    call pltwght(npl,cc(1,1),cc(1,2),foe_data_get_real(foe_obj,"evlow"),foe_data_get_real(foe_obj,"evhigh"),foe_data_get_real(foe_obj,"ef"),foe_data_get_real(foe_obj,"fscale"),temperature)
                  !!    call pltexp(anoise,npl,cc(1,3),foe_data_get_real(foe_obj,"evlow"),foe_data_get_real(foe_obj,"evhigh"))
                  !!end if
                
                
                  if (tmb%linmat%l%nspin==1) then
                      do ipl=1,npl
                          cc(ipl,1)=2.d0*cc(ipl,1)
                          cc(ipl,2)=2.d0*cc(ipl,2)
                          cc(ipl,3)=2.d0*cc(ipl,3)
                          cc_check(ipl,1)=2.d0*cc_check(ipl,1)
                          cc_check(ipl,2)=2.d0*cc_check(ipl,2)
                          cc_check(ipl,3)=2.d0*cc_check(ipl,3)
                      end do
                  end if
                
                
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
        
                  emergency_stop=.false.
                  if (calculate_SHS) then
                      ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                      if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
          !!do i=1,tmb%linmat%l%nvctr
          !!    write(1000+iproc,'(a,2i8,es16.6)') 'ispin, i, val', ispin, i, hamscal_compr(i)
          !!end do
                      call chebyshev_clean(iproc, nproc, npl, cc, &
                           tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nfvctrp, tmb%linmat%l%smmm%isfvctr, foe_obj, &
                           tmb%linmat%l, hamscal_compr, &
                           inv_ovrlp%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), calculate_SHS, &
                           nsize_polynomial, SHS, tmb%linmat%kernel_%matrixp, penalty_ev, chebyshev_polynomials, &
                           emergency_stop)
                  else
                      ! The Chebyshev polynomials are already available
                      if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                           tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nfvctrp, tmb%linmat%l%smmm%isfvctr, &
                          tmb%linmat%l, chebyshev_polynomials, cc, tmb%linmat%kernel_%matrixp)
                  end if 
    
    
    
                 call check_emergency_stop()
                 if (emergency_stop) then
                      eval_bounds_ok(1)=.false.
                      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*1.2d0,ispin)
                      eval_bounds_ok(2)=.false.
                      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                      if (iproc==0) then
                          if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                               (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          call yaml_mapping_close()
                          !call bigdft_utils_flush(unit=6)
                      end if
                      call f_free(cc)
                      call f_free(cc_check)
                      cycle main_loop
                 end if
        
        
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
        
        
                  restart=.false.
        
                  ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                  ! (otherwise this has already been checked in the previous iteration).
                  if (calculate_SHS) then
                      !!call check_eigenvalue_spectrum()
                      call check_eigenvalue_spectrum(nproc, tmb%linmat%l, tmb%linmat%s, tmb%linmat%ovrlp_, ispin, &
                            isshift, 1.2d0, 1.2d0, penalty_ev, anoise, .true., &
                            foe_obj, restart, eval_bounds_ok)
                  end if
        
                  call f_free(cc)
        
                  if (restart) then
                      if(evbounds_shrinked) then
                          ! this shrink was not good, increase the saturation counter
                          call foe_data_set_int(foe_obj,"evboundsshrink_isatur",foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
                      end if
                      call foe_data_set_int(foe_obj,"evbounds_isatur",0)
                      if (iproc==0) then
                          if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                               (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          call yaml_mapping_close()
                          !call bigdft_utils_flush(unit=6)
                      end if
                      call f_free(cc_check)
                      cycle
                  end if
                      
                  ! eigenvalue bounds ok
                  if (calculate_SHS) then
                      call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
                  end if
                
                  call calculate_trace_distributed(tmb%linmat%kernel_%matrixp, sumn)
        
    
                  if (all(eval_bounds_ok) .and. all(bisection_bounds_ok)) then
                      ! Print these informations already now if all entries are true.
                      if (iproc==0) then
                          if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                               (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      end if
                  end if
                  call determine_fermi_level(f, sumn, ef, info)
                  bisection_bounds_ok(1) = fermilevel_get_logical(f,"bisection_bounds_ok(1)")
                  bisection_bounds_ok(2) = fermilevel_get_logical(f,"bisection_bounds_ok(2)")
                  !write(*,*) 'main: efarr', efarr
                  if (info<0) then
                      if (iproc==0) then
                          if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                               (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          call yaml_mapping_close()
                      end if
                      call f_free(cc_check)
                      ! Save the new fermi energy in the foe_obj structure
                      call foe_data_set_real(foe_obj,"ef",ef,ispin)
                      cycle
                  end if
    
                  ! Save the new fermi energy and bisection_shift in the foe_obj structure
                  call foe_data_set_real(foe_obj,"ef",ef,ispin)
                  call foe_data_set_real(foe_obj,"bisection_shift",fermilevel_get_real(f,"bisection_shift"),ispin)
    
                  charge_diff = sumn-foe_data_get_real(foe_obj,"charge",ispin)
        
    
                  ef_old=foe_data_get_real(foe_obj,"ef",ispin)
                  sumn_old=sumn
    
    
        
                  if (iproc==0) then
                      if (foe_verbosity>=1) call yaml_newline()
                      if (foe_verbosity>=1) call yaml_map('iter',it)
                      if (foe_verbosity>=1) call yaml_map('Tr(K)',sumn,fmt='(es16.9)')
                      call yaml_map('charge diff',sumn-foe_data_get_real(foe_obj,"charge",ispin),fmt='(es16.9)')
                  end if
        
                  if (iproc==0) then
                      call yaml_mapping_close()
                      !call bigdft_utils_flush(unit=6)
                  end if
        
                  if (abs(charge_diff)<charge_tolerance) then
                      if (iproc==0) call yaml_sequence_close()
                      ! experimental: calculate a second kernel with a lower
                      ! polynomial degree  and calculate the difference
                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl_check, &
                           tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nfvctrp, tmb%linmat%l%smmm%isfvctr, &
                           tmb%linmat%l, chebyshev_polynomials, cc_check, fermip_check)
                      call f_free(cc_check)
                      diff=0.d0
                      do iorb=1,tmb%linmat%l%smmm%nfvctrp
                          do jorb=1,tmb%linmat%l%nfvctr
                              !SM: need to fix the spin here
                              diff = diff + (tmb%linmat%kernel_%matrixp(jorb,iorb,1)-fermip_check(jorb,iorb))**2
                          end do
                      end do
    
                      if (nproc > 1) then
                          call mpiallred(diff, 1, mpi_sum, bigdft_mpi%mpi_comm)
                      end if
    
                      diff=sqrt(diff)
                      if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')
                      exit
                  end if
    
                  call f_free(cc_check)
        
        
              end do main_loop
        
        
        
    
         call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
              tmb%linmat%kernel_%matrixp(:,1:tmb%linmat%l%smmm%nfvctrp,1), &
              tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr))
    
         call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, fermip_check, fermi_check_compr(1))
    
    
        
        
          ! Calculate S^-1/2 * K * S^-1/2^T
          ! Since S^-1/2 is symmetric, don't use the transpose
          call retransform(tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr))

          !!do i=ilshift+1,ilshift+tmb%linmat%l%nvctr
          !!    write(3000+iproc,'(a,2i8,es16.6)') 'ispin, i, val', ispin, i, tmb%linmat%kernel_%matrix_compr(i)
          !!end do
    
          call retransform(fermi_check_compr)
    
          call calculate_trace_distributed(fermip_check, sumn_check)

          !@NEW ##########################
          sumn = trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                 tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
                 tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)
          sumn_check = trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                       tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
                       fermi_check_compr, ispin)
          !@ENDNEW #######################
        
    
          ! Calculate trace(KH). Since they have the same sparsity pattern and K is
          ! symmetric, this is a simple ddot.
          ebsp=ddot(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(ilshift+1),1 , hamscal_compr, 1)
          ebsp=ebsp/scale_factor+shift_value*sumn
    
          ebs_check=ddot(tmb%linmat%l%nvctr, fermi_check_compr,1 , hamscal_compr, 1)
          ebs_check=ebs_check/scale_factor+shift_value*sumn_check
          diff=abs(ebs_check-ebsp)
          diff=diff/abs(ebsp)
    
          if (iproc==0) then
              call yaml_map('ebs',ebsp)
              call yaml_map('ebs_check',ebs_check)
              call yaml_map('diff',ebs_check-ebsp)
              call yaml_map('relative diff',diff)
          end if
    
          if (foe_data_get_logical(foe_obj,"adjust_FOE_temperature") .and. foe_verbosity>=1) then
              if (diff<5.d-5) then
                  ! can decrease polynomial degree
                  !!call foe_data_set_real(foe_obj,"fscale", 1.25d0*foe_data_get_real(foe_obj,"fscale"))
                  if (iproc==0) call yaml_map('modify fscale','increase')
                  !fscale_new=min(fscale_new,1.25d0*foe_data_get_real(foe_obj,"fscale"))
                  fscale_new=1.25d0*fscale_new
                  degree_sufficient=.true.
              else if (diff>=5.d-5 .and. diff < 1.d-4) then
                  ! polynomial degree seems to be appropriate
                  degree_sufficient=.true.
                  if (iproc==0) call yaml_map('modify fscale','No')
                  !fscale_new=min(fscale_new,foe_data_get_real(foe_obj,"fscale"))
                  fscale_new=fscale_new
              else
                  ! polynomial degree too small, increase and recalculate
                  ! the kernel
                  degree_sufficient=.false.
                  !!call foe_data_set_real(foe_obj,"fscale", 0.5*foe_data_get_real(foe_obj,"fscale"))
                  if (iproc==0) call yaml_map('modify fscale','decrease')
                  !fscale_new=min(fscale_new,0.5d0*foe_data_get_real(foe_obj,"fscale"))
                  fscale_new=0.5d0*fscale_new
              end if
              !if (foe_data_get_real(foe_obj,"fscale")<foe_data_get_real(foe_obj,"fscale_lowerbound")) then
              if (fscale_new<foe_data_get_real(foe_obj,"fscale_lowerbound")) then
                  !call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_lowerbound"))
                  fscale_new=foe_data_get_real(foe_obj,"fscale_lowerbound")
                  if (iproc==0) call yaml_map('fscale reached lower limit; reset to', &
                      foe_data_get_real(foe_obj,"fscale_lowerbound"))
                  reached_limit=.true.
              !else if (foe_data_get_real(foe_obj,"fscale")>foe_data_get_real(foe_obj,"fscale_upperbound")) then
              else if (fscale_new>foe_data_get_real(foe_obj,"fscale_upperbound")) then
                  !call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_upperbound"))
                  fscale_new=foe_data_get_real(foe_obj,"fscale_upperbound")
                  if (iproc==0) call yaml_map('fscale reached upper limit; reset to', &
                      foe_data_get_real(foe_obj,"fscale_upperbound"))
                  reached_limit=.true.
              else
                  reached_limit=.false.
              end if
          end if
        
    
        
      
          ! Purify the kernel
          !tmb%can_use_transposed=.true.
    
          if (.not.purification_quickreturn) then
              if (iproc==0) then
                  call yaml_sequence_open('Final kernel purification')
                  call yaml_newline()
              end if
              overlap_calculated=.true.
              if (itemp==ntemp) then
                  it_shift=20
              else
                  it_shift=1
              end if
              call purify_kernel(iproc, nproc, tmb, overlap_calculated, it_shift, 50, &
                   order_taylor, max_inversion_error, purification_quickreturn, ispin)
              if (iproc==0) then
                  call yaml_sequence_close()
              end if
          end if
        
        
          ! Calculate trace(KS).
          sumn = trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                 tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
                 tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)


          ! Recalculate trace(KH) (needed since the kernel was modified in the above purification). 
          ! If no purification is done, this should not be necessary.
          ! Since K and H have the same sparsity pattern and K is
          ! symmetric, the trace is a simple ddot.
          ebsp=ddot(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(ilshift+1),1 , hamscal_compr, 1)
          ebsp=ebsp/scale_factor+shift_value*sumn
    
    
          if (iproc==0) call yaml_map('trace(KS)',sumn)
    
    
          if (foe_verbosity>=1 .and. iproc==0) then
              call yaml_map('need to repeat with sharper decay (new)',.not.degree_sufficient)
          end if
          if (degree_sufficient) exit temp_loop
          if (reached_limit) then
              if (iproc==0) call yaml_map('limit reached, exit loop',.true.)
              exit temp_loop
          end if
    
    
        
    
      end do temp_loop

      ! Sum up the band structure energy
      ebs = ebs + ebsp

  end do spin_loop


  if (foe_data_get_logical(foe_obj,"adjust_FOE_temperature") .and. foe_verbosity>=1) then
      call foe_data_set_real(foe_obj,"fscale",fscale_new)
   !!   if (diff<5.d-3) then
   !!       ! can decrease polynomial degree
   !!       !!call foe_data_set_real(foe_obj,"fscale", 1.25d0*foe_data_get_real(foe_obj,"fscale"))
   !!       !!if (iproc==0) call yaml_map('modify fscale','increase')
   !!       fscale_new=min(fscale_new,1.25d0*foe_data_get_real(foe_obj,"fscale"))
   !!       degree_sufficient=.true.
   !!   else if (diff>=5.d-3 .and. diff < 1.d-2) then
   !!       ! polynomial degree seems to be appropriate
   !!       degree_sufficient=.true.
   !!       !!if (iproc==0) call yaml_map('modify fscale','No')
   !!       fscale_new=min(fscale_new,foe_data_get_real(foe_obj,"fscale"))
   !!   else
   !!       ! polynomial degree too small, increase and recalculate
   !!       ! the kernel
   !!       degree_sufficient=.false.
   !!       !!call foe_data_set_real(foe_obj,"fscale", 0.5*foe_data_get_real(foe_obj,"fscale"))
   !!       !!if (iproc==0) call yaml_map('modify fscale','decrease')
   !!       fscale_new=min(fscale_new,0.5d0*foe_data_get_real(foe_obj,"fscale"))
   !!   end if
   !!   !!if (foe_data_get_real(foe_obj,"fscale")<foe_data_get_real(foe_obj,"fscale_lowerbound")) then
   !!   !!    call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_lowerbound"))
   !!   !!    if (iproc==0) call yaml_map('fscale reached lower limit; reset to',foe_data_get_real(foe_obj,"fscale_lowerbound"))
   !!   !!    reached_limit=.true.
   !!   !!else if (foe_data_get_real(foe_obj,"fscale")>foe_data_get_real(foe_obj,"fscale_upperbound")) then
   !!   !!    call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_upperbound"))
   !!   !!    if (iproc==0) call yaml_map('fscale reached upper limit; reset to',foe_data_get_real(foe_obj,"fscale_upperbound"))
   !!   !!    reached_limit=.true.
   !!   !!else
   !!   !!    reached_limit=.false.
   !!   !!end if
  end if

  degree_sufficient=.true.


  call f_free_ptr(inv_ovrlp%matrix_compr)
  


!!  scale_factor=1.d0/scale_factor
!!  shift_value=-shift_value
!!
!!
!!  ! Calculate trace(KH). Since they have the same sparsity pattern and K is
!!  ! symmetric, this is a simple ddot.
!!  ebs=ddot(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(ilshift+1),1 , hamscal_compr, 1)
!!  ebs=ebs*scale_factor-shift_value*sumn


  if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


  call f_free(chebyshev_polynomials)
  call f_free(penalty_ev)
  call f_free(hamscal_compr)
  call f_free(fermip_check)
  call f_free(SHS)
  call f_free(fermi_check_compr)

  call timing(iproc, 'FOE_auxiliary ', 'OF')

  call f_release_routine()



      contains

        subroutine overlap_minus_onehalf()
          implicit none
          real(kind=8) :: max_error, mean_error
          integer :: i, j, ii

          call f_routine(id='overlap_minus_onehalf')

          ! Taylor approximation of S^-1/2 up to higher order
          if (imode==DENSE) then
              stop 'overlap_minus_onehalf: DENSE is deprecated'
              tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, &
                                         id='tmb%linmat%ovrlp_%matrix')
              call uncompress_matrix(iproc, tmb%linmat%s, &
                   inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)

              inv_ovrlp%matrix=sparsematrix_malloc_ptr(tmb%linmat%l, &
                                                iaction=DENSE_FULL, id='inv_ovrlp%matrix')
              call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, &
                   imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
                   ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
                   check_accur=.true., max_error=max_error, mean_error=mean_error)
              call compress_matrix(iproc, tmb%linmat%l, inmat=inv_ovrlp%matrix, outmat=inv_ovrlp%matrix_compr)
          end if
          if (imode==SPARSE) then
              call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, &
                   imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
                   ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
                   check_accur=.true., max_error=max_error, mean_error=mean_error)
          end if
          !!ii=0
          !!do i=1,tmb%linmat%l%nspin
          !!    do j=1,tmb%linmat%l%nvctr
          !!        ii=ii+1
          !!        write(2000+iproc*10+i,'(a,3i8,es16.6)') 'i, j, ii, val', i, j, ii, inv_ovrlp%matrix_compr(ii)
          !!    end do
          !!end do
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)
          !if (foe_verbosity>=1 .and. iproc==0) then
          !    call yaml_map('max error of S^-1/2',max_error,fmt='(es9.2)')
          !    call yaml_map('mean error of S^-1/2',mean_error,fmt='(es9.2)')
          !end if


          if (imode==DENSE) then
              call f_free_ptr(inv_ovrlp%matrix)

              call f_free_ptr(tmb%linmat%ovrlp_%matrix)
          end if

          call f_release_routine()
      end subroutine overlap_minus_onehalf



      subroutine retransform(matrix_compr)
          use sparsematrix, only: sequential_acces_matrix_fast, sparsemm
          ! Calling arguments
          real(kind=8),dimension(tmb%linmat%l%nvctr),intent(inout) :: matrix_compr

          ! Local variables
          real(kind=8),dimension(:,:),pointer :: inv_ovrlpp, tempp
          integer,dimension(:,:),pointer :: onedimindices
          real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
          integer,dimension(:),allocatable :: ivectorindex
          integer,dimension(:,:,:),allocatable :: istindexarr
          integer :: nout, nseq

          call f_routine(id='retransform')

          inv_ovrlpp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
          tempp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
          inv_ovrlp_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
          kernel_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
          call sequential_acces_matrix_fast(tmb%linmat%l, matrix_compr, kernel_compr_seq)
          call sequential_acces_matrix_fast(tmb%linmat%l, &
               inv_ovrlp%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), inv_ovrlp_compr_seq)
          call uncompress_matrix_distributed(iproc, tmb%linmat%l, DENSE_MATMUL, &
               inv_ovrlp%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), inv_ovrlpp)

           tempp=0.d0
          call sparsemm(tmb%linmat%l, kernel_compr_seq, inv_ovrlpp, tempp)
          inv_ovrlpp=0.d0
          call sparsemm(tmb%linmat%l, inv_ovrlp_compr_seq, tempp, inv_ovrlpp)

          call to_zero(tmb%linmat%l%nvctr, matrix_compr(1))
          call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, inv_ovrlpp, matrix_compr)

          call f_free_ptr(inv_ovrlpp)
          call f_free_ptr(tempp)
          call f_free(inv_ovrlp_compr_seq)
          call f_free(kernel_compr_seq)

          call f_release_routine()

      end subroutine retransform




      subroutine calculate_trace_distributed(matrixp, trace)
          real(kind=8),dimension(tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp),intent(in) :: matrixp
          real(kind=8),intent(out) :: trace

          call f_routine(id='calculate_trace_distributed')

          trace=0.d0
          if (tmb%linmat%l%smmm%nfvctrp>0) then
              isegstart = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+1)
              isegend = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp) + &
                        tmb%linmat%l%nsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp)-1
              !isegstart=tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr+1)
              !if (tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp<tmb%linmat%l%nfvctr) then
              !    isegend=tmb%linmat%l%istsegline(tmb%linmat%l%isfvctr_par(iproc+1)+1)-1
              !else
              !    isegend=tmb%linmat%l%nseg
              !end if

              !$omp parallel default(private) shared(isegstart, isegend, matrixp, tmb, trace) 
              !$omp do reduction(+:trace)
              do iseg=isegstart,isegend
                  ii=tmb%linmat%l%keyv(iseg)-1
                  do jorb=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/tmb%linmat%l%nfvctr + 1
                      jjorb = jorb - (iiorb-1)*tmb%linmat%l%nfvctr
                      if (jjorb==iiorb) trace = trace + matrixp(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr)
                  end do  
              end do
              !$omp end do
              !$omp end parallel
          end if
    
          if (nproc > 1) then
              call mpiallred(trace, 1, mpi_sum, bigdft_mpi%mpi_comm)
          end if

          call f_release_routine()
      end subroutine calculate_trace_distributed



!!      subroutine determine_new_fermi_level()
!!        implicit none
!!        integer :: info, i
!!        real(kind=8) :: determinant, m, b
!!        real(kind=8),dimension(4,4) :: tmp_matrix
!!        real(kind=8),dimension(4) :: interpol_solution
!!        integer,dimension(4) :: ipiv
!!
!!        ! Shift up the old results.
!!        if (it_solver>4) then
!!            do i=1,4
!!                interpol_matrix(1,i)=interpol_matrix(2,i)
!!                interpol_matrix(2,i)=interpol_matrix(3,i)
!!                interpol_matrix(3,i)=interpol_matrix(4,i)
!!            end do
!!            interpol_vector(1)=interpol_vector(2)
!!            interpol_vector(2)=interpol_vector(3)
!!            interpol_vector(3)=interpol_vector(4)
!!        end if
!!        !LG: if it_solver==0 this index comes out of bounds!
!!        ii=max(min(it_solver,4),1)
!!        interpol_matrix(ii,1)=foe_data_get_real(foe_obj,"ef")**3
!!        interpol_matrix(ii,2)=foe_data_get_real(foe_obj,"ef")**2
!!        interpol_matrix(ii,3)=foe_data_get_real(foe_obj,"ef")
!!        interpol_matrix(ii,4)=1
!!        interpol_vector(ii)=sumn-foe_data_get_real(foe_obj,"charge")
!!    
!!        ! Solve the linear system interpol_matrix*interpol_solution=interpol_vector
!!        if (it_solver>=4) then
!!            do i=1,ii
!!                interpol_solution(i)=interpol_vector(i)
!!                tmp_matrix(i,1)=interpol_matrix(i,1)
!!                tmp_matrix(i,2)=interpol_matrix(i,2)
!!                tmp_matrix(i,3)=interpol_matrix(i,3)
!!                tmp_matrix(i,4)=interpol_matrix(i,4)
!!            end do
!!    
!!            call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
!!            if (info/=0) then
!!               if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
!!            end if
!!    
!!    
!!            call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
!!                 interpol_solution(3), interpol_solution(4), foe_data_get_real(foe_obj,"ef"), ef_interpol)
!!        end if
!!    
!!    
!!    
!!    
!!        ! Calculate the new Fermi energy.
!!        if (foe_verbosity>=1 .and. iproc==0) then
!!            call yaml_newline()
!!            call yaml_mapping_open('Search new eF',flow=.true.)
!!        end if
!!        if (it_solver>=4 .and.  &
!!            abs(sumn-foe_data_get_real(foe_obj,"charge"))<foe_data_get_real(foe_obj,"ef_interpol_chargediff")) then
!!            det=determinant(iproc,4,interpol_matrix)
!!            if (foe_verbosity>=1 .and. iproc==0) then
!!                call yaml_map('det',det,fmt='(es10.3)')
!!                call yaml_map('limit',foe_data_get_real(foe_obj,"ef_interpol_det"),fmt='(es10.3)')
!!            end if
!!            if(abs(det)>foe_data_get_real(foe_obj,"ef_interpol_det")) then
!!                call foe_data_set_real(foe_obj,"ef",ef_interpol)
!!                if (foe_verbosity>=1 .and. iproc==0) call yaml_map('method','cubic interpolation')
!!            else
!!                ! linear interpolation
!!                if (foe_verbosity>=1 .and. iproc==0) call yaml_map('method','linear interpolation')
!!                m = (interpol_vector(4)-interpol_vector(3))/(interpol_matrix(4,3)-interpol_matrix(3,3))
!!                b = interpol_vector(4)-m*interpol_matrix(4,3)
!!                call foe_data_set_real(foe_obj,"ef", -b/m)
!!            end if
!!        else
!!            ! Use mean value of bisection and secant method
!!            ! Secant method solution
!!            call foe_data_set_real(foe_obj,"ef", &
!!                 efarr(2)-(sumnarr(2)-foe_data_get_real(foe_obj,"charge"))*(efarr(2)-efarr(1))/(sumnarr(2)-sumnarr(1)))
!!            ! Add bisection solution
!!            call foe_data_set_real(foe_obj,"ef", foe_data_get_real(foe_obj,"ef") + .5d0*(efarr(1)+efarr(2)))
!!            ! Take the mean value
!!            call foe_data_set_real(foe_obj,"ef", .5d0*foe_data_get_real(foe_obj,"ef"))
!!            if (foe_verbosity>=1 .and. iproc==0) call yaml_map('method','bisection / secant method')
!!        end if
!!        if (foe_verbosity>=1 .and. iproc==0) then
!!            call yaml_mapping_close()
!!            !!call bigdft_utils_flush(unit=6)
!!            !call yaml_newline()
!!        end if
!!
!!      end subroutine determine_new_fermi_level


      !!subroutine check_eigenvalue_spectrum()
      !!  implicit none
      !!  real(kind=8) :: bound_low, bound_up

      !!  call f_routine(id='check_eigenvalue_spectrum')

      !!  ! The penalty function must be smaller than the noise.
      !!  bound_low=0.d0
      !!  bound_up=0.d0
      !!  if (tmb%linmat%l%smmm%nfvctrp>0) then
      !!      isegstart = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+1)
      !!      isegend = tmb%linmat%l%istsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp) + &
      !!                tmb%linmat%l%nsegline(tmb%linmat%l%smmm%isfvctr+tmb%linmat%l%smmm%nfvctrp)-1

      !!      !$omp parallel default(private) &
      !!      !$omp shared(isegstart, isegend, penalty_ev, tmb, bound_low, bound_up, isshift)
      !!      !$omp do reduction(+:bound_low,bound_up)
      !!      do iseg=isegstart,isegend
      !!          ii=tmb%linmat%l%keyv(iseg)-1
      !!          do jorb=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
      !!              ii=ii+1
      !!              iiorb = (jorb-1)/tmb%linmat%l%nfvctr + 1
      !!              jjorb = jorb - (iiorb-1)*tmb%linmat%l%nfvctr
      !!              iismall = matrixindex_in_compressed(tmb%linmat%s, iiorb, jjorb)
      !!              if (iismall>0) then
      !!                  tt=tmb%linmat%ovrlp_%matrix_compr(isshift+iismall)
      !!              else
      !!                  tt=0.d0
      !!              end if
      !!              bound_low = bound_low + penalty_ev(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr,2)*tt
      !!              bound_up = bound_up +penalty_ev(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr,1)*tt
      !!          end do  
      !!      end do
      !!      !$omp end do
      !!      !$omp end parallel
      !!  end if
    
      !!  allredarr(1)=bound_low
      !!  allredarr(2)=bound_up

      !!  if (nproc > 1) then
      !!      call mpiallred(allredarr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
      !!  end if

      !!  allredarr=abs(allredarr) !for some crazy situations this may be negative
      !!  anoise=100.d0*anoise
      !!  if (allredarr(1)>anoise) then
      !!      eval_bounds_ok(1)=.false.
      !!      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*1.2d0,ispin)
      !!      restart=.true.
      !!  else
      !!      eval_bounds_ok(1)=.true.
      !!  end if
      !!  if (allredarr(2)>anoise) then
      !!      eval_bounds_ok(2)=.false.
      !!      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
      !!      restart=.true.
      !!  else
      !!      eval_bounds_ok(2)=.true.
      !!  end if

      !!  call f_release_routine()

      !!end subroutine check_eigenvalue_spectrum


      subroutine scale_and_shift_hamiltonian()
        implicit none

        integer,dimension(2) :: irowcol

        call f_routine(id='scale_and_shift_hamiltonian')

        scale_factor=2.d0/(foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
        shift_value=.5d0*(foe_data_get_real(foe_obj,"evhigh",ispin)+foe_data_get_real(foe_obj,"evlow",ispin))
        !$omp parallel default(none) private(iseg,ii,i,irowcol,iismall_ovrlp,iismall_ham,tt_ovrlp,tt_ham) &
        !$omp shared(tmb,hamscal_compr,scale_factor,shift_value,isshift,imshift)
        !$omp do
        do iseg=1,tmb%linmat%l%nseg
            !do ii=1,tmb%linmat%l%nvctr
            ii=tmb%linmat%l%keyv(iseg)
            do i=tmb%linmat%l%keyg(1,iseg),tmb%linmat%l%keyg(2,iseg)
                irowcol = orb_from_index(tmb%linmat%l,i)
                iismall_ovrlp = matrixindex_in_compressed(tmb%linmat%s, irowcol(1), irowcol(2))
                iismall_ham = matrixindex_in_compressed(tmb%linmat%m, irowcol(1), irowcol(2))
                if (iismall_ovrlp>0) then
                    tt_ovrlp=tmb%linmat%ovrlp_%matrix_compr(isshift+iismall_ovrlp)
                else
                    tt_ovrlp=0.d0
                end if
                if (iismall_ham>0) then
                    tt_ham=tmb%linmat%ham_%matrix_compr(imshift+iismall_ham)
                else
                    tt_ham=0.d0
                end if
                hamscal_compr(ii)=scale_factor*(tt_ham-shift_value*tt_ovrlp)
                ii=ii+1
            end do
        end do
        !$omp end do
        !$omp end parallel

        call f_release_routine()

      end subroutine scale_and_shift_hamiltonian


      subroutine check_emergency_stop()
        implicit none

        call f_routine(id='check_emergency_stop')

        ! Check for an emergency stop, which happens if the kernel explodes, presumably due
        ! to the eigenvalue bounds being too small.
        ! mpi_lor seems not to work on certain systems...
        if (emergency_stop) then
            iflag=1
        else
            iflag=0
        end if
    
        if (nproc > 1) then
            call mpiallred(iflag, 1, mpi_sum, bigdft_mpi%mpi_comm)
        end if
    
        if (iflag>0) then
            emergency_stop=.true.
        else
            emergency_stop=.false.
        end if

        call f_release_routine()

      end subroutine check_emergency_stop

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

  call f_routine(id='chebft')

  if (n>50000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  fac=2.d0/n
  !$omp parallel default(none) shared(bma,bpa,fac,n,tmprtr,cf,fscale,ef,cc) &
  !$omp private(k,y,arg,tt,j)
  !$omp do
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      if (tmprtr.eq.0.d0) then
          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
      else
          cf(k)=1.d0/(1.d0+exp( (arg-ef)*(1.d0/tmprtr) ) )
      end if
  end do
  !$omp end do
  !$omp do
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
  !$omp end do
  !$omp end parallel

  call f_release_routine()

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
  real(kind=8) :: tt, ttt, y, arg, fac, bma, bpa
  real(kind=8),dimension(50000) :: cf

  call f_routine(id='chebft2')

  if (n>50000) stop 'chebft2'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  ! 3 gives broder safety zone than 4
  !ttt=3.0d0*n/(b-a)
  ttt=4.d0*n/(b-a)
  fac=2.d0/n
  !$omp parallel default(none) shared(bma,bpa,ttt,fac,n,cf,b,cc) &
  !$omp private(k,y,arg,tt,j)
  !$omp do
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      cf(k)=exp((arg-b)*ttt)
  end do
  !$omp end do
  !$omp do
  do j=1,n
      tt=0.d0
      do k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
  !$omp end do
  !$omp end parallel

  call f_release_routine()

end subroutine chebft2


! Calculates chebychev expansion of the derivative of Fermi distribution.
subroutine chder(a,b,c,cder,n)
  use dynamic_memory
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(8),dimension(n),intent(in) :: c
  real(8),dimension(n),intent(out) :: cder

  ! Local variables
  integer :: j
  real(kind=8) :: con

  call f_routine(id='chder')

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

  call f_release_routine()

end subroutine chder


!> Determine noise level
subroutine evnoise(npl,cc,evlow,evhigh,anoise)
  use dynamic_memory
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: npl
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),intent(in) :: evlow, evhigh
  real(kind=8),intent(out) :: anoise
  
  ! Local variables
  integer :: i, n
  real(kind=8) :: fact, dist, ddx, cent, tt, x, chebev

  call f_routine(id='evnoise')
  
  fact=1.d0
  dist=(fact*evhigh-fact*evlow)
  ddx=dist/(10*npl)
  cent=.5d0*(fact*evhigh+fact*evlow)
  !!tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  !!do x=ddx,.25d0*dist,ddx
  !!    tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
  !!       & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
  !!end do
  ! Rewritten version of the above loop
  tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  x=ddx
  n=ceiling((0.25d0*dist-ddx)/ddx)
  !$omp parallel default(none) shared(n,ddx,tt,evlow,evhigh,npl,cent,cc) private(i,x)
  !$omp do reduction(max:tt)
  do i=1,n
      x=real(i,kind=8)*ddx
      tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
         & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
      !x=x+ddx
      !if (x>=.25d0*dist) exit
  end do
  !$omp end do
  !$omp end parallel
  !anoise=1.d0*tt
  anoise=20.d0*tt

  call f_release_routine()

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
  double complex :: test
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

  sol_c(1) = -b_c/(3*a_c) &
       - S_c/(3*a_c) &
       - (b_c**2-3*a_c*c_c)/(3*a_c*S_c)
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
      tt=abs(real(sol_c(i),kind=8)-target_solution)
      if (tt<ttmin) then
          ttmin=tt
          solution=real(sol_c(i),kind=8)
      end if
  end do

end subroutine get_roots_of_cubic_polynomial



real(kind=8) function determinant(iproc, n, mat)
    use module_base
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, n
    real(kind=8),dimension(n,n),intent(in) :: mat

    ! Local variables
    integer :: i, info
    integer,dimension(n) :: ipiv
    real(kind=8),dimension(n,n) :: mat_tmp
    real(kind=8) :: sgn

    call vcopy(n**2, mat(1,1), 1, mat_tmp(1,1), 1)

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


subroutine compress_polynomial_vector(iproc, nproc, nsize_polynomial, norb, norbp, isorb, &
           fermi, vector, vector_compressed)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp, isorb
  type(sparse_matrix),intent(in) :: fermi
  real(kind=8),dimension(norb,norbp),intent(in) :: vector
  real(kind=8),dimension(nsize_polynomial),intent(out) :: vector_compressed

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (norbp>0) then
      isegstart = fermi%istsegline(fermi%smmm%isfvctr+1)
      isegend = fermi%istsegline(fermi%smmm%isfvctr+fermi%smmm%nfvctrp) + &
                fermi%nsegline(fermi%smmm%isfvctr+fermi%smmm%nfvctrp)-1
      !!isegstart=fermi%istsegline(fermi%isfvctr+1)
      !!if (isorb+norbp<norb) then
      !!    isegend=fermi%istsegline(fermi%isfvctr_par(iproc+1)+1)-1
      !!else
      !!    isegend=fermi%nseg
      !!end if
      ii=0
      !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermi, vector, vector_compressed)
      !!$omp do
      do iseg=isegstart,isegend
          !ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/norb + 1
              jjorb = jorb - (iiorb-1)*norb
              vector_compressed(ii)=vector(jjorb,iiorb-isorb)
              !write(300,*) 'ii, jjorb, iiorb-isorb', ii, jjorb, iiorb-isorb
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if
end subroutine compress_polynomial_vector



subroutine uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
           norb, norbp, isorb, fermi, vector_compressed, vector)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp, isorb
  type(sparse_matrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial),intent(in) :: vector_compressed
  real(kind=8),dimension(fermi%nfvctr,fermi%smmm%nfvctrp),intent(out) :: vector

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (fermi%smmm%nfvctrp>0) then
      call to_zero(fermi%nfvctr*fermi%smmm%nfvctrp, vector(1,1))
      isegstart = fermi%istsegline(fermi%smmm%isfvctr+1)
      isegend = fermi%istsegline(fermi%smmm%isfvctr+fermi%smmm%nfvctrp) + &
                fermi%nsegline(fermi%smmm%isfvctr+fermi%smmm%nfvctrp)-1
      !!isegstart=fermi%istsegline(fermi%isfvctr+1)
      !!if (fermi%isfvctr+fermi%smmm%nfvctrp<fermi%nfvctr) then
      !!    isegend=fermi%istsegline(fermi%isfvctr_par(iproc+1)+1)-1
      !!else
      !!    isegend=fermi%nseg
      !!end if
      !$omp parallel do default(private) &
      !$omp shared(isegstart, isegend, fermi, vector, vector_compressed)
      do iseg=isegstart,isegend
          ii=fermi%keyv(iseg)-fermi%keyv(isegstart)
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/fermi%nfvctr + 1
              jjorb = jorb - (iiorb-1)*fermi%nfvctr
              vector(jjorb,iiorb-fermi%smmm%isfvctr)=vector_compressed(ii)
              !write(*,*) 'ii, iiorb-fermi%isfvctr, jjorb', ii, iiorb-fermi%isfvctr, jjorb
          end do
      end do
      !$omp end parallel do
  end if
end subroutine uncompress_polynomial_vector


!< Calculates the trace of the matrix product amat*bmat.
!< WARNING: It is mandatory that the sparsity pattern of amat is contained
!< within the sparsity pattern of bmat!
function trace_sparse(iproc, nproc, orbs, asmat, bsmat, amat, bmat, ispin)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix, matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,  nproc, ispin
  type(orbitals_data),intent(in) :: orbs
  type(sparse_matrix),intent(in) :: asmat, bsmat
  real(kind=8),dimension(asmat%nvctr),intent(in) :: amat
  real(kind=8),dimension(bsmat%nvctr),intent(in) :: bmat

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iilarge
  integer :: ierr, iashift, ibshift
  real(kind=8) :: sumn, trace_sparse


  call f_routine(id='trace_sparse')

  iashift = 0!(ispin-1)*asmat%nvctr
  ibshift = 0!(ispin-1)*bsmat%nvctr

  sumn=0.d0
  if (asmat%smmm%nfvctrp>0) then
      isegstart = asmat%istsegline(asmat%smmm%isfvctr+1)
      isegend = asmat%istsegline(asmat%smmm%isfvctr+asmat%smmm%nfvctrp) + &
                asmat%nsegline(asmat%smmm%isfvctr+asmat%smmm%nfvctrp)-1
      !!isegstart=asmat%istsegline(asmat%isfvctr+1)
      !!if (asmat%isfvctr+asmat%nfvctrp<asmat%nfvctr) then
      !!        isegend=asmat%istsegline(asmat%isfvctr_par(iproc+1)+1)-1
      !!else
      !!        isegend=asmat%nseg
      !!end if
      !$omp parallel default(none) &
      !$omp private(iseg, ii, jorb, iiorb, jjorb, iilarge) &
      !$omp shared(isegstart, isegend, bsmat, asmat, amat, bmat, iashift, ibshift, sumn)
      !$omp do reduction(+:sumn)
      do iseg=isegstart,isegend
              ii=iashift+asmat%keyv(iseg)-1
              do jorb=asmat%keyg(1,iseg),asmat%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/asmat%nfvctr + 1
              jjorb = jorb - (iiorb-1)*asmat%nfvctr
              iilarge = ibshift + matrixindex_in_compressed(bsmat, iiorb, jjorb)
              sumn = sumn + amat(ii)*bmat(iilarge)
          end do  
      end do
      !$omp end do
      !$omp end parallel
  end if

  if (nproc > 1) then
      call mpiallred(sumn, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  trace_sparse = sumn

  call f_release_routine()

end function trace_sparse






















! New: chebyshev expansion of the inverse overlap (Inverse Chebyshev Expansion)
subroutine ice(iproc, nproc, norder_polynomial, ovrlp_smat, inv_ovrlp_smat, ex, ovrlp_mat, inv_ovrlp)
  use module_base
  use module_types
  use module_interfaces, except_this_one => ice
  use yaml_output
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, &
                               sparsematrix_malloc0_ptr, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, &
                               matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed, orb_from_index
  use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                      foe_data_set_logical, foe_data_get_logical
  use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                         fermilevel_get_real, fermilevel_get_logical
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, norder_polynomial
  type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
  !type(sparse_matrix),intent(inout) :: ovrlp_smat, inv_ovrlp_smat !for debug inout
  integer :: ex
  type(matrices),intent(in) :: ovrlp_mat
  type(matrices),intent(out) :: inv_ovrlp

  ! Local variables
  integer :: npl, jorb, it, ii, iseg
  integer :: isegstart, isegend, iismall, nsize_polynomial
  integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i
  integer,parameter :: nplx=50000
  real(kind=8),dimension(:,:),allocatable :: cc, chebyshev_polynomials
  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_matrixp
  real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value
  real(kind=8) :: evlow_old, evhigh_old, tt
  real(kind=8) :: tt_ovrlp, tt_ham
  logical :: restart, calculate_SHS, emergency_stop
  real(kind=8),dimension(2) :: allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr, SHS
  logical,dimension(2) :: eval_bounds_ok
  integer,dimension(2) :: irowcol
  integer :: irow, icol, iflag, ispin, isshift, ilshift
  logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
  integer,parameter :: NPL_MIN=5
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_MAX=10.d0
  real(kind=8) :: degree_multiplicator
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2
  integer,parameter :: imode=SPARSE
  type(foe_data) :: foe_obj

  !!real(kind=8),dimension(ovrlp_smat%nfvctr,ovrlp_smat%nfvctr) :: overlap
  !!real(kind=8),dimension(ovrlp_smat%nfvctr) :: eval
  !!integer,parameter :: lwork=100000
  !!real(kind=8),dimension(lwork) :: work
  !!integer :: info

  call f_routine(id='ice')



!@ JUST FOR THE MOMENT.... ########################
     foe_obj%ef = f_malloc0_ptr(ovrlp_smat%nspin,id='(foe_obj%ef)')
     foe_obj%evlow = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evlow')
     foe_obj%evhigh = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evhigh')
     foe_obj%bisection_shift = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%bisection_shift')
     foe_obj%charge = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%charge')
     do ispin=1,ovrlp_smat%nspin
         call foe_data_set_real(foe_obj,"ef",0.d0,ispin)
         call foe_data_set_real(foe_obj,"evlow",0.5d0,ispin)
         call foe_data_set_real(foe_obj,"evhigh",2.5d0,ispin)
         call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,ispin)
         call foe_data_set_real(foe_obj,"charge",0.d0,ispin)
     end do

     call foe_data_set_real(foe_obj,"fscale",1.d-1)
     call foe_data_set_real(foe_obj,"ef_interpol_det",0.d0)
     call foe_data_set_real(foe_obj,"ef_interpol_chargediff",0.d0)
     call foe_data_set_int(foe_obj,"evbounds_isatur",0)
     call foe_data_set_int(foe_obj,"evboundsshrink_isatur",0)
     call foe_data_set_int(foe_obj,"evbounds_nsatur",10)
     call foe_data_set_int(foe_obj,"evboundsshrink_nsatur",10)
     call foe_data_set_real(foe_obj,"fscale_lowerbound",1.d-2)
     call foe_data_set_real(foe_obj,"fscale_upperbound",0.d0)
     call foe_data_set_logical(foe_obj,"adjust_FOE_temperature",.false.)
!@ ################################################




  call timing(iproc, 'FOE_auxiliary ', 'ON')



  penalty_ev = f_malloc((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,2/),id='penalty_ev')
  SHS = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='SHS')


  hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='hamscal_compr')

    
  ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
  nsize_polynomial=0
  if (inv_ovrlp_smat%smmm%nfvctrp>0) then
      isegstart = inv_ovrlp_smat%istsegline(inv_ovrlp_smat%smmm%isfvctr+1)
      isegend = inv_ovrlp_smat%istsegline(inv_ovrlp_smat%smmm%isfvctr+inv_ovrlp_smat%smmm%nfvctrp) + &
                inv_ovrlp_smat%nsegline(inv_ovrlp_smat%smmm%isfvctr+inv_ovrlp_smat%smmm%nfvctrp)-1
      !$omp parallel default(private) shared(isegstart, isegend, inv_ovrlp_smat, nsize_polynomial)
      !$omp do reduction(+:nsize_polynomial)
      do iseg=isegstart,isegend
          do jorb=inv_ovrlp_smat%keyg(1,iseg),inv_ovrlp_smat%keyg(2,iseg)
              nsize_polynomial=nsize_polynomial+1
          end do
      end do
      !$omp end do
      !$omp end parallel
  end if
  
  
  ! Fake allocation, will be modified later
  chebyshev_polynomials = f_malloc((/nsize_polynomial,1/),id='chebyshev_polynomials')


  inv_ovrlp_matrixp = sparsematrix_malloc0_ptr(inv_ovrlp_smat, &
                           iaction=DENSE_MATMUL, id='inv_ovrlp_matrixp')


      spin_loop: do ispin=1,ovrlp_smat%nspin

          degree_multiplicator = real(norder_polynomial,kind=8)/ &
                                 (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
          degree_multiplicator = min(degree_multiplicator,DEGREE_MULTIPLICATOR_MAX)

          isshift=(ispin-1)*ovrlp_smat%nvctr
          ilshift=(ispin-1)*inv_ovrlp_smat%nvctr

          evlow_old=1.d100
          evhigh_old=-1.d100
          
    
        
              !!calculate_SHS=.true.
        
              if (inv_ovrlp_smat%smmm%nfvctrp>0) then
                  call to_zero(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_matrixp(1,1))
              end if
        
        
              it=0
              eval_bounds_ok=.false.
              !!bisection_bounds_ok=.false.
              main_loop: do 
                  
                  it=it+1
        
                  ! Scale the Hamiltonian such that all eigenvalues are in the intervall [0:1]
                  if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                      shift_value=.5d0*(foe_data_get_real(foe_obj,"evhigh",ispin)+foe_data_get_real(foe_obj,"evlow",ispin))
                      scale_factor=2.d0/(foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
                      !$omp parallel default(none) private(iseg,i,ii,irowcol,iismall_ovrlp,iismall_ham,tt_ovrlp,tt_ham) &
                      !$omp shared(inv_ovrlp_smat,ovrlp_smat,hamscal_compr,ovrlp_mat,scale_factor,shift_value,isshift)
                      !$omp do
                      do iseg=1,inv_ovrlp_smat%nseg
                          ii = inv_ovrlp_smat%keyv(iseg)
                          do i=inv_ovrlp_smat%keyg(1,iseg),inv_ovrlp_smat%keyg(2,iseg)
                              irowcol = orb_from_index(inv_ovrlp_smat,i)
                              iismall_ovrlp = matrixindex_in_compressed(ovrlp_smat, irowcol(1), irowcol(2))
                              if (iismall_ovrlp>0) then
                                  if (irowcol(1)==irowcol(2)) then
                                      tt_ovrlp = 1.d0
                                  else
                                      tt_ovrlp = 0.d0
                                  end if
                                  tt_ham = ovrlp_mat%matrix_compr(isshift+iismall_ovrlp)
                              else
                                  tt_ovrlp = 0.d0
                                  tt_ham = 0.d0
                              end if
                              hamscal_compr(ii) = scale_factor*(tt_ham-shift_value*tt_ovrlp)
                              ii = ii + 1
                          end do
                      end do
                      !$omp end do
                      !$omp end parallel
                      calculate_SHS=.true.
                  else
                      calculate_SHS=.false.
                  end if
                  evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
                  evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)
    
    
                  !call uncompress_matrix(iproc,ovrlp_smat,ovrlp_mat%matrix_compr,overlap)
                  !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
                  !if (iproc==0) write(*,*) 'ovrlp_mat%matrix_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
                  !call uncompress_matrix(iproc,inv_ovrlp_smat,hamscal_compr,overlap)
                  !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
                  !if (iproc==0) write(*,*) 'hamscal_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
        
        
                  ! Determine the degree of the polynomial
                  npl=nint(degree_multiplicator* &
                       (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)))
                  npl=max(npl,NPL_MIN)
                  npl_boundaries = nint(degree_multiplicator* &
                      (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)) &
                          /foe_data_get_real(foe_obj,"fscale_lowerbound")) ! max polynomial degree for given eigenvalue boundaries
                  if (npl>npl_boundaries) then
                      npl=npl_boundaries
                      if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
                      if (iproc==0) write(*,*) 'STOP SINCE THIS WILL CREATE PROBLEMS WITH NPL_CHECK'
                      stop
                  end if
                  if (npl>nplx) stop 'npl>nplx'
        
                  ! Array the holds the Chebyshev polynomials. Needs to be recalculated
                  ! every time the Hamiltonian has been modified.
                  if (calculate_SHS) then
                      call f_free(chebyshev_polynomials)
                      chebyshev_polynomials = f_malloc((/nsize_polynomial,npl/),id='chebyshev_polynomials')
                  end if
                  if (iproc==0) then
                      call yaml_newline()
                      call yaml_mapping_open('ICE')
                      call yaml_map('eval bounds',&
                           (/foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                      call yaml_map('mult.',degree_multiplicator,fmt='(f5.2)')
                      call yaml_map('pol. deg.',npl)
                      call yaml_mapping_close()
                  end if
    
        
                  cc = f_malloc((/npl,3/),id='cc')
        
                  !!if (foe_data_get_real(foe_obj,"evlow")>=0.d0) then
                  !!    stop 'ERROR: lowest eigenvalue must be negative'
                  !!end if
                  if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      stop 'ERROR: highest eigenvalue must be positive'
                  end if
        
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
                  call timing(iproc, 'chebyshev_coef', 'ON')
        
                  call cheb_exp(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,1), ex)
                  call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), cc(1,1), cc(1,2), npl)
                  call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,3))
                  call evnoise(npl, cc(1,3), foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
    
                  call timing(iproc, 'chebyshev_coef', 'OF')
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
                
                
                
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
        
                  emergency_stop=.false.
                  if (calculate_SHS) then
                      ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                      !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
                      call chebyshev_clean(iproc, nproc, npl, cc, &
                           inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
                           inv_ovrlp_smat%smmm%isfvctr, foe_obj, &
                           inv_ovrlp_smat, hamscal_compr, &
                           inv_ovrlp%matrix_compr, .false., &
                           nsize_polynomial, SHS, inv_ovrlp_matrixp, penalty_ev, chebyshev_polynomials, &
                           emergency_stop)
                  else
                      ! The Chebyshev polynomials are already available
                      !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                           inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
                           inv_ovrlp_smat%smmm%isfvctr, &
                           inv_ovrlp_smat, chebyshev_polynomials, cc, inv_ovrlp_matrixp)
                  end if 
    
    
    
                 ! Check for an emergency stop, which happens if the kernel explodes, presumably due
                 ! to the eigenvalue bounds being too small.
                 ! mpi_lor seems not to work on certain systems...
                 if (emergency_stop) then
                     iflag=1
                 else
                     iflag=0
                 end if
    
                 if (nproc > 1) then
                     call mpiallred(iflag, 1, mpi_sum, bigdft_mpi%mpi_comm)
                 end if
    
                 if (iflag>0) then
                     emergency_stop=.true.
                 else
                     emergency_stop=.false.
                 end if
                 if (emergency_stop) then
                      eval_bounds_ok(1)=.false.
                      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
                      eval_bounds_ok(2)=.false.
                      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                      call f_free(cc)
                      cycle main_loop
                 end if
        
        
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
        
        
                  restart=.false.
        
                  ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                  ! (otherwise this has already been checked in the previous iteration).
                  if (calculate_SHS) then
                      !call check_eigenvalue_spectrum()
                      call check_eigenvalue_spectrum(nproc, inv_ovrlp_smat, ovrlp_smat, ovrlp_mat, 1, &
                           0, 1.2d0, 1.d0/1.2d0, penalty_ev, anoise, .false., &
                           foe_obj, restart, eval_bounds_ok)
                  end if
        
                  call f_free(cc)
        
                  if (restart) then
                      if(evbounds_shrinked) then
                          ! this shrink was not good, increase the saturation counter
                          call foe_data_set_int(foe_obj,"evboundsshrink_isatur",foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
                      end if
                      call foe_data_set_int(foe_obj,"evbounds_isatur",0)
                      cycle
                  end if
                      
                  ! eigenvalue bounds ok
                  if (calculate_SHS) then
                      call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
                  end if
                
    
                  exit
        
        
              end do main_loop
        
        
    
          call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, inv_ovrlp_matrixp, &
               inv_ovrlp%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr))
    

      end do spin_loop

  call f_free_ptr(inv_ovrlp_matrixp)
  call f_free(chebyshev_polynomials)
  call f_free(penalty_ev)
  call f_free(hamscal_compr)
  call f_free(SHS)

  call f_free_ptr(foe_obj%ef)
  call f_free_ptr(foe_obj%evlow)
  call f_free_ptr(foe_obj%evhigh)
  call f_free_ptr(foe_obj%bisection_shift)
  call f_free_ptr(foe_obj%charge)

  call timing(iproc, 'FOE_auxiliary ', 'OF')

  call f_release_routine()


      !!contains


      !!subroutine check_eigenvalue_spectrum()
      !!  implicit none
      !!  real(kind=8) :: bound_low, bound_up

      !!  ! The penalty function must be smaller than the noise.
      !!  bound_low=0.d0
      !!  bound_up=0.d0
      !!  if (inv_ovrlp_smat%smmm%nfvctrp>0) then
      !!      isegstart = inv_ovrlp_smat%istsegline(inv_ovrlp_smat%smmm%isfvctr+1)
      !!      isegend = inv_ovrlp_smat%istsegline(inv_ovrlp_smat%smmm%isfvctr+inv_ovrlp_smat%smmm%nfvctrp) + &
      !!                inv_ovrlp_smat%nsegline(inv_ovrlp_smat%smmm%isfvctr+inv_ovrlp_smat%smmm%nfvctrp)-1
      !!      !$omp parallel default(none) &
      !!      !$omp private(iseg, ii, jorb, irow, icol, iismall, tt) &
      !!      !$omp shared(isegstart, isegend, inv_ovrlp_smat, penalty_ev, bound_low, bound_up)
      !!      !$omp do reduction(+:bound_low,bound_up)
      !!      do iseg=isegstart,isegend
      !!          ii=inv_ovrlp_smat%keyv(iseg)-1
      !!          do jorb=inv_ovrlp_smat%keyg(1,iseg),inv_ovrlp_smat%keyg(2,iseg)
      !!              ii=ii+1
      !!              irow = (jorb-1)/inv_ovrlp_smat%nfvctr + 1
      !!              icol = jorb - (irow-1)*inv_ovrlp_smat%nfvctr
      !!              iismall = matrixindex_in_compressed(inv_ovrlp_smat, irow, icol)
      !!              if (iismall>0) then
      !!                  if (irow==jorb) then
      !!                      tt=1.d0
      !!                  else
      !!                      tt=0.d0
      !!                  end if
      !!              else
      !!                  tt=0.d0
      !!              end if
      !!              bound_low = bound_low + penalty_ev(icol,irow-inv_ovrlp_smat%smmm%isfvctr,2)*tt
      !!              bound_up = bound_up +penalty_ev(icol,irow-inv_ovrlp_smat%smmm%isfvctr,1)*tt
      !!          end do  
      !!      end do
      !!      !$omp end do
      !!      !$omp end parallel
      !!  end if
    
      !!  allredarr(1)=bound_low
      !!  allredarr(2)=bound_up

      !!  if (nproc > 1) then
      !!      call mpiallred(allredarr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
      !!  end if

      !!  allredarr=abs(allredarr) !for some crazy situations this may be negative
      !!  anoise=100.d0*anoise
      !!  !write(*,*) 'allredarr, anoise', allredarr, anoise
      !!  if (allredarr(1)>anoise) then
      !!      eval_bounds_ok(1)=.false.
      !!      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",1)/1.2d0,1)
      !!      restart=.true.
      !!  else
      !!      eval_bounds_ok(1)=.true.
      !!  end if
      !!  if (allredarr(2)>anoise) then
      !!      eval_bounds_ok(2)=.false.
      !!      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",1)*1.2d0,1)
      !!      restart=.true.
      !!  else
      !!      eval_bounds_ok(2)=.true.
      !!  end if

      !!end subroutine check_eigenvalue_spectrum

end subroutine ice



! Calculates chebychev expansion of x**ex, where ex is any value (typically -1, -1/2, 1/2)
! Taken from numerical receipes: press et al
subroutine cheb_exp(A,B,N,cc,ex)
  use module_base, pi => pi_param
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: A, B
  integer,intent(in) :: n, ex
  real(8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8) :: bma, bpa, y, arg, fac, tt, erfcc
  real(kind=8),dimension(50000) :: cf
  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  call f_routine(id='chebft')

  if (n>50000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  fac=2.d0/n
  !$omp parallel default(none) shared(bma,bpa,fac,n,cf,cc,ex) &
  !$omp private(k,y,arg,j,tt)
  !$omp do
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      !cf(k)=arg**ex
      select case(ex)
      case (-2)
          !ex=-0.5d0
          cf(k)=1.d0/sqrt(arg)
      case (2)
          !ex=0.5d0
          cf(k)=sqrt(arg)
      case (1)
          !ex=-1.d0
          cf(k)=1.d0/arg
      case default
          stop 'wrong power'
      end select
  end do
  !$omp end do
  !$omp do
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
  !$omp end do
  !$omp end parallel

  call f_release_routine()

end subroutine cheb_exp



subroutine check_eigenvalue_spectrum(nproc, smat_l, smat_s, mat, ispin, isshift, &
           factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
           foe_obj, restart, eval_bounds_ok)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use foe_base, only: foe_data, foe_data_set_real, foe_data_get_real
  implicit none

  ! Calling arguments
  type(sparse_matrix),intent(in) :: smat_l, smat_s
  type(matrices),intent(in) :: mat
  integer,intent(in) :: nproc, ispin, isshift
  real(kind=8),intent(in) :: factor_high, factor_low, anoise
  real(kind=8),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
  logical,intent(in) :: trace_with_overlap
  type(foe_data),intent(inout) :: foe_obj
  logical,intent(out) :: restart
  logical,dimension(2),intent(out) :: eval_bounds_ok

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, irow, icol, iismall
  real(kind=8) :: bound_low, bound_up, tt, noise
  real(kind=8),dimension(2) :: allredarr

  call f_routine(id='check_eigenvalue_spectrum')

  ! The penalty function must be smaller than the noise.
  bound_low=0.d0
  bound_up=0.d0
  if (smat_l%smmm%nfvctrp>0) then
      isegstart = smat_l%istsegline(smat_l%smmm%isfvctr+1)
      isegend = smat_l%istsegline(smat_l%smmm%isfvctr+smat_l%smmm%nfvctrp) + &
                smat_l%nsegline(smat_l%smmm%isfvctr+smat_l%smmm%nfvctrp)-1
      !$omp parallel default(none) &
      !$omp private(iseg, ii, jorb, irow, icol, iismall, tt) &
      !$omp shared(isegstart, isegend, smat_l, smat_s, mat, penalty_ev) &
      !$omp shared(bound_low, bound_up, isshift, trace_with_overlap) 
      !$omp do reduction(+:bound_low,bound_up)
      do iseg=isegstart,isegend
          ii=smat_l%keyv(iseg)-1
          do jorb=smat_l%keyg(1,iseg),smat_l%keyg(2,iseg)
              ii=ii+1
              irow = (jorb-1)/smat_l%nfvctr + 1
              icol = jorb - (irow-1)*smat_l%nfvctr
              iismall = matrixindex_in_compressed(smat_s, irow, icol)
              if (iismall>0) then
                  if (trace_with_overlap) then
                      ! Take the trace of the product matrix times overlap
                      tt=mat%matrix_compr(isshift+iismall)
                  else
                      ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
                      if (irow==jorb) then
                          tt=1.d0
                      else
                          tt=0.d0
                      end if
                  end if
              else
                  tt=0.d0
              end if
              bound_low = bound_low + penalty_ev(icol,irow-smat_l%smmm%isfvctr,2)*tt
              bound_up = bound_up +penalty_ev(icol,irow-smat_l%smmm%isfvctr,1)*tt
          end do  
      end do
      !$omp end do
      !$omp end parallel
  end if

  allredarr(1)=bound_low
  allredarr(2)=bound_up

  if (nproc > 1) then
      call mpiallred(allredarr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  allredarr=abs(allredarr) !for some crazy situations this may be negative
  noise=100.d0*anoise
  !write(*,*) 'allredarr, anoise', allredarr, anoise
  if (allredarr(1)>noise) then
      eval_bounds_ok(1)=.false.
      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
      restart=.true.
  else
      eval_bounds_ok(1)=.true.
  end if
  if (allredarr(2)>noise) then
      eval_bounds_ok(2)=.false.
      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*factor_high,ispin)
      restart=.true.
  else
      eval_bounds_ok(2)=.true.
  end if

  call f_release_routine()

end subroutine check_eigenvalue_spectrum
