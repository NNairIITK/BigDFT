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
           ebs, itout, it_scc, order_taylor, max_inversion_error, purification_quickreturn, &
           calculate_minusonehalf, foe_verbosity, &
           accuracy_level, tmb, foe_obj)
  use module_base
  use module_types
  use module_interfaces, except_this_one => foe
  use yaml_output
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, SPARSE_TASKGROUP, &
                               matrices
  use sparsematrix_init, only: matrixindex_in_compressed, get_line_and_column
  use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed, &
                          uncompress_matrix_distributed, orb_from_index, &
                          transform_sparsity_pattern, compress_matrix_distributed_new2
  use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                      foe_data_get_logical
  use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                         fermilevel_get_real, fermilevel_get_logical
  use chebyshev, only: chebyshev_clean, chebyshev_fast
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc,itout,it_scc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  real(kind=8),intent(in) :: tmprtr
  real(kind=8),intent(out) :: ebs
  logical,intent(in) :: purification_quickreturn
  logical,intent(in) :: calculate_minusonehalf
  integer,intent(in) :: foe_verbosity
  integer,intent(in) :: accuracy_level
  type(DFT_wavefunction),intent(inout) :: tmb
  type(foe_data),intent(inout) :: foe_obj

  ! Local variables
  integer :: npl, jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
  integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
  integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
  integer,parameter :: nplx=50000
  real(kind=8),dimension(:,:,:),allocatable :: cc, cc_check
  real(kind=8),dimension(:,:),allocatable :: chebyshev_polynomials, fermip_check
  real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value, sumn, sumn_check, charge_diff, ef_interpol, ddot
  real(kind=8) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
  real(kind=8) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
  logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible, emergency_stop
  real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr, fermi_check_compr
  real(kind=8),dimension(4,4) :: interpol_matrix
  real(kind=8),dimension(4) :: interpol_vector
  real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion
  logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
  real(kind=8) :: trace_sparse, temp_multiplicator, ebs_check, ef, ebsp
  integer :: irow, icol, itemp, iflag,info, ispin, isshift, imshift, ilshift, ilshift2, i, j, itg, ncount, istl, ists
  logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
  real(kind=8),parameter :: FSCALE_LOWER_LIMIT=5.d-3
  real(kind=8),parameter :: FSCALE_UPPER_LIMIT=5.d-2
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
  real(kind=8),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
  real(kind=8),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
  real(kind=8),parameter :: CHECK_RATIO=1.25d0
  integer,parameter :: NPL_MIN=100
  !!type(matrices) :: inv_ovrlp
  integer,parameter :: NTEMP_ACCURATE=4
  integer,parameter :: NTEMP_FAST=1
  real(kind=8) :: degree_multiplicator
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2
  integer,parameter :: imode=SPARSE
  type(fermi_aux) :: f
  real(kind=8),dimension(2) :: temparr
  real(kind=8),dimension(:,:),allocatable :: penalty_ev_new
  real(kind=8),dimension(:),allocatable :: fermi_new, fermi_check_new, fermi_small_new
  integer :: iline, icolumn, icalc
  


  call f_routine(id='foe')


  if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')

  if (accuracy_level/=FOE_ACCURATE .and. accuracy_level/=FOE_FAST) then
      stop 'wrong value of accuracy_level'
  end if


  call timing(iproc, 'FOE_auxiliary ', 'ON')


  evbounds_shrinked=.false.


  !!penalty_ev = f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp,2/),id='penalty_ev')
  !!fermip_check = f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp/),id='fermip_check')
  fermi_check_compr = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_TASKGROUP, id='fermi_check_compr')

  penalty_ev_new = f_malloc((/tmb%linmat%l%smmm%nvctrp,2/),id='penalty_ev_new')
  fermi_check_new = f_malloc((/tmb%linmat%l%smmm%nvctrp_mm/),id='fermip_check_new')
  fermi_new = f_malloc((/tmb%linmat%l%smmm%nvctrp/),id='fermi_new')
  fermi_small_new = f_malloc((/tmb%linmat%l%smmm%nvctrp_mm/),id='fermi_small_new')


  call timing(iproc, 'FOE_auxiliary ', 'OF')
  if (calculate_minusonehalf) then
      if (iproc==0) call yaml_map('S^-1/2','recalculate')
      call overlap_minus_onehalf() ! has internal timer
  else
      if (iproc==0) call yaml_map('S^-1/2','from memory')
  end if
  call timing(iproc, 'FOE_auxiliary ', 'ON')


  hamscal_compr = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

    
  ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
  nsize_polynomial = tmb%linmat%l%smmm%nvctrp_mm
  
  
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

      isshift=(ispin-1)*tmb%linmat%s%nvctrp_tg
      imshift=(ispin-1)*tmb%linmat%m%nvctrp_tg
      ilshift=(ispin-1)*tmb%linmat%l%nvctrp_tg
      ilshift2=(ispin-1)*tmb%linmat%l%nvctrp_tg

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
        
              !if (tmb%linmat%l%smmm%nfvctrp>0) then
              !    call f_zero(tmb%linmat%l%nfvctr*tmb%linmat%l%smmm%nfvctrp*tmb%linmat%l%nspin,tmb%linmat%kernel_%matrixp(1,1,1))
              !end if
        
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
              
        
                  ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
                  if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                      !!call scale_and_shift_hamiltonian()
                      call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, tmb%linmat%l, &
                           tmb%linmat%m, tmb%linmat%ham_, imshift, &
                           smat2=tmb%linmat%s, mat2=tmb%linmat%ovrlp_, i2shift=isshift, &
                           matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
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
        
        
                  cc = f_malloc((/npl,3,1/),id='cc')
                  cc_check = f_malloc((/npl,3,1/),id='cc_check')
        
                  if (foe_data_get_real(foe_obj,"evlow",ispin)>=0.d0) then
                      stop 'ERROR: lowest eigenvalue must be negative'
                  end if
                  if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      stop 'ERROR: highest eigenvalue must be positive'
                  end if
        
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
                  call timing(iproc, 'chebyshev_coef', 'ON')
        
                  call chebft(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,1,1), &
                       foe_data_get_real(foe_obj,"ef",ispin), fscale, tmprtr)
                  call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), cc(1,1,1), cc(1,2,1), npl)
                  call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,3,1))
                  call evnoise(npl, cc(1,3,1), foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
    
                  call chebft(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,1,1), &
                       foe_data_get_real(foe_obj,"ef",ispin), fscale_check, tmprtr)
                  call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), &
                       cc_check(1,1,1), cc_check(1,2,1), npl_check)
                  call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,3,1))
        
                  call timing(iproc, 'chebyshev_coef', 'OF')
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
                
                  !!if (iproc==0) then
                  !!    call pltwght(npl,cc(1,1),cc(1,2),foe_data_get_real(foe_obj,"evlow"),foe_data_get_real(foe_obj,"evhigh"),foe_data_get_real(foe_obj,"ef"),foe_data_get_real(foe_obj,"fscale"),temperature)
                  !!    call pltexp(anoise,npl,cc(1,3),foe_data_get_real(foe_obj,"evlow"),foe_data_get_real(foe_obj,"evhigh"))
                  !!end if
                
                
                  !!write(1000+iproc,*) 'foe_data_get_real(foe_obj,"evlow",ispin)',foe_data_get_real(foe_obj,"evlow",ispin)
                  !!write(1000+iproc,*) 'foe_data_get_real(foe_obj,"evhigh",ispin)', foe_data_get_real(foe_obj,"evhigh",ispin)
                  !!write(1000+iproc,*) 'npl_check', npl_check
                  !!write(1000+iproc,*) 'foe_data_get_real(foe_obj,"ef",ispin)',foe_data_get_real(foe_obj,"ef",ispin)
                  !!write(1000+iproc,*) 'fscale_check',fscale_check
                  !!write(1000+iproc,*) 'tmprtr',tmprtr

                  if (tmb%linmat%l%nspin==1) then
                      do ipl=1,npl
                          cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                          cc(ipl,2,1)=2.d0*cc(ipl,2,1)
                          cc(ipl,3,1)=2.d0*cc(ipl,3,1)
                          cc_check(ipl,1,1)=2.d0*cc_check(ipl,1,1)
                          cc_check(ipl,2,1)=2.d0*cc_check(ipl,2,1)
                          cc_check(ipl,3,1)=2.d0*cc_check(ipl,3,1)
                      end do
                  end if
                
                
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
        
                  emergency_stop=.false.
                  if (calculate_SHS) then
                      ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                      if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
                      call chebyshev_clean(iproc, nproc, npl, cc, &
                           tmb%linmat%l, hamscal_compr, &
                           tmb%linmat%ovrlppowers_(2)%matrix_compr(ilshift2+1:), calculate_SHS, &
                           nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                           emergency_stop)
                      call transform_sparsity_pattern(tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nvctrp_mm, tmb%linmat%l%smmm%isvctr_mm, &
                           tmb%linmat%l%nseg, tmb%linmat%l%keyv, tmb%linmat%l%keyg, &
                           tmb%linmat%l%smmm%nvctrp, tmb%linmat%l%smmm%isvctr, &
                           tmb%linmat%l%smmm%nseg, tmb%linmat%l%smmm%keyv, tmb%linmat%l%smmm%keyg, &
                           fermi_new, fermi_small_new)


                      !!do i=1,tmb%linmat%l%smmm%nvctrp
                      !!    ii = tmb%linmat%l%smmm%isvctr + i
                      !!    call get_line_and_column(ii, tmb%linmat%l%smmm%nseg, tmb%linmat%l%smmm%keyv, tmb%linmat%l%smmm%keyg, iline, icolumn)
                      !!!!    tmb%linmat%kernel_%matrixp(icolumn,iline-tmb%linmat%l%smmm%isfvctr,1) = fermi_new(i)
                      !!    penalty_ev(icolumn,iline-tmb%linmat%l%smmm%isfvctr,1) = penalty_ev_new(i,1)
                      !!    penalty_ev(icolumn,iline-tmb%linmat%l%smmm%isfvctr,2) = penalty_ev_new(i,2)
                      !!end do

                  else
                      ! The Chebyshev polynomials are already available
                      if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                           tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nfvctrp, &
                          tmb%linmat%l, chebyshev_polynomials, 1, cc, fermi_small_new)
                      !!call calculate_trace_distributed_new(fermi_new, sumn)
                      !!write(*,*) 'trace debug', sumn

                      !!call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                      !!     tmb%linmat%l, fermi_small_new, tmb%linmat%kernel_%matrixp(:,:,1))
                      !!!!call calculate_trace_distributed(tmb%linmat%kernel_%matrixp, sumn)
                      !!write(*,'(a,2es16.8)') 'sum(fermi_new), sum(tmb%linmat%kernel_%matrix(:,:,1)', sum(abs(fermi_new)), sum(abs(tmb%linmat%kernel_%matrixp(:,:,1)))
                  end if 
    
    
    
                 !call check_emergency_stop()
                 !!call check_emergency_stop(nproc,emergency_stop)
                 !!if (emergency_stop) then
                 !!     eval_bounds_ok(1)=.false.
                 !!     call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*1.2d0,ispin)
                 !!     eval_bounds_ok(2)=.false.
                 !!     call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                 !!     if (iproc==0) then
                 !!         if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                 !!              (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                 !!         call yaml_mapping_close()
                 !!         !call bigdft_utils_flush(unit=6)
                 !!     end if
                 !!     call f_free(cc)
                 !!     call f_free(cc_check)
                 !!     cycle main_loop
                 !!end if
        
        
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
        
        
                  restart=.false.
        
                  ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                  ! (otherwise this has already been checked in the previous iteration).
                  if (calculate_SHS) then
                      !!call check_eigenvalue_spectrum()
                      !!call check_eigenvalue_spectrum(nproc, tmb%linmat%l, tmb%linmat%s, tmb%linmat%ovrlp_, ispin, &
                      !!      isshift, 1.2d0, 1.2d0, penalty_ev, anoise, .true., emergency_stop, &
                      !!      foe_obj, restart, eval_bounds_ok)
                      call check_eigenvalue_spectrum_new(nproc, tmb%linmat%l, tmb%linmat%s, tmb%linmat%ovrlp_, ispin, &
                            isshift, 1.2d0, 1.2d0, penalty_ev_new, anoise, .true., emergency_stop, &
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
                
                  !call calculate_trace_distributed(tmb%linmat%kernel_%matrixp, sumn)
                  call calculate_trace_distributed_new(fermi_small_new, sumn)
                  write(*,*) 'sumn',sumn
        
    
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
                           tmb%linmat%l%nfvctr, tmb%linmat%l%smmm%nfvctrp, &
                           tmb%linmat%l, chebyshev_polynomials, 1, cc_check, fermi_check_new)
                      !!call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                      !!     tmb%linmat%l, fermi_check_new, fermip_check)
                      call f_free(cc_check)
                      diff=0.d0
                      !do iorb=1,tmb%linmat%l%smmm%nfvctrp
                      !    do jorb=1,tmb%linmat%l%nfvctr
                      !        !SM: need to fix the spin here
                      !        diff = diff + (tmb%linmat%kernel_%matrixp(jorb,iorb,1)-fermip_check(jorb,iorb))**2
                      !    end do
                      !end do
                      do i=1,tmb%linmat%l%smmm%nvctrp_mm
                          diff = diff + (fermi_small_new(i)-fermi_check_new(i))**2
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
        
        
        
    
         !!call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
         !!     tmb%linmat%kernel_%matrixp(:,1:tmb%linmat%l%smmm%nfvctrp,1), &
         !!     tmb%linmat%kernel_%matrix_compr(ilshift+1:))
         call compress_matrix_distributed_new2(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
              fermi_small_new, &
              tmb%linmat%kernel_%matrix_compr(ilshift+1:))
         !!write(*,*) 'sum(tmb%linmat%kernel_%matrix_compr(ilshift+1:))',sum(tmb%linmat%kernel_%matrix_compr(ilshift+1:))
         !!tmb%linmat%kernel_%matrix_compr(ilshift+1:) = fermi_small_new
    
         !!call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
         !!     fermip_check, fermi_check_compr(1))
         call compress_matrix_distributed_new2(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
              fermi_check_new, fermi_check_compr(1))
         !!fermi_check_compr = fermi_check_new
    
    
        
        
          ! Calculate S^-1/2 * K * S^-1/2^T
          ! Since S^-1/2 is symmetric, don't use the transpose
          call retransform(tmb%linmat%kernel_%matrix_compr(ilshift+1:))

          !!do i=ilshift+1,ilshift+tmb%linmat%l%nvctr
          !!    write(3000+iproc,'(a,2i8,es16.6)') 'ispin, i, val', ispin, i, tmb%linmat%kernel_%matrix_compr(i)
          !!end do
    
          call retransform(fermi_check_compr)
    
          !call calculate_trace_distributed(fermip_check, sumn_check)
          call calculate_trace_distributed_new(fermi_check_new, sumn_check)

          !@NEW ##########################
          sumn = trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                 tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
                 tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)
          sumn_check = trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                       tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
                       fermi_check_compr, ispin)
          !@ENDNEW #######################
        
    
          ! Calculate trace(KH). Since they have the same sparsity pattern and K is
          ! symmetric, this is a simple ddot.
          ncount = tmb%linmat%l%smmm%istartend_mm_dj(2) - tmb%linmat%l%smmm%istartend_mm_dj(1) + 1
          istl = tmb%linmat%l%smmm%istartend_mm_dj(1)-tmb%linmat%l%isvctrp_tg
          ebsp = ddot(ncount, tmb%linmat%kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)

          ncount = tmb%linmat%l%smmm%istartend_mm_dj(2) - tmb%linmat%l%smmm%istartend_mm_dj(1) + 1
          istl = tmb%linmat%l%smmm%istartend_mm_dj(1)
          ebs_check = ddot(ncount, fermi_check_compr(istl-tmb%linmat%l%isvctrp_tg), 1, &
                      hamscal_compr(istl-tmb%linmat%l%isvctrp_tg), 1)

          temparr(1) = ebsp
          temparr(2) = ebs_check
          if (nproc>1) then
              call mpiallred(temparr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
          end if
          ebsp = temparr(1)
          ebs_check = temparr(2)


          ebsp=ebsp/scale_factor+shift_value*sumn
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
                 tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
                 tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)


          ! Recalculate trace(KH) (needed since the kernel was modified in the above purification). 
          ! If no purification is done, this should not be necessary.
          ! Since K and H have the same sparsity pattern and K is
          ! symmetric, the trace is a simple ddot.
          ncount = tmb%linmat%l%smmm%istartend_mm_dj(2) - tmb%linmat%l%smmm%istartend_mm_dj(1) + 1
          istl = tmb%linmat%l%smmm%istartend_mm_dj(1) - tmb%linmat%l%isvctrp_tg
          ebsp = ddot(ncount, tmb%linmat%kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)
          if (nproc>1) then
              call mpiallred(ebsp, 1, mpi_sum, bigdft_mpi%mpi_comm)
          end if
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
  end if

  degree_sufficient=.true.

  if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


  call f_free(chebyshev_polynomials)
  !!call f_free(penalty_ev)
  call f_free(hamscal_compr)
  !!call f_free(fermip_check)
  call f_free(fermi_check_compr)

  call f_free(penalty_ev_new)
  call f_free(fermi_check_new)
  call f_free(fermi_new)
  call f_free(fermi_small_new)

  call timing(iproc, 'FOE_auxiliary ', 'OF')

  call f_release_routine()



      contains

        subroutine overlap_minus_onehalf()
          use sparsematrix_base, only: sparsematrix_malloc, SPARSE_FULL
          use sparsematrix, only: extract_taskgroup_inplace
          implicit none
          real(kind=8) :: max_error, mean_error
          integer :: i, j, ii
          real(kind=8),dimension(:),allocatable :: tmparr

          call f_routine(id='overlap_minus_onehalf')

          ! Taylor approximation of S^-1/2 up to higher order
          if (imode==DENSE) then
              stop 'overlap_minus_onehalf: DENSE is deprecated'
              !!tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, &
              !!                           id='tmb%linmat%ovrlp_%matrix')
              !!call uncompress_matrix(iproc, tmb%linmat%s, &
              !!     inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)

              !!inv_ovrlp%matrix=sparsematrix_malloc_ptr(tmb%linmat%l, &
              !!                                  iaction=DENSE_FULL, id='inv_ovrlp%matrix')
              !!call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, &
              !!     imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
              !!     ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, &
              !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
              !!call compress_matrix(iproc, tmb%linmat%l, inmat=inv_ovrlp%matrix, outmat=inv_ovrlp%matrix_compr)
          end if
          if (imode==SPARSE) then
              !tmparr = sparsematrix_malloc(tmb%linmat%s,iaction=SPARSE_FULL,id='tmparr')
              !call vcopy(tmb%linmat%s%nvctr*tmb%linmat%s%nspin, tmb%linmat%ovrlp_%matrix_compr(1), 1, tmparr(1), 1)
              !call extract_taskgroup_inplace(tmb%linmat%s, tmb%linmat%ovrlp_)
              call overlapPowerGeneral(iproc, nproc, order_taylor, 1, (/-2/), -1, &
                   imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
                   ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(2), &
                   check_accur=.true., max_error=max_error, mean_error=mean_error)
              !call vcopy(tmb%linmat%s%nvctr*tmb%linmat%s%nspin, tmparr(1), 1, tmb%linmat%ovrlp_%matrix_compr(1), 1)
              !call f_free(tmparr)
          end if
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)

          call f_release_routine()
      end subroutine overlap_minus_onehalf



      subroutine retransform(matrix_compr)
          use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, sparsemm, &
                                  uncompress_matrix_distributed, compress_matrix_distributed, uncompress_matrix_distributed2
          ! Calling arguments
          real(kind=8),dimension(tmb%linmat%l%nvctrp_tg),intent(inout) :: matrix_compr

          ! Local variables
          real(kind=8),dimension(:,:),pointer :: inv_ovrlpp, tempp
          integer,dimension(:,:),pointer :: onedimindices
          real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
          integer,dimension(:,:,:),allocatable :: istindexarr
          integer :: nout, nseq

          call f_routine(id='retransform')

          inv_ovrlpp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
          tempp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
          inv_ovrlp_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
          kernel_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
          call sequential_acces_matrix_fast2(tmb%linmat%l, matrix_compr, kernel_compr_seq)
          call sequential_acces_matrix_fast2(tmb%linmat%l, &
               tmb%linmat%ovrlppowers_(2)%matrix_compr(ilshift2+1:), inv_ovrlp_compr_seq)
          call uncompress_matrix_distributed2(iproc, tmb%linmat%l, DENSE_MATMUL, &
               tmb%linmat%ovrlppowers_(2)%matrix_compr(ilshift2+1:), inv_ovrlpp)

           tempp=0.d0
          call sparsemm(tmb%linmat%l, kernel_compr_seq, inv_ovrlpp, tempp)
          inv_ovrlpp=0.d0
          call sparsemm(tmb%linmat%l, inv_ovrlp_compr_seq, tempp, inv_ovrlpp)

          call f_zero(matrix_compr)
          call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
               inv_ovrlpp, matrix_compr)

          call f_free_ptr(inv_ovrlpp)
          call f_free_ptr(tempp)
          call f_free(inv_ovrlp_compr_seq)
          call f_free(kernel_compr_seq)

          call f_release_routine()

      end subroutine retransform




      !!subroutine calculate_trace_distributed(matrixp, trace)
      !!    real(kind=8),dimension(tmb%linmat%l%nfvctr,tmb%linmat%l%smmm%nfvctrp),intent(in) :: matrixp
      !!    real(kind=8),intent(out) :: trace

      !!    call f_routine(id='calculate_trace_distributed')


      !!    trace=0.d0
      !!    if (tmb%linmat%l%smmm%nfvctrp>0) then
      !!        !$omp parallel default(private) shared(matrixp, tmb, trace) 
      !!        !$omp do reduction(+:trace)
      !!        do iseg=tmb%linmat%l%smmm%isseg,tmb%linmat%l%smmm%ieseg
      !!            ii=tmb%linmat%l%keyv(iseg)-1
      !!            ! A segment is always on one line, therefore no double loop
      !!            do jorb=tmb%linmat%l%keyg(1,1,iseg),tmb%linmat%l%keyg(2,1,iseg)
      !!                ii=ii+1
      !!                if (ii<tmb%linmat%l%smmm%isvctr_mm+1) cycle
      !!                if (ii>tmb%linmat%l%smmm%isvctr_mm+tmb%linmat%l%smmm%nvctrp_mm) exit
      !!                iiorb = tmb%linmat%l%keyg(1,2,iseg)
      !!                jjorb = jorb
      !!                if (jjorb==iiorb) then
      !!                    !write(900,*) iiorb, matrixp(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr)
      !!                    trace = trace + matrixp(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr)
      !!                end if
      !!            end do  
      !!        end do
      !!        !$omp end do
      !!        !$omp end parallel
      !!    end if

      !!    if (nproc > 1) then
      !!        call mpiallred(trace, 1, mpi_sum, bigdft_mpi%mpi_comm)
      !!    end if

      !!    call f_release_routine()
      !!end subroutine calculate_trace_distributed


      subroutine calculate_trace_distributed_new(matrixp, trace)
          real(kind=8),dimension(tmb%linmat%l%smmm%nvctrp_mm),intent(in) :: matrixp
          real(kind=8),intent(out) :: trace
          integer :: i, ii

          call f_routine(id='calculate_trace_distributed_new')

          !!trace=0.d0
          !!if (tmb%linmat%l%smmm%nfvctrp>0) then
          !!    !$omp parallel default(private) shared(matrixp, tmb, trace) 
          !!    !$omp do reduction(+:trace)
          !!    do iseg=tmb%linmat%l%smmm%isseg,tmb%linmat%l%smmm%ieseg
          !!        ii=tmb%linmat%l%keyv(iseg)-1
          !!        ! A segment is always on one line, therefore no double loop
          !!        do jorb=tmb%linmat%l%keyg(1,1,iseg),tmb%linmat%l%keyg(2,1,iseg)
          !!            ii=ii+1
          !!            if (ii<tmb%linmat%l%smmm%isvctr_mm+1) cycle
          !!            if (ii>tmb%linmat%l%smmm%isvctr_mm+tmb%linmat%l%smmm%nvctrp_mm) exit
          !!            iiorb = tmb%linmat%l%keyg(1,2,iseg)
          !!            jjorb = jorb
          !!            if (jjorb==iiorb) trace = trace + matrixp(jjorb,iiorb-tmb%linmat%l%smmm%isfvctr)
          !!        end do  
          !!    end do
          !!    !$omp end do
          !!    !$omp end parallel
          !!end if

          trace = 0.d0
          do i=1,tmb%linmat%l%smmm%nvctrp_mm
              ii = tmb%linmat%l%smmm%isvctr_mm + i
              call get_line_and_column(ii, tmb%linmat%l%nseg, tmb%linmat%l%keyv, tmb%linmat%l%keyg, iline, icolumn)
              if (iline==icolumn) then
                  !write(901,*) iiorb, matrixp(i)
                  trace = trace + matrixp(i)
              end if
          end do

          if (nproc > 1) then
              call mpiallred(trace, 1, mpi_sum, bigdft_mpi%mpi_comm)
          end if

          call f_release_routine()
      end subroutine calculate_trace_distributed_new


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
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iel

  call f_routine(id='compress_polynomial_vector')

  if (norbp>0) then
      ii=0
      !!$omp parallel default(private) shared(fermi, vector, vector_compressed)
      !!$omp do
      !do iseg=isegstart,isegend
      do iseg=fermi%smmm%isseg,fermi%smmm%ieseg
          iel = fermi%keyv(iseg) - 1
          ! A segment is always on one line, therefore no double loop
          do jorb=fermi%keyg(1,1,iseg),fermi%keyg(2,1,iseg)
              iel = iel + 1
              if (iel<fermi%smmm%isvctr_mm+1) cycle
              if (iel>fermi%smmm%isvctr_mm+fermi%smmm%nvctrp_mm) exit
              ii=ii+1
              iiorb = fermi%keyg(1,2,iseg)
              jjorb = jorb
              vector_compressed(ii)=vector(jjorb,iiorb-isorb)
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if

  if (ii/=fermi%smmm%nvctrp_mm) then
      write(*,*) 'ii, fermi%nvctrp, size(vector_compressed)', ii, fermi%smmm%nvctrp_mm, size(vector_compressed)
      stop 'compress_polynomial_vector: ii/=fermi%nvctrp'
  end if

  call f_release_routine()

end subroutine compress_polynomial_vector




subroutine compress_polynomial_vector_new(iproc, nproc, nsize_polynomial, norb, norbp, isorb, &
           fermi, vector_compr, vector_compressed)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: get_line_and_column
  use sparsematrix, only: transform_sparsity_pattern
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp, isorb
  type(sparse_matrix),intent(in) :: fermi
  real(kind=8),dimension(fermi%smmm%nvctrp),intent(in) :: vector_compr
  real(kind=8),dimension(nsize_polynomial),intent(out) :: vector_compressed

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iel, i, iline, icolumn
  real(kind=8),dimension(:,:),allocatable :: vector

  call f_routine(id='compress_polynomial_vector')

  call transform_sparsity_pattern(fermi%nfvctr, fermi%smmm%nvctrp_mm, fermi%smmm%isvctr_mm, &
       fermi%nseg, fermi%keyv, fermi%keyg, &
       fermi%smmm%nvctrp, fermi%smmm%isvctr, fermi%smmm%nseg, fermi%smmm%keyv, fermi%smmm%keyg, &
       vector_compr, vector_compressed)


  !!vector = f_malloc((/norb,norbp/),id='vector')

  !!   do i=1,fermi%smmm%nvctrp
  !!       ii = fermi%smmm%isvctr + i
  !!       call get_line_and_column(ii, fermi%smmm%nseg, fermi%smmm%keyv, fermi%smmm%keyg, iline, icolumn)
  !!       vector(icolumn,iline-fermi%smmm%isfvctr) = vector_compr(i)
  !!   end do


  !!if (norbp>0) then
  !!    ii=0
  !!    !!$omp parallel default(private) shared(fermi, vector, vector_compressed)
  !!    !!$omp do
  !!    !do iseg=isegstart,isegend
  !!    do iseg=fermi%smmm%isseg,fermi%smmm%ieseg
  !!        iel = fermi%keyv(iseg) - 1
  !!        ! A segment is always on one line, therefore no double loop
  !!        do jorb=fermi%keyg(1,1,iseg),fermi%keyg(2,1,iseg)
  !!            iel = iel + 1
  !!            if (iel<fermi%smmm%isvctr_mm+1) cycle
  !!            if (iel>fermi%smmm%isvctr_mm+fermi%smmm%nvctrp_mm) exit
  !!            ii=ii+1
  !!            iiorb = fermi%keyg(1,2,iseg)
  !!            jjorb = jorb
  !!            vector_compressed(ii)=vector(jjorb,iiorb-isorb)
  !!        end do
  !!    end do
  !!    !!$omp end do
  !!    !!$omp end parallel
  !!end if

  !!if (ii/=fermi%smmm%nvctrp_mm) then
  !!    write(*,*) 'ii, fermi%nvctrp, size(vector_compressed)', ii, fermi%smmm%nvctrp_mm, size(vector_compressed)
  !!    stop 'compress_polynomial_vector: ii/=fermi%nvctrp'
  !!end if

  !!call f_free(vector)

  call f_release_routine()

end subroutine compress_polynomial_vector_new



subroutine uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
           fermi, vector_compressed, vector)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsize_polynomial
  type(sparse_matrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial),intent(in) :: vector_compressed
  real(kind=8),dimension(fermi%nfvctr,fermi%smmm%nfvctrp),intent(out) :: vector

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iel


  if (fermi%smmm%nfvctrp>0) then
      call f_zero(vector)
      !!$omp parallel do default(private) &
      !!$omp shared(isegstart, isegend, fermi, vector, vector_compressed)
      ii = 0
      do iseg=fermi%smmm%isseg,fermi%smmm%ieseg
          iel = fermi%keyv(iseg) - 1
          !ii=fermi%keyv(iseg)-fermi%keyv(isegstart)
          ! A segment is always on one line, therefore no double loop
          do jorb=fermi%keyg(1,1,iseg),fermi%keyg(2,1,iseg)
              iel = iel + 1
              if (iel<fermi%smmm%isvctr_mm+1) cycle
              if (iel>fermi%smmm%isvctr_mm+fermi%smmm%nvctrp_mm) exit
              ii=ii+1
              iiorb = fermi%keyg(1,2,iseg)
              jjorb = jorb
              vector(jjorb,iiorb-fermi%smmm%isfvctr)=vector_compressed(ii)
              !write(*,*) 'ii, iiorb-fermi%isfvctr, jjorb', ii, iiorb-fermi%isfvctr, jjorb
          end do
      end do
      !!$omp end parallel do
  end if

  if (ii/=nsize_polynomial) stop 'ERROR uncompress_polynomial_vector: ii/=nsize_polynomial'

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
  real(kind=8),dimension(asmat%nvctrp_tg),intent(in) :: amat
  real(kind=8),dimension(bsmat%nvctrp_tg),intent(in) :: bmat

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iilarge
  integer :: ierr, iashift, ibshift, iel
  real(kind=8) :: sumn, trace_sparse


  call f_routine(id='trace_sparse')

  iashift = 0!(ispin-1)*asmat%nvctr
  ibshift = 0!(ispin-1)*bsmat%nvctr

  sumn=0.d0
  if (asmat%smmm%nfvctrp>0) then
      !$omp parallel default(none) &
      !$omp private(iseg, ii, jorb, iiorb, jjorb, iilarge, iel) &
      !$omp shared(bsmat, asmat, amat, bmat, iashift, ibshift, sumn)
      !$omp do reduction(+:sumn)
      !do iseg=isegstart,isegend
      do iseg=asmat%smmm%isseg,asmat%smmm%ieseg
          iel = asmat%keyv(iseg) - 1
          ii=iashift+asmat%keyv(iseg)-1
          ! A segment is always on one line, therefore no double loop
          do jorb=asmat%keyg(1,1,iseg),asmat%keyg(2,1,iseg)
              iel = iel + 1
              if (iel<asmat%smmm%isvctr_mm+1) cycle
              if (iel>asmat%smmm%isvctr_mm+asmat%smmm%nvctrp_mm) exit
              ii=ii+1
              iiorb = asmat%keyg(1,2,iseg)
              jjorb = jorb
              iilarge = ibshift + matrixindex_in_compressed(bsmat, iiorb, jjorb)
              sumn = sumn + amat(ii-asmat%isvctrp_tg)*bmat(iilarge-bsmat%isvctrp_tg)
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
subroutine ice(iproc, nproc, norder_polynomial, ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp)
  use module_base
  use module_types
  use module_interfaces, except_this_one_A => ice
  use yaml_output
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, &
                               sparsematrix_malloc0_ptr, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, SPARSE_TASKGROUP, &
                               matrices
  use sparsematrix_init, only: matrixindex_in_compressed, get_line_and_column
  use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed, orb_from_index, &
                          compress_matrix_distributed_new2, transform_sparsity_pattern
  use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                      foe_data_set_logical, foe_data_get_logical
  use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                         fermilevel_get_real, fermilevel_get_logical
  use chebyshev, only: chebyshev_clean, chebyshev_fast
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, norder_polynomial, ncalc
  type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
  integer,dimension(ncalc) :: ex
  type(matrices),intent(in) :: ovrlp_mat
  type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp

  ! Local variables
  integer :: npl, jorb, it, ii, iseg
  integer :: isegstart, isegend, iismall, nsize_polynomial
  integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i
  integer,parameter :: nplx=50000
  real(kind=8),dimension(:,:),allocatable :: chebyshev_polynomials
  real(kind=8),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
  real(kind=8),dimension(:,:,:),allocatable :: cc, penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value
  real(kind=8) :: evlow_old, evhigh_old, tt
  real(kind=8) :: tt_ovrlp, tt_ham
  logical :: restart, calculate_SHS, emergency_stop
  real(kind=8),dimension(2) :: allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr
  logical,dimension(2) :: eval_bounds_ok
  integer,dimension(2) :: irowcol
  integer :: irow, icol, iflag, ispin, isshift, ilshift, ilshift2
  logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
  integer,parameter :: NPL_MIN=5
  real(kind=8),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
  real(kind=8) :: degree_multiplicator
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2
  integer,parameter :: imode=SPARSE
  type(foe_data) :: foe_obj
  real(kind=8),dimension(:),allocatable :: eval, work
  real(kind=8),dimension(:,:),allocatable :: tempmat
  integer :: lwork, info, j, icalc, iline, icolumn
  real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_new
  real(kind=8),dimension(:,:),allocatable :: penalty_ev_new
  real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_small_new

  !!real(kind=8),dimension(ovrlp_smat%nfvctr,ovrlp_smat%nfvctr) :: overlap
  !!real(kind=8),dimension(ovrlp_smat%nfvctr) :: eval
  !!integer,parameter :: lwork=100000
  !!real(kind=8),dimension(lwork) :: work
  !!integer :: info

  call f_routine(id='ice')


  penalty_ev_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,2/),id='penalty_ev_new')
  inv_ovrlp_matrixp_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,ncalc/),id='inv_ovrlp_matrixp_new')
  inv_ovrlp_matrixp_small_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp_mm,ncalc/),id='inv_ovrlp_matrixp_small_new')


!@ JUST FOR THE MOMENT.... ########################
     foe_obj%ef = f_malloc0_ptr(ovrlp_smat%nspin,id='(foe_obj%ef)')
     foe_obj%evlow = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evlow')
     foe_obj%evhigh = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evhigh')
     foe_obj%bisection_shift = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%bisection_shift')
     foe_obj%charge = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%charge')
     do ispin=1,ovrlp_smat%nspin
         call foe_data_set_real(foe_obj,"ef",0.d0,ispin)
         call foe_data_set_real(foe_obj,"evlow",0.5d0,ispin)
         call foe_data_set_real(foe_obj,"evhigh",1.5d0,ispin)
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


  evbounds_shrinked = .false.

  !!!@ TEMPORARY: eigenvalues of  the overlap matrix ###################
  !!tempmat = f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/),id='tempmat')
  !!do iseg=1,ovrlp_smat%nseg
  !!    ii=ovrlp_smat%keyv(iseg)
  !!    do i=ovrlp_smat%keyg(1,1,iseg),ovrlp_smat%keyg(2,1,iseg)
  !!        tempmat(i,ovrlp_smat%keyg(1,2,iseg)) = ovrlp_mat%matrix_compr(ii)
  !!        ii = ii + 1
  !!    end do
  !!end do
  !!!!if (iproc==0) then
  !!!!    do i=1,ovrlp_smat%nfvctr
  !!!!        do j=1,ovrlp_smat%nfvctr
  !!!!            write(*,'(a,2i6,es17.8)') 'i,j,val',i,j,tempmat(j,i)
  !!!!        end do
  !!!!    end do
  !!!!end if
  !!eval = f_malloc(ovrlp_smat%nfvctr,id='eval')
  !!lwork=100*ovrlp_smat%nfvctr
  !!work = f_malloc(lwork,id='work')
  !!call dsyev('n','l', ovrlp_smat%nfvctr, tempmat, ovrlp_smat%nfvctr, eval, work, lwork, info)
  !!!if (iproc==0) write(*,*) 'eval',eval
  !!if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(ovrlp_smat%nfvctr)/),fmt='(es16.6)')

  !!call f_free(tempmat)
  !!call f_free(eval)
  !!call f_free(work)

  !!!@ END TEMPORARY: eigenvalues of  the overlap matrix ###############


  call timing(iproc, 'FOE_auxiliary ', 'ON')



  !!penalty_ev = f_malloc((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,2/),id='penalty_ev')


  hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

    
  ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
  nsize_polynomial = inv_ovrlp_smat%smmm%nvctrp_mm
  
  
  ! Fake allocation, will be modified later
  chebyshev_polynomials = f_malloc((/nsize_polynomial,1/),id='chebyshev_polynomials')


  !inv_ovrlp_matrixp = sparsematrix_malloc0_ptr(inv_ovrlp_smat, &
  !                         iaction=DENSE_MATMUL, id='inv_ovrlp_matrixp')
  !!inv_ovrlp_matrixp = f_malloc_ptr((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,ncalc/),&
  !!                                  id='inv_ovrlp_matrixp')


      spin_loop: do ispin=1,ovrlp_smat%nspin

          degree_multiplicator = real(norder_polynomial,kind=8)/ &
                                 (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
          degree_multiplicator = min(degree_multiplicator,DEGREE_MULTIPLICATOR_MAX)

          isshift=(ispin-1)*ovrlp_smat%nvctr
          ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
          ilshift2=(ispin-1)*inv_ovrlp_smat%nvctr

          evlow_old=1.d100
          evhigh_old=-1.d100
          
    
        
              !!calculate_SHS=.true.
        
          !if (inv_ovrlp_smat%smmm%nfvctrp>0) then !LG: this conditional seems decorrelated
          !call f_zero(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp*ncalc, inv_ovrlp_matrixp(1,1,1))
          !end if
          !!    call f_zero(inv_ovrlp_matrixp)
              
        
              it=0
              eval_bounds_ok=.false.
              !!bisection_bounds_ok=.false.
              main_loop: do 
                  
                  it=it+1
        
                  ! Scale the Hamiltonian such that all eigenvalues are in the intervall [0:1]
                  if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                      !!call scale_and_shift_matrix()
                      call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
                           ovrlp_smat, ovrlp_mat, isshift, &
                           matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
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
        
                  ! Array that holds the Chebyshev polynomials. Needs to be recalculated
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
    
        
                  cc = f_malloc((/npl,3,ncalc/),id='cc')
        
                  !!if (foe_data_get_real(foe_obj,"evlow")>=0.d0) then
                  !!    stop 'ERROR: lowest eigenvalue must be negative'
                  !!end if
                  if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      stop 'ERROR: highest eigenvalue must be positive'
                  end if
        
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
                  call timing(iproc, 'chebyshev_coef', 'ON')
        
                  do icalc=1,ncalc
                      call cheb_exp(foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,1,icalc), ex(icalc))
                      call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), cc(1,1,icalc), cc(1,2,icalc), npl)
                      call chebft2(foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,3,icalc))
                      call evnoise(npl, cc(1,3,icalc), foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
                  end do
    
                  call timing(iproc, 'chebyshev_coef', 'OF')
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
                
                
                
                  call timing(iproc, 'FOE_auxiliary ', 'OF')
        
                  emergency_stop=.false.
                  if (calculate_SHS) then
                      ! Passing inv_ovrlp(1)%matrix_compr as it will not be
                      ! used, to be improved...
                      call chebyshev_clean(iproc, nproc, npl, cc, &
                           inv_ovrlp_smat, hamscal_compr, &
                           inv_ovrlp(1)%matrix_compr(ilshift2+1:), .false., &
                           nsize_polynomial, ncalc, inv_ovrlp_matrixp_new, penalty_ev_new, chebyshev_polynomials, &
                           emergency_stop)
                      do icalc=1,ncalc
                          call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                               inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                               inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                               inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                               inv_ovrlp_matrixp_new(1,icalc), inv_ovrlp_matrixp_small_new(1,icalc))
                      end do

                       !write(*,'(a,i5,2es24.8)') 'iproc, sum(inv_ovrlp_matrixp(:,:,1:2)', (sum(inv_ovrlp_matrixp(:,:,icalc)),icalc=1,ncalc)
                      !!do i=1,inv_ovrlp_smat%smmm%nvctrp
                      !!    ii = inv_ovrlp_smat%smmm%isvctr + i
                      !!    call get_line_and_column(ii, inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, iline, icolumn)
                      !!    do icalc=1,ncalc
                      !!        inv_ovrlp_matrixp(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,icalc) = inv_ovrlp_matrixp_new(i,icalc)
                      !!    end do
                      !!    !!penalty_ev(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,1) = penalty_ev_new(i,1)
                      !!    !!penalty_ev(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,2) = penalty_ev_new(i,2)
                      !!end do
                  else
                      ! The Chebyshev polynomials are already available
                      !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                           inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
                           inv_ovrlp_smat, chebyshev_polynomials, ncalc, cc, inv_ovrlp_matrixp_new)
                      do icalc=1,ncalc
                          write(*,*) 'sum(inv_ovrlp_matrixp_new(:,icalc))',sum(inv_ovrlp_matrixp_new(:,icalc))
                      end do
                      !!do icalc=1,ncalc
                      !!    call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                      !!         inv_ovrlp_smat, inv_ovrlp_matrixp_new, inv_ovrlp_matrixp(:,:,icalc))
                      !!end do
                  end if 
    
    
    
                 !!! Check for an emergency stop, which happens if the kernel explodes, presumably due
                 !!! to the eigenvalue bounds being too small.
                 !!call check_emergency_stop(nproc,emergency_stop)
                 !!if (emergency_stop) then
                 !!     eval_bounds_ok(1)=.false.
                 !!     call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
                 !!     eval_bounds_ok(2)=.false.
                 !!     call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                 !!     call f_free(cc)
                 !!     cycle main_loop
                 !!end if
        
        
                  call timing(iproc, 'FOE_auxiliary ', 'ON')
        
        
                  restart=.false.
        
                  ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                  ! (otherwise this has already been checked in the previous iteration).
                  if (calculate_SHS) then
                      !call check_eigenvalue_spectrum()
                      !!call check_eigenvalue_spectrum(nproc, inv_ovrlp_smat, ovrlp_smat, ovrlp_mat, 1, &
                      !!     0, 1.2d0, 1.d0/1.2d0, penalty_ev, anoise, .false., emergency_stop, &
                      !!     foe_obj, restart, eval_bounds_ok)
                      call check_eigenvalue_spectrum_new(nproc, inv_ovrlp_smat, ovrlp_smat, ovrlp_mat, 1, &
                           0, 1.2d0, 1.d0/1.2d0, penalty_ev_new, anoise, .false., emergency_stop, &
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
        
        
    
          do icalc=1,ncalc
              !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, inv_ovrlp_matrixp(1:,1:,icalc), &
              !!     inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
              call compress_matrix_distributed_new2(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, inv_ovrlp_matrixp_small_new(:,icalc), &
                   inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
          end do
    

      end do spin_loop

  call f_free_ptr(inv_ovrlp_matrixp)
  call f_free(inv_ovrlp_matrixp_small_new)
  call f_free(inv_ovrlp_matrixp_new)
  call f_free(chebyshev_polynomials)
  !!call f_free(penalty_ev)
  call f_free(penalty_ev_new)
  call f_free(hamscal_compr)

  call f_free_ptr(foe_obj%ef)
  call f_free_ptr(foe_obj%evlow)
  call f_free_ptr(foe_obj%evhigh)
  call f_free_ptr(foe_obj%bisection_shift)
  call f_free_ptr(foe_obj%charge)

  call timing(iproc, 'FOE_auxiliary ', 'OF')

  call f_release_routine()


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



!!subroutine check_eigenvalue_spectrum(nproc, smat_l, smat_s, mat, ispin, isshift, &
!!           factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
!!           emergency_stop, foe_obj, restart, eval_bounds_ok)
!!  use module_base
!!  use sparsematrix_base, only: sparse_matrix, matrices
!!  use sparsematrix_init, only: matrixindex_in_compressed
!!  use foe_base, only: foe_data, foe_data_set_real, foe_data_get_real
!!  use yaml_output
!!  implicit none
!!
!!  ! Calling arguments
!!  type(sparse_matrix),intent(in) :: smat_l, smat_s
!!  type(matrices),intent(in) :: mat
!!  integer,intent(in) :: nproc, ispin, isshift
!!  real(kind=8),intent(in) :: factor_high, factor_low, anoise
!!  real(kind=8),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
!!  logical,intent(in) :: trace_with_overlap, emergency_stop
!!  type(foe_data),intent(inout) :: foe_obj
!!  logical,intent(out) :: restart
!!  logical,dimension(2),intent(out) :: eval_bounds_ok
!!
!!  ! Local variables
!!  integer :: isegstart, isegend, iseg, ii, jorb, irow, icol, iismall, iel
!!  real(kind=8) :: bound_low, bound_up, tt, noise
!!  real(kind=8),dimension(2) :: allredarr
!!
!!  call f_routine(id='check_eigenvalue_spectrum')
!!
!!  if (.not.emergency_stop) then
!!      ! The penalty function must be smaller than the noise.
!!      bound_low=0.d0
!!      bound_up=0.d0
!!      if (smat_l%smmm%nfvctrp>0) then
!!          !$omp parallel default(none) &
!!          !$omp private(iseg, ii, jorb, irow, icol, iismall, tt, iel) &
!!          !$omp shared(isegstart, isegend, smat_l, smat_s, mat, penalty_ev) &
!!          !$omp shared(bound_low, bound_up, isshift, trace_with_overlap) 
!!          !$omp do reduction(+:bound_low,bound_up)
!!          !!do iseg=isegstart,isegend
!!          do iseg=smat_l%smmm%isseg,smat_l%smmm%ieseg
!!              iel = smat_l%keyv(iseg) - 1
!!              ii=smat_l%keyv(iseg)-1
!!              ! A segment is always on one line, therefore no double loop
!!              do jorb=smat_l%keyg(1,1,iseg),smat_l%keyg(2,1,iseg)
!!                  iel = iel + 1
!!                  if (iel<smat_l%smmm%isvctr_mm+1) cycle
!!                  if (iel>smat_l%smmm%isvctr_mm+smat_l%smmm%nvctrp_mm) exit
!!                  ii=ii+1
!!                  irow = smat_l%keyg(1,2,iseg)
!!                  icol = jorb
!!                  iismall = matrixindex_in_compressed(smat_s, irow, icol)
!!                  if (iismall>0) then
!!                      if (trace_with_overlap) then
!!                          ! Take the trace of the product matrix times overlap
!!                          tt=mat%matrix_compr(isshift+iismall-smat_s%isvctrp_tg)
!!                      else
!!                          ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
!!                          if (irow==icol) then
!!                              tt=1.d0
!!                          else
!!                              tt=0.d0
!!                          end if
!!                      end if
!!                  else
!!                      tt=0.d0
!!                  end if
!!                  bound_low = bound_low + penalty_ev(icol,irow-smat_l%smmm%isfvctr,2)*tt
!!                  bound_up = bound_up +penalty_ev(icol,irow-smat_l%smmm%isfvctr,1)*tt
!!              end do  
!!          end do
!!          !$omp end do
!!          !$omp end parallel
!!      end if
!!  else
!!      ! This means that the Chebyshev expansion exploded, so take a very large
!!      ! value for the error function such that eigenvalue bounds will be enlarged
!!      bound_low = 1.d10
!!      bound_up = 1.d10
!!  end if
!!
!!  allredarr(1)=bound_low
!!  allredarr(2)=bound_up
!!
!!  if (nproc > 1) then
!!      call mpiallred(allredarr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
!!  end if
!!
!!
!!  allredarr=abs(allredarr) !for some crazy situations this may be negative
!!  noise=1000.d0*anoise
!!
!!  if (bigdft_mpi%iproc==0) then
!!      call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
!!  end if
!!  !write(*,*) 'allredarr, anoise', allredarr, anoise
!!  if (allredarr(1)>noise) then
!!      eval_bounds_ok(1)=.false.
!!      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
!!      restart=.true.
!!      !!if (bigdft_mpi%iproc==0) then
!!      !!    call yaml_map('adjust lower bound',.true.)
!!      !!end if
!!  else
!!      eval_bounds_ok(1)=.true.
!!  end if
!!  if (allredarr(2)>noise) then
!!      eval_bounds_ok(2)=.false.
!!      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*factor_high,ispin)
!!      restart=.true.
!!      !!if (bigdft_mpi%iproc==0) then
!!      !!    call yaml_map('adjust upper bound',.true.)
!!      !!end if
!!  else
!!      eval_bounds_ok(2)=.true.
!!  end if
!!
!!  call f_release_routine()
!!
!!end subroutine check_eigenvalue_spectrum




subroutine check_eigenvalue_spectrum_new(nproc, smat_l, smat_s, mat, ispin, isshift, &
           factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
           emergency_stop, foe_obj, restart, eval_bounds_ok)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  use sparsematrix_init, only: matrixindex_in_compressed, get_line_and_column
  use foe_base, only: foe_data, foe_data_set_real, foe_data_get_real
  use yaml_output
  implicit none

  ! Calling arguments
  type(sparse_matrix),intent(in) :: smat_l, smat_s
  type(matrices),intent(in) :: mat
  integer,intent(in) :: nproc, ispin, isshift
  real(kind=8),intent(in) :: factor_high, factor_low, anoise
  !real(kind=8),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
  real(kind=8),dimension(smat_l%smmm%nvctrp,2),intent(in) :: penalty_ev
  logical,intent(in) :: trace_with_overlap, emergency_stop
  type(foe_data),intent(inout) :: foe_obj
  logical,intent(out) :: restart
  logical,dimension(2),intent(out) :: eval_bounds_ok

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, irow, icol, iismall, iel, i, iline, icolumn
  real(kind=8) :: bound_low, bound_up, tt, noise
  real(kind=8),dimension(2) :: allredarr

  call f_routine(id='check_eigenvalue_spectrum_new')

  if (.not.emergency_stop) then
      ! The penalty function must be smaller than the noise.
      bound_low=0.d0
      bound_up=0.d0
      !!if (smat_l%smmm%nfvctrp>0) then
      !!    !$omp parallel default(none) &
      !!    !$omp private(iseg, ii, jorb, irow, icol, iismall, tt, iel) &
      !!    !$omp shared(isegstart, isegend, smat_l, smat_s, mat, penalty_ev) &
      !!    !$omp shared(bound_low, bound_up, isshift, trace_with_overlap) 
      !!    !$omp do reduction(+:bound_low,bound_up)
      !!    !!do iseg=isegstart,isegend
      !!    do iseg=smat_l%smmm%isseg,smat_l%smmm%ieseg
      !!        iel = smat_l%keyv(iseg) - 1
      !!        ii=smat_l%keyv(iseg)-1
      !!        ! A segment is always on one line, therefore no double loop
      !!        do jorb=smat_l%keyg(1,1,iseg),smat_l%keyg(2,1,iseg)
      !!            iel = iel + 1
      !!            if (iel<smat_l%smmm%isvctr_mm+1) cycle
      !!            if (iel>smat_l%smmm%isvctr_mm+smat_l%smmm%nvctrp_mm) exit
      !!            ii=ii+1
      !!            irow = smat_l%keyg(1,2,iseg)
      !!            icol = jorb
      !!            iismall = matrixindex_in_compressed(smat_s, irow, icol)
      !!            if (iismall>0) then
      !!                if (trace_with_overlap) then
      !!                    ! Take the trace of the product matrix times overlap
      !!                    tt=mat%matrix_compr(isshift+iismall-smat_s%isvctrp_tg)
      !!                else
      !!                    ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
      !!                    if (irow==icol) then
      !!                        tt=1.d0
      !!                    else
      !!                        tt=0.d0
      !!                    end if
      !!                end if
      !!            else
      !!                tt=0.d0
      !!            end if
      !!            bound_low = bound_low + penalty_ev(icol,irow-smat_l%smmm%isfvctr,2)*tt
      !!            bound_up = bound_up +penalty_ev(icol,irow-smat_l%smmm%isfvctr,1)*tt
      !!        end do  
      !!    end do
      !!    !$omp end do
      !!    !$omp end parallel
      !!end if


      do i=1,smat_l%smmm%nvctrp
          ii = smat_l%smmm%isvctr + i
          call get_line_and_column(ii, smat_l%smmm%nseg, smat_l%smmm%keyv, smat_l%smmm%keyg, iline, icolumn)
          !!iismall = matrixindex_in_compressed_fn(icolumn, iline, &
          !!          smat_s%nfvctr, smat_l%smmm%nseg_mm, smat_l%smmm%keyv_mm, smat_l%smmm%keyg_mm)
          iismall = matrixindex_in_compressed(smat_s, icolumn, iline)
          if (iismall>0) then
              if (trace_with_overlap) then
                  ! Take the trace of the product matrix times overlap
                  tt=mat%matrix_compr(isshift+iismall-smat_s%isvctrp_tg)
              else
                  ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
                  if (irow==icol) then
                      tt=1.d0
                  else
                      tt=0.d0
                  end if
              end if
          else
              tt=0.d0
          end if
          !!bound_low = bound_low + penalty_ev(icol,irow-smat_l%smmm%isfvctr,2)*tt
          !!bound_up = bound_up +penalty_ev(icol,irow-smat_l%smmm%isfvctr,1)*tt
          bound_low = bound_low + penalty_ev(i,2)*tt
          bound_up = bound_up +penalty_ev(i,1)*tt
      end do
  else
      ! This means that the Chebyshev expansion exploded, so take a very large
      ! value for the error function such that eigenvalue bounds will be enlarged
      bound_low = 1.d10
      bound_up = 1.d10
  end if

  allredarr(1)=bound_low
  allredarr(2)=bound_up

  if (nproc > 1) then
      call mpiallred(allredarr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
  end if


  allredarr=abs(allredarr) !for some crazy situations this may be negative
  noise=1000.d0*anoise

  if (bigdft_mpi%iproc==0) then
      call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
  end if
  !write(*,*) 'allredarr, anoise', allredarr, anoise
  if (allredarr(1)>noise) then
      eval_bounds_ok(1)=.false.
      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
      restart=.true.
      !!if (bigdft_mpi%iproc==0) then
      !!    call yaml_map('adjust lower bound',.true.)
      !!end if
  else
      eval_bounds_ok(1)=.true.
  end if
  if (allredarr(2)>noise) then
      eval_bounds_ok(2)=.false.
      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*factor_high,ispin)
      restart=.true.
      !!if (bigdft_mpi%iproc==0) then
      !!    call yaml_map('adjust upper bound',.true.)
      !!end if
  else
      eval_bounds_ok(2)=.true.
  end if

  call f_release_routine()

end subroutine check_eigenvalue_spectrum_new





subroutine check_emergency_stop(nproc,emergency_stop)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in) :: nproc
  logical,intent(inout) :: emergency_stop

  ! Local variables
  integer :: iflag

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


subroutine scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
           smat1, mat1, i1shift, smat2, mat2, i2shift, &
           matscal_compr, scale_factor, shift_value)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  use foe_base, only: foe_data, foe_data_get_real
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: orb_from_index
  implicit none
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ispin, i1shift
  type(foe_data),intent(in) :: foe_obj
  type(sparse_matrix),intent(in) :: smatl, smat1
  type(matrices),intent(in) :: mat1
  type(sparse_matrix),intent(in),optional :: smat2
  type(matrices),intent(in),optional :: mat2
  integer,intent(in),optional :: i2shift
  real(kind=8),dimension(smatl%nvctrp_tg),intent(out) :: matscal_compr
  real(kind=8),intent(out) :: scale_factor, shift_value

  ! Local variables
  integer :: iseg, ii, i, ii1, ii2, isegstart, isegend, ierr
  integer :: itaskgroup, iitaskgroup, j
  integer,dimension(2) :: irowcol
  real(kind=8) :: tt1, tt2
  logical :: with_overlap
  real(kind=8),dimension(:),pointer :: matscal_compr_local
  integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
  integer,parameter :: comm_strategy=GET
  integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX

  

  call f_routine(id='scale_and_shift_matrix')
  call timing(iproc,'foe_aux_mcpy  ','ON')

  ! smat2 and mat2 must be present at the same time
  if (all((/present(smat2),present(mat2),present(i2shift)/))) then
      with_overlap = .true.
  else
      if (any((/present(smat2),present(mat2),present(i2shift)/))) then
          stop 'smat2, mat2 and i2shift must be present at the same time'
      end if
      with_overlap = .false.
  end if

  scale_factor=2.d0/(foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
  shift_value=.5d0*(foe_data_get_real(foe_obj,"evhigh",ispin)+foe_data_get_real(foe_obj,"evlow",ispin))

  if (data_strategy==GLOBAL_MATRIX) then
      stop 'scale_and_shift_matrix: data_strategy=GLOBAL_MATRIX is deprecated'
      !!isegstart = smatl%istsegline(smatl%isfvctr+1)
      !!isegend = smatl%istsegline(smatl%isfvctr+smatl%nfvctrp) + &
      !!          smatl%nsegline(smatl%isfvctr+smatl%nfvctrp)-1
      !!if (nproc>1) then
      !!    matscal_compr_local = f_malloc_ptr(smatl%nvctrp,id='matscal_compr_local')
      !!else
      !!    matscal_compr_local => matscal_compr
      !!end if
      !!!$omp parallel default(none) private(iseg,ii,i,irowcol,ii2,ii1,tt2,tt1) &
      !!!$omp shared(matscal_compr_local,scale_factor,shift_value,i2shift,i1shift,smatl,smat1,smat2,mat1,mat2,with_overlap) &
      !!!$omp shared(isegstart,isegend)
      !!!$omp do
      !!do iseg=isegstart,isegend
      !!    ii=smatl%keyv(iseg)
      !!    do i=smatl%keyg(1,iseg),smatl%keyg(2,iseg)
      !!        irowcol = orb_from_index(smatl,i)
      !!        ii1 = matrixindex_in_compressed(smat1, irowcol(1), irowcol(2))
      !!        if (ii1>0) then
      !!            tt1=mat1%matrix_compr(i1shift+ii1)
      !!        else
      !!            tt1=0.d0
      !!        end if
      !!        if (with_overlap) then
      !!            ii2 = matrixindex_in_compressed(smat2, irowcol(1), irowcol(2))
      !!            if (ii2>0) then
      !!                tt2=mat2%matrix_compr(i2shift+ii2)
      !!            else
      !!                tt2=0.d0
      !!            end if
      !!        else
      !!            if (irowcol(1)==irowcol(2)) then
      !!                tt2 = 1.d0
      !!            else
      !!                tt2 = 0.d0
      !!            end if
      !!        end if
      !!        !write(*,*) 'ii, tt1, tt2', ii, tt1, tt2
      !!        matscal_compr_local(ii-smatl%isvctr)=scale_factor*(tt1-shift_value*tt2)
      !!        ii=ii+1
      !!    end do
      !!end do
      !!!$omp end do
      !!!$omp end parallel

      !!call timing(iproc,'foe_aux_mcpy  ','OF')
      !!call timing(iproc,'foe_aux_comm  ','ON')
      !!if (nproc>1) then
      !!    !!call mpi_allgatherv(matscal_compr_local(1), smatl%nvctrp, mpi_double_precision, &
      !!    !!     matscal_compr(1), smatl%nvctr_par, smatl%isvctr_par, mpi_double_precision, &
      !!    !!     bigdft_mpi%mpi_comm, ierr)
      !!    if (comm_strategy==ALLGATHERV) then
      !!        call mpi_allgatherv(matscal_compr_local(1), smatl%nvctrp, mpi_double_precision, &
      !!             matscal_compr(1), smatl%nvctr_par, smatl%isvctr_par, mpi_double_precision, &
      !!             bigdft_mpi%mpi_comm, ierr)
      !!        call f_free_ptr(matscal_compr_local)
      !!    else if (comm_strategy==GET) then
      !!        !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, smatl%nvctrp, matscal_compr_local, &
      !!        !!     smatl%nvctr_par, smatl%isvctr_par, smatl%nvctr, matscal_compr)
      !!        call mpi_get_to_allgatherv(matscal_compr_local(1), smatl%nvctrp, matscal_compr(1), &
      !!             smatl%nvctr_par, smatl%isvctr_par, bigdft_mpi%mpi_comm)
      !!    else
      !!        stop 'scale_and_shift_matrix: wrong communication strategy'
      !!    end if
      !!    call f_free_ptr(matscal_compr_local)
      !!end if
      !!call timing(iproc,'foe_aux_comm  ','OF')

  else if (data_strategy==SUBMATRIX) then
      !$omp parallel default(none) private(ii,i,j,ii2,ii1,tt2,tt1,iseg) &
      !$omp shared(matscal_compr,scale_factor,shift_value,i2shift,i1shift,smatl,smat1,smat2,mat1,mat2,with_overlap)
      !$omp do
      do iseg=smatl%smmm%istartendseg_mm(1),smatl%smmm%istartendseg_mm(2)
          !if (smatl%keyv(min(iseg+1,smatl%nseg))<smatl%smmm%istartend_mm(1)) cycle
          !if (smatl%keyv(iseg)>smatl%smmm%istartend_mm(2)) exit
          ! A segment is always on one line, therefore no double loop
          j = smatl%keyg(1,2,iseg)
          do i=smatl%keyg(1,1,iseg),smatl%keyg(2,1,iseg) !this is too much, but for the moment ok 
              ii1 = matrixindex_in_compressed(smat1, i, j)
              if (ii1>0) then
                  tt1=mat1%matrix_compr(i1shift+ii1-smat1%isvctrp_tg)
              else
                  tt1=0.d0
              end if
              if (with_overlap) then
                  ii2 = matrixindex_in_compressed(smat2, i, j)
                  if (ii2>0) then
                      tt2=mat2%matrix_compr(i2shift+ii2-smat2%isvctrp_tg)
                  else
                      tt2=0.d0
                  end if
              else
                  if (i==j) then
                      tt2 = 1.d0
                  else
                      tt2 = 0.d0
                  end if
              end if
              ii=matrixindex_in_compressed(smatl, i, j)
              !write(*,*) 'i, ii, tt1, tt2', i, ii, tt1, tt2
              matscal_compr(ii-smatl%isvctrp_tg)=scale_factor*(tt1-shift_value*tt2)
          end do
      end do
      !$omp end do
      !$omp end parallel
      call timing(iproc,'foe_aux_mcpy  ','OF')
  else
      stop 'scale_and_shift_matrix: wrong data strategy'
  end if

  !!do i=1,smatl%nvctr
  !!    write(500+iproc,*) 'i, matscal_compr(i)', i, matscal_compr(i)
  !!end do
  call f_release_routine()

end subroutine scale_and_shift_matrix
