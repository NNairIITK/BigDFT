module ice
  implicit none

  private

  !> Public routines
  public :: inverse_chebyshev_expansion
  public :: inverse_chebyshev_expansion_new

  contains

    ! New: chebyshev expansion of the inverse overlap (Inverse Chebyshev Expansion)
    subroutine inverse_chebyshev_expansion(iproc, nproc, norder_polynomial, &
               ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp, &
               npl_auto)
      use module_base
      use yaml_output
      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, &
                                   sparsematrix_malloc0_ptr, assignment(=), &
                                   SPARSE_TASKGROUP, SPARSE_MATMUL_SMALL, &
                                   matrices, sparse_matrix, matrices_null, deallocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed_wrapper, &
                              transform_sparsity_pattern
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_set_logical, foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: scale_and_shift_matrix, &
                            evnoise, check_eigenvalue_spectrum_new, get_chebyshev_expansion_coefficients, &
                            get_chebyshev_polynomials
      use module_func
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norder_polynomial, ncalc
      type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
      real(kind=8),dimension(ncalc),intent(in) :: ex
      type(matrices),intent(in) :: ovrlp_mat
      type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp
      logical,intent(in),optional :: npl_auto
    
      ! Local variables
      integer :: npl, jorb, it, ii, iseg
      integer :: isegstart, isegend, iismall, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i, ipl
      integer,parameter :: nplx=50000
      real(kind=8),dimension(:,:),pointer :: chebyshev_polynomials
      real(kind=8),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
      real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=8),dimension(:,:,:),pointer :: cc
      real(kind=8) :: anoise, scale_factor, shift_value
      real(kind=8) :: evlow_old, evhigh_old, tt
      real(kind=8) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=8) :: tt_ovrlp, tt_ham, eval_multiplicator, eval_multiplicator_total
      logical :: restart, calculate_SHS
      logical,dimension(2) :: emergency_stop
      real(kind=8),dimension(2) :: allredarr
      real(kind=8),dimension(:),allocatable :: hamscal_compr
      logical,dimension(2) :: eval_bounds_ok
      integer,dimension(2) :: irowcol
      integer :: irow, icol, iflag, ispin, isshift, ilshift, ilshift2
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit, npl_auto_
      integer,parameter :: NPL_MIN=5
      real(kind=8),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
      real(kind=8) :: degree_multiplicator
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(foe_data) :: foe_obj
      real(kind=8),dimension(:),allocatable :: eval, work, x_max_error, max_error, mean_error
      real(kind=8),dimension(:,:),allocatable :: tempmat
      integer :: lwork, info, j, icalc, iline, icolumn
      real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_new
      real(kind=8),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_small_new
      type(matrices) :: ovrlp_scaled
      character(len=3),parameter :: old='old'
      character(len=3),parameter :: new='new'
      character(len=3) :: mode=old
    
      !!real(kind=8),dimension(ovrlp_smat%nfvctr,ovrlp_smat%nfvctr) :: overlap
      !!real(kind=8),dimension(ovrlp_smat%nfvctr) :: eval
      !!integer,parameter :: lwork=100000
      !!real(kind=8),dimension(lwork) :: work
      !!integer :: info
    
      call f_routine(id='inverse_chebyshev_expansion')

      if (present(npl_auto)) then
          npl_auto_ = npl_auto
      else
          npl_auto_ = .false.
      end if
    
    
      penalty_ev_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,2/),id='penalty_ev_new')
      inv_ovrlp_matrixp_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp,1),ncalc/),id='inv_ovrlp_matrixp_new')
      inv_ovrlp_matrixp_small_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp_mm,1),ncalc/),id='inv_ovrlp_matrixp_small_new')
    
    
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
    !@ ################################################
    
    
      evbounds_shrinked = .false.
    
      !@ TEMPORARY: eigenvalues of  the overlap matrix ###################
      call get_minmax_eigenvalues(iproc, ovrlp_smat, ovrlp_mat)

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
      !@ END TEMPORARY: eigenvalues of  the overlap matrix ###############
    
    
      call timing(iproc, 'FOE_auxiliary ', 'ON')
    
    
    
      !!penalty_ev = f_malloc((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,2/),id='penalty_ev')
    
    
      if (npl_auto_) then
          ovrlp_scaled = matrices_null()
          ovrlp_scaled%matrix_compr = sparsematrix_malloc_ptr(ovrlp_smat, &
              iaction=SPARSE_TASKGROUP, id='ovrlp_scaled%matrix_compr')
          call f_memcpy(src=ovrlp_mat%matrix_compr,dest=ovrlp_scaled%matrix_compr)
      end if
      hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')
    
        
      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = inv_ovrlp_smat%smmm%nvctrp_mm
      
      
      ! Fake allocation, will be modified later
      chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,1/),id='chebyshev_polynomials')
    
    
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
              
              eval_multiplicator = 1.d0
              eval_multiplicator_total = 1.d0
        
            
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
                          if (npl_auto_) then
                              call dscal(size(ovrlp_scaled%matrix_compr), eval_multiplicator, ovrlp_scaled%matrix_compr(1), 1)
                              eval_multiplicator_total = eval_multiplicator_total*eval_multiplicator
                              call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
                                   ovrlp_smat, ovrlp_scaled, isshift, &
                                   matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
                              if (iproc==0) then
                                  write(*,*) 'eval_multiplicator, eval_multiplicator_total', &
                                              eval_multiplicator, eval_multiplicator_total
                              end if
                          else
                              call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
                                   ovrlp_smat, ovrlp_mat, isshift, &
                                   matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
                          end if
                          calculate_SHS=.true.
                      else
                          calculate_SHS=.false.
                      end if
                      !!do i=1,size(ovrlp_mat%matrix_compr)
                      !!    write(900+iproc,*) i, ovrlp_mat%matrix_compr(i)
                      !!end do
                      !!do i=1,size(hamscal_compr)
                      !!    write(950+iproc,*) i, hamscal_compr(i)
                      !!end do
                      evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
                      evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)
        
        
                      !call uncompress_matrix(iproc,ovrlp_smat,ovrlp_mat%matrix_compr,overlap)
                      !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
                      !if (iproc==0) write(*,*) 'ovrlp_mat%matrix_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
                      !call uncompress_matrix(iproc,inv_ovrlp_smat,hamscal_compr,overlap)
                      !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
                      !if (iproc==0) write(*,*) 'hamscal_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
            
            
                      ! Determine the degree of the polynomial
                      if (.not. npl_auto_) then
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
                      else
                          call get_poynomial_degree(iproc, nproc, ispin, ncalc, ex, foe_obj, 5, 100, 1.d-10, &
                               npl, cc, anoise)
                      end if

                      ! Array that holds the Chebyshev polynomials. Needs to be recalculated
                      ! every time the Hamiltonian has been modified.
                      if (iproc==0) then
                          call yaml_newline()
                          call yaml_mapping_open('ICE')
                          call yaml_map('eval bounds',&
                               (/foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          call yaml_map('mult.',degree_multiplicator,fmt='(f5.2)')
                          call yaml_map('pol. deg.',npl)
                      end if
                      if (calculate_SHS) then
                          call f_free_ptr(chebyshev_polynomials)
                          chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,npl/),id='chebyshev_polynomials')
                      end if
        
            
                      if (.not. npl_auto_) then
                          cc = f_malloc_ptr((/npl,3,ncalc/),id='cc')
            
                          !!if (foe_data_get_real(foe_obj,"evlow")>=0.d0) then
                          !!    stop 'ERROR: lowest eigenvalue must be negative'
                          !!end if
                          if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                              stop 'ERROR: highest eigenvalue must be positive'
                          end if
            
                          call timing(iproc, 'FOE_auxiliary ', 'OF')
                          call timing(iproc, 'chebyshev_coef', 'ON')
            
                          max_error = f_malloc(ncalc,id='max_error')
                          x_max_error = f_malloc(ncalc,id='x_max_error')
                          mean_error = f_malloc(ncalc,id='mean_error')
                          do icalc=1,ncalc
                              call func_set(FUNCTION_POLYNOMIAL, powerx=ex(icalc))
                              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                                   foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1:,1:,icalc), &
                                   x_max_error(icalc), max_error(icalc), mean_error(icalc))
                              !!##call cheb_exp(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                              !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1:,1:,icalc:), ex(icalc), &
                              !!##     x_max_error(icalc), max_error(icalc), mean_error(icalc))
                              !call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                              !     foe_data_get_real(foe_obj,"evhigh",ispin), cc(1:,1:,icalc:), cc(1:,2:,icalc:), npl)
                              call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                                   muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                                   foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1:,2:,icalc), &
                                   x_max_error_fake, max_error_fake, mean_error_fake)
                              do ipl=1,npl
                                 cc(ipl,3,1) = -cc(ipl,2,1)
                              end do
                              !!##call chebyshev_coefficients_penalyfunction(foe_data_get_real(foe_obj,"evlow",ispin), &
                              !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1:,2:,icalc:), max_error_fake)
                              call evnoise(npl, cc(1:,2:,icalc:), foe_data_get_real(foe_obj,"evlow",ispin), &
                                   foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
                          end do
                          if (iproc==0) then
                              call yaml_mapping_open('accuracy (x, max err, mean err)')
                              do icalc=1,ncalc
                                  call yaml_map('Operation '//trim(yaml_toa(icalc)), &
                                      (/x_max_error(icalc),max_error(icalc),mean_error(icalc)/),fmt='(es9.2)')
                              end do
                              call yaml_mapping_close()
                          end if
                          call f_free(mean_error)
                          call f_free(max_error)
                          call f_free(x_max_error)
        
                          call timing(iproc, 'chebyshev_coef', 'OF')
                          call timing(iproc, 'FOE_auxiliary ', 'ON')
                      end if
                      call yaml_mapping_close()
                    
                      !!do j=1,npl
                      !!    write(*,*) 'in main: j, cc(j,1,1), cc(j,2,1)', j, cc(j,1,1), cc(j,2,1)
                      !!end do
                    
                    
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
            
                      emergency_stop=.false.
                      if (calculate_SHS) then
                          ! Passing inv_ovrlp(1)%matrix_compr as it will not be
                          ! used, to be improved...
                          call chebyshev_clean(iproc, nproc, npl, cc, &
                               inv_ovrlp_smat, hamscal_compr, &
                               .false., &
                               nsize_polynomial, ncalc, inv_ovrlp_matrixp_new, penalty_ev_new, chebyshev_polynomials, &
                               emergency_stop)
                           !@NEW#####################################################################################
                           if (mode==new) then
                               call f_free_ptr(chebyshev_polynomials)
                               call get_chebyshev_polynomials(iproc, nproc, 1, 2, npl, ovrlp_smat, inv_ovrlp_smat, &     
                                                         ovrlp_scaled, foe_obj, chebyshev_polynomials, eval_bounds_ok)
                               call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                                    inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
                                    inv_ovrlp_smat, chebyshev_polynomials, ncalc, cc, inv_ovrlp_matrixp_new)
                               penalty_ev_new = 0.d0
                           end if
                           !@END NEW#################################################################################
                           !!do i=1,size(inv_ovrlp_matrixp_new,1)
                           !!    write(400+iproc,*) i, inv_ovrlp_matrixp_new(i,1)
                           !!end do
                          if (inv_ovrlp_smat%smmm%nvctrp>0) then
                              do icalc=1,ncalc
                                  call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                       inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                       inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                       inv_ovrlp_smat%smmm%line_and_column_mm, &
                                       inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                       inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                       inv_ovrlp_smat%smmm%istsegline, 'large_to_small', &
                                       inv_ovrlp_matrixp_small_new(1,icalc), inv_ovrlp_matrixp_new(1,icalc))
                                 !!do i=1,size(inv_ovrlp_matrixp_small_new,1)
                                 !!    write(410+iproc,*) i, inv_ovrlp_matrixp_small_new(i,icalc)
                                 !!end do
                              end do
                          end if
    
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
                          if (mode==old) then
                              call check_eigenvalue_spectrum_new(nproc, inv_ovrlp_smat, ispin, &
                                   0, 1.2d0, 1.d0/1.2d0, penalty_ev_new, anoise, .false., emergency_stop, &
                                   foe_obj, restart, eval_bounds_ok, eval_multiplicator)
                          else if (mode==new) then
                              if (.not.eval_bounds_ok(1)) then
                                  ! lower bound too large
                                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
                                  restart=.true.
                                  eval_multiplicator = 2.0d0
                              else if (.not.eval_bounds_ok(2)) then
                                  ! upper bound too small
                                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                                  restart=.true.
                                  eval_multiplicator = 1.d0/2.0d0
                              end if
                          end if
                      end if
            
                      call f_free_ptr(cc)
            
                      if (restart) then
                          if(evbounds_shrinked) then
                              ! this shrink was not good, increase the saturation counter
                              call foe_data_set_int(foe_obj,"evboundsshrink_isatur", &
                                   foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
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
            
            
        
              !if (inv_ovrlp_smat%smmm%nvctrp>0) then
                  do icalc=1,ncalc
                      !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, inv_ovrlp_matrixp(1:,1:,icalc), &
                      !!     inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
                      call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, &
                           SPARSE_MATMUL_SMALL, inv_ovrlp_matrixp_small_new(:,icalc), &
                           inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
                  end do
              !end if

              if (npl_auto_) then
                  do icalc=1,ncalc
                      call dscal(inv_ovrlp_smat%nvctrp_tg, 1.d0/eval_multiplicator**ex(icalc), &
                           inv_ovrlp(icalc)%matrix_compr(ilshift2+1), 1)
                  end do
              end if
        
    
          end do spin_loop
    
      !call f_free_ptr(inv_ovrlp_matrixp)
      call f_free(inv_ovrlp_matrixp_small_new)
      call f_free(inv_ovrlp_matrixp_new)
      call f_free_ptr(chebyshev_polynomials)
      !!call f_free(penalty_ev)
      call f_free(penalty_ev_new)
      call f_free(hamscal_compr)
      if (npl_auto_) then
          call deallocate_matrices(ovrlp_scaled)
      end if
    
      call f_free_ptr(foe_obj%ef)
      call f_free_ptr(foe_obj%evlow)
      call f_free_ptr(foe_obj%evhigh)
      call f_free_ptr(foe_obj%bisection_shift)
      call f_free_ptr(foe_obj%charge)
    
      call timing(iproc, 'FOE_auxiliary ', 'OF')
    
      call f_release_routine()
    
    
    end subroutine inverse_chebyshev_expansion



    ! Determine the polynomial degree which yields the desired precision
    subroutine get_poynomial_degree(iproc, nproc, ispin, ncalc, ex, foe_obj, &
               npl_min, npl_max, max_polynomial_degree, npl, cc, anoise)
      use module_base
      use foe_base, only: foe_data, foe_data_get_real
      use foe_common, only: evnoise, get_chebyshev_expansion_coefficients
      use yaml_output
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ispin, ncalc
      integer,intent(in) :: npl_min, npl_max
      real(kind=8),dimension(ncalc),intent(in) :: ex
      type(foe_data),intent(in) :: foe_obj
      real(kind=8),intent(in) :: max_polynomial_degree
      integer,intent(out) :: npl
      real(kind=8),dimension(:,:,:),pointer,intent(inout) :: cc
      real(kind=8),intent(out) :: anoise

      ! Local variables
      integer :: ipl, icalc, j, jpl
      logical :: error_ok, found_degree
      real(kind=8),dimension(:),allocatable :: max_error, x_max_error, mean_error
      real(kind=8),dimension(:,:,:),allocatable :: cc_trial
      real(kind=8) :: x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction

      call f_routine(id='get_poynomial_degree')

      max_error = f_malloc(ncalc,id='max_error')
      x_max_error = f_malloc(ncalc,id='x_max_error')
      mean_error = f_malloc(ncalc,id='mean_error')

      if (npl_min<3) then
          call f_err_throw('npl_min must be at least 3')
      end if
      if (npl_min>npl_max) then
          call f_err_throw('npl_min must be smaller or equal than npl_max')
      end if

      if (iproc==0) then
          call yaml_sequence_open('Determine polynomial degree')
      end if

      cc_trial = f_malloc((/npl_max,3,ncalc/),id='cc_trial')

      found_degree = .false.
      degree_loop: do ipl=npl_min,npl_max
          
          if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
              stop 'ERROR: highest eigenvalue must be positive'
          end if
          
          call timing(iproc, 'FOE_auxiliary ', 'OF')
          call timing(iproc, 'chebyshev_coef', 'ON')
          
          do icalc=1,ncalc
              call func_set(FUNCTION_POLYNOMIAL, powerx=ex(icalc))
              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                   foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,1,icalc), &
                   x_max_error(icalc), max_error(icalc), mean_error(icalc))
              call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                   muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                   foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,2,icalc), &
                   x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction)
              do jpl=1,ipl
                  cc_trial(jpl,3,icalc) = -cc_trial(jpl,2,icalc)
              end do
          end do

          call timing(iproc, 'chebyshev_coef', 'OF')
          call timing(iproc, 'FOE_auxiliary ', 'ON')

          if (iproc==0) then
              call yaml_mapping_open(flow=.true.)
              call yaml_map('ipl',ipl)
              do icalc=1,ncalc
                  call yaml_map('Operation '//trim(yaml_toa(icalc)), &
                      (/x_max_error(icalc),max_error(icalc),mean_error(icalc),max_error_penaltyfunction/),fmt='(es9.2)')
              end do
              call yaml_mapping_close()
          end if

          error_ok = .true.
          do icalc=1,ncalc
              if (max_error(icalc)>max_polynomial_degree) then
                  error_ok = .false.
                  exit
              end if
          end do
          if (error_ok) then
              if (max_error_penaltyfunction>1.d-2) then
                  error_ok = .false.
              end if
          end if
          if (error_ok) then
              npl = ipl
              found_degree = .true.
              exit degree_loop
          end if


      end do degree_loop

      if (.not.found_degree) then
          call yaml_warning('Not possible to reach desired accuracy, using highest available polynomial degree')
          npl = npl_max
      end if

      if (iproc==0) then
          call yaml_sequence_close()
      end if

      cc = f_malloc_ptr((/npl,3,ncalc/),id='cc')
      do icalc=1,ncalc
          do j=1,3
              do ipl=1,npl
                  cc(ipl,j,icalc)=cc_trial(ipl,j,icalc)
              end do
          end do
      end do
      call f_free(cc_trial)
      call f_free(mean_error)
      call f_free(max_error)
      call f_free(x_max_error)

      call f_release_routine

    end subroutine get_poynomial_degree


    subroutine get_minmax_eigenvalues(iproc, ovrlp_smat, ovrlp_mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(in) :: ovrlp_smat
      type(matrices),intent(in) :: ovrlp_mat

      ! Local variables
      integer :: iseg, ii, i, lwork, info
      real(kind=8),dimension(:,:),allocatable :: tempmat
      real(kind=8),dimension(:),allocatable :: eval, work

      call f_routine(id='get_minmax_eigenvalues')

      tempmat = f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/),id='tempmat')
      do iseg=1,ovrlp_smat%nseg
          ii=ovrlp_smat%keyv(iseg)
          do i=ovrlp_smat%keyg(1,1,iseg),ovrlp_smat%keyg(2,1,iseg)
              tempmat(i,ovrlp_smat%keyg(1,2,iseg)) = ovrlp_mat%matrix_compr(ii)
              ii = ii + 1
          end do
      end do
      !!if (iproc==0) then
      !!    do i=1,ovrlp_smat%nfvctr
      !!        do j=1,ovrlp_smat%nfvctr
      !!            write(*,'(a,2i6,es17.8)') 'i,j,val',i,j,tempmat(j,i)
      !!        end do
      !!    end do
      !!end if
      eval = f_malloc(ovrlp_smat%nfvctr,id='eval')
      lwork=100*ovrlp_smat%nfvctr
      work = f_malloc(lwork,id='work')
      call dsyev('n','l', ovrlp_smat%nfvctr, tempmat, ovrlp_smat%nfvctr, eval, work, lwork, info)
      !if (iproc==0) write(*,*) 'eval',eval
      if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(ovrlp_smat%nfvctr)/),fmt='(es16.6)')
    
      call f_free(tempmat)
      call f_free(eval)
      call f_free(work)

      call f_release_routine()
    
    end subroutine get_minmax_eigenvalues


    ! New: chebyshev expansion of the inverse overlap (Inverse Chebyshev Expansion)
    subroutine inverse_chebyshev_expansion_new(iproc, nproc, norder_polynomial, &
               ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp, &
               npl_auto)
      use module_base
      use yaml_output
      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, &
                                   sparsematrix_malloc0_ptr, assignment(=), &
                                   SPARSE_TASKGROUP, SPARSE_MATMUL_SMALL, &
                                   matrices, sparse_matrix, matrices_null, deallocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed_wrapper, &
                              transform_sparsity_pattern
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_set_logical, foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: scale_and_shift_matrix, &
                            evnoise, check_eigenvalue_spectrum_new, get_chebyshev_expansion_coefficients, &
                            get_chebyshev_polynomials
      use module_func
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norder_polynomial, ncalc
      type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
      real(kind=8),dimension(ncalc),intent(in) :: ex
      type(matrices),intent(in) :: ovrlp_mat
      type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp
      logical,intent(in),optional :: npl_auto
    
      ! Local variables
      integer :: npl, jorb, it, ii, iseg
      integer :: isegstart, isegend, iismall, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i, ipl
      integer,parameter :: nplx=50000
      real(kind=8),dimension(:,:),pointer :: chebyshev_polynomials
      real(kind=8),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
      real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=8),dimension(:,:,:),pointer :: cc
      real(kind=8) :: anoise, scale_factor, shift_value
      real(kind=8) :: evlow_old, evhigh_old, tt
      real(kind=8) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=8) :: tt_ovrlp, tt_ham, eval_multiplicator, eval_multiplicator_total
      logical :: restart, calculate_SHS
      logical,dimension(2) :: emergency_stop
      real(kind=8),dimension(2) :: allredarr
      real(kind=8),dimension(:),allocatable :: hamscal_compr
      logical,dimension(2) :: eval_bounds_ok
      integer,dimension(2) :: irowcol
      integer :: irow, icol, iflag, ispin, isshift, ilshift, ilshift2
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit, npl_auto_
      integer,parameter :: NPL_MIN=5
      real(kind=8),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
      real(kind=8) :: degree_multiplicator
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(foe_data) :: foe_obj
      real(kind=8),dimension(:),allocatable :: eval, work, x_max_error, max_error, mean_error
      real(kind=8),dimension(:,:),allocatable :: tempmat
      integer :: lwork, info, j, icalc, iline, icolumn
      real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_new
      real(kind=8),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=8),dimension(:,:),allocatable :: inv_ovrlp_matrixp_small_new
      type(matrices) :: ovrlp_scaled
      character(len=3),parameter :: old='old'
      character(len=3),parameter :: new='new'
      character(len=3) :: mode=old

      call f_routine(id='inverse_chebyshev_expansion')

      inv_ovrlp_matrixp_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp,1),ncalc/),id='inv_ovrlp_matrixp_new')
      inv_ovrlp_matrixp_small_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp_mm,1),ncalc/),id='inv_ovrlp_matrixp_small_new')
      hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

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
           call foe_data_set_real(foe_obj,"evlow_min",0.d0)
           call foe_data_set_real(foe_obj,"evhigh_max",200.d0)
      !@ ################################################

      !@ TEMPORARY: eigenvalues of  the overlap matrix ###################
      call get_minmax_eigenvalues(iproc, ovrlp_smat, ovrlp_mat)

      ovrlp_scaled = matrices_null()
      ovrlp_scaled%matrix_compr = sparsematrix_malloc_ptr(ovrlp_smat, &
          iaction=SPARSE_TASKGROUP, id='ovrlp_scaled%matrix_compr')
      call f_memcpy(src=ovrlp_mat%matrix_compr,dest=ovrlp_scaled%matrix_compr)

      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = inv_ovrlp_smat%smmm%nvctrp_mm
      
      ! Fake allocation, will be modified later
      chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,1/),id='chebyshev_polynomials')

      spin_loop: do ispin=1,ovrlp_smat%nspin

          isshift=(ispin-1)*ovrlp_smat%nvctr
          ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
          ilshift2=(ispin-1)*inv_ovrlp_smat%nvctr

          eval_multiplicator = 1.d0
          eval_multiplicator_total = 1.d0

          bounds_loop: do
              call dscal(size(ovrlp_scaled%matrix_compr), eval_multiplicator, ovrlp_scaled%matrix_compr(1), 1)
              eval_multiplicator_total = eval_multiplicator_total*eval_multiplicator
              call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
                   ovrlp_smat, ovrlp_scaled, isshift, &
                   matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
              if (iproc==0) then
                  write(*,*) 'eval_multiplicator, eval_multiplicator_total', &
                              eval_multiplicator, eval_multiplicator_total
              end if
              call get_poynomial_degree(iproc, nproc, ispin, ncalc, ex, foe_obj, 5, 100, 1.d-10, &
                   npl, cc, anoise)
              call f_free_ptr(chebyshev_polynomials)
              call get_chebyshev_polynomials(iproc, nproc, 1, 2, npl, ovrlp_smat, inv_ovrlp_smat, &     
                                        ovrlp_scaled, foe_obj, chebyshev_polynomials, eval_bounds_ok)
              if (all(eval_bounds_ok)) then
                  exit bounds_loop
              else
                  if (.not.eval_bounds_ok(1)) then
                      ! lower bound too large
                      call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
                      eval_multiplicator = 2.0d0
                  else if (.not.eval_bounds_ok(2)) then
                      ! upper bound too small
                      call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                      eval_multiplicator = 1.d0/2.0d0
                  end if
              end if
              call f_free_ptr(cc)
              write(*,*) 'eval_bounds_ok',eval_bounds_ok
              write(*,*) 'evlow, evhigh',foe_data_get_real(foe_obj,"evlow",ispin), foe_data_get_real(foe_obj,"evhigh",ispin)
          end do bounds_loop
          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
               inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
               inv_ovrlp_smat, chebyshev_polynomials, ncalc, cc, inv_ovrlp_matrixp_new)
          ! TEST ##################################################
          call foe_data_set_real(foe_obj,"ef",1.d0,ispin)
          call foe_data_set_real(foe_obj,"charge",10.d0,ispin)
          call find_fermi_level(iproc, nproc, npl, chebyshev_polynomials, &
               2, 'test', inv_ovrlp_smat, foe_obj, inv_ovrlp(1))
          ! END TEST ##############################################
          if (inv_ovrlp_smat%smmm%nvctrp>0) then
              do icalc=1,ncalc
                  call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                       inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                       inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                       inv_ovrlp_smat%smmm%line_and_column_mm, &
                       inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                       inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                       inv_ovrlp_smat%smmm%istsegline, 'large_to_small', &
                       inv_ovrlp_matrixp_small_new(1,icalc), inv_ovrlp_matrixp_new(1,icalc))
              end do
          end if
          do icalc=1,ncalc
              call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, &
                   SPARSE_MATMUL_SMALL, inv_ovrlp_matrixp_small_new(:,icalc), &
                   inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
              call dscal(inv_ovrlp_smat%nvctrp_tg, 1.d0/eval_multiplicator**ex(icalc), &
                   inv_ovrlp(icalc)%matrix_compr(ilshift2+1), 1)
          end do

      end do spin_loop

      call f_free(inv_ovrlp_matrixp_small_new)
      call f_free(inv_ovrlp_matrixp_new)
      call f_free_ptr(chebyshev_polynomials)
      call f_free(hamscal_compr)
      call deallocate_matrices(ovrlp_scaled)
      call f_free_ptr(cc)
    
      call f_free_ptr(foe_obj%ef)
      call f_free_ptr(foe_obj%evlow)
      call f_free_ptr(foe_obj%evhigh)
      call f_free_ptr(foe_obj%bisection_shift)
      call f_free_ptr(foe_obj%charge)


      call f_release_routine()

    end subroutine inverse_chebyshev_expansion_new


    !!subroutine adjust_eval_bounds
    !!  implicit none

    !!  ! Calling arguments
    !!  logical,dimension(2),intent(in) :: eval_bounds_ok
    !!  type(foe_data),intent(inout) :: foe_obj

    !!  if (.not. eval_bounds_ok(1)) then
    !!      ! Lower bounds too large

    !!end subroutine adjust_eval_bounds



    subroutine find_fermi_level(iproc, nproc, npl, chebyshev_polynomials, &
               foe_verbosity, label, smatl, foe_obj, kernel_)
      use module_base
      use yaml_output
      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                                   SPARSE_FULL, SPARSE_MATMUL_SMALL, &
                                   SPARSE_MATMUL_LARGE, SPARSEMM_SEQ, SPARSE_TASKGROUP, &
                                   matrices, sparse_matrix
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, compress_matrix_distributed_wrapper, &
                              trace_sparse
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: scale_and_shift_matrix, evnoise, &
                            check_eigenvalue_spectrum_new, retransform_ext, get_chebyshev_expansion_coefficients
      use module_func
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npl
      type(sparse_matrix),intent(in) :: smatl
      real(kind=8),dimension(smatl%smmm%nvctrp_mm,npl) :: chebyshev_polynomials
      integer,intent(in) :: foe_verbosity
      character(len=*),intent(in) :: label
      type(foe_data),intent(inout) :: foe_obj
      type(matrices),intent(inout) :: kernel_
    
      ! Local variables
      integer :: jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
      integer,parameter :: nplx=50000
      real(kind=8),dimension(:,:,:),allocatable :: cc, cc_check
      real(kind=8),dimension(:,:),allocatable :: fermip_check
      real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=8) :: anoise, scale_factor, shift_value, sumn, sumn_check, charge_diff, ef_interpol, ddot
      real(kind=8) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
      real(kind=8) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=8) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible
      logical,dimension(2) :: emergency_stop
      real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
      real(kind=8),dimension(:),allocatable :: hamscal_compr, fermi_check_compr
      real(kind=8),dimension(4,4) :: interpol_matrix
      real(kind=8),dimension(4) :: interpol_vector
      real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
      real(kind=8) :: temp_multiplicator, ebs_check, ef, ebsp
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
      real(kind=8) :: degree_multiplicator, x_max_error, max_error, x_max_error_check, max_error_check
      real(kind=8) :: mean_error, mean_error_check
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      real(kind=8),dimension(2) :: temparr
      real(kind=8),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=8),dimension(:),allocatable :: fermi_new, fermi_check_new, fermi_small_new
      integer :: iline, icolumn, icalc
      
    
    
      call f_routine(id='find_fermi_level')

    
      if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')
    
    
      call timing(iproc, 'FOE_auxiliary ', 'ON')
    
    
      evbounds_shrinked=.false.
    
    
      !!penalty_ev = f_malloc((/smatl%nfvctr,smatl%smmm%nfvctrp,2/),id='penalty_ev')
      !!fermip_check = f_malloc((/smatl%nfvctr,smatl%smmm%nfvctrp/),id='fermip_check')
    
      !penalty_ev_new = f_malloc((/smatl%smmm%nvctrp,2/),id='penalty_ev_new')
      !fermi_new = f_malloc((/smatl%smmm%nvctrp/),id='fermi_new')
      fermi_small_new = f_malloc(max(smatl%smmm%nvctrp_mm,1),id='fermi_small_new')
    
    
      !!call timing(iproc, 'FOE_auxiliary ', 'OF')
      !!if (calculate_minusonehalf) then
      !!    if (iproc==0) call yaml_map('S^-1/2','recalculate')
      !!    call overlap_minus_onehalf() ! has internal timer
      !!else
      !!    if (iproc==0) call yaml_map('S^-1/2','from memory')
      !!end if
      !!call timing(iproc, 'FOE_auxiliary ', 'ON')
    
    
      hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')
    
        
      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm
      
      
      !!! Fake allocation, will be modified later
      !!chebyshev_polynomials = f_malloc((/nsize_polynomial,1/),id='chebyshev_polynomials')
    
    
      !!! try to decrease the eigenvalue spectrum a bit
      !!if (foe_data_get_int(foe_obj,"evbounds_isatur")>foe_data_get_int(foe_obj,"evbounds_nsatur") .and. &
      !!    foe_data_get_int(foe_obj,"evboundsshrink_isatur")<=foe_data_get_int(foe_obj,"evboundsshrink_nsatur")) then
      !!    do ispin=1,smatl%nspin
      !!        call foe_data_set_real(foe_obj,"evlow",0.9d0*foe_data_get_real(foe_obj,"evlow",ispin),ispin)
      !!        call foe_data_set_real(foe_obj,"evhigh",0.9d0*foe_data_get_real(foe_obj,"evhigh",ispin),ispin)
      !!    end do
      !!    evbounds_shrinked=.true.
      !!else
      !!    evbounds_shrinked=.false.
      !!end if
    
      ntemp = NTEMP_ACCURATE
      degree_multiplicator = DEGREE_MULTIPLICATOR_ACCURATE
      temp_multiplicator = TEMP_MULTIPLICATOR_ACCURATE
    
      fscale_new=1.d100
    
    
      spin_loop: do ispin=1,smatl%nspin
    
          !isshift=(ispin-1)*smats%nvctrp_tg
          !imshift=(ispin-1)*smatm%nvctrp_tg
          ilshift=(ispin-1)*smatl%nvctrp_tg
          ilshift2=(ispin-1)*smatl%nvctrp_tg

          !call get_minmax_eigenvalues(iproc, smatm, ham_, imshift, smats, ovrlp_, isshift)
    
          degree_sufficient=.true.
    
          fscale_new = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")
    
          !temp_loop: do itemp=1,ntemp
    
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
                  call foe_data_set_real(foe_obj, &
                       "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",ispin),1.d-4), &
                       ispin)
            
                  efarr(1)=foe_data_get_real(foe_obj,"ef",ispin)-foe_data_get_real(foe_obj,"bisection_shift",ispin)
                  efarr(2)=foe_data_get_real(foe_obj,"ef",ispin)+foe_data_get_real(foe_obj,"bisection_shift",ispin)
    
                  sumnarr(1)=0.d0
                  sumnarr(2)=1.d100
                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
                  call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)
            
                  !adjust_lower_bound=.true.
                  !adjust_upper_bound=.true.
            
                  !calculate_SHS=.true.
            
                  !if (smatl%smmm%nfvctrp>0) then
                  !    call f_zero(smatl%nfvctr*smatl%smmm%nfvctrp*smatl%nspin,kernel_%matrixp(1,1,1))
                  !end if
            
                  if (iproc==0) then
                      !call yaml_sequence(advance='no')
                      if (foe_verbosity>=1) then
                          !!call yaml_sequence_open('FOE to determine density kernel',label=&
                          !!     'it_foe'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//&
                          !!     trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)')))//'-'//&
                          !!     trim(adjustl(yaml_toa(itemp,fmt='(i2.2)')))//'-'//&
                          !!     trim(adjustl(yaml_toa(ispin,fmt='(i2.2)'))))
                          call yaml_sequence_open('FOE to determine density kernel',&
                               label='it_foe'//trim(label)//'-'//&
                               trim(adjustl(yaml_toa(itemp,fmt='(i2.2)')))//'-'//&
                               trim(adjustl(yaml_toa(ispin,fmt='(i2.2)'))))
                      else
                          call yaml_sequence_open('FOE to determine density kernel')
                          if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='-')
                      end if
                  end if
            
            
            
                  it=0
                  eval_bounds_ok=.true.
                  bisection_bounds_ok=.false.
                  main_loop: do 
                      
                      it=it+1
            
                      if (iproc==0) then
                          call yaml_newline()
                          call yaml_sequence(advance='no')
                          call yaml_mapping_open(flow=.true.)
                          if (foe_verbosity>=1) call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
                      end if
                  
            
                      !! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
                      !if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      !    foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                      !    !!call scale_and_shift_hamiltonian()
                      !    call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
                      !         smatm, ham_, imshift, &
                      !         smat2=smats, mat2=ovrlp_, i2shift=isshift, &
                      !         matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
                      !    calculate_SHS=.true.
                      !else
                      !    calculate_SHS=.false.
                      !end if
                      !evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
                      !evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)
            
            
                      !!! Determine the degree of the polynomial
                      !!npl=nint(degree_multiplicator* &
                      !!    (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))/fscale)
                      !!npl=max(npl,NPL_MIN)
                      !!!npl_check = nint(degree_multiplicator*(foe_data_get_real(foe_obj,"evhigh")-foe_data_get_real(foe_obj,"evlow"))/fscale_check)
                      !!!npl_check = max(npl_check,nint(real(npl,kind=8)/CHECK_RATIO)) ! this is necessary if npl was set to the minimal value
                      !!npl_check = nint(real(npl,kind=8)/CHECK_RATIO)
                      !!npl_boundaries = nint(degree_multiplicator* &
                      !!    (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)) &
                      !!        /foe_data_get_real(foe_obj,"fscale_lowerbound")) ! max polynomial degree for given eigenvalue boundaries
                      !!if (npl>npl_boundaries) then
                      !!    npl=npl_boundaries
                      !!    if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
                      !!    if (iproc==0) write(*,*) 'STOP SINCE THIS WILL CREATE PROBLEMS WITH NPL_CHECK'
                      !!    stop
                      !!end if
                      !!if (npl>nplx) stop 'npl>nplx'
            
                      !!! Array the holds the Chebyshev polynomials. Needs to be recalculated
                      !!! every time the Hamiltonian has been modified.
                      !!if (calculate_SHS) then
                      !!    call f_free(chebyshev_polynomials)
                      !!    chebyshev_polynomials = f_malloc((/nsize_polynomial,npl/),id='chebyshev_polynomials')
                      !!end if
            
                      !if (foe_verbosity>=1 .and. iproc==0) then
                      if (iproc==0) then
                          if (foe_verbosity>=1) then
                              call yaml_map('bisec/eval bounds',&
                                   (/fermilevel_get_real(f,"efarr(1)"),fermilevel_get_real(f,"efarr(2)"),&
                                   foe_data_get_real(foe_obj,"evlow",ispin), &
                                   foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          else
                              call yaml_map('eval bounds',&
                                   (/foe_data_get_real(foe_obj,"evlow",ispin), &
                                   foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          end if
                          call yaml_map('pol deg',npl,fmt='(i0)')
                          if (foe_verbosity>=1) call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es16.9)')
                      end if
            
            
                      cc = f_malloc((/npl,3,1/),id='cc')
                      !!cc_check = f_malloc((/npl_check,3,1/),id='cc_check')
            
                      !!if (foe_data_get_real(foe_obj,"evlow",ispin)>=0.d0) then
                      !!    call f_err_throw('Lowest eigenvalue must be negative')
                      !!end if
                      !!if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      !!    call f_err_throw('Highest eigenvalue must be positive')
                      !!end if
            
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call timing(iproc, 'chebyshev_coef', 'ON')
            
                      !if (foe_data_get_real(foe_obj,"tmprtr")/=0.d0) call f_err_throw('tmprtr must be zero')
                      call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale)
                      call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,1), &
                           x_max_error, max_error, mean_error)
                      !!##call chebft(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,1,1), &
                      !!##     foe_data_get_real(foe_obj,"ef",ispin), fscale, foe_data_get_real(foe_obj,"tmprtr"), &
                      !!##     x_max_error, max_error, mean_error)
                      !!call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!     foe_data_get_real(foe_obj,"evhigh",ispin), cc(1,1,1), cc(1,2,1), npl)
                      call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                           muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,2,1), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      do ipl=1,npl
                         cc(ipl,3,1) = -cc(ipl,2,1)
                      end do
                      !!##call chebyshev_coefficients_penalyfunction(foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1,2,1), max_error_fake)
                      call evnoise(npl, cc(1,2,1), foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
        
                      !!call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale_check)
                      !!call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!     foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, func, cc_check(1,1,1), &
                      !!     x_max_error_check, max_error_check, mean_error_check)
                      !!##call chebft(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,1,1), &
                      !!##     foe_data_get_real(foe_obj,"ef",ispin), fscale_check, foe_data_get_real(foe_obj,"tmprtr"), &
                      !!##     x_max_error_check, max_error_check, mean_error_check)
                      !call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
                      !     foe_data_get_real(foe_obj,"evhigh",ispin), &
                      !     cc_check(1,1,1), cc_check(1,2,1), npl_check)
                      !!call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                      !!     muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      !!call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!     foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, func, cc_check(1,2,1), &
                      !!     x_max_error_check, max_error_check, mean_error_check)
                      !!do ipl=1,npl_check
                      !!   cc_check(ipl,3,1) = -cc_check(ipl,2,1)
                      !!end do
                      !!##call chebyshev_coefficients_penalyfunction(foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl_check, cc_check(1,2,1), max_error_fake)

                      if (iproc==0 .and. foe_verbosity>=1) then
                          call yaml_newline()
                          call yaml_mapping_open('accuracy (x, max err, mean err)')
                          call yaml_map('main',(/x_max_error,max_error,mean_error/),fmt='(es9.2)')
                          call yaml_map('check',(/x_max_error_check,max_error_check,max_error/),fmt='(es9.2)')
                          call yaml_mapping_close()
                          call yaml_newline()
                      end if
            
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
    
                      if (smatl%nspin==1) then
                          do ipl=1,npl
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,2,1)=2.d0*cc(ipl,2,1)
                              cc(ipl,3,1)=2.d0*cc(ipl,3,1)
                          end do
                          !!do ipl=1,npl_check
                          !!    cc_check(ipl,1,1)=2.d0*cc_check(ipl,1,1)
                          !!    cc_check(ipl,2,1)=2.d0*cc_check(ipl,2,1)
                          !!    cc_check(ipl,3,1)=2.d0*cc_check(ipl,3,1)
                          !!end do
                      end if
                    
                    
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
            
                      emergency_stop=.false.
                      !!if (calculate_SHS) then
                      !!    ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                      !!    if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
                      !!    !!write(*,*) 'calling chebyshev_clean, iproc', iproc
                      !!    call chebyshev_clean(iproc, nproc, npl, cc, &
                      !!         smatl, hamscal_compr, &
                      !!         calculate_SHS, &
                      !!         nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                      !!         emergency_stop, invovrlp_compr=ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:))
                      !!     
                      !!    !!write(*,*) 'before mpi_barrier, iproc', iproc
                      !!    !!call mpi_barrier(bigdft_mpi%mpi_comm, ipl)
                      !!    !!write(*,*) 'after chebyshev_clean, iproc', iproc
                      !!    call transform_sparsity_pattern(smatl%nfvctr, &
                      !!         smatl%smmm%nvctrp_mm, smatl%smmm%isvctr_mm, &
                      !!         smatl%nseg, smatl%keyv, smatl%keyg, smatl%smmm%line_and_column_mm, &
                      !!         smatl%smmm%nvctrp, smatl%smmm%isvctr, &
                      !!         smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, &
                      !!         smatl%smmm%istsegline, 'large_to_small', fermi_small_new, fermi_new)
                      !!    !!write(*,*) 'after transform_sparsity_pattern, iproc', iproc
    
    
                      !!    !!do i=1,smatl%smmm%nvctrp
                      !!    !!    ii = smatl%smmm%isvctr + i
                      !!    !!    call get_line_and_column(ii, smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, iline, icolumn)
                      !!    !!!!    kernel_%matrixp(icolumn,iline-smatl%smmm%isfvctr,1) = fermi_new(i)
                      !!    !!    penalty_ev(icolumn,iline-smatl%smmm%isfvctr,1) = penalty_ev_new(i,1)
                      !!    !!    penalty_ev(icolumn,iline-smatl%smmm%isfvctr,2) = penalty_ev_new(i,2)
                      !!    !!end do
    
                      !!else
                          ! The Chebyshev polynomials are already available
                          if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                               smatl%nfvctr, smatl%smmm%nfvctrp, &
                              smatl, chebyshev_polynomials, 1, cc, fermi_small_new)
                          !!call calculate_trace_distributed_new(fermi_new, sumn)
                          !!write(*,*) 'trace debug', sumn
    
                          !!call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                          !!     smatl, fermi_small_new, kernel_%matrixp(:,:,1))
                          !!!!call calculate_trace_distributed(kernel_%matrixp, sumn)
                          !!write(*,'(a,2es16.8)') 'sum(fermi_new), sum(kernel_%matrix(:,:,1)', sum(abs(fermi_new)), sum(abs(kernel_%matrixp(:,:,1)))
                      !!end if 
        
        
        
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
                      !!if (calculate_SHS) then
                      !!    !!call check_eigenvalue_spectrum()
                      !!    !!call check_eigenvalue_spectrum(nproc, smatl, smats, ovrlp_, ispin, &
                      !!    !!      isshift, 1.2d0, 1.2d0, penalty_ev, anoise, .true., emergency_stop, &
                      !!    !!      foe_obj, restart, eval_bounds_ok)
                      !!    call check_eigenvalue_spectrum_new(nproc, smatl, ispin, &
                      !!          isshift, 1.2d0, 1.2d0, penalty_ev_new, anoise, .true., emergency_stop, &
                      !!          foe_obj, restart, eval_bounds_ok, smat_s=smats, mat=ovrlp_)
                      !!end if
            
                      call f_free(cc)
            
                      !!if (restart) then
                      !!    if(evbounds_shrinked) then
                      !!        ! this shrink was not good, increase the saturation counter
                      !!        call foe_data_set_int(foe_obj,"evboundsshrink_isatur", &
                      !!             foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
                      !!    end if
                      !!    call foe_data_set_int(foe_obj,"evbounds_isatur",0)
                      !!    if (iproc==0) then
                      !!        if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                      !!             (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      !!        call yaml_mapping_close()
                      !!        !call bigdft_utils_flush(unit=6)
                      !!    end if
                      !!    call f_free(cc_check)
                      !!    cycle
                      !!end if
                          
                      !!! eigenvalue bounds ok
                      !!if (calculate_SHS) then
                      !!    call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
                      !!end if
                    
                      !call calculate_trace_distributed(kernel_%matrixp, sumn)
                      call calculate_trace_distributed_new(iproc, nproc, smatl, fermi_small_new, sumn)
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
                          !!call f_free(cc_check)
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
                          !call chebyshev_fast(iproc, nproc, nsize_polynomial, npl_check, &
                          !     smatl%nfvctr, smatl%smmm%nfvctrp, &
                          !     smatl, chebyshev_polynomials, 1, cc_check, fermi_check_new)
                          !!call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                          !!     smatl, fermi_check_new, fermip_check)
                          !call f_free(cc_check)
                          diff=0.d0
                          !do iorb=1,smatl%smmm%nfvctrp
                          !    do jorb=1,smatl%nfvctr
                          !        !SM: need to fix the spin here
                          !        diff = diff + (kernel_%matrixp(jorb,iorb,1)-fermip_check(jorb,iorb))**2
                          !    end do
                          !end do
                          !do i=1,smatl%smmm%nvctrp_mm
                          !    diff = diff + (fermi_small_new(i)-fermi_check_new(i))**2
                          !end do
        
                          if (nproc > 1) then
                              call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
                          end if
        
                          diff=sqrt(diff)
                          if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')
                          exit
                      end if
        
                      !call f_free(cc_check)
            
            
                  end do main_loop
            
            
             write(*,*) 'before compress_matrix_distributed_wrapper'
        
             call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                  fermi_small_new, &
                  kernel_%matrix_compr(ilshift+1:))

             write(*,*) 'after compress_matrix_distributed_wrapper'
             !!call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
             !!     fermi_check_new, fermi_check_compr)
        
        
            
            
              !!! Calculate S^-1/2 * K * S^-1/2^T
              !!! Since S^-1/2 is symmetric, don't use the transpose
              !!istl = smatl%smmm%istartend_mm_dj(1)-smatl%isvctrp_tg
              !!!write(*,*) 'before kernel_%matrix_compr(ilshift+istl)',iproc, kernel_%matrix_compr(ilshift+istl)
              !!call retransform_ext(iproc, nproc, smatl, &
              !!     ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:), kernel_%matrix_compr(ilshift+1:))
              !!!write(*,*) 'after kernel_%matrix_compr(ilshift+istl)',iproc, kernel_%matrix_compr(ilshift+istl)
    
        
              !!call retransform_ext(iproc, nproc, smatl, &
              !!     ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:), fermi_check_compr)
        
              !!call calculate_trace_distributed_new(fermi_check_new, sumn_check)
    
              !@NEW ##########################
              !!sumn = trace_sparse(iproc, nproc, smats, smatl, &
              !!       ovrlp_%matrix_compr(isshift+1:), &
              !!       kernel_%matrix_compr(ilshift+1:), ispin)
              !!sumn_check = trace_sparse(iproc, nproc, smats, smatl, &
              !!             ovrlp_%matrix_compr(isshift+1:), &
              !!             fermi_check_compr, ispin)
              !@ENDNEW #######################
            
        
            !!  ! Calculate trace(KH). Since they have the same sparsity pattern and K is
            !!  ! symmetric, this is a simple ddot.
            !!  !write(*,*) 'iproc, smatl%smmm%istartend_mm_dj', iproc, smatl%smmm%istartend_mm_dj
            !!  ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
            !!  istl = smatl%smmm%istartend_mm_dj(1)-smatl%isvctrp_tg
            !!  !write(*,*) 'ddot kernel_%matrix_compr(ilshift+istl)', &
            !!  !    iproc, kernel_%matrix_compr(ilshift+istl), ilshift+istl, hamscal_compr(istl)
            !!  ebsp = ddot(ncount, kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)
            !!  !write(*,*) 'iproc, ncount, ebsp', iproc, ncount, ebsp
            !!  !!write(*,'(a,3i8,3es16.8)') 'iproc, ncount, istl, sum(k), sum(h), ebsp', &
            !!  !!    iproc, ncount, istl, sum(kernel_%matrix_compr(ilshift+istl:)), sum(hamscal_compr(istl:)), ebsp
    
            !!  ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
            !!  istl = smatl%smmm%istartend_mm_dj(1)
            !!  ebs_check = ddot(ncount, fermi_check_compr(istl-smatl%isvctrp_tg), 1, &
            !!              hamscal_compr(istl-smatl%isvctrp_tg), 1)
    
            !!  temparr(1) = ebsp
            !!  temparr(2) = ebs_check
            !!  if (nproc>1) then
            !!      call mpiallred(temparr, mpi_sum, comm=bigdft_mpi%mpi_comm)
            !!  end if
            !!  ebsp = temparr(1)
            !!  ebs_check = temparr(2)
    
            !!  !!write(*,'(a,i6,4es16.8)') 'iproc, ebsp, scale_factor, shift_value, sumn', iproc, ebsp, scale_factor, shift_value, sumn
            !!  ebsp=ebsp/scale_factor+shift_value*sumn
            !!  ebs_check=ebs_check/scale_factor+shift_value*sumn_check
            !!  diff=abs(ebs_check-ebsp)
            !!  diff=diff/abs(ebsp)
        
            !!  if (iproc==0) then
            !!      call yaml_map('ebs',ebsp)
            !!      call yaml_map('ebs_check',ebs_check)
            !!      call yaml_map('diff',ebs_check-ebsp)
            !!      call yaml_map('relative diff',diff)
            !!  end if
        
            !!  if (diff<5.d-5) then
            !!      ! can decrease polynomial degree
            !!      !!call foe_data_set_real(foe_obj,"fscale", 1.25d0*foe_data_get_real(foe_obj,"fscale"))
            !!      if (iproc==0) call yaml_map('modify fscale','increase')
            !!      !fscale_new=min(fscale_new,1.25d0*foe_data_get_real(foe_obj,"fscale"))
            !!      fscale_new=1.25d0*fscale_new
            !!      degree_sufficient=.true.
            !!  else if (diff>=5.d-5 .and. diff < 1.d-4) then
            !!      ! polynomial degree seems to be appropriate
            !!      degree_sufficient=.true.
            !!      if (iproc==0) call yaml_map('modify fscale','No')
            !!      !fscale_new=min(fscale_new,foe_data_get_real(foe_obj,"fscale"))
            !!      fscale_new=fscale_new
            !!  else
            !!      ! polynomial degree too small, increase and recalculate
            !!      ! the kernel
            !!      degree_sufficient=.false.
            !!      !!call foe_data_set_real(foe_obj,"fscale", 0.5*foe_data_get_real(foe_obj,"fscale"))
            !!      if (iproc==0) call yaml_map('modify fscale','decrease')
            !!      !fscale_new=min(fscale_new,0.5d0*foe_data_get_real(foe_obj,"fscale"))
            !!      fscale_new=0.5d0*fscale_new
            !!  end if
            !!  !if (foe_data_get_real(foe_obj,"fscale")<foe_data_get_real(foe_obj,"fscale_lowerbound")) then
            !!  if (fscale_new<foe_data_get_real(foe_obj,"fscale_lowerbound")) then
            !!      !call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_lowerbound"))
            !!      fscale_new=foe_data_get_real(foe_obj,"fscale_lowerbound")
            !!      if (iproc==0) call yaml_map('fscale reached lower limit; reset to', &
            !!          foe_data_get_real(foe_obj,"fscale_lowerbound"))
            !!      reached_limit=.true.
            !!  !else if (foe_data_get_real(foe_obj,"fscale")>foe_data_get_real(foe_obj,"fscale_upperbound")) then
            !!  else if (fscale_new>foe_data_get_real(foe_obj,"fscale_upperbound")) then
            !!      !call foe_data_set_real(foe_obj,"fscale",foe_data_get_real(foe_obj,"fscale_upperbound"))
            !!      fscale_new=foe_data_get_real(foe_obj,"fscale_upperbound")
            !!      if (iproc==0) call yaml_map('fscale reached upper limit; reset to', &
            !!          foe_data_get_real(foe_obj,"fscale_upperbound"))
            !!      reached_limit=.true.
            !!  else
            !!      reached_limit=.false.
            !!  end if
            
        
            
          
              ! Purify the kernel
              !tmb%can_use_transposed=.true.
        
              !!if (.not.purification_quickreturn) then
              !!    call f_err_throw('calling purify_kernel from FOE is deprecated',err_name='BIGDFT_RUNTIME_ERROR')
              !!    !if (iproc==0) then
              !!    !    call yaml_sequence_open('Final kernel purification')
              !!    !    call yaml_newline()
              !!    !end if
              !!    !overlap_calculated=.true.
              !!    !if (itemp==ntemp) then
              !!    !    it_shift=20
              !!    !else
              !!    !    it_shift=1
              !!    !end if
              !!    !call purify_kernel(iproc, nproc, tmb, overlap_calculated, it_shift, 50, &
              !!    !     order_taylor, max_inversion_error, purification_quickreturn, ispin)
              !!    !if (iproc==0) then
              !!    !    call yaml_sequence_close()
              !!    !end if
              !!end if
            
            
           !!   ! Calculate trace(KS).
           !!   sumn = trace_sparse(iproc, nproc, smats, smatl, &
           !!          ovrlp_%matrix_compr(isshift+1:), &
           !!          kernel_%matrix_compr(ilshift+1:), ispin)
    
    
           !!   ! Recalculate trace(KH) (needed since the kernel was modified in the above purification). 
           !!   ! If no purification is done, this should not be necessary.
           !!   ! Since K and H have the same sparsity pattern and K is
           !!   ! symmetric, the trace is a simple ddot.
           !!   ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
           !!   istl = smatl%smmm%istartend_mm_dj(1) - smatl%isvctrp_tg
           !!   ebsp = ddot(ncount, kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)
           !!   if (nproc>1) then
           !!       call mpiallred(ebsp, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
           !!   end if
           !!   ebsp=ebsp/scale_factor+shift_value*sumn
        
        
           !!   if (iproc==0) call yaml_map('trace(KS)',sumn)
        
        
              !if (iproc==0) then
              !    call yaml_map('need to repeat with sharper decay (new)',.not.degree_sufficient)
              !end if
              !if (degree_sufficient) exit temp_loop
              !if (reached_limit) then
              !    if (iproc==0) call yaml_map('limit reached, exit loop',.true.)
              !    exit temp_loop
              !end if
        
        
            
        
          !end do temp_loop
    
          ! Sum up the band structure energy
       !!   ebs = ebs + ebsp
    
      end do spin_loop
    
    
      call foe_data_set_real(foe_obj,"fscale",fscale_new)
    
      degree_sufficient=.true.
    
      if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')
    
    
      !!call f_free(chebyshev_polynomials)
      !!call f_free(penalty_ev)
      !!call f_free(hamscal_compr)
      !!call f_free(fermip_check)
      !!call f_free(fermi_check_compr)
    
      !call f_free(penalty_ev_new)
      !call f_free(fermi_check_new)
      !call f_free(fermi_new)
      call f_free(fermi_small_new)
    
      call timing(iproc, 'FOE_auxiliary ', 'OF')
    
      call f_release_routine()
    
    
    
  !!        contains
  !!  
  !!          subroutine overlap_minus_onehalf()
  !!            use sparsematrix_base, only: sparsematrix_malloc, SPARSE_FULL
  !!            use sparsematrix, only: extract_taskgroup_inplace
  !!            use matrix_operations, only: overlapPowerGeneral, check_taylor_order
  !!            implicit none
  !!            real(kind=8) :: max_error, mean_error
  !!            integer :: i, j, ii
  !!            real(kind=8),dimension(:),allocatable :: tmparr
  !!  
  !!            call f_routine(id='overlap_minus_onehalf')
  !!  
  !!            ! Taylor approximation of S^-1/2 up to higher order
  !!            if (imode==DENSE) then
  !!                stop 'overlap_minus_onehalf: DENSE is deprecated'
  !!            end if
  !!            if (imode==SPARSE) then
  !!                call overlapPowerGeneral(iproc, nproc, order_taylor, 1, (/-2/), -1, &
  !!                     imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
  !!                     ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_minus_one_half_, &
  !!                     check_accur=.true., max_error=max_error, mean_error=mean_error)
  !!            end if
  !!            call check_taylor_order(mean_error, max_inversion_error, order_taylor)
  !!  
  !!            call f_release_routine()
  !!        end subroutine overlap_minus_onehalf
  !!  
  !!  
  !!  
  !!        subroutine retransform(matrix_compr)
  !!            use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
  !!                                    compress_matrix_distributed_wrapper, &
  !!                                    sparsemm_new
  !!            ! Calling arguments
  !!            real(kind=8),dimension(smatl%nvctrp_tg),intent(inout) :: matrix_compr
  !!  
  !!            ! Local variables
  !!            real(kind=8),dimension(:,:),pointer :: inv_ovrlpp, tempp
  !!            real(kind=8),dimension(:),pointer :: inv_ovrlpp_new, tempp_new
  !!            real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
  !!            integer,dimension(:,:,:),allocatable :: istindexarr
  !!            integer :: nout, nseq
  !!  
  !!            call f_routine(id='retransform')
  !!  
  !!            !!inv_ovrlpp = sparsematrix_malloc_ptr(smatl, iaction=DENSE_MATMUL, id='inv_ovrlpp')
  !!            inv_ovrlpp_new = f_malloc_ptr(smatl%smmm%nvctrp, id='inv_ovrlpp_new')
  !!            !!tempp = sparsematrix_malloc_ptr(smatl, iaction=DENSE_MATMUL, id='tmpp')
  !!            tempp_new = f_malloc_ptr(smatl%smmm%nvctrp, id='tempp_new')
  !!            inv_ovrlp_compr_seq = sparsematrix_malloc(smatl, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
  !!            kernel_compr_seq = sparsematrix_malloc(smatl, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
  !!            call sequential_acces_matrix_fast2(smatl, matrix_compr, kernel_compr_seq)
  !!            call sequential_acces_matrix_fast2(smatl, &
  !!                 ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:), inv_ovrlp_compr_seq)
  !!            !!call uncompress_matrix_distributed2(iproc, smatl, DENSE_MATMUL, &
  !!            !!     ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:), inv_ovrlpp)
  !!            !! write(*,*) 'sum(matrix_compr) 0', iproc, sum(ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:))
  !!            !!  write(*,*) 'smatl%nvctrp, smatl%smmm%nvctrp_mm', smatl%nvctrp, smatl%smmm%nvctrp_mm
  !!            !!  write(*,*) 'smatl%isvctr, smatl%smmm%isvctr_mm', smatl%isvctr, smatl%smmm%isvctr_mm
  !!            call transform_sparsity_pattern(smatl%nfvctr, smatl%smmm%nvctrp_mm, smatl%smmm%isvctr_mm, &
  !!                 smatl%nseg, smatl%keyv, smatl%keyg, smatl%smmm%line_and_column_mm, &
  !!                 smatl%smmm%nvctrp, smatl%smmm%isvctr, &
  !!                 smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, &
  !!                 smatl%smmm%istsegline, 'small_to_large', &
  !!                 ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+smatl%smmm%isvctr_mm-smatl%isvctrp_tg+1:), &
  !!                 inv_ovrlpp_new)
  !!            !!  write(*,*) 'sum(matrix_compr) 1', iproc, sum(ovrlp_minus_one_half_(1)%matrix_compr(ilshift2+1:))
  !!            !!  write(*,*) 'sum(inv_ovrlpp_new) 1', iproc, sum(inv_ovrlpp_new)
  !!            !!  write(*,*) 'sum(inv_ovrlpp) 1', iproc, sum(inv_ovrlpp)
  !!  
  !!  
  !!              !!!!          call transform_sparsity_pattern(smatl%nfvctr, smatl%smmm%nvctrp_mm, smatl%smmm%isvctr_mm, &
  !!              !!!!               smatl%nseg, smatl%keyv, smatl%keyg, &
  !!              !!!!               smatl%smmm%nvctrp, smatl%smmm%isvctr, &
  !!              !!!!               smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, &
  !!              !!!!               fermi_new, fermi_small_new)
  !!  
  !!  
  !!  
  !!            !!call f_zero(tempp)
  !!            !!call sparsemm(smatl, kernel_compr_seq, inv_ovrlpp, tempp)
  !!            !!write(*,*) 'sum(tempp) 2',iproc, sum(tempp)
  !!            call sparsemm_new(smatl, kernel_compr_seq, inv_ovrlpp_new, tempp_new)
  !!            !!write(*,*) 'sum(tempp_new) 2',iproc, sum(tempp_new)
  !!            !!inv_ovrlpp=0.d0
  !!            !!call sparsemm(smatl, inv_ovrlp_compr_seq, tempp, inv_ovrlpp)
  !!            call sparsemm_new(smatl, inv_ovrlp_compr_seq, tempp_new, inv_ovrlpp_new)
  !!             !!write(*,*) 'sum(inv_ovrlpp) 3', iproc, sum(inv_ovrlpp)
  !!             !!write(*,*) 'sum(inv_ovrlpp_new) 3', iproc, sum(inv_ovrlpp_new)
  !!  
  !!            !!call f_zero(matrix_compr)
  !!            !!call compress_matrix_distributed(iproc, nproc, smatl, DENSE_MATMUL, &
  !!            !!     inv_ovrlpp, matrix_compr)
  !!            !!  write(*,*) 'sum(matrix_compr) old', iproc, sum(matrix_compr)
  !!            call f_zero(matrix_compr)
  !!            call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_LARGE, &
  !!                 inv_ovrlpp_new, matrix_compr)
  !!              !!write(*,*) 'sum(matrix_compr) new', iproc, sum(matrix_compr)
  !!  
  !!            !!call f_free_ptr(inv_ovrlpp)
  !!            !!call f_free_ptr(tempp)
  !!            call f_free_ptr(inv_ovrlpp_new)
  !!            call f_free_ptr(tempp_new)
  !!            call f_free(inv_ovrlp_compr_seq)
  !!            call f_free(kernel_compr_seq)
  !!  
  !!            call f_release_routine()
  !!  
  !!        end subroutine retransform
  !!  
  !!  
  !!        subroutine calculate_trace_distributed_new(matrixp, trace)
  !!            real(kind=8),dimension(smatl%smmm%nvctrp_mm),intent(in) :: matrixp
  !!            real(kind=8),intent(out) :: trace
  !!            integer :: i, ii
  !!  
  !!            call f_routine(id='calculate_trace_distributed_new')
  !!  
  !!            trace = 0.d0
  !!            !$omp parallel default(none) &
  !!            !$omp shared(trace, smatl, matrixp) &
  !!            !$omp private(i, iline, icolumn)
  !!            !$omp do reduction(+:trace)
  !!            do i=1,smatl%smmm%nvctrp_mm
  !!                iline = smatl%smmm%line_and_column_mm(1,i)
  !!                icolumn = smatl%smmm%line_and_column_mm(2,i)
  !!                if (iline==icolumn) then
  !!                    trace = trace + matrixp(i)
  !!                end if
  !!            end do
  !!            !$omp end do
  !!            !$omp end parallel
  !!  
  !!            if (nproc > 1) then
  !!                call mpiallred(trace, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
  !!            end if
  !!  
  !!            call f_release_routine()
  !!        end subroutine calculate_trace_distributed_new

    
    
    end subroutine find_fermi_level


    subroutine calculate_trace_distributed_new(iproc, nproc, smatl, matrixp, trace)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smatl
      real(kind=8),dimension(smatl%smmm%nvctrp_mm),intent(in) :: matrixp
      real(kind=8),intent(out) :: trace
      integer :: i, ii, iline, icolumn

      call f_routine(id='calculate_trace_distributed_new')

      trace = 0.d0
      !$omp parallel default(none) &
      !$omp shared(trace, smatl, matrixp) &
      !$omp private(i, iline, icolumn)
      !$omp do reduction(+:trace)
      do i=1,smatl%smmm%nvctrp_mm
          iline = smatl%smmm%line_and_column_mm(1,i)
          icolumn = smatl%smmm%line_and_column_mm(2,i)
          if (iline==icolumn) then
              trace = trace + matrixp(i)
          end if
      end do
      !$omp end do
      !$omp end parallel

      if (nproc > 1) then
          call mpiallred(trace, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      call f_release_routine()
    end subroutine calculate_trace_distributed_new


end module ice
