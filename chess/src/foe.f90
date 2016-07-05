module foe
  use sparsematrix_base
  implicit none

  private

  !> Public routines
  public :: fermi_operator_expansion_new
  public :: get_selected_eigenvalues

  contains



    subroutine fermi_operator_expansion_new(iproc, nproc, comm, &
               ebs, &
               calculate_minusonehalf, foe_verbosity, &
               smats, smatm, smatl, ham_, ovrlp_, ovrlp_minus_one_half_, kernel_, foe_obj, ice_obj, &
               symmetrize_kernel, calculate_energy_density_kernel, energy_kernel_)
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, compress_matrix_distributed_wrapper, &
                              trace_sparse_matrix_product, symmetrize_matrix, max_asymmetry_of_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: evnoise, &
                            retransform_ext, get_chebyshev_expansion_coefficients, &
                            find_fermi_level, get_polynomial_degree, &
                            calculate_trace_distributed_new, get_bounds_and_polynomials
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      real(kind=mp),intent(out) :: ebs
      logical,intent(in) :: calculate_minusonehalf, symmetrize_kernel, calculate_energy_density_kernel
      integer,intent(in) :: foe_verbosity
      type(sparse_matrix),intent(in) :: smats, smatm, smatl
      type(matrices),intent(in) :: ham_, ovrlp_
      type(matrices),dimension(1),intent(inout) :: ovrlp_minus_one_half_
      type(matrices),intent(inout) :: kernel_
      type(foe_data),intent(inout) :: foe_obj, ice_obj
      type(matrices),intent(inout),optional :: energy_kernel_

      ! Local variables
      integer :: npl, jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
      integer,parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),pointer :: cc
      real(kind=mp),dimension(:,:,:),allocatable :: cc_check
      real(kind=mp),dimension(:,:,:),pointer :: chebyshev_polynomials
      real(kind=mp),dimension(:,:),allocatable :: fermip_check
      real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp) :: anoise, scale_factor, shift_value, sumn, sumn_check, charge_diff, ef_interpol, ddot
      real(kind=mp) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=mp) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new, fscale_newx, asymm_K
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible
      logical,dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: efarr, sumnarr, allredarr
      real(kind=mp),dimension(:),allocatable :: hamscal_compr, fermi_check_compr, kernel_tmp
      real(kind=mp),dimension(4,4) :: interpol_matrix
      real(kind=mp),dimension(4) :: interpol_vector
      real(kind=mp),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
      real(kind=mp) :: temp_multiplicator, ebs_check, ef, ebsp
      integer :: irow, icol, itemp, iflag,info, ispin, isshift, imshift, ilshift, i, j, itg, ncount, istl, ists
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
      real(kind=mp),parameter :: FSCALE_LOWER_LIMIT=5.d-3
      real(kind=mp),parameter :: FSCALE_UPPER_LIMIT=5.d-2
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
      real(kind=mp),parameter :: CHECK_RATIO=1.25d0
      !integer,parameter :: NPL_MIN=100
      !!type(matrices) :: inv_ovrlp
      integer,parameter :: NTEMP_ACCURATE=4
      integer,parameter :: NTEMP_FAST=1
      real(kind=mp) :: degree_multiplicator, betax, ebsp_allspins
      real(kind=mp),dimension(1) :: x_max_error, max_error, x_max_error_check, max_error_check, mean_error, mean_error_check
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      real(kind=mp),dimension(2) :: temparr
      real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:),allocatable :: fermi_new, fermi_check_new, fermi_small_new
      integer :: iline, icolumn, icalc, npl_min
      real(kind=mp),dimension(:),allocatable :: ham_large
      real(kind=mp),dimension(2) :: fscale_ispin
      real(kind=mp),dimension(1) :: ef_arr, fscale_arr
      real(mp) :: ebs_check_allspins
      real(mp),dimension(:),allocatable :: sumn_allspins
      !!integer,parameter :: NPL_MAX = 10000
      !!integer,parameter :: NPL_STRIDE = 10
      integer :: npl_max, npl_stride



      call f_routine(id='fermi_operator_expansion_new')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')


      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      npl_min = foe_data_get_int(foe_obj,"npl_min")
      npl_max = foe_data_get_int(foe_obj,"npl_max")
      npl_stride = foe_data_get_int(foe_obj,"npl_stride")
      betax = foe_data_get_real(foe_obj,"betax")

      evbounds_shrinked=.false.


      !!penalty_ev = f_malloc((/smatl%nfvctr,smatl%smmm%nfvctrp,2/),id='penalty_ev')
      !!fermip_check = f_malloc((/smatl%nfvctr,smatl%smmm%nfvctrp/),id='fermip_check')
      fermi_check_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='fermi_check_compr')
      kernel_tmp = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='kernel_tmp')
      sumn_allspins = f_malloc(smatl%nspin,id='sumn_allspins')

      fermi_check_new = f_malloc(max(smatl%smmm%nvctrp_mm,1),id='fermip_check_new')
      !!fermi_new = f_malloc((/smatl%smmm%nvctrp/),id='fermi_new')
!      fermi_small_new = f_malloc(max(smatl%smmm%nvctrp_mm,1),id='fermi_small_new')


      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')
      if (iproc==0) call yaml_mapping_open('S^-1/2')
      if (calculate_minusonehalf) then
          !if (iproc==0) call yaml_map('S^-1/2','recalculate')
          !!call overlap_minus_onehalf() ! has internal timer
          if (iproc==0) call yaml_map('Can take from memory',.false.)
          call overlap_minus_onehalf(iproc, nproc, comm, smats, smatl, ovrlp_, ovrlp_minus_one_half_, &
               ice_obj=ice_obj) !has internal timer
      else
          if (iproc==0) call yaml_map('Can take from memory',.true.)
      end if
      if (iproc==0) call yaml_mapping_close()
      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')



      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm


      !! Fake allocation, will be modified later
      !chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,1/),id='chebyshev_polynomials')


      ! try to decrease the eigenvalue spectrum a bit
      if (foe_data_get_int(foe_obj,"evbounds_isatur")>foe_data_get_int(foe_obj,"evbounds_nsatur") .and. &
          foe_data_get_int(foe_obj,"evboundsshrink_isatur")<=foe_data_get_int(foe_obj,"evboundsshrink_nsatur")) then
          do ispin=1,smatl%nspin
              call foe_data_set_real(foe_obj,"evlow",0.9d0*foe_data_get_real(foe_obj,"evlow",ispin),ispin)
              call foe_data_set_real(foe_obj,"evhigh",0.9d0*foe_data_get_real(foe_obj,"evhigh",ispin),ispin)
          end do
          evbounds_shrinked=.true.
      else
          evbounds_shrinked=.false.
      end if

      !write(*,*) 'evlow, evhigh', foe_data_get_real(foe_obj,"evlow",1), foe_data_get_real(foe_obj,"evhigh",1)

      ntemp = NTEMP_ACCURATE
      degree_multiplicator = DEGREE_MULTIPLICATOR_ACCURATE
      temp_multiplicator = TEMP_MULTIPLICATOR_ACCURATE

      fscale_new=1.d100

      ebs=0.d0

       hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

      fscale_newx = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")

      if (iproc==0) then
          call yaml_sequence_open('Kernel calculation')
      end if

      if (smatl%nspin>2) then
          call f_err_throw('smatl%nspin>2')
      end if
      fscale_ispin(1:2) = huge(fscale_ispin(1:2))

      !spin_loop: do ispin=1,smatl%nspin


          fscale_new = fscale_newx
          call foe_data_set_real(foe_obj,"fscale",fscale_new)

      !!    isshift=(ispin-1)*smats%nvctrp_tg
      !!    imshift=(ispin-1)*smatm%nvctrp_tg
      !!    ilshift=(ispin-1)*smatl%nvctrp_tg

      !!do ispin=1,smatl%nspin
      !!    isshift=(ispin-1)*smats%nvctrp_tg
      !!    imshift=(ispin-1)*smatm%nvctrp_tg
      !!    call get_minmax_eigenvalues(iproc, smatm, ham_, imshift, smats, ovrlp_, isshift)
      !!end do

          degree_sufficient=.true.


          npl_min = 10

          temp_loop: do itemp=1,ntemp

              if (iproc==0) then
                  call yaml_sequence(advance='no')
                  call yaml_comment('ispin:'//trim(yaml_toa(1))//', itemp:'//trim(yaml_toa(itemp)),hfill='-')
                  call yaml_newline()
              end if

              fscale = fscale_new
              fscale = max(fscale,FSCALE_LOWER_LIMIT)
              fscale = min(fscale,FSCALE_UPPER_LIMIT)
              fscale_check = CHECK_RATIO*fscale

              evlow_old=1.d100
              evhigh_old=-1.d100

              if (iproc==0) then
                  call yaml_map('decay length of error function',fscale,fmt='(es10.3)')
                  !call yaml_map('decay length multiplicator',temp_multiplicator,fmt='(es10.3)')
                  !call yaml_map('polynomial degree multiplicator',degree_multiplicator,fmt='(es10.3)')
              end if


                  ! Don't let this value become too small.
                  call foe_data_set_real(foe_obj, &
                       "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",1),1.d-4), &
                       1)


                  sumnarr(1)=0.d0
                  sumnarr(2)=1.d100
                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",1), foe_data_get_real(foe_obj,"ef",1), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",1), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)

                  ! Use kernel_%matrix_compr as workarray to save memory
                  efarr(1) = foe_data_get_real(foe_obj,"ef",1)
                  fscale_arr(1) = foe_data_get_real(foe_obj,"fscale",1)
                  call get_bounds_and_polynomials(iproc, nproc, comm, 2, 1, NPL_MAX, NPL_STRIDE, betax, &
                       1, FUNCTION_ERRORFUNCTION, .false., 1.2_mp, 1.2_mp, foe_verbosity, &
                       smatm, smatl, ham_, foe_obj, npl_min, kernel_%matrix_compr(ilshift+1:), &
                       chebyshev_polynomials, npl, scale_factor, shift_value, hamscal_compr, &
                       smats=smats, ovrlp_=ovrlp_, ovrlp_minus_one_half_=ovrlp_minus_one_half_(1), &
                       efarr=efarr, fscale_arr=fscale_arr, max_errorx=max_error)

                  if (iproc==0) then
                      call yaml_mapping_open('summary',flow=.true.)
                      call yaml_map('npl',npl)
                      call yaml_map('bounds', &
                           (/foe_data_get_real(ice_obj,"evlow",1),foe_data_get_real(ice_obj,"evhigh",1)/),fmt='(f6.2)')
                      call yaml_map('exp accur',max_error,fmt='(es8.2)')
                      call yaml_mapping_close()
                  end if

                  call find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
                       foe_verbosity, 'test', smatl, 1, foe_obj, kernel_)

                  npl_check = nint(real(npl,kind=mp)/CHECK_RATIO)
                  cc_check = f_malloc0((/npl_check,1,3/),id='cc_check')
                  call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",1), fscalex=fscale_check)
                  call get_chebyshev_expansion_coefficients(iproc, nproc, comm, &
                       foe_data_get_real(foe_obj,"evlow",1), &
                       foe_data_get_real(foe_obj,"evhigh",1), npl_check, func, cc_check(1,1,1), &
                       x_max_error_check(1), max_error_check(1), mean_error_check(1))
                  if (smatl%nspin==1) then
                      do ipl=1,npl_check
                          cc_check(ipl,1,1)=2.d0*cc_check(ipl,1,1)
                          cc_check(ipl,1,2)=2.d0*cc_check(ipl,1,2)
                          cc_check(ipl,1,3)=2.d0*cc_check(ipl,1,3)
                      end do
                  end if

                  ebsp_allspins = 0.0_mp
                  ebs_check_allspins = 0.0_mp

                  spin_loop: do ispin=1,smatl%nspin


                      fscale_new = fscale_newx
                      call foe_data_set_real(foe_obj,"fscale",fscale_new)

                      isshift=(ispin-1)*smats%nvctrp_tg
                      imshift=(ispin-1)*smatm%nvctrp_tg
                      ilshift=(ispin-1)*smatl%nvctrp_tg

                      !write(*,*) 'start spinloop', trace_sparse_matrix_product(iproc, nproc, comm, smats, smatl, &
                      !       ovrlp_%matrix_compr(isshift+1:), &
                      !       kernel_%matrix_compr(ilshift+1:))

                      call chebyshev_fast(iproc, nproc, nsize_polynomial, npl_check, &
                           smatl%nfvctr, smatl%smmm%nfvctrp, &
                           smatl, chebyshev_polynomials(:,:,ispin), 1, cc_check, fermi_check_new)
                      !!call f_free(cc)

                      call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                           fermi_check_new, fermi_check_compr(ilshift+1:))
                      ! Calculate S^-1/2 * K * S^-1/2^T
                      ! Since S^-1/2 is symmetric, don't use the transpose
                      istl = smatl%smmm%istartend_mm_dj(1)-smatl%isvctrp_tg
                      !write(*,*) 'before kernel_%matrix_compr(ilshift+istl)',iproc, kernel_%matrix_compr(ilshift+istl)
                      !!write(*,*) 'BEFORE RET, ispin, sum(K)', ispin, sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))
                      !!write(*,*) 'BEFORE RET, ispin, sum(S-)', &
                      !!    ispin, sum(ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:ilshift+smatl%nvctr))
                      call retransform_ext(iproc, nproc, smatl, &
                           ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), kernel_%matrix_compr(ilshift+1:))
                      !write(*,*) 'after kernel_%matrix_compr(ilshift+istl)',iproc, kernel_%matrix_compr(ilshift+istl)
                      !!write(*,*) 'AFTER RET, ispin, sum(K)', ispin, sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))

                      call retransform_ext(iproc, nproc, smatl, &
                           ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), fermi_check_compr(ilshift+1:))
                      ! Explicitly symmetrize the kernel, use fermi_check_compr as temporary array
                      call max_asymmetry_of_matrix(iproc, nproc, comm, &
                           smatl, kernel_%matrix_compr, asymm_K, ispinx=ispin)
                      !!write(*,*) 'BEFORE SYM, ispin, sum(K)', ispin, sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))
                      if (symmetrize_kernel) then
                          call f_memcpy(src=kernel_%matrix_compr, dest=kernel_tmp)
                          call symmetrize_matrix(smatl, 'plus', kernel_tmp, kernel_%matrix_compr, ispinx=ispin)
                          call f_memcpy(src=fermi_check_compr, dest=kernel_tmp)
                          call symmetrize_matrix(smatl, 'plus', kernel_tmp, fermi_check_compr, ispinx=ispin)
                      end if
                      !!write(*,*) 'AFTER SYM, ispin, sum(K)', ispin, sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))

                      call calculate_trace_distributed_new(iproc, nproc, comm, smatl, fermi_check_new, sumn_check)

                      !@NEW ##########################
                      sumn = trace_sparse_matrix_product(iproc, nproc, comm, smats, smatl, &
                             ovrlp_%matrix_compr(isshift+1:), &
                             kernel_%matrix_compr(ilshift+1:))
                      !!write(*,*) 'ispin, sumn, sum(K)', ispin, sumn, sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))
                      sumn_check = trace_sparse_matrix_product(iproc, nproc, comm, smats, smatl, &
                                   ovrlp_%matrix_compr(isshift+1:), &
                                   fermi_check_compr(ilshift+1:))
                      !write(*,*) 'sumn, sumn_check', sumn, sumn_check
                      !@ENDNEW #######################


                      ! Calculate trace(KH). Since they have the same sparsity pattern and K is
                      ! symmetric, this is a simple ddot.
                      !write(*,*) 'iproc, smatl%smmm%istartend_mm_dj', iproc, smatl%smmm%istartend_mm_dj
                      ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
                      istl = smatl%smmm%istartend_mm_dj(1)-smatl%isvctrp_tg
                      !write(*,*) 'ddot kernel_%matrix_compr(ilshift+istl)', &
                      !    iproc, kernel_%matrix_compr(ilshift+istl), ilshift+istl, hamscal_compr(istl)
                      !ebsp = ddot(ncount, kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)
                      ebsp = ddot(ncount, kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(ilshift+istl), 1)
                      !!write(*,*) 'ispin, sum(kernel_%matrix_compr), sum(hamscal_compr)', &
                      !!    sum(kernel_%matrix_compr), sum(hamscal_compr)
                      !write(*,*) 'ebsp',ebsp
                      !write(*,*) 'iproc, ncount, ebsp', iproc, ncount, ebsp
                      !!write(*,'(a,3i8,3es16.8)') 'iproc, ncount, istl, sum(k), sum(h), ebsp', &
                      !!    iproc, ncount, istl, sum(kernel_%matrix_compr(ilshift+istl:)), sum(hamscal_compr(istl:)), ebsp

                      ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
                      istl = smatl%smmm%istartend_mm_dj(1) - smatl%isvctrp_tg
                      ebs_check = ddot(ncount, fermi_check_compr(ilshift+istl), 1, &
                                  hamscal_compr(istl), 1)
                      !write(*,*) 'ebs_check',ebs_check

                      temparr(1) = ebsp
                      temparr(2) = ebs_check
                      if (nproc>1) then
                          call mpiallred(temparr, mpi_sum, comm=comm)
                      end if
                      ebsp = temparr(1)
                      ebs_check = temparr(2)

                      !write(*,'(a,i6,5es16.8)') 'iproc, ebsp, scale_factor, shift_value, sumn, sum(hamscal_compr)', &
                      !    iproc, ebsp, scale_factor, shift_value, sumn, sum(hamscal_compr)
                      ebsp=ebsp/scale_factor+shift_value*sumn
                      ebs_check=ebs_check/scale_factor+shift_value*sumn_check
                      !write(*,*) 'ispin, ebsp, ebs_check', ispin, ebsp, ebs_check
                      diff=abs(ebs_check-ebsp)
                      diff=diff/abs(ebsp)

                      !write(*,*) 'ispin, ebsp', ispin, ebsp
                      ebsp_allspins = ebsp_allspins + ebsp
                      ebs_check_allspins = ebs_check_allspins + ebs_check

                      ! Calculate trace(KS).
                      sumn = trace_sparse_matrix_product(iproc, nproc, comm, smats, smatl, &
                             ovrlp_%matrix_compr(isshift+1:), &
                             kernel_%matrix_compr(ilshift+1:))
                      sumn_allspins(ispin) = sumn

                  end do spin_loop

                  call f_free(cc_check)


                  if (iproc==0) then
                      call yaml_map('Asymmetry of kernel',asymm_K,fmt='(es8.2)')
                      call yaml_map('symmetrize_kernel',symmetrize_kernel)
                      call yaml_map('EBS',ebsp_allspins,fmt='(es19.12)')
                      call yaml_map('EBS higher temperature',ebs_check_allspins,fmt='(es19.12)')
                      call yaml_map('difference',ebs_check_allspins-ebsp_allspins,fmt='(es19.12)')
                      call yaml_map('relative difference',diff,fmt='(es19.12)')
                      if (smatl%nspin==1) then
                          call yaml_map('trace(KS)',sumn_allspins(1))
                      else
                          call yaml_map('trace(KS)',sumn_allspins)
                      end if
                  end if

                  if (diff<5.d-5) then
                      ! can decrease polynomial degree
                      !!call foe_data_set_real(foe_obj,"fscale", 1.25d0*foe_data_get_real(foe_obj,"fscale"))
                      if (iproc==0) call yaml_map('modify error function decay length','increase')
                      !fscale_new=min(fscale_new,1.25d0*foe_data_get_real(foe_obj,"fscale"))
                      fscale_new=1.25d0*fscale_new
                      degree_sufficient=.true.
                  else if (diff>=5.d-5 .and. diff < 1.d-4) then
                      ! polynomial degree seems to be appropriate
                      degree_sufficient=.true.
                      if (iproc==0) call yaml_map('modify error function decay length','No')
                      !fscale_new=min(fscale_new,foe_data_get_real(foe_obj,"fscale"))
                      fscale_new=fscale_new
                  else
                      ! polynomial degree too small, increase and recalculate
                      ! the kernel
                      degree_sufficient=.false.
                      !!call foe_data_set_real(foe_obj,"fscale", 0.5*foe_data_get_real(foe_obj,"fscale"))
                      if (iproc==0) call yaml_map('modify error function decay length','decrease')
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

                  call foe_data_set_real(foe_obj,"fscale",fscale_new)


!!$                  !mpimaxdiff might be used (however fermi_small_new is not initialised there
!!$                  diff=0.d0
!!$                  do i=1,smatl%smmm%nvctrp_mm
!!$                      diff = diff + (fermi_small_new(i)-fermi_check_new(i))**2
!!$                  end do
!!$
!!$                  if (nproc > 1) then
!!$                      call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$                  end if
!!$
!!$                  diff=sqrt(diff)
!!$                  !if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')




!!!                  ! Recalculate trace(KH) (needed since the kernel was modified in the above purification).
!!!                  ! If no purification is done, this should not be necessary.
!!!                  ! Since K and H have the same sparsity pattern and K is
!!!                  ! symmetric, the trace is a simple ddot.
!!!                  ncount = smatl%smmm%istartend_mm_dj(2) - smatl%smmm%istartend_mm_dj(1) + 1
!!!                  istl = smatl%smmm%istartend_mm_dj(1) - smatl%isvctrp_tg
!!!                  ebsp = ddot(ncount, kernel_%matrix_compr(ilshift+istl), 1, hamscal_compr(istl), 1)
!!!                  if (nproc>1) then
!!!                      call mpiallred(ebsp, 1, mpi_sum, comm=comm)
!!!                  end if
!!!                  ebsp=ebsp/scale_factor+shift_value*sumn




                  if (iproc==0) then
                      call yaml_map('need to repeat with sharper decay (new)',.not.degree_sufficient)
                  end if
                  if (degree_sufficient) then
                      !!call f_free_ptr(chebyshev_polynomials)
                      exit temp_loop
                  end if
                  if (reached_limit) then
                      if (iproc==0) call yaml_map('limit reached, exit loop',.true.)
                      !!call f_free_ptr(chebyshev_polynomials)
                      exit temp_loop
                  end if

                  call f_free_ptr(chebyshev_polynomials)

                  !# END NEW ################################################

          end do temp_loop


          ! Sum up the band structure energy
          ebs = ebs + ebsp_allspins

          fscale_ispin(1) = fscale_new


          if (calculate_energy_density_kernel) then
              if (.not.present(energy_kernel_)) then
                  call f_err_throw('energy_kernel_ not present',err_name='SPARSEMATRIX_RUNTIME_ERROR')
              end if
              cc_check = f_malloc0((/npl,1,3/),id='cc_check')
              call func_set(FUNCTION_XTIMESERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",1), fscalex=fscale)
              call get_chebyshev_expansion_coefficients(iproc, nproc, comm, &
                   foe_data_get_real(foe_obj,"evlow",1), &
                   foe_data_get_real(foe_obj,"evhigh",1), npl, func, cc_check(1,1,1), &
                   x_max_error_check(1), max_error_check(1), mean_error_check(1))
              if (smatl%nspin==1) then
                  do ipl=1,npl
                      cc_check(ipl,1,1)=2.d0*cc_check(ipl,1,1)
                      cc_check(ipl,1,2)=2.d0*cc_check(ipl,1,2)
                      cc_check(ipl,1,3)=2.d0*cc_check(ipl,1,3)
                  end do
              end if
              call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                   smatl%nfvctr, smatl%smmm%nfvctrp, &
                   smatl, chebyshev_polynomials, 1, cc_check, fermi_check_new)
              !!call f_free(cc)
              call f_free(cc_check)

              call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                   fermi_check_new, energy_kernel_%matrix_compr(ilshift+1:))
              ! Calculate S^-1/2 * K * S^-1/2^T
              ! Since S^-1/2 is symmetric, don't use the transpose
              istl = smatl%smmm%istartend_mm_dj(1)-smatl%isvctrp_tg

              call retransform_ext(iproc, nproc, smatl, &
                   ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), energy_kernel_%matrix_compr(ilshift+1:))

              !ebsp=ebsp/scale_factor+shift_value*sumn
              !energy_kernel_%matrix_compr = energy_kernel_%matrix_compr/scale_factor + shift_value*ovrlp_%matrix_compr
          end if

          call f_free_ptr(chebyshev_polynomials)

      !end do spin_loop

      if (iproc==0) then
          call yaml_sequence_close()
      end if

      !!! This always takes the value for ispin=2... should be improved
      call foe_data_set_real(foe_obj,"fscale",minval(fscale_ispin))

      degree_sufficient=.true.

    
      if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


      call f_free(sumn_allspins)
      call f_free(hamscal_compr)
      call f_free(fermi_check_compr)
      call f_free(kernel_tmp)
      call f_free(fermi_check_new)
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      call f_release_routine()



        !!!  contains

        !!!    subroutine overlap_minus_onehalf()
        !!!      use ice, only: inverse_chebyshev_expansion_new
        !!!      implicit none
        !!!      integer :: i, j, ii
        !!!      real(kind=mp) :: max_error, mean_error
        !!!      real(kind=mp), dimension(1) :: ex
        !!!      real(kind=mp),dimension(:),allocatable :: tmparr
    
        !!!      call f_routine(id='overlap_minus_onehalf')
    
        !!!      ! Can't use the wrapper, since it is at a higher level in the hierarchy (to be improved)
        !!!      ex=-0.5d0
        !!!      call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
        !!!           ovrlp_smat=smats, inv_ovrlp_smat=smatl, ncalc=1, ex=ex, &
        !!!           ovrlp_mat=ovrlp_, inv_ovrlp=ovrlp_minus_one_half_, &
        !!!           npl_auto=.true., ice_objx=ice_obj)


        !!!      call f_release_routine()
        !!!  end subroutine overlap_minus_onehalf



    end subroutine fermi_operator_expansion_new



    subroutine get_minmax_eigenvalues(iproc, ham_smat, ham_mat, imshift, ovrlp_smat, ovrlp_mat, isshift)
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, imshift, isshift
      type(sparse_matrix),intent(in) :: ham_smat, ovrlp_smat
      type(matrices),intent(in) :: ham_mat, ovrlp_mat

      ! Local variables
      integer :: iseg, ii, i, lwork, info, ieval
      real(kind=mp),dimension(:,:,:),allocatable :: tempmat
      real(kind=mp),dimension(:),allocatable :: eval, work
      !!real(mp) :: tt5, tt7

      call f_routine(id='get_minmax_eigenvalues')

      if (ham_smat%nfvctr/=ovrlp_smat%nfvctr) call f_err_throw('ham_smat&nfvctr/=ovrlp_smat%nfvctr')

      tempmat = f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr,2/),id='tempmat')
      do iseg=1,ham_smat%nseg
          ii=ham_smat%keyv(iseg)
          do i=ham_smat%keyg(1,1,iseg),ham_smat%keyg(2,1,iseg)
              tempmat(i,ham_smat%keyg(1,2,iseg),1) = ham_mat%matrix_compr(imshift+ii)
              ii = ii + 1
          end do
      end do
      do iseg=1,ovrlp_smat%nseg
          ii=ovrlp_smat%keyv(iseg)
          do i=ovrlp_smat%keyg(1,1,iseg),ovrlp_smat%keyg(2,1,iseg)
              tempmat(i,ovrlp_smat%keyg(1,2,iseg),2) = ovrlp_mat%matrix_compr(isshift+ii)
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
      call sygv(1, 'n','l', ovrlp_smat%nfvctr, tempmat(1,1,1), ovrlp_smat%nfvctr, tempmat(1,1,2), ovrlp_smat%nfvctr, &
           eval(1), work(1), lwork, info)
      !!if (iproc==0) then
      !!    tt5 = 0.d0
      !!    tt7 = 0.d0
      !!    do ieval=1,ovrlp_smat%nfvctr
      !!        write(*,*) 'ieval',ieval,eval(ieval)
      !!        if (ieval<=5) tt5 = tt5 + eval(ieval)
      !!        if (ieval<=7) tt7 = tt7 + eval(ieval)
      !!    end do 
      !!    write(*,*) 'SUM of evals up to 5', tt5
      !!    write(*,*) 'SUM of evals up to 7', tt7
      !!end if
      if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(ovrlp_smat%nfvctr)/),fmt='(es16.6)')

      call f_free(tempmat)
      call f_free(eval)
      call f_free(work)

      call f_release_routine()
    
    end subroutine get_minmax_eigenvalues



    subroutine get_selected_eigenvalues(iproc, nproc, comm, calculate_minusonehalf, foe_verbosity, &
               iev_min, iev_max, fscale, &
               smats, smatm, smatl, ham_, ovrlp_, ovrlp_minus_one_half_, eval)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: symmetrize_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical, foe_data_null, foe_data_deallocate
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level
      use foe_common, only: retransform_ext, find_fermi_level, get_bounds_and_polynomials, init_foe
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, iev_min, iev_max
      real(mp),intent(in) :: fscale
      logical,intent(in) :: calculate_minusonehalf
      integer,intent(in) :: foe_verbosity
      type(sparse_matrix),intent(in) :: smats, smatm, smatl
      type(matrices),intent(in) :: ham_, ovrlp_
      type(matrices),dimension(1),intent(inout) :: ovrlp_minus_one_half_
      real(mp),dimension(iev_min:iev_max),intent(out) :: eval

      ! Local variables
      integer :: iev, i, ispin, ilshift, npl, npl_min, ind
      real(mp) :: dq, q, scale_factor, shift_value, factor, betax
      real(mp),dimension(:),allocatable :: charges
      type(matrices) :: kernel
      real(mp),dimension(1),parameter :: EF = 0.0_mp
      !real(mp),dimension(1),parameter :: FSCALE = 2.e-2_mp
      type(foe_data) :: foe_obj
      type(fermi_aux) :: f
      integer,parameter :: NPL_MAX = 10000
      integer,parameter :: NPL_STRIDE = 100
      real(mp),dimension(:,:,:),pointer :: chebyshev_polynomials
      real(mp),dimension(:),allocatable :: hamscal_compr
      type(f_progress_bar) :: bar

      call f_routine(id='get_selected_eigenvalues')


      kernel = matrices_null()
      kernel%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='kernel%matrix_compr')

      betax = foe_data_get_real(foe_obj,"betax")

      ! the occupation numbers...
      if (smatl%nspin==1) then
          factor = sqrt(0.5_mp)
      else
          factor = 1.0_mp
      end if

      ilshift = 0
      foe_obj = foe_data_null()
      charges = f_malloc(smatl%nspin,id='charges')
      do ispin=1,smatl%nspin
          charges(ispin) = real(iev_min,kind=mp)/real(smatl%nspin,kind=mp)
      end do
      call init_foe(iproc, nproc, smatl%nspin, charges, foe_obj, 0.0_mp, fscale=fscale)
      call f_free(charges)

      hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

      !if (iproc==0) call yaml_map('S^-1/2','recalculate')
      !!call overlap_minus_onehalf() ! has internal timer
      call overlap_minus_onehalf(iproc, nproc, comm, smats, smatl, ovrlp_, ovrlp_minus_one_half_, &
          verbosity=0) !has internal timer

      ! Use kernel_%matrix_compr as workarray to save memory
      npl_min = 10
      ispin = 1 !hack
      call get_bounds_and_polynomials(iproc, nproc, comm, 2, ispin, NPL_MAX, NPL_STRIDE, betax, &
           1, FUNCTION_ERRORFUNCTION, .false., 2.2_mp, 2.2_mp, 0, &
           smatm, smatl, ham_, foe_obj, npl_min, kernel%matrix_compr(ilshift+1:), &
           chebyshev_polynomials, npl, scale_factor, shift_value, hamscal_compr, &
           smats=smats, ovrlp_=ovrlp_, ovrlp_minus_one_half_=ovrlp_minus_one_half_(1), &
           efarr=EF, fscale_arr=(/fscale/))

     ! To determine the HOMO/LUMO, subtract/add one electrom for closed shell
      ! systems of one half electron for open shell systems.
      if (smatl%nspin==1) then
          dq = 1.d0
      else if (smatl%nspin==2) then
          dq = 0.5d0
      end if



      if (iproc==0) then
          call yaml_map('decay length of error function',fscale,fmt='(es8.2)')
          call yaml_map('degree of polynomial',npl)
          bar=f_progress_bar_new(nstep=iev_max-iev_min+1)
      end if
      do iev=iev_min,iev_max
          ! Calculate the 'lower' kernel
          do ispin=1,smatl%nspin
              q = real(iev*2,kind=mp)/real(smatl%nspin,kind=mp)-dq
              call foe_data_set_real(foe_obj,"charge",q,ispin)
          end do
          ispin = 1 !hack
          !!call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
          !!     foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
          !!     foe_data_get_real(foe_obj,"ef_interpol_det"), 0) !foe_verbosity)
          call find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
               0, 'test', smatl, ispin, foe_obj, kernel)
          eval(iev) = foe_data_get_real(foe_obj,"ef",ispin)
          !call retransform_ext(iproc, nproc, smatl, &
          !     ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), kernel(1)%matrix_compr(ilshift+1:))
      
          !!! Calculate the 'lower' kernel
          !!do ispin=1,smatl%nspin
          !!    q = real(iev*2,kind=mp)/real(smatl%nspin,kind=mp)-2.0_mp*dq
          !!    call foe_data_set_real(foe_obj,"charge",q,ispin)
          !!end do
          !!ispin = 1 !hack
          !!call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
          !!     foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
          !!     foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
          !!call find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
          !!     2, 'test', smatl, ispin, foe_obj, kernel(1))
          !!call retransform_ext(iproc, nproc, smatl, &
          !!     ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), kernel(1)%matrix_compr(ilshift+1:))

          !!! Calculate the 'upper' kernel
          !!do ispin=1,smatl%nspin
          !!    q = real(iev*2,kind=mp)/real(smatl%nspin,kind=mp)+0.0_mp*dq
          !!    call foe_data_set_real(foe_obj,"charge",q,ispin)
          !!end do
          !!ispin = 1 !hack
          !!call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
          !!     foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
          !!     foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
          !!call find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
          !!     2, 'test', smatl, ispin, foe_obj, kernel(2))
          !!call retransform_ext(iproc, nproc, smatl, &
          !!     ovrlp_minus_one_half_(1)%matrix_compr(ilshift+1:), kernel(2)%matrix_compr(ilshift+1:))

          !!! Calculate the square root of the difference, which is the eigenvector
          !!do i=1,smatl%nfvctr
          !!     ind = matrixindex_in_compressed(smatl, i, i)
          !!     write(*,*) 'i, evec', i, factor*sqrt(kernel(2)%matrix_compr(ind) - kernel(1)%matrix_compr(ind))
          !!end do

          if (iproc==0) call dump_progress_bar(bar,step=iev-iev_min+1)

       end do

       call f_free_ptr(chebyshev_polynomials)
       call deallocate_matrices(kernel)
       call f_free(hamscal_compr)
       call foe_data_deallocate(foe_obj)
       
      call f_release_routine()

     !!  contains
     !!       subroutine overlap_minus_onehalf()
     !!         use ice, only: inverse_chebyshev_expansion_new
     !!         implicit none
     !!         integer :: i, j, ii
     !!         real(kind=mp) :: max_error, mean_error
     !!         real(kind=mp), dimension(1) :: ex
     !!         real(kind=mp),dimension(:),allocatable :: tmparr
    
     !!         call f_routine(id='overlap_minus_onehalf')
    
     !!         ! Can't use the wrapper, since it is at a higher level in the hierarchy (to be improved)
     !!         ex=-0.5d0
     !!         call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
     !!              ovrlp_smat=smats, inv_ovrlp_smat=smatl, ncalc=1, ex=ex, &
     !!              ovrlp_mat=ovrlp_, inv_ovrlp=ovrlp_minus_one_half_, &
     !!              verbosity=0)

    
     !!         call f_release_routine()
     !!     end subroutine overlap_minus_onehalf

    end subroutine get_selected_eigenvalues


    subroutine overlap_minus_onehalf(iproc, nproc, comm, smats, smatl, ovrlp_, ovrlp_minus_one_half_, &
               verbosity, ice_obj)
      use foe_base, only: foe_data
      use ice, only: inverse_chebyshev_expansion_new
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smats, smatl
      type(matrices),intent(in) :: ovrlp_
      type(matrices),dimension(1),intent(out) :: ovrlp_minus_one_half_
      integer,intent(in),optional :: verbosity
      type(foe_data),intent(inout),optional :: ice_obj
      ! Local variables
      integer :: verbosity_
      real(mp),dimension(1) :: ex
    
      call f_routine(id='overlap_minus_onehalf')

      verbosity_ = 1
      if (present(verbosity)) verbosity_ = verbosity
    
      ! Can't use the wrapper, since it is at a higher level in the hierarchy (to be improved)
      ex=-0.5d0
      if (present(ice_obj)) then
          call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
               ovrlp_smat=smats, inv_ovrlp_smat=smatl, ncalc=1, ex=ex, &
               ovrlp_mat=ovrlp_, inv_ovrlp=ovrlp_minus_one_half_, &
               verbosity=verbosity_, ice_objx=ice_obj)
      else
          call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
               ovrlp_smat=smats, inv_ovrlp_smat=smatl, ncalc=1, ex=ex, &
               ovrlp_mat=ovrlp_, inv_ovrlp=ovrlp_minus_one_half_, &
               verbosity=verbosity_)
      end if
    
      call f_release_routine()
    end subroutine overlap_minus_onehalf


end module foe
