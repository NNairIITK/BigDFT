module module_func
  use sparsematrix_base
  !use module_base, only: safe_exp
  implicit none

  private

  ! Shared variables within the modules
  integer :: ifunc
  real(kind=mp) :: power, ef, fscale, beta, mua, mub

  ! Public routines
  public :: func_set
  public :: func

  ! Public parameters
  integer,parameter,public :: FUNCTION_POLYNOMIAL = 101
  integer,parameter,public :: FUNCTION_ERRORFUNCTION = 102
  integer,parameter,public :: FUNCTION_EXPONENTIAL = 103

  contains

    subroutine func_set(ifuncx, powerx, efx, fscalex, betax, muax, mubx)
      implicit none
      integer,intent(in) :: ifuncx
      real(kind=mp),intent(in),optional :: powerx, efx, fscalex, betax, muax, mubx

      call f_routine(id='func_set')

      select case (ifuncx)
      case(FUNCTION_POLYNOMIAL)
          ifunc = FUNCTION_POLYNOMIAL
          if (.not.present(powerx)) call f_err_throw("'powerx' not present")
          power = powerx
      case(FUNCTION_ERRORFUNCTION)
          ifunc = FUNCTION_ERRORFUNCTION
          if (.not.present(efx)) call f_err_throw("'efx' not present")
          if (.not.present(fscalex)) call f_err_throw("'fscalex' not present")
          ef = efx
          fscale = fscalex
      case(FUNCTION_EXPONENTIAL)
          ifunc = FUNCTION_EXPONENTIAL
          if (.not.present(betax)) call f_err_throw("'betax' not present")
          if (.not.present(muax)) call f_err_throw("'muax' not present")
          if (.not.present(mubx)) call f_err_throw("'mubx' not present")
          beta = betax
          mua = muax
          mub = mubx
      case default
          call f_err_throw("wrong value of 'ifuncx'")
      end select

      call f_release_routine()

    end subroutine func_set

    function func(x)
      implicit none
      real(kind=mp),intent(in) :: x
      real(kind=mp) :: func
      select case (ifunc)
      case(FUNCTION_POLYNOMIAL)
          func = x**power
      case(FUNCTION_ERRORFUNCTION)
          func = 0.5d0*erfcc((x-ef)*(1.d0/fscale))
      case(FUNCTION_EXPONENTIAL)
          !func = safe_exp(beta*(x-mu))
          func = safe_exp(beta*(x-mua)) - safe_exp(-beta*(x-mub))
      case default
          call f_err_throw("wrong value of 'ifunc'")
      end select
    end function func



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
      erfcc=t*safe_exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
            & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
            & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.D0-erfcc

    end function erfcc

end module module_func


module foe_common
  use foe_base
  use sparsematrix_base
  implicit none

  private

  !> Public routines
  public :: get_chebyshev_expansion_coefficients
  !public :: chebft
  !public :: chebyshev_coefficients_penalyfunction
  public :: chder
  public :: evnoise
  !public :: erfcc
  !public :: chebev
  public :: pltwght
  public :: pltexp
  !public :: check_eigenvalue_spectrum_new
  !public :: scale_and_shift_matrix
  public :: retransform_ext
  !!public :: cheb_exp
  public :: init_foe
  !public :: get_chebyshev_polynomials
  public :: find_fermi_level
  public :: get_polynomial_degree
  public :: calculate_trace_distributed_new
  public :: get_bounds_and_polynomials


  contains


    subroutine get_chebyshev_expansion_coefficients(iproc, nproc, comm, A, B, N, func, cc, x_max_error,max_error,mean_error)
      use yaml_output
      use sparsematrix_init, only: distribute_on_tasks
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, n
      real(kind=mp),intent(in) :: A, B
      real(kind=mp),external :: func
      real(8),dimension(n),intent(out) :: cc
      real(kind=mp),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: k, j, is, np, ii, jj
      real(kind=mp) :: bma, bpa, y, arg, fac, tt, one_over_n
      real(kind=mp),dimension(:),allocatable :: cf

      call f_routine(id='get_chebyshev_expansion_coefficients')

      ! MPI parallelization... maybe only worth for large n?
      !call chebyshev_coefficients_init_parallelization(iproc, nproc, comm, n, np, is)
      ! Initialize the parallelization.
      call distribute_on_tasks(n, iproc, nproc, np, is)

      call chebyshev_coefficients_calculate(n, a, b, np, is, func, cc)

      call chebyshev_coefficients_communicate(comm, n, cc)

      call accuracy_of_chebyshev_expansion(iproc, nproc, comm, n, cc, A,B, &
           1.d-3, func, x_max_error, max_error, mean_error)

      call f_release_routine()

    end subroutine get_chebyshev_expansion_coefficients

    !!! Calculates chebychev expansion of fermi distribution.
    !!! Taken from numerical receipes: press et al
    !!subroutine chebft(iproc,nproc,A,B,N,cc,ef,fscale,tmprtr,x_max_error,max_error,mean_error)
    !!  use module_base
    !!  use module_func
    !!  use yaml_output
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  real(kind=mp),intent(in) :: A, B, ef, fscale, tmprtr
    !!  integer,intent(in) :: iproc, nproc, n
    !!  real(8),dimension(n),intent(out) :: cc
    !!  real(kind=mp),intent(out) :: x_max_error, max_error, mean_error
    !!
    !!  ! Local variables
    !!  integer :: k, j, is, np, ii, jj
    !!  real(kind=mp) :: bma, bpa, y, arg, fac, tt
    !!  real(kind=mp),dimension(50000) :: cf
    !!  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!
    !!  call f_routine(id='chebft')

    !!  if (tmprtr/=0.d0) call f_err_throw('tmprtr should be zero for the moment')

    !!  ! MPI parallelization... maybe only worth for large n?
    !!  ii = n/nproc
    !!  np = ii
    !!  is = iproc*ii
    !!  ii = n - nproc*ii
    !!  if (iproc<ii) then
    !!      np = np + 1
    !!  end if
    !!  is = is + min(iproc,ii)

    !!  !write(*,*) 'iproc, nproc, is, np, n', iproc, nproc, is, np, n
    !!  call f_zero(cc)
    !!
    !!  if (n>50000) stop 'chebft'
    !!  bma=0.5d0*(b-a)
    !!  bpa=0.5d0*(b+a)
    !!  fac=2.d0/n
    !!  !$omp parallel default(none) shared(bma,bpa,fac,n,tmprtr,cf,fscale,ef,cc,is,np,tt) &
    !!  !$omp private(k,y,arg,j,jj)
    !!  !$omp do
    !!  do k=1,n
    !!      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!      arg=y*bma+bpa
    !!      if (tmprtr.eq.0.d0) then
    !!          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
    !!      else
    !!          cf(k)=1.d0/(1.d0+safe_exp( (arg-ef)*(1.d0/tmprtr) ) )
    !!      end if
    !!  end do
    !!  !$omp end do
    !!  !$omp end parallel
    !!  do j=1,np
    !!      jj = j + is
    !!      tt=0.d0
    !!      !$omp parallel do default(none) shared(n,cf,jj) private(k) reduction(+:tt)
    !!      do  k=1,n
    !!          tt=tt+cf(k)*cos((pi*(jj-1))*((k-0.5d0)*(1.d0/n)))
    !!      end do
    !!      !$omp end parallel do
    !!      cc(jj)=fac*tt
    !!  end do

    !!  call mpiallred(cc, mpi_sum, comm=bigdft_mpi%mpi_comm)

    !!  call func_set(FUNCTION_ERRORFUNCTION, efx=ef, fscalex=fscale)
    !!  call accuracy_of_chebyshev_expansion(n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
    !!  !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    !!
    !!  call f_release_routine()
    !!
    !!end subroutine chebft



    !!!! Calculates chebychev expansion of fermi distribution.
    !!!! Taken from numerical receipes: press et al
    !!!subroutine chebyshev_coefficients_penalyfunction(a,b,n,cc,max_error)
    !!!  use module_base
    !!!  use module_func
    !!!  implicit none
    !!!
    !!!  ! Calling arguments
    !!!  real(kind=mp),intent(in) :: a, b
    !!!  integer,intent(in) :: n
    !!!  real(kind=mp),dimension(n,2),intent(out) :: cc
    !!!  real(kind=mp),intent(out) :: max_error
    !!!
    !!!  ! Local variables
    !!!  integer :: k, j
    !!!  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!!  real(kind=mp) :: tt1, tt2, ttt, y, arg, fac, bma, bpa, x_max, max_err, mean_err
    !!!  real(kind=mp),dimension(50000) :: cf
    !!!
    !!!  call f_routine(id='chebyshev_coefficients_penalyfunction')
    !!!
    !!!  if (n>50000) stop 'chebyshev_coefficients_penalyfunction'
    !!!  bma=0.5d0*(b-a)
    !!!  bpa=0.5d0*(b+a)
    !!!  ! 3 gives broder safety zone than 4
    !!!  !ttt=3.0d0*n/(b-a)
    !!!  !ttt=4.d0*n/(b-a)
    !!!  ttt=40.d0
    !!!  fac=2.d0/n
    !!!  !$omp parallel default(none) shared(bma,bpa,ttt,fac,n,cf,a,b,cc) &
    !!!  !$omp private(k,y,arg,tt1,tt2,j)
    !!!  !$omp do
    !!!  do k=1,n
    !!!      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!!      arg=y*bma+bpa
    !!!      !write(*,*) 'arg, safe_exp(-(arg-a)*ttt)', arg, safe_exp(-(arg-a)*ttt)
    !!!      cf(k)= safe_exp(-(arg-a)*ttt)-safe_exp((arg-b)*ttt)
    !!!      !cf(k,2)=-safe_exp(-(arg-a)*ttt)+safe_exp((arg-b)*ttt)
    !!!      !cf(k,1)= safe_exp(-(arg-a)*ttt)
    !!!      !cf(k,2)= safe_exp((arg-b)*ttt)
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp do
    !!!  do j=1,n
    !!!      tt1=0.d0
    !!!      tt2=0.d0
    !!!      do k=1,n
    !!!          tt1=tt1+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!!          !tt2=tt2+cf(k,2)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!!      end do
    !!!      cc(j,1)=fac*tt1
    !!!      !cc(j,2)=fac*tt2
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp do
    !!!  do j=1,n
    !!!      cc(j,2) = -cc(j,1)
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp end parallel
    !!!
    !!!  !!do j=1,n
    !!!  !!    write(*,*) 'j, cc(j,1), cc(j,2)', j, cc(j,1), cc(j,2)
    !!!  !!end do
    !!!  call func_set(FUNCTION_EXPONENTIAL, betax=-ttt, muax=a, mubx=b)
    !!!  call accuracy_of_chebyshev_expansion(n, cc(:,1), (/A,B/), 1.d-3, func, x_max, max_err, mean_err)
    !!!  max_error = max_err
    !!!  call func_set(FUNCTION_EXPONENTIAL, betax=ttt, muax=b, mubx=a)
    !!!  call accuracy_of_chebyshev_expansion(n, cc(:,2), (/A,B/), 1.d-3, func, x_max, max_err, mean_err)
    !!!  max_error = max(max_error,max_err)
    !!!  call f_release_routine()
    !!!
    !!!end subroutine chebyshev_coefficients_penalyfunction

    ! Calculates chebychev expansion of the derivative of Fermi distribution.
    subroutine chder(a,b,c,cder,n)
      use dynamic_memory
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: a, b
      integer,intent(in) :: n
      real(8),dimension(n),intent(in) :: c
      real(8),dimension(n),intent(out) :: cder

      ! Local variables
      integer :: j
      real(kind=mp) :: con

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
      real(kind=mp),dimension(npl),intent(in) :: cc
      real(kind=mp),intent(in) :: evlow, evhigh
      real(kind=mp),intent(out) :: anoise

      ! Local variables
      integer :: i, n
      real(kind=mp) :: fact, dist, ddx, cent, tt, x

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
          x=real(i,kind=mp)*ddx
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





    !> Evaluates chebychev expansion
    function chebev(a,b,m,x,cc)
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: a, b, x
      integer,intent(in) :: m
      real(kind=mp),dimension(m),intent(in) :: cc
      real(kind=mp) :: chebev

      ! Local variables
      integer :: j
      real(kind=mp) :: d, dd, y, sv

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
              real(kind=mp),dimension(npl),intent(in) :: cc, cder
              real(kind=mp),intent(in) :: evlow, evhigh, ef, fscale, tmprtr

              ! Local variables
              integer :: ic
              real(kind=mp) :: ddx, x, tt, err

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
            real(kind=mp),dimension(npl),intent(in) :: cc
            real(kind=mp),intent(in) :: anoise, evlow, evhigh

            ! Local variables
            integer :: ic
            real(kind=mp) :: fact, ddx, tt, x

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



    subroutine check_eigenvalue_spectrum_new(iproc, nproc, comm, smat_l, ispin, isshift, &
               factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
               emergency_stop, foe_obj, restart, eval_bounds_ok, &
               verbosity, eval_multiplicator, smat_s, mat)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: transform_sparse_matrix
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin, isshift
      type(sparse_matrix),intent(in) :: smat_l
      real(kind=mp),intent(in) :: factor_high, factor_low, anoise
      !real(kind=mp),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
      real(kind=mp),dimension(smat_l%smmm%nvctrp),intent(in) :: penalty_ev
      logical,intent(in) :: trace_with_overlap
      logical,dimension(2) :: emergency_stop
      type(foe_data),intent(inout) :: foe_obj
      logical,intent(inout) :: restart
      logical,dimension(2),intent(out) :: eval_bounds_ok
      integer,intent(in),optional :: verbosity
      real(kind=mp),intent(inout),optional :: eval_multiplicator
      type(sparse_matrix),intent(in),optional :: smat_s
      type(matrices),intent(in),optional :: mat

      ! Local variables
      integer :: isegstart, isegend, iseg, ii, jorb, irow, icol, iismall, iel, i, iline, icolumn, ibound, verbosity_
      integer :: ishift
      real(mp),dimension(:),allocatable :: mat_large
      real(kind=mp) :: tt, noise, penalty
      !real(kind=mp),dimension(1) :: allredarr

      call f_routine(id='check_eigenvalue_spectrum_new')

      if (.not.smat_l%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      if (present(verbosity)) then
          verbosity_ = verbosity
      else
          verbosity_ = 1
      end if

      if (trace_with_overlap) then
          if (.not.present(smat_s) .or. .not.present(mat)) then
              call f_err_throw('not all required optional arguments are present')
          end if
          mat_large = sparsematrix_malloc(smat_l, iaction=sparse_taskgroup, id='mat_large')
          call transform_sparse_matrix(iproc, smat_s, smat_l, sparse_taskgroup, 'small_to_large', &
               smat_in=mat%matrix_compr, lmat_out=mat_large)
      end if

      penalty=0.d0
      ishift = (ispin-1)*smat_l%nvctrp_tg
      !$omp parallel if (smat_l%smmm%nvctrp>1000) &
      !$omp default(none) &
      !$omp shared(smat_l, penalty, trace_with_overlap, mat_large, ishift, penalty_ev) &
      !$omp private(i, ii, iline, icolumn, iismall, tt)
      !$omp do schedule(static) reduction(+:penalty)
      do i=1,smat_l%smmm%nvctrp
          ii = smat_l%smmm%isvctr + i
          iline = smat_l%smmm%line_and_column(1,i)
          icolumn = smat_l%smmm%line_and_column(2,i)
          if (trace_with_overlap) then
              ! Take the trace of the product matrix times overlap
              tt = mat_large(ishift+i)
          else
              ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
              if (iline==icolumn) then
                  tt=1.d0
              else
                  tt=0.d0
              end if
          end if
          penalty = penalty + penalty_ev(i)*tt
      end do
      !$omp end do
      !$omp end parallel

      if (trace_with_overlap) then
          call f_free(mat_large)
      end if

      ! Divide the traces by the matrix dimension, to make them size independent
      penalty = penalty/real(smat_l%nfvctr,kind=mp)

      !!if (nproc > 1) then
      !!    call mpiallred(penalty, mpi_sum, comm=comm)
      !!end if
      call penalty_communicate(nproc, comm, penalty)


      !noise=10.d0*anoise
      noise = 1.d-5

      if (iproc==0 .and. verbosity_>0) then
          !call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
          !call yaml_map('pnlty',(/allredarr(1),allredarr(2)/),fmt='(es8.1)')
          call yaml_map('penalty',penalty,fmt='(es8.1)')
      end if

      eval_bounds_ok(1) = .true.
      eval_bounds_ok(2) = .true.
      if (any((/abs(penalty)>noise/))) then
          if (all((/abs(penalty)>noise/))) then
              if (penalty>0.d0) then
                  ! lower bound too large
                  eval_bounds_ok(1)=.false.
                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
                  restart=.true.
                  if (present(eval_multiplicator)) then
                      !eval_multiplicator = eval_multiplicator*2.0d0
                      eval_multiplicator = 2.0d0
                  end if
              else if (penalty<0.d0) then
                  ! upper bound too small
                  eval_bounds_ok(2)=.false.
                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*factor_high,ispin)
                  restart=.true.
                  if (present(eval_multiplicator)) then
                      !eval_multiplicator = eval_multiplicator/2.0d0
                      eval_multiplicator = 1.d0/2.0d0
                  end if
              else
                  call f_err_throw('The errors should have opposite signs')
              end if
          else
              call f_err_throw('The errors should have the same magnitude')
          end if
      end if

      call f_release_routine()

    end subroutine check_eigenvalue_spectrum_new


    subroutine scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
               smat1, mat1, i1shift, smat2, mat2, i2shift, &
               matscal_compr, scale_factor, shift_value)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: transform_sparse_matrix
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ispin, i1shift
      type(foe_data),intent(in) :: foe_obj
      type(sparse_matrix),intent(in) :: smatl, smat1
      type(matrices),intent(in) :: mat1
      type(sparse_matrix),intent(in),optional :: smat2
      type(matrices),intent(in),optional :: mat2
      integer,intent(in),optional :: i2shift
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(out) :: matscal_compr
      real(kind=mp),intent(out) :: scale_factor, shift_value

      ! Local variables
      integer :: iseg, ii, i, ii1, ii2, isegstart, isegend, ierr, jj
      integer :: itaskgroup, iitaskgroup, j, ishift
      integer,dimension(2) :: irowcol
      real(kind=mp) :: tt1, tt2
      logical :: with_overlap
      real(kind=mp),dimension(:),pointer :: matscal_compr_local
      real(kind=mp),dimension(:),allocatable :: mat1_large, mat2_large
      integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: comm_strategy=GET
      integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX


      call f_routine(id='scale_and_shift_matrix')
      !call timing(iproc,'foe_aux_mcpy  ','ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      call f_zero(matscal_compr)

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
      else if (data_strategy==SUBMATRIX) then

          ishift = (ispin-1)*smatl%nvctrp_tg

          ! Transform all matrices to the large sparsity pattern.
          ! Takes some memory, but probably faster than the old way...
          mat1_large = sparsematrix_malloc(smatl, iaction=sparse_taskgroup, id='mat1_large')
          call transform_sparse_matrix(iproc, smat1, smatl, sparse_taskgroup, 'small_to_large', &
               smat_in=mat1%matrix_compr, lmat_out=mat1_large)
          if (with_overlap) then
              mat2_large = sparsematrix_malloc(smatl, iaction=sparse_taskgroup, id='mat2_large')
              call transform_sparse_matrix(iproc, smat2, smatl, sparse_taskgroup, 'small_to_large', &
                   smat_in=mat2%matrix_compr, lmat_out=mat2_large)
          end if

          !$omp parallel if (smatl%smmm%istartendseg_mm(2)-smatl%smmm%istartendseg_mm(1)>100) &
          !$omp default(none) &
          !$omp shared(smatl, mat1_large, mat2_large, with_overlap, matscal_compr, scale_factor, shift_value) &
          !$omp shared(ishift) &
          !$omp private(iseg, j, ii, i, jj, tt1, tt2)
          !$omp do schedule(guided)
          do iseg=smatl%smmm%istartendseg_mm(1),smatl%smmm%istartendseg_mm(2)
              ! A segment is always on one line, therefore no double loop
              j = smatl%keyg(1,2,iseg)
              ii=smatl%keyv(iseg)-1
              do i=smatl%keyg(1,1,iseg),smatl%keyg(2,1,iseg)
                  ii = ii + 1
                  jj = ii-smatl%isvctrp_tg
                  tt1=mat1_large(ishift+jj)
                  if (with_overlap) then
                      tt2 = mat2_large(ishift+jj)
                  else
                      if (i==j) then
                          tt2 = 1.d0
                      else
                          tt2 = 0.d0
                      end if
                  end if
                  !!write(200+iproc,*) 'jj, tt1, tt2', jj, tt1, tt2
                  matscal_compr(jj)=scale_factor*(tt1-shift_value*tt2)
              end do
          end do
          !$omp end do
          !$omp end parallel

          call f_free(mat1_large)
          if (with_overlap) then
              call f_free(mat2_large)
          end if


          !!!write(*,*) 'smatl%smmm%istartendseg_mm',smatl%smmm%istartendseg_mm
          !!!$omp parallel default(none) private(ii,i,j,ii2,ii1,tt2,tt1,iseg) &
          !!!$omp shared(matscal_compr,scale_factor,shift_value,i2shift,i1shift,smatl,smat1,smat2,mat1,mat2,with_overlap)
          !!!$omp do
          !!do iseg=smatl%smmm%istartendseg_mm(1),smatl%smmm%istartendseg_mm(2)
          !!    !if (smatl%keyv(min(iseg+1,smatl%nseg))<smatl%smmm%istartend_mm(1)) cycle
          !!    !if (smatl%keyv(iseg)>smatl%smmm%istartend_mm(2)) exit
          !!    ! A segment is always on one line, therefore no double loop
          !!    j = smatl%keyg(1,2,iseg)
          !!    ii=smatl%keyv(iseg)-1
          !!    do i=smatl%keyg(1,1,iseg),smatl%keyg(2,1,iseg)
          !!        ii = ii + 1
          !!        ii1 = matrixindex_in_compressed(smat1, i, j)
          !!        if (ii1>0) then
          !!            tt1=mat1%matrix_compr(i1shift+ii1-smat1%isvctrp_tg)
          !!        else
          !!            tt1=0.d0
          !!        end if
          !!        if (with_overlap) then
          !!            ii2 = matrixindex_in_compressed(smat2, i, j)
          !!            if (ii2>0) then
          !!                tt2=mat2%matrix_compr(i2shift+ii2-smat2%isvctrp_tg)
          !!            else
          !!                tt2=0.d0
          !!            end if
          !!        else
          !!            if (i==j) then
          !!                tt2 = 1.d0
          !!            else
          !!                tt2 = 0.d0
          !!            end if
          !!        end if
          !!        !ii=matrixindex_in_compressed(smatl, i, j)
          !!        !write(*,*) 'i, ii, ii1, tt1, tt2', i, ii, ii1, tt1, tt2, i1shift, smat1%isvctrp_tg, i1shift+ii1-smat1%isvctrp_tg
          !!        write(300+iproc,*) 'jj, tt1, tt2', ii-smatl%isvctrp_tg, tt1, tt2
          !!        matscal_compr(ii-smatl%isvctrp_tg)=scale_factor*(tt1-shift_value*tt2)
          !!    end do
          !!end do
          !!!$omp end do
          !!!$omp end parallel
          !!!call timing(iproc,'foe_aux_mcpy  ','OF')
          call f_timing(TCAT_CME_AUXILIARY,'OF')
      else
          stop 'scale_and_shift_matrix: wrong data strategy'
      end if

      call f_release_routine()

    end subroutine scale_and_shift_matrix


    subroutine retransform_ext(iproc, nproc, smat, inv_ovrlp, kernel)
        use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
                                compress_matrix_distributed_wrapper, &
                                sparsemm_new, transform_sparsity_pattern
        implicit none
        ! Calling arguments
        integer,intent(in) :: iproc, nproc
        type(sparse_matrix),intent(in) :: smat
        real(kind=mp),dimension(smat%nvctrp_tg),intent(inout) :: inv_ovrlp
        real(kind=mp),dimension(smat%nvctrp_tg),intent(inout) :: kernel


        ! Local variables
        real(kind=mp),dimension(:),pointer :: inv_ovrlpp_new, tempp_new
        real(kind=mp),dimension(:),allocatable :: mat_compr_seq

        call f_routine(id='retransform_ext')

        if (.not.smat%smatmul_initialized) then
            call f_err_throw('sparse matrix multiplication not initialized', &
                 err_name='SPARSEMATRIX_RUNTIME_ERROR')
        end if
    
        inv_ovrlpp_new = f_malloc_ptr(smat%smmm%nvctrp, id='inv_ovrlpp_new')
        tempp_new = f_malloc_ptr(smat%smmm%nvctrp, id='tempp_new')
        mat_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='mat_compr_seq')
        !!kernel_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='kernel_compr_seq')
        if (smat%smmm%nvctrp_mm>0) then !to avoid an out of bounds error
            call transform_sparsity_pattern(iproc, smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
                 smat%nseg, smat%keyv, smat%keyg, smat%smmm%line_and_column_mm, &
                 smat%smmm%nvctrp, smat%smmm%isvctr, &
                 smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
                 'small_to_large', &
                 matrix_s_in=inv_ovrlp(smat%smmm%isvctr_mm-smat%isvctrp_tg+1), matrix_l_out=inv_ovrlpp_new)
        end if
        call sequential_acces_matrix_fast2(smat, kernel, mat_compr_seq)
        call sparsemm_new(iproc, smat, mat_compr_seq, inv_ovrlpp_new, tempp_new)
        call sequential_acces_matrix_fast2(smat, inv_ovrlp, mat_compr_seq)
        call sparsemm_new(iproc, smat, mat_compr_seq, tempp_new, inv_ovrlpp_new)
        call f_zero(kernel)
        call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_LARGE, &
             inv_ovrlpp_new, kernel)
        call f_free_ptr(inv_ovrlpp_new)
        call f_free_ptr(tempp_new)
        call f_free(mat_compr_seq)
        !!call f_free(kernel_compr_seq)

        call f_release_routine()

    end subroutine retransform_ext



    !!##! Calculates chebychev expansion of x**ex, where ex is any value (typically -1, -1/2, 1/2)
    !!##! Taken from numerical receipes: press et al
    !!##subroutine cheb_exp(iproc, nproc, A,B,N,cc,ex,x_max_error,max_error,mean_error)
    !!##  use module_base
    !!##  use module_func
    !!##  use yaml_output
    !!##  implicit none
    !!##
    !!##  ! Calling arguments
    !!##  integer,intent(in) :: iproc, nproc
    !!##  real(kind=mp),intent(in) :: A, B
    !!##  integer,intent(in) :: n
    !!##  real(kind=mp),intent(in) :: ex
    !!##  real(8),dimension(n),intent(out) :: cc
    !!##  real(kind=mp),intent(out) :: x_max_error, max_error,mean_error
    !!##
    !!##  ! Local variables
    !!##  integer :: k, j
    !!##  real(kind=mp) :: bma, bpa, y, arg, fac, tt
    !!##  real(kind=mp),dimension(50000) :: cf
    !!##  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!##
    !!##  call f_routine(id='chebft')

    !!##
    !!##  if (n>50000) stop 'chebft'
    !!##  bma=0.5d0*(b-a)
    !!##  bpa=0.5d0*(b+a)
    !!##  fac=2.d0/n
    !!##  !$omp parallel default(none) shared(bma,bpa,fac,n,cf,cc,ex) &
    !!##  !$omp private(k,y,arg,j,tt)
    !!##  !$omp do
    !!##  do k=1,n
    !!##      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!##      arg=y*bma+bpa
    !!##      cf(k)=arg**ex
    !!##  end do
    !!##  !$omp end do
    !!##  !$omp do
    !!##  do j=1,n
    !!##      tt=0.d0
    !!##      do  k=1,n
    !!##          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!##      end do
    !!##      cc(j)=fac*tt
    !!##  end do
    !!##  !$omp end do
    !!##  !$omp end parallel

    !!##  call func_set(FUNCTION_POLYNOMIAL, powerx=ex)
    !!##  call accuracy_of_chebyshev_expansion(iproc, nproc, n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
    !!##  !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    !!##
    !!##  call f_release_routine()
    !!##
    !!##end subroutine cheb_exp


    subroutine init_foe(iproc, nproc, nspin, charge, foe_obj, tmprtr, evbounds_nsatur, evboundsshrink_nsatur, &
               evlow, evhigh, fscale, ef_interpol_det, ef_interpol_chargediff, &
               fscale_lowerbound, fscale_upperbound, eval_multiplicator)
      use foe_base, only: foe_data, foe_data_set_int, foe_data_set_real, foe_data_set_logical, foe_data_get_real, foe_data_null
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin
      real(kind=mp),dimension(nspin),intent(in) :: charge
      type(foe_data),intent(out) :: foe_obj
      integer,intent(in),optional :: evbounds_nsatur
      integer,intent(in),optional :: evboundsshrink_nsatur
      real(kind=mp),intent(in),optional :: evlow
      real(kind=mp),intent(in),optional :: evhigh
      real(kind=mp),intent(in),optional :: fscale
      real(kind=mp),intent(in),optional :: ef_interpol_det
      real(kind=mp),intent(in),optional :: ef_interpol_chargediff
      real(kind=mp),intent(in),optional :: fscale_lowerbound
      real(kind=mp),intent(in),optional :: fscale_upperbound
      real(kind=mp),intent(in),optional :: tmprtr
      real(kind=mp),intent(in),optional :: eval_multiplicator

      ! Local variables
      character(len=*), parameter :: subname='init_foe'
      integer :: ispin
      integer :: evbounds_nsatur_
      integer :: evboundsshrink_nsatur_
      real(kind=mp) :: evlow_
      real(kind=mp) :: evhigh_
      real(kind=mp) :: fscale_
      real(kind=mp) :: ef_interpol_det_
      real(kind=mp) :: ef_interpol_chargediff_
      real(kind=mp) :: fscale_lowerbound_
      real(kind=mp) :: fscale_upperbound_
      real(kind=mp) :: tmprtr_
      real(kind=mp) :: eval_multiplicator_

      !call timing(iproc,'init_matrCompr','ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      ! Define the default values... Is there a way to get them from input_variables_definition.yaml?
      evbounds_nsatur_ = 3
      evboundsshrink_nsatur_ =4
      evlow_ = -0.5_mp
      evhigh_ = 0.5_mp
      fscale_ = 2.d-2
      ef_interpol_det_ = 1.d-12
      ef_interpol_chargediff_ = 1.0_mp
      fscale_lowerbound_ = 5.d-3
      fscale_upperbound_ = 5.d-2
      tmprtr_ = 0.0_mp
      eval_multiplicator_ = 1.0_mp

      if (present(evbounds_nsatur)) evbounds_nsatur_ = evbounds_nsatur
      if (present(evboundsshrink_nsatur)) evboundsshrink_nsatur_ = evboundsshrink_nsatur
      if (present(evlow)) evlow_ = evlow
      if (present(evhigh)) evhigh_ = evhigh
      if (present(fscale)) fscale_ = fscale
      if (present(ef_interpol_det)) ef_interpol_det_ = ef_interpol_det
      if (present(ef_interpol_chargediff)) ef_interpol_chargediff_ = ef_interpol_chargediff
      if (present(fscale_lowerbound)) fscale_lowerbound_ = fscale_lowerbound
      if (present(fscale_upperbound)) fscale_upperbound_ = fscale_upperbound
      if (present(tmprtr)) tmprtr_ = tmprtr
      if (present(eval_multiplicator)) eval_multiplicator_ = eval_multiplicator
    
      foe_obj = foe_data_null()

      call foe_data_set_real(foe_obj,"fscale",fscale_)
      call foe_data_set_real(foe_obj,"ef_interpol_det",ef_interpol_det_)
      call foe_data_set_real(foe_obj,"ef_interpol_chargediff",ef_interpol_chargediff_)
      call foe_data_set_int(foe_obj,"evbounds_isatur",0)
      call foe_data_set_int(foe_obj,"evboundsshrink_isatur",0)
      call foe_data_set_int(foe_obj,"evbounds_nsatur",evbounds_nsatur_)
      call foe_data_set_int(foe_obj,"evboundsshrink_nsatur",evboundsshrink_nsatur_)
      call foe_data_set_real(foe_obj,"fscale_lowerbound",fscale_lowerbound_)
      call foe_data_set_real(foe_obj,"fscale_upperbound",fscale_upperbound_)
      call foe_data_set_real(foe_obj,"tmprtr",tmprtr_)

      foe_obj%charge = f_malloc0_ptr(nspin,id='foe_obj%charge')
      foe_obj%ef = f_malloc0_ptr(nspin,id='(foe_obj%ef)')
      foe_obj%evlow = f_malloc0_ptr(nspin,id='foe_obj%evlow')
      foe_obj%evhigh = f_malloc0_ptr(nspin,id='foe_obj%evhigh')
      foe_obj%bisection_shift = f_malloc0_ptr(nspin,id='foe_obj%bisection_shift')
      foe_obj%eval_multiplicator = f_malloc0_ptr(nspin,id='foe_obj%eval_multiplicator')
      do ispin=1,nspin
          call foe_data_set_real(foe_obj,"charge",charge(ispin),ispin)
          call foe_data_set_real(foe_obj,"ef",0.d0,ispin)
          call foe_data_set_real(foe_obj,"evhigh",evhigh_,ispin)
          call foe_data_set_real(foe_obj,"evlow",evlow_,ispin)
          call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,ispin)
          call foe_data_set_real(foe_obj,"eval_multiplicator",eval_multiplicator_,ispin)
      end do

      !call timing(iproc,'init_matrCompr','OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')


    end subroutine init_foe


    subroutine accuracy_of_chebyshev_expansion(iproc, nproc, comm, npl, coeff, bound_lower, bound_upper, &
               h, func, x_max_error, max_error, mean_error)
      use sparsematrix_init, only: distribute_on_tasks
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, npl
      real(kind=mp),dimension(npl),intent(in) :: coeff
      real(kind=mp),intent(in) :: bound_lower, bound_upper
      real(kind=mp),intent(in) :: h
      real(kind=mp),external :: func
      real(kind=mp),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: isx, iex, i, ipl, n, iimin, iimax, ii, is, np, jproc, ithread, nthread
      real(kind=mp) :: x, xx, val_chebyshev, val_function, xxm1, xxm2, xxx, sigma, tau, error, val_chebyshev1
      real(kind=mp),dimension(:,:),allocatable :: max_errors
      real(kind=mp),dimension(:),allocatable :: max_error_arr, x_max_error_arr
      !$ integer :: omp_get_thread_num, omp_get_max_threads

      call f_routine(id='accuracy_of_chebyshev_expansion')

      sigma = 2.d0/(bound_upper-bound_lower)
      tau = (bound_lower+bound_upper)/2.d0

      isx = ceiling(bound_lower/h)
      iex = floor(bound_upper/h)
      n = iex - isx + 1

      ! MPI parallelization... maybe only worth for large n?
      !!ii = n/nproc
      !!np = ii
      !!is = iproc*ii
      !!ii = n - nproc*ii
      !!if (iproc<ii) then
      !!    np = np + 1
      !!end if
      !!is = is + min(iproc,ii)
      call distribute_on_tasks(n, iproc, nproc, np, is)

      is = is + isx - 1 !shift of the starting index
      !check
      ii = np
      call mpiallred(ii, 1, mpi_sum, comm=comm)
      if (ii/=n) then
          call f_err_throw('wrong partition of n')
      end if
      iimin = 1 + is
      call mpiallred(iimin, 1, mpi_min, comm=comm)
      if (iimin/=isx) then
          call f_err_throw('wrong starting index')
      end if
      iimax = np + is
      call mpiallred(iimax, 1, mpi_max, comm=comm)
      if (iimax/=iex) then
          call f_err_throw('wrong ending index')
      end if

      nthread = 1
      !$ nthread = omp_get_max_threads()
      max_error_arr = f_malloc(0.to.nthread,id='max_error_arr')
      x_max_error_arr = f_malloc(0.to.nthread,id='x_max_error_arr')

      max_error_arr(:) = 0.d0
      mean_error = 0.d0
      val_chebyshev1 = 0.5d0*coeff(1)

      ithread = 0
      !$omp parallel if (np>1 .and. np*npl>1000) &
      !$omp default(none) &
      !$omp shared(np, is, h, sigma, tau, val_chebyshev1, coeff, npl, mean_error, max_error_arr, x_max_error_arr) &
      !$omp private(i, ii, x, xx, val_chebyshev, xxm2, xxm1, ipl, xxx, val_function, error) &
      !$omp firstprivate(ithread)
      !$ ithread = omp_get_thread_num()
      !$omp do reduction(+: mean_error)
      do i=1,np
          ii = i + is
          x = real(ii,kind=mp)*h
          xx = sigma*(x-tau)
          val_chebyshev = val_chebyshev1 + coeff(2)*xx
          xxm2 = 1.d0
          xxm1 = xx
          xx=2.d0*xx
          do ipl=3,npl
              xxx = xx*xxm1 - xxm2
              val_chebyshev = val_chebyshev + coeff(ipl)*xxx
              xxm2 = xxm1
              xxm1 = xxx
          end do
          val_function = func(x)
          error = abs(val_chebyshev-val_function)
          if (error>max_error_arr(ithread)) then
              max_error_arr(ithread) = error
              x_max_error_arr(ithread) = x
          end if
          mean_error = mean_error + error
          !if (abs(bounds(1)-0.15d0)<1.d-1 .and. abs(bounds(2)-30.0d0)<1.d-1 .and. npl==100) then
          !    write(*,*) 'x, val_chebyshev, val_function, max_error', x, val_chebyshev, val_function, max_error
          !end if
          !write(*,*) 'x, val_chebyshev, exp(x-bounds(2))', x, val_chebyshev, exp(x-bounds(2))
      end do
      !$omp end do
      !$omp end parallel
      !write(*,*) 'max_error',max_error
      mean_error = mean_error/real(iex-isx+1,kind=mp)

      ! Get the maximum among the OpenMP threads
      max_error = 0.d0
      do ithread=0,nthread
          if (max_error_arr(ithread)>max_error) then
              max_error = max_error_arr(ithread)
              x_max_error = x_max_error_arr(ithread)
          end if
      end do
      call f_free(max_error_arr)
      call f_free(x_max_error_arr)

      ! Communicate the results... for the maximum in an array since also the position is required
      call mpiallred(mean_error, 1, mpi_sum, comm=comm)
      max_errors = f_malloc0((/1.to.2,0.to.nproc-1/),id='max_errors')
      max_errors(1,iproc) = max_error
      max_errors(2,iproc) = x_max_error
      call mpiallred(max_errors, mpi_sum, comm=comm)
      max_error = 0.d0
      do jproc=0,nproc-1
          if (max_errors(1,jproc)>max_error) then
              max_error = max_errors(1,jproc)
              x_max_error = max_errors(2,jproc)
          end if
      end do
      call f_free(max_errors)

      call f_release_routine()

    end subroutine accuracy_of_chebyshev_expansion


    !!pure function x_power(x, power)
    !!  implicit none
    !!  real(kind=mp),intent(in) :: x
    !!  real(kind=mp),intent(in) :: power
    !!  real(kind=mp) :: x_power
    !!  x_power = x**power
    !!end function x_power


    subroutine get_chebyshev_polynomials(iproc, nproc, comm, itype, foe_verbosity, npl, smatm, smatl, &
               ham_, workarr_compr, foe_obj, chebyshev_polynomials, ispin, eval_bounds_ok, &
               hamscal_compr, scale_factor, shift_value, smats, ovrlp_, ovrlp_minus_one_half)
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, compress_matrix_distributed_wrapper, &
                              trace_sparse!, max_asymmetry_of_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, itype, ispin
      integer,intent(in) :: foe_verbosity
      integer,intent(in) :: npl
      type(sparse_matrix),intent(in) :: smatm, smatl
      type(matrices),intent(in) :: ham_
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(inout) :: workarr_compr
      type(foe_data),intent(inout) :: foe_obj
      real(kind=mp),dimension(:,:),pointer,intent(inout) :: chebyshev_polynomials
      logical,dimension(2),intent(out) :: eval_bounds_ok
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(out) :: hamscal_compr
      real(kind=mp),intent(out) :: scale_factor, shift_value
      type(sparse_matrix),intent(in),optional :: smats
      type(matrices),intent(in),optional :: ovrlp_
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(in),optional :: ovrlp_minus_one_half

      ! Local variables
      integer :: jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
      integer,parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),allocatable :: cc, cc_check
      real(kind=mp),dimension(:,:),allocatable ::  fermip_check
      real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp) :: anoise, sumn, sumn_check, charge_diff, ef_interpol, ddot
      real(kind=mp) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=mp) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible, with_overlap
      logical,dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: efarr, sumnarr, allredarr
      real(kind=mp),dimension(:),allocatable :: fermi_check_compr
      real(kind=mp),dimension(4,4) :: interpol_matrix
      real(kind=mp),dimension(4) :: interpol_vector
      real(kind=mp),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: bisection_bounds_ok
      real(kind=mp) :: temp_multiplicator, ebs_check, ef, ebsp
      integer :: irow, icol, itemp, iflag,info, i, j, itg, ncount, istl, ists, isshift, imshift
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
      real(kind=mp),parameter :: FSCALE_LOWER_LIMIT=5.d-3
      real(kind=mp),parameter :: FSCALE_UPPER_LIMIT=5.d-2
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
      real(kind=mp),parameter :: CHECK_RATIO=1.25d0
      integer,parameter :: NPL_MIN=100
      !!type(matrices) :: inv_ovrlp
      integer,parameter :: NTEMP_ACCURATE=4
      integer,parameter :: NTEMP_FAST=1
      real(kind=mp) :: degree_multiplicator, x_max_error, max_error, x_max_error_check, max_error_check
      real(kind=mp) :: mean_error, mean_error_check
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      real(kind=mp),dimension(2) :: temparr
      real(kind=mp),dimension(:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:),allocatable :: fermi_new, fermi_check_new, fermi_small_new
      integer :: iline, icolumn, icalc



      call f_routine(id='get_chebyshev_polynomials')

      !if (iproc==0) call yaml_comment('get Chebyshev polynomials',hfill='~')


      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      ! Check the arguments
      select case (itype)
      case (1) !standard eigenvalue problem, i.e. the overlap matrix is the identity and is not required
          with_overlap = .false.
      case (2) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(smats)) call f_err_throw('smats not present')
          if (.not.present(ovrlp_)) call f_err_throw('ovrlp_ not present')
          if (.not.present(ovrlp_minus_one_half)) call f_err_throw('ovrlp_minus_one_half not present')
          isshift = (ispin-1)*smats%nvctrp_tg
          with_overlap = .true.
      case default
          call f_err_throw('wrong value for itype')
      end select

      imshift = (ispin-1)*smatm%nvctrp_tg

      evbounds_shrinked=.false.



      penalty_ev_new = f_malloc((/smatl%smmm%nvctrp/),id='penalty_ev_new')
      fermi_new = f_malloc((/smatl%smmm%nvctrp/),id='fermi_new')


      !hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')


      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm


      ! Fake allocation, will be modified later
      chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,1/),id='chebyshev_polynomials')


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

      !!fscale_new=1.d100




    !ispin = 1




          !!fscale_new = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")


              !fscale = fscale_new
              !fscale = max(fscale,FSCALE_LOWER_LIMIT)
              !fscale = min(fscale,FSCALE_UPPER_LIMIT)

              evlow_old=1.d100
              evhigh_old=-1.d100

              !if (iproc==0) then
              !    call yaml_map('decay length of error function',fscale,fmt='(es10.3)')
              !    call yaml_map('decay length multiplicator',temp_multiplicator,fmt='(es10.3)')
              !    call yaml_map('polynomial degree multiplicator',degree_multiplicator,fmt='(es10.3)')
              !end if


                  ! Don't let this value become too small.
                  call foe_data_set_real(foe_obj, &
                       "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",ispin),1.d-4), &
                       ispin)

                  efarr(1)=foe_data_get_real(foe_obj,"ef",ispin)-foe_data_get_real(foe_obj,"bisection_shift",ispin)
                  efarr(2)=foe_data_get_real(foe_obj,"ef",ispin)+foe_data_get_real(foe_obj,"bisection_shift",ispin)

                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
                  !call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)

                  adjust_lower_bound=.true.
                  adjust_upper_bound=.true.




                  it=0
                  eval_bounds_ok=.false.
                  bisection_bounds_ok=.false.
                  !main_loop: do

                      it=it+1

                      !if (iproc==0) then
                      !    call yaml_newline()
                      !    call yaml_sequence(advance='no')
                      !    call yaml_mapping_open(flow=.true.)
                      !    if (foe_verbosity>=1) call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
                      !end if


                      ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
                      !if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
                      !    foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
                          !!call scale_and_shift_hamiltonian()
                      if (with_overlap) then
                          call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
                               smatm, ham_, i1shift=imshift, &
                               smat2=smats, mat2=ovrlp_, i2shift=isshift, &
                               matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
                      else
                          call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
                               smatm, ham_, i1shift=imshift, &
                               matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
                      end if
                      !!write(*,*) 'scale_factor, shift_value, sum(hamscal_compr), sum(ham_%matrix_compr)', &
                      !!            scale_factor, shift_value, sum(hamscal_compr), sum(ham_%matrix_compr)
                          calculate_SHS=.true.
                      !else
                      !    calculate_SHS=.false.
                      !end if
                      evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
                      evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)


                      ! Determine the degree of the polynomial
                      !!npl=nint(degree_multiplicator* &
                      !!    (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))/fscale)
                      !!npl=max(npl,NPL_MIN)
                      !!npl_boundaries = nint(degree_multiplicator* &
                      !!    (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)) &
                      !!        /foe_data_get_real(foe_obj,"fscale_lowerbound")) ! max polynomial degree for given eigenvalue boundaries
                      !!if (npl>npl_boundaries) then
                      !!    npl=npl_boundaries
                      !!    if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
                      !!    if (iproc==0) write(*,*) 'STOP SINCE THIS WILL CREATE PROBLEMS WITH NPL_CHECK'
                      !!    stop
                      !!end if
                      if (npl>nplx) stop 'npl>nplx'

                      ! Array the holds the Chebyshev polynomials. Needs to be recalculated
                      ! every time the Hamiltonian has been modified.
                      if (calculate_SHS) then
                          call f_free_ptr(chebyshev_polynomials)
                          chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,npl/),id='chebyshev_polynomials')
                      end if

                      !!if (iproc==0) then
                      !!    if (foe_verbosity>=1) then
                      !!        call yaml_map('eval bounds',&
                      !!             (/foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!             foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                      !!    else
                      !!        call yaml_map('eval bounds',&
                      !!             (/foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!             foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                      !!    end if
                      !!    call yaml_map('pol deg',npl,fmt='(i0)')
                      !!    if (foe_verbosity>=1) call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es16.9)')
                      !!end if


                      cc = f_malloc((/npl,2,1/),id='cc')

                      !!if (foe_data_get_real(foe_obj,"evlow",ispin)>=0.d0) then
                      !!    call f_err_throw('Lowest eigenvalue must be negative')
                      !!end if
                      !!if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      !!    call f_err_throw('Highest eigenvalue must be positive')
                      !!end if

                      !call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call f_timing(TCAT_CME_AUXILIARY,'OF')
                      !call timing(iproc, 'chebyshev_coef', 'ON')
                      call f_timing(TCAT_CME_COEFFICIENTS,'ON')

                      !!if (foe_data_get_real(foe_obj,"tmprtr")/=0.d0) call f_err_throw('tmprtr must be zero')
                      !!call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale)
                      !!call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!     foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,1), &
                      !!     x_max_error, max_error, mean_error)
                      cc = 0.d0
                      call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                           muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,2,1), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      !!do ipl=1,npl
                      !!   cc(ipl,3,1) = -cc(ipl,2,1)
                      !!end do
                      call evnoise(npl, cc(1,2,1), foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
                      !write(*,*) 'ef', foe_data_get_real(foe_obj,"ef",ispin)


                      !if (iproc==0 .and. foe_verbosity>=1) then
                      !    call yaml_newline()
                      !    call yaml_mapping_open('accuracy (x, max err, mean err)')
                      !    call yaml_map('main',(/x_max_error,max_error,mean_error/),fmt='(es9.2)')
                      !    call yaml_map('check',(/x_max_error_check,max_error_check,max_error/),fmt='(es9.2)')
                      !    call yaml_mapping_close()
                      !    call yaml_newline()
                      !end if

                      !call timing(iproc, 'chebyshev_coef', 'OF')
                      call f_timing(TCAT_CME_COEFFICIENTS,'OF')
                      !call timing(iproc, 'FOE_auxiliary ', 'ON')
                      call f_timing(TCAT_CME_AUXILIARY,'ON')


                      if (smatl%nspin==1) then
                          do ipl=1,npl
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,2,1)=2.d0*cc(ipl,2,1)
                              !!cc(ipl,3,1)=2.d0*cc(ipl,3,1)
                          end do
                      end if


                      !call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call f_timing(TCAT_CME_AUXILIARY,'OF')

                      emergency_stop=.false.
                          ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
                          if (with_overlap) then
                              !!call max_asymmetry_of_matrix(iproc, nproc, comm, &
                              !!     smatl, hamscal_compr, tt)
                              !!if (iproc==0) call yaml_map('max assymetry of Hscal',tt)
                              !!call max_asymmetry_of_matrix(iproc, nproc, comm, &
                              !!     smatl, ovrlp_minus_one_half, tt)
                              !!if (iproc==0) call yaml_map('max assymetry of inv',tt)
                              call chebyshev_clean(iproc, nproc, npl, cc, &
                                   smatl, hamscal_compr, &
                                   .true., workarr_compr, &
                                   nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                                   emergency_stop, invovrlp_compr=ovrlp_minus_one_half)
                          else
                              call chebyshev_clean(iproc, nproc, npl, cc, &
                                   smatl, hamscal_compr, &
                                   .false., workarr_compr, &
                                   nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                                   emergency_stop)
                          end if
                          !write(*,*) 'sum(hamscal_compr)',sum(hamscal_compr)
                          !write(*,*) 'npl, sum(cc(:,1,1)), sum(chebyshev_polynomials)', &
                          !            npl, sum(cc(:,1,1)), sum(chebyshev_polynomials)



                          !!call transform_sparsity_pattern(smatl%nfvctr, &
                          !!     smatl%smmm%nvctrp_mm, smatl%smmm%isvctr_mm, &
                          !!     smatl%nseg, smatl%keyv, smatl%keyg, smatl%smmm%line_and_column_mm, &
                          !!     smatl%smmm%nvctrp, smatl%smmm%isvctr, &
                          !!     smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, &
                          !!     smatl%smmm%istsegline, 'large_to_small', fermi_small_new, fermi_new)

                      !call timing(iproc, 'FOE_auxiliary ', 'ON')
                      call f_timing(TCAT_CME_AUXILIARY,'ON')


                      restart=.false.

                      ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                      ! (otherwise this has already been checked in the previous iteration).
                      call check_eigenvalue_spectrum_new(iproc, nproc, comm, smatl, ispin, &
                            0, 1.0d0, 1.0d0, penalty_ev_new, anoise, .false., emergency_stop, &
                            foe_obj, restart, eval_bounds_ok, foe_verbosity)

                      call f_free(cc)

                      if (restart) then
                          if(evbounds_shrinked) then
                              ! this shrink was not good, increase the saturation counter
                              call foe_data_set_int(foe_obj,"evboundsshrink_isatur", &
                                   foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
                          end if
                          call foe_data_set_int(foe_obj,"evbounds_isatur",0)
                          !if (iproc==0) then
                          !    if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                          !         (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          !    call yaml_mapping_close()
                          !    !call bigdft_utils_flush(unit=6)
                          !end if
                          !cycle
                      end if

                      ! eigenvalue bounds ok
                      if (calculate_SHS) then
                          call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
                      end if

                      !!!call calculate_trace_distributed(kernel_%matrixp, sumn)
                      !!call calculate_trace_distributed_new(fermi_small_new, sumn)


                      !if (all(eval_bounds_ok) .and. all(bisection_bounds_ok)) then
                      !    ! Print these informations already now if all entries are true.
                      !    if (iproc==0) then
                      !        if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                      !             (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      !    end if
                      !end if

                      ! If we are here exit the main loop
                      !call yaml_mapping_close()
                      !exit main_loop


                  !end do main_loop



      !if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


      !!call f_free(chebyshev_polynomials)
      !call f_free(hamscal_compr)

      call f_free(penalty_ev_new)
      call f_free(fermi_new)

      !write(*,*) 'end get_chebyshev_polynomials: ef', foe_data_get_real(foe_obj,"ef",ispin)

      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      call f_release_routine()


    end subroutine get_chebyshev_polynomials


    subroutine find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
               foe_verbosity, label, smatl, ispin, foe_obj, kernel_)
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, compress_matrix_distributed_wrapper, &
                              trace_sparse, max_asymmetry_of_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      !!use foe_common, only: scale_and_shift_matrix, evnoise, &
      !!                      check_eigenvalue_spectrum_new, retransform_ext, get_chebyshev_expansion_coefficients
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, npl, ispin
      type(sparse_matrix),intent(in) :: smatl
      real(kind=mp),dimension(smatl%smmm%nvctrp_mm,npl),intent(in) :: chebyshev_polynomials
      integer,intent(in) :: foe_verbosity
      character(len=*),intent(in) :: label
      type(foe_data),intent(inout) :: foe_obj
      type(matrices),intent(inout) :: kernel_

      ! Local variables
      integer :: jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries, ilshift
      integer,parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),allocatable :: cc, cc_check
      real(kind=mp),dimension(:,:),allocatable :: fermip_check
      real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp) :: anoise, scale_factor, shift_value, sumn, sumn_check, charge_diff, ef_interpol, ddot
      real(kind=mp) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=mp) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible
      logical,dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: efarr, sumnarr, allredarr
      real(kind=mp),dimension(:),allocatable :: hamscal_compr, fermi_check_compr
      real(kind=mp),dimension(4,4) :: interpol_matrix
      real(kind=mp),dimension(4) :: interpol_vector
      real(kind=mp),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
      real(kind=mp) :: temp_multiplicator, ebs_check, ef, ebsp
      integer :: irow, icol, itemp, iflag,info, isshift, imshift, ilshift2, i, j, itg, ncount, istl, ists
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit
      real(kind=mp),parameter :: FSCALE_LOWER_LIMIT=5.d-3
      real(kind=mp),parameter :: FSCALE_UPPER_LIMIT=5.d-2
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
      real(kind=mp),parameter :: CHECK_RATIO=1.25d0
      integer,parameter :: NPL_MIN=100
      !!type(matrices) :: inv_ovrlp
      integer,parameter :: NTEMP_ACCURATE=4
      integer,parameter :: NTEMP_FAST=1
      real(kind=mp) :: degree_multiplicator, x_max_error, max_error, x_max_error_check, max_error_check
      real(kind=mp) :: mean_error, mean_error_check
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      real(kind=mp),dimension(2) :: temparr
      real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:),allocatable :: fermi_new, fermi_check_new, fermi_small_new
      integer :: iline, icolumn, icalc



      call f_routine(id='find_fermi_level')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      !if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')


      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      evbounds_shrinked=.false.


      fermi_small_new = f_malloc(max(smatl%smmm%nvctrp_mm,1),id='fermi_small_new')




      !hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')


      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm



      ntemp = NTEMP_ACCURATE
      degree_multiplicator = DEGREE_MULTIPLICATOR_ACCURATE
      temp_multiplicator = TEMP_MULTIPLICATOR_ACCURATE

      fscale_new=1.d100


      !spin_loop: do ispin=1,smatl%nspin

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

              evlow_old=1.d100
              evhigh_old=-1.d100


                  ! Don't let this value become too small.
                  call foe_data_set_real(foe_obj, &
                       "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",ispin),1.d-4), &
                       ispin)

                  efarr(1)=foe_data_get_real(foe_obj,"ef",ispin)-foe_data_get_real(foe_obj,"bisection_shift",ispin)
                  efarr(2)=foe_data_get_real(foe_obj,"ef",ispin)+foe_data_get_real(foe_obj,"bisection_shift",ispin)
                  !write(*,*) 'ef, efarr', foe_data_get_real(foe_obj,"ef",ispin), efarr

                  sumnarr(1)=0.d0
                  sumnarr(2)=1.d100
                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), verbosity=foe_verbosity)
                  call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)


                  if (iproc==0 .and. foe_verbosity>0) then
                      !if (foe_verbosity>=1) then
                      !    call yaml_sequence_open('FOE to determine density kernel',&
                      !         label='it_foe'//trim(label)//'-'//&
                      !         trim(adjustl(yaml_toa(itemp,fmt='(i2.2)')))//'-'//&
                      !         trim(adjustl(yaml_toa(ispin,fmt='(i2.2)'))))
                      !else
                      !    call yaml_sequence_open('FOE to determine density kernel')
                      !    if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='-')
                      !end if
                      call yaml_sequence_open('determine Fermi energy')
                  end if



                  it=0
                  eval_bounds_ok=.true.
                  bisection_bounds_ok=.false.
                  main_loop: do

                      it=it+1

                      if (iproc==0 .and. foe_verbosity>0) then
                          call yaml_newline()
                          call yaml_sequence(advance='no')
                          call yaml_mapping_open(flow=.true.)
                          if (foe_verbosity>=1) call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
                      end if

                      if (iproc==0) then
                          !if (foe_verbosity>=1) then
                          !    call yaml_map('bisec/eval bounds',&
                          !         (/fermilevel_get_real(f,"efarr(1)"),fermilevel_get_real(f,"efarr(2)"),&
                          !         foe_data_get_real(foe_obj,"evlow",ispin), &
                          !         foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          !else
                          !    call yaml_map('eval bounds',&
                          !         (/foe_data_get_real(foe_obj,"evlow",ispin), &
                          !         foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          !end if
                          !call yaml_map('pol deg',npl,fmt='(i0)')
                          !if (foe_verbosity>=1) call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es16.9)')
                      end if


                      cc = f_malloc((/npl,1,3/),id='cc')

                      !call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call f_timing(TCAT_CME_AUXILIARY,'OF')
                      !call timing(iproc, 'chebyshev_coef', 'ON')
                      call f_timing(TCAT_CME_COEFFICIENTS,'ON')

                      call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale)
                      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,1), &
                           x_max_error, max_error, mean_error)
                      call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                           muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,2), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      do ipl=1,npl
                         cc(ipl,1,3) = -cc(ipl,1,2)
                      end do
                      call evnoise(npl, cc(1,1,2), foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), anoise)


                      !if (iproc==0 .and. foe_verbosity>=1) then
                      !    call yaml_newline()
                      !    call yaml_mapping_open('accuracy (x, max err, mean err)')
                      !    call yaml_map('main',(/x_max_error,max_error,mean_error/),fmt='(es9.2)')
                      !    !call yaml_map('check',(/x_max_error_check,max_error_check,max_error/),fmt='(es9.2)')
                      !    call yaml_mapping_close()
                      !    call yaml_newline()
                      !end if

                      !call timing(iproc, 'chebyshev_coef', 'OF')
                      call f_timing(TCAT_CME_COEFFICIENTS,'OF')
                      !call timing(iproc, 'FOE_auxiliary ', 'ON')
                      call f_timing(TCAT_CME_AUXILIARY,'ON')


                      if (smatl%nspin==1) then
                          do ipl=1,npl
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,1,2)=2.d0*cc(ipl,1,2)
                              cc(ipl,1,3)=2.d0*cc(ipl,1,3)
                          end do
                      end if


                      !call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call f_timing(TCAT_CME_AUXILIARY,'OF')

                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                               smatl%nfvctr, smatl%smmm%nfvctrp, &
                              smatl, chebyshev_polynomials, 1, cc, fermi_small_new)


                      !call timing(iproc, 'FOE_auxiliary ', 'ON')
                      call f_timing(TCAT_CME_AUXILIARY,'ON')


                      call f_free(cc)


                      call calculate_trace_distributed_new(iproc, nproc, comm, smatl, fermi_small_new, sumn)
                      !write(*,*) 'sumn',sumn


                      if (all(eval_bounds_ok) .and. all(bisection_bounds_ok)) then
                          ! Print these informations already now if all entries are true.
                          if (iproc==0) then
                              !!if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                              !!     (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          end if
                      end if

                      if (iproc==0 .and. foe_verbosity>0) then
                          !call yaml_newline()
                          !call yaml_map('iter',it)
                          call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es13.6)')
                          !call yaml_map('bisec bounds ok',&
                          !     (/bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          call yaml_map('Tr(K)',sumn,fmt='(es14.7)')
                          call yaml_map('D Tr(K)',sumn-foe_data_get_real(foe_obj,"charge",ispin),fmt='(es9.2)')
                      end if
                      call determine_fermi_level(iproc, f, sumn, ef, info)
                      bisection_bounds_ok(1) = fermilevel_get_logical(f,"bisection_bounds_ok(1)")
                      bisection_bounds_ok(2) = fermilevel_get_logical(f,"bisection_bounds_ok(2)")

                      charge_diff = sumn-foe_data_get_real(foe_obj,"charge",ispin)

                      ! If the charge difference is smaller than the threshold, there is no need to cycle even though we
                      ! are in principle still looking for the bisection bounds.
                      if (info<0 .and. abs(charge_diff)>=charge_tolerance) then
                          if (iproc==0 .and. foe_verbosity>0) then
                              !if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                              !     (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
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



                      ef_old=foe_data_get_real(foe_obj,"ef",ispin)
                      sumn_old=sumn




                      if (iproc==0 .and. foe_verbosity>0) then
                          call yaml_mapping_close()
                      end if

                      if (abs(charge_diff)<charge_tolerance) then
                          if (iproc==0 .and. foe_verbosity>0) call yaml_sequence_close()
                          diff=0.d0

                          if (nproc > 1) then
                              call mpiallred(diff, 1, mpi_sum, comm=comm)
                          end if

                          diff=sqrt(diff)
                          !if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')
                          exit
                      end if

                  end do main_loop

                  if (iproc==0) then
                      call yaml_mapping_open('summary',flow=.true.)
                      call yaml_map('nit',it)
                      call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es13.6)')
                      call yaml_map('Tr(K)',sumn,fmt='(es14.7)')
                      call yaml_map('D Tr(K)',sumn-foe_data_get_real(foe_obj,"charge",ispin),fmt='(es9.2)')
                      call yaml_mapping_close()
                  end if
            
             call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                  fermi_small_new, &
                  kernel_%matrix_compr(ilshift+1:))

      !end do spin_loop



      degree_sufficient=.true.

      !if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


      call f_free(fermi_small_new)

      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      call f_release_routine()




    end subroutine find_fermi_level


    subroutine calculate_trace_distributed_new(iproc, nproc, comm, smatl, matrixp, trace)
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smatl
      real(kind=mp),dimension(smatl%smmm%nvctrp_mm),intent(in) :: matrixp
      real(kind=mp),intent(out) :: trace
      integer :: i, ii, iline, icolumn

      call f_routine(id='calculate_trace_distributed_new')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

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
          call mpiallred(trace, 1, mpi_sum, comm=comm)
      end if

      call f_release_routine()
    end subroutine calculate_trace_distributed_new


    !> Determine the polynomial degree which yields the desired precision
    subroutine get_polynomial_degree(iproc, nproc, comm, ispin, ncalc, fun, foe_obj, &
               npl_min, npl_max, npl_stride, max_polynomial_degree, verbosity, npl, cc, &
               max_error, x_max_error, mean_error, anoise, &
               ex, ef, fscale)
      use foe_base, only: foe_data, foe_data_get_real
      use yaml_output
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin, ncalc, fun, verbosity
      integer,intent(in) :: npl_min, npl_max, npl_stride
      type(foe_data),intent(in) :: foe_obj
      real(kind=mp),intent(in) :: max_polynomial_degree
      integer,intent(out) :: npl
      real(kind=mp),dimension(:,:,:),pointer,intent(inout) :: cc
      real(kind=mp),dimension(ncalc),intent(out) :: max_error, x_max_error, mean_error
      real(kind=mp),intent(out) :: anoise
      real(kind=mp),dimension(ncalc),intent(in),optional :: ex, ef, fscale

      ! Local variables
      integer :: ipl, icalc, j, jpl
      logical :: error_ok, found_degree
      real(kind=mp),dimension(:,:,:),allocatable :: cc_trial
      real(kind=mp) :: x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction

      call f_routine(id='get_polynomial_degree')

      ! Check the arguments
      select case (fun)
      case (FUNCTION_POLYNOMIAL)
          if (.not. present(ex)) call f_err_throw("arguments 'ex' is not present")
      case (FUNCTION_ERRORFUNCTION)
          if (.not. present(ef)) call f_err_throw("arguments 'ef' is not present")
          if (.not. present(fscale)) call f_err_throw("arguments 'fscale' is not present")
          !write(*,*) 'iproc, ef, fscale, evlow, evhigh', &
          !    iproc, ef, fscale, foe_data_get_real(foe_obj,"evlow",ispin), foe_data_get_real(foe_obj,"evhigh",ispin)
      case default
          call f_err_throw("wrong value of argument 'fun'")
      end select

      !max_error = f_malloc(ncalc,id='max_error')
      !x_max_error = f_malloc(ncalc,id='x_max_error')
      !mean_error = f_malloc(ncalc,id='mean_error')

      if (npl_min<3) then
          call f_err_throw('npl_min must be at least 3')
      end if
      if (npl_min>npl_max) then
          call f_err_throw('npl_min must be smaller or equal than npl_max')
      end if

      if (iproc==0 .and. verbosity>0) then
          call yaml_sequence_open('Determine polynomial degree')
      end if

      cc_trial = f_malloc0((/npl_max,ncalc,3/),id='cc_trial')

      found_degree = .false.
      degree_loop: do ipl=npl_min,npl_max,npl_stride

          if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
              stop 'ERROR: highest eigenvalue must be positive'
          end if

          !call timing(iproc, 'FOE_auxiliary ', 'OF')
          call f_timing(TCAT_CME_AUXILIARY,'OF')
          !call timing(iproc, 'chebyshev_coef', 'ON')
          call f_timing(TCAT_CME_COEFFICIENTS,'ON')

          do icalc=1,ncalc
              select case (fun)
              case (FUNCTION_POLYNOMIAL)
                  call func_set(FUNCTION_POLYNOMIAL, powerx=ex(icalc))
              case (FUNCTION_ERRORFUNCTION)
                  call func_set(FUNCTION_ERRORFUNCTION, efx=ef(icalc), fscalex=fscale(icalc))
              end select
              call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                   foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,icalc,1), &
                   x_max_error(icalc), max_error(icalc), mean_error(icalc))
              !write(*,*) 'icalc, sum(cc_trial(:,icalc,1))', icalc, sum(cc_trial(:,icalc,1)), ex(icalc)
          end do

          !call timing(iproc, 'chebyshev_coef', 'OF')
          call f_timing(TCAT_CME_COEFFICIENTS,'OF')
          !call timing(iproc, 'FOE_auxiliary ', 'ON')
          call f_timing(TCAT_CME_AUXILIARY,'ON')

          if (iproc==0 .and. verbosity>0) then
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
              do icalc=1,ncalc
                  call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                       muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                  call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,icalc,2), &
                       x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction)
                  do jpl=1,ipl
                      cc_trial(jpl,icalc,3) = -cc_trial(jpl,icalc,2)
                  end do
                  if (max_error_penaltyfunction>1.d-2) then
                      error_ok = .false.
                  end if
              end do
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

      if (iproc==0 .and. verbosity>0) then
          call yaml_sequence_close()
      end if

      cc = f_malloc_ptr((/npl,ncalc,3/),id='cc')
      do j=1,3
          do icalc=1,ncalc
              do ipl=1,npl
                  cc(ipl,icalc,j)=cc_trial(ipl,icalc,j)
                  !write(*,*) 'icalc, ipl, cc(ipl,icalc,1)', icalc, ipl, cc(ipl,icalc,1)
              end do
          end do
      end do
      call f_free(cc_trial)
      !call f_free(mean_error)
      !call f_free(max_error)
      !call f_free(x_max_error)

      call f_release_routine

    end subroutine get_polynomial_degree


    subroutine chebyshev_coefficients_init_parallelization(iproc, nproc, comm, n, np, is)
      implicit none
      ! Caling arguments
      integer,intent(in) :: iproc, nproc, comm, n
      integer,intent(out) :: np, is

      ! Local variables
      integer :: ii

      call f_routine(id='chebyshev_coefficients_init_parallelization')

      ii = n/nproc
      np = ii
      is = iproc*ii
      ii = n - nproc*ii
      if (iproc<ii) then
          np = np + 1
      end if
      is = is + min(iproc,ii)
      !check
      ii = np
      call mpiallred(ii, 1, mpi_sum, comm=comm)
      if (ii/=n) then
          call f_err_throw('wrong partition of n')
      end if

      call f_release_routine()

    end subroutine chebyshev_coefficients_init_parallelization


    subroutine chebyshev_coefficients_calculate(n, a, b, np, is, func, cc)
      implicit none

      ! Calling arguments
      integer,intent(in) :: n, np, is
      real(kind=mp),intent(in) :: a, b
      real(kind=mp),external :: func
      real(kind=mp),dimension(n),intent(out) :: cc

      ! Local variables
      integer :: k, j, ii, jj
      real(kind=mp) :: bma, bpa, y, arg, fac, tt, one_over_n
      real(kind=mp),dimension(:),allocatable :: cf

      call f_routine(id='chebyshev_coefficients_calculate')


      call f_zero(cc)
      cf = f_malloc0(n,id='cf')

      bma=0.5d0*(b-a)
      bpa=0.5d0*(b+a)
      fac=2.d0/real(n,kind=mp)
      one_over_n = 1.d0/real(n,kind=mp)
      !$omp parallel default(none) shared(bma,bpa,fac,n,cf,cc,is,np,tt,one_over_n) &
      !$omp private(k,y,arg,j,jj)
      !$omp do
      do k=1,n
          y=cos(pi*(real(k,kind=mp)-0.5d0)*(one_over_n))
          arg=y*bma+bpa
          cf(k)=func(arg)
      end do
      !$omp end do
      !$omp end parallel

      do j=1,np
          jj = j + is
          tt=0.d0
          !$omp parallel do default(none) shared(n,cf,jj,one_over_n) private(k) reduction(+:tt)
          do  k=1,n
              tt=tt+cf(k)*cos((pi*real(jj-1,kind=mp))*((real(k,kind=mp)-0.5d0)*(one_over_n)))
          end do
          !$omp end parallel do
          cc(jj)=fac*tt
      end do
      call f_free(cf)

      call f_release_routine()

    end subroutine chebyshev_coefficients_calculate


    ! This routine is basically just here to get the profiling...
    subroutine chebyshev_coefficients_communicate(comm, n, cc)
      implicit none

      ! Calling arguments
      integer,intent(in) :: comm, n
      real(kind=mp),dimension(n),intent(inout) :: cc

      call f_routine(id='chebyshev_coefficients_communicate')

      call mpiallred(cc, mpi_sum, comm=comm)

      call f_release_routine()

    end subroutine chebyshev_coefficients_communicate


    ! This routine is basically just here to get the profiling...
    subroutine penalty_communicate(nproc, comm, penalty)
      implicit none

      ! Calling arguments
      integer,intent(in) :: nproc, comm
      real(mp),intent(inout) :: penalty

      call f_routine(id='penalty_communicate')

      if (nproc > 1) then
          call mpiallred(penalty, 1, mpi_sum, comm=comm)
      end if

      call f_release_routine()

    end subroutine penalty_communicate



    subroutine get_bounds_and_polynomials(iproc, nproc, comm, itype, ispin, npl_max, npl_stride, ncalc, func_name, &
               do_scaling, bounds_factor_low, bounds_factor_up, foe_verbosity, &
               smatm, smatl, ham_, foe_obj, npl_min, workarr_compr, chebyshev_polynomials, &
               npl, scale_factor, shift_value, hamscal_compr, &
               smats, ovrlp_, ovrlp_minus_one_half_, efarr, fscale_arr, ex, &
               scaling_factor_low, scaling_factor_up, eval_multiplicator, eval_multiplicator_total, cc, max_errorx)
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, itype, ispin, npl_max, npl_stride, ncalc, func_name, foe_verbosity
      type(sparse_matrix),intent(in) :: smatm, smatl
      type(matrices),intent(in) :: ham_
      logical,intent(in) :: do_scaling
      real(mp),intent(in),optional :: bounds_factor_low, bounds_factor_up
      type(foe_data),intent(inout) :: foe_obj
      integer,intent(inout) :: npl_min
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(inout) :: workarr_compr
      real(mp),dimension(:,:),pointer,intent(inout) :: chebyshev_polynomials
      integer,intent(out) :: npl
      real(mp),intent(out) :: scale_factor, shift_value
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(out) :: hamscal_compr
      type(sparse_matrix),intent(in),optional :: smats
      type(matrices),intent(in),optional :: ovrlp_, ovrlp_minus_one_half_
      real(kind=mp),dimension(ncalc),intent(in),optional :: efarr
      real(kind=mp),dimension(ncalc),intent(in),optional :: fscale_arr
      real(kind=mp),dimension(ncalc),intent(in),optional :: ex
      real(mp),intent(in),optional :: scaling_factor_low, scaling_factor_up
      real(mp),intent(inout),optional :: eval_multiplicator, eval_multiplicator_total
      real(kind=mp),dimension(:,:,:),pointer,intent(out),optional :: cc
      real(mp),dimension(ncalc),intent(out),optional :: max_errorx

      ! Local variables
      integer :: ilshift
      real(mp),dimension(:),allocatable :: max_error, x_max_error, mean_error
      real(mp) :: anoise
      real(kind=mp),dimension(:,:,:),pointer :: cc_
      logical,dimension(2) :: eval_bounds_ok
      type(matrices) :: ham_scaled

      call f_routine(id='get_bounds_and_polynomials')

      ! Check the arguments
      select case (itype)
      case (1) !standard eigenvalue problem, i.e. the overlap matrix is the identity and is not required
      case (2) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(smats)) call f_err_throw('smats not present')
          if (.not.present(ovrlp_)) call f_err_throw('ovrlp_ not present')
          if (.not.present(ovrlp_minus_one_half_)) call f_err_throw('ovrlp_minus_one_half_ not present')
      case default
          call f_err_throw('wrong value for itype')
      end select

      select case (func_name)
      case (FUNCTION_ERRORFUNCTION) 
          if (.not.present(efarr)) call f_err_throw('efarr not present')
          if (.not.present(fscale_arr)) call f_err_throw('fscale_arr not present')
      case (FUNCTION_POLYNOMIAL) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(ex)) call f_err_throw('ex not present')
      case default
          call f_err_throw('wrong value for func_name')
      end select

      if (do_scaling) then
          if (.not.present(scaling_factor_low)) call f_err_throw('scaling_factor_low not present')
          if (.not.present(scaling_factor_up)) call f_err_throw('scaling_factor_up not present')
          if (.not.present(eval_multiplicator)) call f_err_throw('eval_multiplicator not present')
          if (.not.present(eval_multiplicator_total)) call f_err_throw('eval_multiplicator_total not present')
      end if

      ilshift = (ispin-1)*smatl%nvctrp_tg
      max_error = f_malloc(ncalc,id='max_error')
      x_max_error = f_malloc(ncalc,id='x_max_error')
      mean_error = f_malloc(ncalc,id='mean_error')

      ham_scaled = matrices_null()
      if (do_scaling) then
          ham_scaled%matrix_compr = sparsematrix_malloc_ptr(smatm, &
              iaction=SPARSE_TASKGROUP, id='ham_scaled%matrix_compr')
          call f_memcpy(src=ham_%matrix_compr,dest=ham_scaled%matrix_compr)
      else
          ham_scaled%matrix_compr => ham_%matrix_compr
      end if


      if (iproc==0 .and. foe_verbosity>0) then
          call yaml_sequence_open('determine eigenvalue bounds')
      end if
      bounds_loop: do
          !efarr(1) = foe_data_get_real(foe_obj,"ef",ispin)
          !fscale_arr(1) = foe_data_get_real(foe_obj,"fscale",ispin)
          if (do_scaling) then
              call dscal(size(ham_scaled%matrix_compr), eval_multiplicator, ham_scaled%matrix_compr(1), 1)
              eval_multiplicator_total = eval_multiplicator_total*eval_multiplicator
          end if

          if (func_name==FUNCTION_ERRORFUNCTION) then
              call get_polynomial_degree(iproc, nproc, comm, ispin, ncalc, FUNCTION_ERRORFUNCTION, foe_obj, &
                   npl_min, npl_max, npl_stride, 1.d-5, 0, npl, cc_, &
                   max_error, x_max_error, mean_error, anoise, &
                   ef=efarr, fscale=fscale_arr)
          else if (func_name==FUNCTION_POLYNOMIAL) then
              call get_polynomial_degree(iproc, nproc, comm, ispin, ncalc, FUNCTION_POLYNOMIAL, foe_obj, &
                   npl_min, npl_max, npl_stride, 1.d-8, 0, npl, cc_, &
                   max_error, x_max_error, mean_error, anoise, &
                   ex=ex)
          end if
          npl_min = npl !to be used to speed up the search for npl in a following iteration in case the temperature must be lowered
          if (iproc==0 .and. foe_verbosity>0) then
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('npl',npl)
              if (do_scaling) call yaml_map('scale',eval_multiplicator_total,fmt='(es9.2)')
              call yaml_map('bounds', &
                   (/foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f7.3)')
          end if

          if (itype==2) then
              call get_chebyshev_polynomials(iproc, nproc, comm, &
                   itype, foe_verbosity, npl, smatm, smatl, &
                   ham_scaled, workarr_compr, foe_obj, &
                   chebyshev_polynomials, ispin, eval_bounds_ok, hamscal_compr, &
                   scale_factor, shift_value, &
                   smats=smats, ovrlp_=ovrlp_, &
                   ovrlp_minus_one_half=ovrlp_minus_one_half_%matrix_compr(ilshift+1:))
          else if (itype==1) then
              call get_chebyshev_polynomials(iproc, nproc, comm, &
                   itype, foe_verbosity, npl, smatm, smatl, &
                   ham_scaled, workarr_compr, foe_obj, &
                   chebyshev_polynomials, ispin, eval_bounds_ok, hamscal_compr, &
                   scale_factor, shift_value)
          end if
          if (iproc==0 .and. foe_verbosity>0) then
              call yaml_map('ok',eval_bounds_ok)
              call yaml_map('exp accur',max_error,fmt='(es8.2)')
              call yaml_mapping_close()
          end if
          if (all(eval_bounds_ok)) then
              exit bounds_loop
          else
              if (.not.eval_bounds_ok(1)) then
                  ! lower bound not ok
                  !!call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*1.2d0,ispin)
                  !!eval_multiplicator = 2.0d0
                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*bounds_factor_low,ispin)
                  if (do_scaling) then
                      eval_multiplicator = scaling_factor_low
                  end if
              else if (.not.eval_bounds_ok(2)) then
                  ! upper bound not ok
                  !!call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                  !!eval_multiplicator = 1.d0/2.0d0
                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*bounds_factor_up,ispin)
                  if (do_scaling) then
                      eval_multiplicator = scaling_factor_up
                  end if
              end if
          end if
          call f_free_ptr(cc_)
          call f_free_ptr(chebyshev_polynomials)
      end do bounds_loop
      if (iproc==0 .and. foe_verbosity>0) then
          call yaml_sequence_close()
      end if

      if (do_scaling) then
          call deallocate_matrices(ham_scaled)
      end if
      if (present(cc)) then
          !f_malloc((/npl,ncalc,3/),id='cc')
           cc = f_malloc_ptr((/npl,ncalc,3/), id='cc')
           call f_memcpy(src=cc_, dest=cc)
      end if
      call f_free_ptr(cc_)
      if (present(max_errorx)) then
          call f_memcpy(src=max_error, dest=max_errorx)
      end if
      call f_free(max_error)
      call f_free(x_max_error)
      call f_free(mean_error)

      call f_release_routine()

    end subroutine get_bounds_and_polynomials

end module foe_common
