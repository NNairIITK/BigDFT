module module_func
  use module_base
  implicit none

  private

  ! Shared variables within the modules
  integer :: ifunc
  real(kind=8) :: power, ef, fscale, beta, mua, mub

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
      real(kind=8),intent(in),optional :: powerx, efx, fscalex, betax, muax, mubx

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
      real(kind=8),intent(in) :: x
      real(kind=8) :: func
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
  use module_defs, only: uninitialized
  use module_base
  use foe_base
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
  public :: check_eigenvalue_spectrum_new
  public :: scale_and_shift_matrix
  public :: retransform_ext
  !!public :: cheb_exp
  public :: init_foe
  public :: get_chebyshev_polynomials
  public :: find_fermi_level


  contains


    subroutine get_chebyshev_expansion_coefficients(iproc, nproc, A, B, N, func, cc, x_max_error,max_error,mean_error)
      use module_base
      use yaml_output
      implicit none
      
      ! Calling arguments
      real(kind=8),intent(in) :: A, B
      integer,intent(in) :: iproc, nproc, n
      real(kind=8),external :: func
      real(8),dimension(n),intent(out) :: cc
      real(kind=8),intent(out) :: x_max_error, max_error, mean_error
    
      ! Local variables
      integer :: k, j, is, np, ii, jj
      real(kind=8) :: bma, bpa, y, arg, fac, tt, one_over_n
      real(kind=8),dimension(:),allocatable :: cf
    
      call f_routine(id='get_chebyshev_expansion_coefficients')

      ! MPI parallelization... maybe only worth for large n?
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
      call mpiallred(ii, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      if (ii/=n) then
          call f_err_throw('wrong partition of n')
      end if
      

      call f_zero(cc)
      cf = f_malloc0(n,id='cf')
    
      bma=0.5d0*(b-a)
      bpa=0.5d0*(b+a)
      fac=2.d0/real(n,kind=8)
      one_over_n = 1.d0/real(n,kind=8)
      !$omp parallel default(none) shared(bma,bpa,fac,n,cf,cc,is,np,tt,one_over_n) &
      !$omp private(k,y,arg,j,jj)
      !$omp do
      do k=1,n
          y=cos(pi*(real(k,kind=8)-0.5d0)*(one_over_n))
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
              tt=tt+cf(k)*cos((pi*real(jj-1,kind=8))*((real(k,kind=8)-0.5d0)*(one_over_n)))
          end do
          !$omp end parallel do
          cc(jj)=fac*tt
      end do

      call mpiallred(cc, mpi_sum, comm=bigdft_mpi%mpi_comm)

      call f_free(cf)
      call accuracy_of_chebyshev_expansion(iproc, nproc, n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
    
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
    !!  real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
    !!  integer,intent(in) :: iproc, nproc, n
    !!  real(8),dimension(n),intent(out) :: cc
    !!  real(kind=8),intent(out) :: x_max_error, max_error, mean_error
    !!
    !!  ! Local variables
    !!  integer :: k, j, is, np, ii, jj
    !!  real(kind=8) :: bma, bpa, y, arg, fac, tt
    !!  real(kind=8),dimension(50000) :: cf
    !!  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
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
    !!!  real(kind=8),intent(in) :: a, b
    !!!  integer,intent(in) :: n
    !!!  real(kind=8),dimension(n,2),intent(out) :: cc
    !!!  real(kind=8),intent(out) :: max_error
    !!!
    !!!  ! Local variables
    !!!  integer :: k, j
    !!!  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
    !!!  real(kind=8) :: tt1, tt2, ttt, y, arg, fac, bma, bpa, x_max, max_err, mean_err
    !!!  real(kind=8),dimension(50000) :: cf
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
      real(kind=8) :: fact, dist, ddx, cent, tt, x
    
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
              real(kind=8) :: ddx, x, tt, err
    
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
            real(kind=8) :: fact, ddx, tt, x
    
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
    
    

    subroutine check_eigenvalue_spectrum_new(nproc, smat_l, ispin, isshift, &
               factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
               emergency_stop, foe_obj, restart, eval_bounds_ok, &
               eval_multiplicator, smat_s, mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat_l
      integer,intent(in) :: nproc, ispin, isshift
      real(kind=8),intent(in) :: factor_high, factor_low, anoise
      !real(kind=8),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
      real(kind=8),dimension(smat_l%smmm%nvctrp,2),intent(in) :: penalty_ev
      logical,intent(in) :: trace_with_overlap
      logical,dimension(2) :: emergency_stop
      type(foe_data),intent(inout) :: foe_obj
      logical,intent(inout) :: restart
      logical,dimension(2),intent(out) :: eval_bounds_ok
      real(kind=8),intent(inout),optional :: eval_multiplicator
      type(sparse_matrix),intent(in),optional :: smat_s
      type(matrices),intent(in),optional :: mat
    
      ! Local variables
      integer :: isegstart, isegend, iseg, ii, jorb, irow, icol, iismall, iel, i, iline, icolumn, ibound
      real(kind=8) :: bound_low, bound_up, tt, noise
      real(kind=8),dimension(2) :: allredarr
    
      call f_routine(id='check_eigenvalue_spectrum_new')

      if (trace_with_overlap) then
          if (.not.present(smat_s) .or. .not.present(mat)) then
              call f_err_throw('not all required optional arguments are present')
          end if
      end if
    
      bound_low=0.d0
      bound_up=0.d0
      do ibound=1,2
          !if (.not.emergency_stop(ibound)) then
              ! The penalty function must be smaller than the noise.
    
              !$omp parallel default(none) &
              !$omp shared(ibound,bound_low, bound_up, smat_l, smat_s, trace_with_overlap, mat, isshift, penalty_ev) &
              !$omp private(i, ii, iline, icolumn, iismall, tt)
              !$omp do reduction(+:bound_low, bound_up)
              do i=1,smat_l%smmm%nvctrp
                  ii = smat_l%smmm%isvctr + i
                  iline = smat_l%smmm%line_and_column(1,i)
                  icolumn = smat_l%smmm%line_and_column(2,i)
                  if (trace_with_overlap) then
                      ! Take the trace of the product matrix times overlap
                      iismall = matrixindex_in_compressed(smat_s, icolumn, iline)
                      if (iismall>0) then
                          tt=mat%matrix_compr(isshift+iismall-smat_s%isvctrp_tg)
                      else
                          tt=0.d0
                      end if
                  else
                      ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
                      if (iline==icolumn) then
                          tt=1.d0
                      else
                          tt=0.d0
                      end if
                  end if
                  ! This should be improved...
                  if (ibound==1) bound_low = bound_low + penalty_ev(i,1)*tt
                  if (ibound==2) bound_up = bound_up + penalty_ev(i,2)*tt
              end do
              !$omp end do
              !$omp end parallel
          !else
          !    ! This means that the Chebyshev expansion exploded, so take a very large
          !    ! value for the error function such that eigenvalue bounds will be enlarged
          !    if (ibound==1) bound_low = 1.d10
          !    if (ibound==2) bound_up = 1.d10
          !    !!if (emergency_stop(1)) then
          !    !!    ! Need to enlarge boundary
          !    !!    bound_low = 1.d10
          !    !!else
          !    !!    ! No need to enlarge boundary
          !    !!    bound_low = 0.d0
          !    !!end if
          !    !!if (emergency_stop(2)) then
          !    !!    ! Need to enlarge boundary
          !    !!    bound_up = 1.d10
          !    !!else
          !    !!    ! No need to enlarge boundary
          !    !!    bound_up = 0.d0
          !    !!end if
          !end if
      end do
    
      allredarr(1)=bound_low
      allredarr(2)=bound_up
    
      if (nproc > 1) then
          call mpiallred(allredarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
    
      !allredarr=abs(allredarr) !for some crazy situations this may be negative
      !noise=1000.d0*anoise
      noise=10.d0*anoise
      noise = 1.d-1
    
      if (bigdft_mpi%iproc==0) then
          !call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
          !call yaml_map('pnlty',(/allredarr(1),allredarr(2)/),fmt='(es8.1)')
          call yaml_map('penalty',allredarr(1),fmt='(es8.1)')
      end if

      eval_bounds_ok(1) = .true.
      eval_bounds_ok(2) = .true.
      if (any((/abs(allredarr(1))>noise,abs(allredarr(2))>noise/))) then
          if (all((/abs(allredarr(1))>noise,abs(allredarr(2))>noise/))) then
              if (allredarr(1)>0.d0 .and. allredarr(2)<0.d0) then
                  ! lower bound too large
                  eval_bounds_ok(1)=.false.
                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
                  restart=.true.
                  if (present(eval_multiplicator)) then
                      !eval_multiplicator = eval_multiplicator*2.0d0
                      eval_multiplicator = 2.0d0
                  end if
              else if (allredarr(1)<0.d0 .and. allredarr(2)>0.d0) then
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

      !!write(*,*) 'allredarr, anoise', allredarr, anoise
      !if (allredarr(1)>noise) then
      !    eval_bounds_ok(1)=.false.
      !    call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*factor_low,ispin)
      !    restart=.true.
      !    !!if (bigdft_mpi%iproc==0) then
      !    !!    call yaml_map('adjust lower bound',.true.)
      !    !!end if
      !else
      !    eval_bounds_ok(1)=.true.
      !end if
      !if (allredarr(2)>noise) then
      !    eval_bounds_ok(2)=.false.
      !    call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*factor_high,ispin)
      !    restart=.true.
      !    !!if (bigdft_mpi%iproc==0) then
      !    !!    call yaml_map('adjust upper bound',.true.)
      !    !!end if
      !else
      !    eval_bounds_ok(2)=.true.
      !end if
    
      call f_release_routine()
    
    end subroutine check_eigenvalue_spectrum_new


    subroutine scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
               smat1, mat1, i1shift, smat2, mat2, i2shift, &
               matscal_compr, scale_factor, shift_value)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
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
          !write(*,*) 'smatl%smmm%istartendseg_mm',smatl%smmm%istartendseg_mm
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
                  !write(*,*) 'i, ii, ii1, tt1, tt2', i, ii, ii1, tt1, tt2, i1shift, smat1%isvctrp_tg, i1shift+ii1-smat1%isvctrp_tg
                  matscal_compr(ii-smatl%isvctrp_tg)=scale_factor*(tt1-shift_value*tt2)
              end do
          end do
          !$omp end do
          !$omp end parallel
          call timing(iproc,'foe_aux_mcpy  ','OF')
      else
          stop 'scale_and_shift_matrix: wrong data strategy'
      end if
    
      call f_release_routine()
    
    end subroutine scale_and_shift_matrix


    subroutine retransform_ext(iproc, nproc, smat, inv_ovrlp, kernel)
        use module_base
        use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), &
                                     SPARSEMM_SEQ, SPARSE_MATMUL_LARGE
        use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
                                compress_matrix_distributed_wrapper, &
                                sparsemm_new, transform_sparsity_pattern
        implicit none
        ! Calling arguments
        integer,intent(in) :: iproc, nproc
        type(sparse_matrix),intent(in) :: smat
        real(kind=8),dimension(smat%nvctrp_tg),intent(inout) :: inv_ovrlp
        real(kind=8),dimension(smat%nvctrp_tg),intent(inout) :: kernel
        
    
        ! Local variables
        real(kind=8),dimension(:),pointer :: inv_ovrlpp_new, tempp_new
        real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
    
        call f_routine(id='retransform_ext')
    
        inv_ovrlpp_new = f_malloc_ptr(smat%smmm%nvctrp, id='inv_ovrlpp_new')
        tempp_new = f_malloc_ptr(smat%smmm%nvctrp, id='tempp_new')
        inv_ovrlp_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
        kernel_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
        call sequential_acces_matrix_fast2(smat, kernel, kernel_compr_seq)
        call sequential_acces_matrix_fast2(smat, &
             inv_ovrlp, inv_ovrlp_compr_seq)
        if (smat%smmm%nvctrp_mm>0) then !to avoid an out of bounds error
            call transform_sparsity_pattern(smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
                 smat%nseg, smat%keyv, smat%keyg, smat%smmm%line_and_column_mm, &
                 smat%smmm%nvctrp, smat%smmm%isvctr, &
                 smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
                 'small_to_large', inv_ovrlp(smat%smmm%isvctr_mm-smat%isvctrp_tg+1), inv_ovrlpp_new)
        end if
        call sparsemm_new(smat, kernel_compr_seq, inv_ovrlpp_new, tempp_new)
        call sparsemm_new(smat, inv_ovrlp_compr_seq, tempp_new, inv_ovrlpp_new)
        call f_zero(kernel)
        call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_LARGE, &
             inv_ovrlpp_new, kernel)
        call f_free_ptr(inv_ovrlpp_new)
        call f_free_ptr(tempp_new)
        call f_free(inv_ovrlp_compr_seq)
        call f_free(kernel_compr_seq)
    
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
    !!##  real(kind=8),intent(in) :: A, B
    !!##  integer,intent(in) :: n
    !!##  real(kind=8),intent(in) :: ex
    !!##  real(8),dimension(n),intent(out) :: cc
    !!##  real(kind=8),intent(out) :: x_max_error, max_error,mean_error
    !!##
    !!##  ! Local variables
    !!##  integer :: k, j
    !!##  real(kind=8) :: bma, bpa, y, arg, fac, tt
    !!##  real(kind=8),dimension(50000) :: cf
    !!##  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
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


    subroutine init_foe(iproc, nproc, nspin, charge, tmprtr, evbounds_nsatur, evboundsshrink_nsatur, &
               evlow, evhigh, fscale, ef_interpol_det, ef_interpol_chargediff, &
               fscale_lowerbound, fscale_upperbound, foe_obj)
      use module_base
      use foe_base, only: foe_data, foe_data_set_int, foe_data_set_real, foe_data_set_logical, foe_data_get_real, foe_data_null
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, nspin, evbounds_nsatur, evboundsshrink_nsatur
      real(kind=8),dimension(nspin),intent(in) :: charge
      real(kind=8),intent(in) :: evlow, evhigh, fscale, ef_interpol_det
      real(kind=8),intent(in) :: ef_interpol_chargediff, fscale_lowerbound, fscale_upperbound
      real(kind=8),intent(in) :: tmprtr
      type(foe_data), intent(out) :: foe_obj
      
      ! Local variables
      character(len=*), parameter :: subname='init_foe'
      integer :: iorb, ispin
      real(kind=8) :: incr
    
      call timing(iproc,'init_matrCompr','ON')
    
      foe_obj = foe_data_null()
    
      foe_obj%ef = f_malloc0_ptr(nspin,id='(foe_obj%ef)')
      call foe_data_set_real(foe_obj,"ef",0.d0,1)
      if (nspin==2) then
          call foe_data_set_real(foe_obj,"ef",0.d0,2)
      end if
      foe_obj%evlow = f_malloc0_ptr(nspin,id='foe_obj%evlow')
      call foe_data_set_real(foe_obj,"evlow",evlow,1)
      if (nspin==2) then
          call foe_data_set_real(foe_obj,"evlow",evlow,2)
      end if
      foe_obj%evhigh = f_malloc0_ptr(nspin,id='foe_obj%evhigh')
      call foe_data_set_real(foe_obj,"evhigh",evhigh,1)
      if (nspin==2) then
          call foe_data_set_real(foe_obj,"evhigh",evhigh,2)
      end if
      foe_obj%bisection_shift = f_malloc0_ptr(nspin,id='foe_obj%bisection_shift')
      call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,1)
      if (nspin==2) then
          call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,2)
      end if
      call foe_data_set_real(foe_obj,"fscale",fscale)
      call foe_data_set_real(foe_obj,"ef_interpol_det",ef_interpol_det)
      call foe_data_set_real(foe_obj,"ef_interpol_chargediff",ef_interpol_chargediff)
      foe_obj%charge = f_malloc0_ptr(nspin,id='foe_obj%charge')
      call foe_data_set_real(foe_obj,"charge",0.d0,1)
      !!do iorb=1,orbs_KS%norbu
      !!    call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",1)+orbs_KS%occup(iorb),1)
      !!end do
      !!if (nspin==2) then
      !!    call foe_data_set_real(foe_obj,"charge",0.d0,2)
      !!    do iorb=orbs_KS%norbu+1,orbs_KS%norb
      !!         call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",2)+orbs_KS%occup(iorb),2)
      !!    end do
      !!end if
      do ispin=1,nspin
          call foe_data_set_real(foe_obj,"charge",charge(ispin),ispin)
      end do
      call foe_data_set_int(foe_obj,"evbounds_isatur",0)
      call foe_data_set_int(foe_obj,"evboundsshrink_isatur",0)
      call foe_data_set_int(foe_obj,"evbounds_nsatur",evbounds_nsatur)
      call foe_data_set_int(foe_obj,"evboundsshrink_nsatur",evboundsshrink_nsatur)
      call foe_data_set_real(foe_obj,"fscale_lowerbound",fscale_lowerbound)
      call foe_data_set_real(foe_obj,"fscale_upperbound",fscale_upperbound)
      call foe_data_set_real(foe_obj,"tmprtr",tmprtr)
    
      call timing(iproc,'init_matrCompr','OF')
    
    
    end subroutine init_foe


    subroutine accuracy_of_chebyshev_expansion(iproc, nproc, npl, coeff, bounds, h, func, x_max_error, max_error, mean_error)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npl
      real(kind=8),dimension(npl),intent(in) :: coeff
      real(kind=8),dimension(2),intent(in) :: bounds
      real(kind=8),intent(in) :: h
      real(kind=8),external :: func
      real(kind=8),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: isx, iex, i, ipl, n, iimin, iimax, ii, is, np, jproc
      real(kind=8) :: x, xx, val_chebyshev, val_function, xxm1, xxm2, xxx, sigma, tau, error
      real(kind=8),dimension(:,:),allocatable :: max_errors

      call f_routine(id='accuracy_of_chebyshev_expansion')

      sigma = 2.d0/(bounds(2)-bounds(1))
      tau = (bounds(1)+bounds(2))/2.d0

      isx = nint(bounds(1)/h)
      iex = nint(bounds(2)/h)
      n = iex - isx + 1

      ! MPI parallelization... maybe only worth for large n?
      ii = n/nproc
      np = ii
      is = iproc*ii
      ii = n - nproc*ii
      if (iproc<ii) then
          np = np + 1
      end if
      is = is + min(iproc,ii)
      is = is + isx - 1 !shift of the starting index
      !check
      ii = np
      call mpiallred(ii, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      if (ii/=n) then
          call f_err_throw('wrong partition of n')
      end if
      iimin = 1 + is
      call mpiallred(iimin, 1, mpi_min, comm=bigdft_mpi%mpi_comm)
      if (iimin/=isx) then
          call f_err_throw('wrong starting index')
      end if
      iimax = np + is
      call mpiallred(iimax, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      if (iimax/=iex) then
          call f_err_throw('wrong ending index')
      end if


      max_error = 0.d0
      mean_error = 0.d0
      !do i=isx,iex
      do i=1,np
          ii = i + is
          x = real(ii,kind=8)*h
          val_chebyshev = 0.5d0*coeff(1)*1.d0
          xx = sigma*(x-tau)
          val_chebyshev = val_chebyshev + coeff(2)*xx
          xxm2 = 1.d0
          xxm1 = xx
          do ipl=3,npl
              xx = sigma*(x-tau)
              xxx = 2.d0*xx*xxm1 - xxm2
              val_chebyshev = val_chebyshev + coeff(ipl)*xxx
              xxm2 = xxm1
              xxm1 = xxx
          end do
          val_function = func(x)
          error = abs(val_chebyshev-val_function)
          if (error>max_error) then
              max_error = error
              x_max_error = x
          end if
          mean_error = mean_error + error
          !if (abs(bounds(1)-0.15d0)<1.d-1 .and. abs(bounds(2)-30.0d0)<1.d-1 .and. npl==100) then
          !    write(*,*) 'x, val_chebyshev, val_function, max_error', x, val_chebyshev, val_function, max_error
          !end if
          !write(*,*) 'x, val_chebyshev, exp(x-bounds(2))', x, val_chebyshev, exp(x-bounds(2))
      end do
      !write(*,*) 'max_error',max_error
      mean_error = mean_error/real(iex-isx+1,kind=8)

      ! Communicate the results... for the maximum in an array since also the position is required
      call mpiallred(mean_error, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      max_errors = f_malloc0((/1.to.2,0.to.nproc-1/),id='max_errors')
      max_errors(1,iproc) = max_error
      max_errors(2,iproc) = x_max_error
      call mpiallred(max_errors, mpi_sum, comm=bigdft_mpi%mpi_comm)
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
    !!  real(kind=8),intent(in) :: x
    !!  real(kind=8),intent(in) :: power
    !!  real(kind=8) :: x_power
    !!  x_power = x**power
    !!end function x_power


    subroutine get_chebyshev_polynomials(iproc, nproc, itype, foe_verbosity, npl, smatm, smatl, &
               ham_, foe_obj, chebyshev_polynomials, ispin, eval_bounds_ok, &
               hamscal_compr, scale_factor, shift_value, smats, ovrlp_, ovrlp_minus_one_half)
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
      use module_func
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, itype, ispin
      integer,intent(in) :: foe_verbosity
      integer,intent(in) :: npl
      type(sparse_matrix),intent(in) :: smatm, smatl
      type(matrices),intent(in) :: ham_
      type(foe_data),intent(inout) :: foe_obj
      real(kind=8),dimension(:,:),pointer,intent(inout) :: chebyshev_polynomials
      logical,dimension(2),intent(out) :: eval_bounds_ok
      real(kind=8),dimension(smatl%nvctrp_tg),intent(out) :: hamscal_compr
      real(kind=8),intent(out) :: scale_factor, shift_value
      type(sparse_matrix),intent(in),optional :: smats
      type(matrices),intent(in),optional :: ovrlp_
      real(kind=8),dimension(smatl%nvctrp_tg),intent(in),optional :: ovrlp_minus_one_half
    
      ! Local variables
      integer :: jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries
      integer,parameter :: nplx=50000
      real(kind=8),dimension(:,:,:),allocatable :: cc, cc_check
      real(kind=8),dimension(:,:),allocatable ::  fermip_check
      real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=8) :: anoise, sumn, sumn_check, charge_diff, ef_interpol, ddot
      real(kind=8) :: evlow_old, evhigh_old, det, determinant, sumn_old, ef_old, tt
      real(kind=8) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=8) :: fscale, tt_ovrlp, tt_ham, diff, fscale_check, fscale_new
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible, with_overlap
      logical,dimension(2) :: emergency_stop
      real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
      real(kind=8),dimension(:),allocatable :: fermi_check_compr
      real(kind=8),dimension(4,4) :: interpol_matrix
      real(kind=8),dimension(4) :: interpol_vector
      real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: bisection_bounds_ok
      real(kind=8) :: temp_multiplicator, ebs_check, ef, ebsp
      integer :: irow, icol, itemp, iflag,info, i, j, itg, ncount, istl, ists, isshift, imshift
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
      
    
    
      call f_routine(id='get_chebyshev_polynomials')
    
      !if (iproc==0) call yaml_comment('get Chebyshev polynomials',hfill='~')
    
    
      call timing(iproc, 'FOE_auxiliary ', 'ON')

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
    
    
    
      penalty_ev_new = f_malloc((/smatl%smmm%nvctrp,2/),id='penalty_ev_new')
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
    
      fscale_new=1.d100


    
    
    !ispin = 1
    

    
    
          fscale_new = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")
    
    
              fscale = fscale_new
              fscale = max(fscale,FSCALE_LOWER_LIMIT)
              fscale = min(fscale,FSCALE_UPPER_LIMIT)
        
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
            
            
                      cc = f_malloc((/npl,3,1/),id='cc')
            
                      !!if (foe_data_get_real(foe_obj,"evlow",ispin)>=0.d0) then
                      !!    call f_err_throw('Lowest eigenvalue must be negative')
                      !!end if
                      !!if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
                      !!    call f_err_throw('Highest eigenvalue must be positive')
                      !!end if
            
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call timing(iproc, 'chebyshev_coef', 'ON')
            
                      !!if (foe_data_get_real(foe_obj,"tmprtr")/=0.d0) call f_err_throw('tmprtr must be zero')
                      !!call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale)
                      !!call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                      !!     foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,1), &
                      !!     x_max_error, max_error, mean_error)
                      cc = 0.d0
                      call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                           muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,2,1), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      do ipl=1,npl
                         cc(ipl,3,1) = -cc(ipl,2,1)
                      end do
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
            
                      call timing(iproc, 'chebyshev_coef', 'OF')
                      call timing(iproc, 'FOE_auxiliary ', 'ON')
                    
    
                      if (smatl%nspin==1) then
                          do ipl=1,npl
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,2,1)=2.d0*cc(ipl,2,1)
                              cc(ipl,3,1)=2.d0*cc(ipl,3,1)
                          end do
                      end if
                    
                    
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
            
                      emergency_stop=.false.
                          ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','recalculated')
                          if (with_overlap) then
                              call chebyshev_clean(iproc, nproc, npl, cc, &
                                   smatl, hamscal_compr, &
                                   .true., &
                                   nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                                   emergency_stop, invovrlp_compr=ovrlp_minus_one_half)
                          else
                              call chebyshev_clean(iproc, nproc, npl, cc, &
                                   smatl, hamscal_compr, &
                                   .false., &
                                   nsize_polynomial, 1, fermi_new, penalty_ev_new, chebyshev_polynomials, &
                                   emergency_stop)
                          end if
                           
                          !!call transform_sparsity_pattern(smatl%nfvctr, &
                          !!     smatl%smmm%nvctrp_mm, smatl%smmm%isvctr_mm, &
                          !!     smatl%nseg, smatl%keyv, smatl%keyg, smatl%smmm%line_and_column_mm, &
                          !!     smatl%smmm%nvctrp, smatl%smmm%isvctr, &
                          !!     smatl%smmm%nseg, smatl%smmm%keyv, smatl%smmm%keyg, &
                          !!     smatl%smmm%istsegline, 'large_to_small', fermi_small_new, fermi_new)
            
                      call timing(iproc, 'FOE_auxiliary ', 'ON')
            
            
                      restart=.false.
            
                      ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
                      ! (otherwise this has already been checked in the previous iteration).
                      call check_eigenvalue_spectrum_new(nproc, smatl, ispin, &
                            0, 1.0d0, 1.0d0, penalty_ev_new, anoise, .false., emergency_stop, &
                            foe_obj, restart, eval_bounds_ok)
            
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
    
      call timing(iproc, 'FOE_auxiliary ', 'OF')
    
      call f_release_routine()
    
    
    end subroutine get_chebyshev_polynomials


    subroutine find_fermi_level(iproc, nproc, npl, chebyshev_polynomials, &
               foe_verbosity, label, smatl, ispin, foe_obj, kernel_)
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
      !!use foe_common, only: scale_and_shift_matrix, evnoise, &
      !!                      check_eigenvalue_spectrum_new, retransform_ext, get_chebyshev_expansion_coefficients
      use module_func
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npl, ispin
      type(sparse_matrix),intent(in) :: smatl
      real(kind=8),dimension(smatl%smmm%nvctrp_mm,npl) :: chebyshev_polynomials
      integer,intent(in) :: foe_verbosity
      character(len=*),intent(in) :: label
      type(foe_data),intent(inout) :: foe_obj
      type(matrices),intent(inout) :: kernel_
    
      ! Local variables
      integer :: jorb, ipl, it, ii, iiorb, jjorb, iseg, iorb
      integer :: isegstart, isegend, iismall, iilarge, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, ntemp, it_shift, npl_check, npl_boundaries, ilshift
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
      integer :: irow, icol, itemp, iflag,info, isshift, imshift, ilshift2, i, j, itg, ncount, istl, ists
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
                  !write(*,*) 'ef, efarr', foe_data_get_real(foe_obj,"ef",ispin), efarr
    
                  sumnarr(1)=0.d0
                  sumnarr(2)=1.d100
                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",ispin), foe_data_get_real(foe_obj,"ef",ispin), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",ispin), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
                  call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)
            
            
                  if (iproc==0) then
                      if (foe_verbosity>=1) then
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
            
            
                      cc = f_malloc((/npl,3,1/),id='cc')
            
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call timing(iproc, 'chebyshev_coef', 'ON')
            
                      call func_set(FUNCTION_ERRORFUNCTION, efx=foe_data_get_real(foe_obj,"ef",ispin), fscalex=fscale)
                      !!write(*,*) 'evlow, evhigh, ef, fscale', &
                      !!     foe_data_get_real(foe_obj,"evlow",ispin), foe_data_get_real(foe_obj,"evhigh",ispin), &
                      !!     foe_data_get_real(foe_obj,"ef",ispin), fscale
                      call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,1,1), &
                           x_max_error, max_error, mean_error)
                      call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
                           muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1,2,1), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      do ipl=1,npl
                         cc(ipl,3,1) = -cc(ipl,2,1)
                      end do
                      call evnoise(npl, cc(1,2,1), foe_data_get_real(foe_obj,"evlow",ispin), &
                           foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
        

                      !if (iproc==0 .and. foe_verbosity>=1) then
                      !    call yaml_newline()
                      !    call yaml_mapping_open('accuracy (x, max err, mean err)')
                      !    call yaml_map('main',(/x_max_error,max_error,mean_error/),fmt='(es9.2)')
                      !    !call yaml_map('check',(/x_max_error_check,max_error_check,max_error/),fmt='(es9.2)')
                      !    call yaml_mapping_close()
                      !    call yaml_newline()
                      !end if
            
                      call timing(iproc, 'chebyshev_coef', 'OF')
                      call timing(iproc, 'FOE_auxiliary ', 'ON')
                    
    
                      if (smatl%nspin==1) then
                          do ipl=1,npl
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,2,1)=2.d0*cc(ipl,2,1)
                              cc(ipl,3,1)=2.d0*cc(ipl,3,1)
                          end do
                      end if
                    
                    
                      call timing(iproc, 'FOE_auxiliary ', 'OF')
            
                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
                          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                               smatl%nfvctr, smatl%smmm%nfvctrp, &
                              smatl, chebyshev_polynomials, 1, cc, fermi_small_new)
            
            
                      call timing(iproc, 'FOE_auxiliary ', 'ON')
            
            
                      call f_free(cc)
            
                    
                      call calculate_trace_distributed_new(iproc, nproc, smatl, fermi_small_new, sumn)
                      !write(*,*) 'sumn',sumn
            
                      call yaml_map('bisec bounds ok',&
                           (/bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                      call yaml_map('ef',ef)
        
                      if (all(eval_bounds_ok) .and. all(bisection_bounds_ok)) then
                          ! Print these informations already now if all entries are true.
                          if (iproc==0) then
                              !!if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                              !!     (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          end if
                      end if
                      call determine_fermi_level(f, sumn, ef, info)
                      bisection_bounds_ok(1) = fermilevel_get_logical(f,"bisection_bounds_ok(1)")
                      bisection_bounds_ok(2) = fermilevel_get_logical(f,"bisection_bounds_ok(2)")
                      if (info<0) then
                          if (iproc==0) then
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
                      end if
            
                      if (abs(charge_diff)<charge_tolerance) then
                          if (iproc==0) call yaml_sequence_close()
                          diff=0.d0
        
                          if (nproc > 1) then
                              call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
                          end if
        
                          diff=sqrt(diff)
                          if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')
                          exit
                      end if
        
                  end do main_loop
            
             call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                  fermi_small_new, &
                  kernel_%matrix_compr(ilshift+1:))

      !end do spin_loop
    
    
      !call foe_data_set_real(foe_obj,"fscale",fscale_new)
    
      degree_sufficient=.true.
    
      !if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')
    
    
      call f_free(fermi_small_new)
    
      call timing(iproc, 'FOE_auxiliary ', 'OF')
    
      call f_release_routine()
    

    
    
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


end module foe_common
