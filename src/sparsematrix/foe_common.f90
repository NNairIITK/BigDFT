module module_func
  use module_base
  implicit none

  private

  ! Shared variables within the modules
  integer :: ifunc
  real(kind=8) :: power, ef, fscale

  ! Public routines
  public :: func_set
  public :: func

  ! Public parameters
  integer,parameter,public :: FUNCTION_POLYNOMIAL = 101
  integer,parameter,public :: FUNCTION_ERRORFUNCTION = 102

  contains

    subroutine func_set(ifuncx, powerx, efx, fscalex)
      implicit none
      integer,intent(in) :: ifuncx
      real(kind=8),intent(in),optional :: powerx, efx, fscalex
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
      case default
          call f_err_throw("wrong value of 'ifuncx'")
      end select
    end subroutine func_set

    function func(x)
      implicit none
      real(kind=8),intent(in) :: x
      real(kind=8) :: func
      select case (ifunc)
      case(FUNCTION_POLYNOMIAL)
          func = x**power
      case(FUNCTION_ERRORFUNCTION)
          func = 0.5d0*erfc((x-ef)*(1.d0/fscale))
      case default
          call f_err_throw("wrong value of 'ifunc'")
      end select
    end function func

end module module_func


module foe_common
  use module_defs, only: uninitialized
  use module_base
  use foe_base
  implicit none

  private

  !> Public routines
  public :: chebft
  public :: chebft2
  public :: chder
  public :: evnoise
  !public :: erfcc
  !public :: chebev
  public :: pltwght
  public :: pltexp
  public :: get_roots_of_cubic_polynomial
  public :: determinant
  public :: check_eigenvalue_spectrum_new
  public :: scale_and_shift_matrix
  public :: retransform_ext
  public :: cheb_exp
  public :: init_foe


  contains


    ! Calculates chebychev expansion of fermi distribution.
    ! Taken from numerical receipes: press et al
    subroutine chebft(A,B,N,cc,ef,fscale,tmprtr,x_max_error,max_error,mean_error)
      use module_base
      use module_func
      use yaml_output
      implicit none
      
      ! Calling arguments
      real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
      integer,intent(in) :: n
      real(8),dimension(n),intent(out) :: cc
      real(kind=8),intent(out) :: x_max_error, max_error, mean_error
    
      ! Local variables
      integer :: k, j
      real(kind=8) :: bma, bpa, y, arg, fac, tt
      real(kind=8),dimension(50000) :: cf
      !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
    
      call f_routine(id='chebft')

      if (tmprtr/=0.d0) call f_err_throw('tmprtr should be zero for the moment')
    
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
              cf(k)=1.d0/(1.d0+safe_exp( (arg-ef)*(1.d0/tmprtr) ) )
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

      call func_set(FUNCTION_ERRORFUNCTION, efx=ef, fscalex=fscale)
      call accuracy_of_chebyshev_expansion(n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
      !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    
      call f_release_routine()
    
    end subroutine chebft
    
    
    
    ! Calculates chebychev expansion of fermi distribution.
    ! Taken from numerical receipes: press et al
    subroutine chebft2(a,b,n,cc)
      use module_base
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
          cf(k)=safe_exp((arg-b)*ttt)
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


    subroutine check_eigenvalue_spectrum_new(nproc, smat_l, smat_s, mat, ispin, isshift, &
               factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
               emergency_stop, foe_obj, restart, eval_bounds_ok)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
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
      logical,intent(inout) :: restart
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
    
          !$omp parallel default(none) &
          !$omp shared(bound_low, bound_up, smat_l, smat_s, trace_with_overlap, mat, isshift, penalty_ev) &
          !$omp private(i, ii, iline, icolumn, iismall, tt)
          !$omp do reduction(+:bound_low, bound_up)
          do i=1,smat_l%smmm%nvctrp
              ii = smat_l%smmm%isvctr + i
              iline = smat_l%smmm%line_and_column(1,i)
              icolumn = smat_l%smmm%line_and_column(2,i)
              iismall = matrixindex_in_compressed(smat_s, icolumn, iline)
              if (iismall>0) then
                  if (trace_with_overlap) then
                      ! Take the trace of the product matrix times overlap
                      tt=mat%matrix_compr(isshift+iismall-smat_s%isvctrp_tg)
                  else
                      ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
                      if (iline==icolumn) then
                          tt=1.d0
                      else
                          tt=0.d0
                      end if
                  end if
              else
                  tt=0.d0
              end if
              bound_low = bound_low + penalty_ev(i,2)*tt
              bound_up = bound_up +penalty_ev(i,1)*tt
          end do
          !$omp end do
          !$omp end parallel
      else
          ! This means that the Chebyshev expansion exploded, so take a very large
          ! value for the error function such that eigenvalue bounds will be enlarged
          bound_low = 1.d10
          bound_up = 1.d10
      end if
    
      allredarr(1)=bound_low
      allredarr(2)=bound_up
    
      if (nproc > 1) then
          call mpiallred(allredarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
    
      allredarr=abs(allredarr) !for some crazy situations this may be negative
      noise=1000.d0*anoise
    
      !if (bigdft_mpi%iproc==0) then
      !    call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
      !end if
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



    ! Calculates chebychev expansion of x**ex, where ex is any value (typically -1, -1/2, 1/2)
    ! Taken from numerical receipes: press et al
    subroutine cheb_exp(A,B,N,cc,ex,x_max_error,max_error,mean_error)
      use module_base
      use module_func
      use yaml_output
      implicit none
      
      ! Calling arguments
      real(kind=8),intent(in) :: A, B
      integer,intent(in) :: n
      real(kind=8),intent(in) :: ex
      real(8),dimension(n),intent(out) :: cc
      real(kind=8),intent(out) :: x_max_error, max_error,mean_error
    
      ! Local variables
      integer :: k, j
      real(kind=8) :: bma, bpa, y, arg, fac, tt
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
          cf(k)=arg**ex
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

      call func_set(FUNCTION_POLYNOMIAL, powerx=ex)
      call accuracy_of_chebyshev_expansion(n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
      !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    
      call f_release_routine()
    
    end subroutine cheb_exp


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


    subroutine accuracy_of_chebyshev_expansion(npl, coeff, bounds, h, func, x_max_error, max_error, mean_error)
      implicit none

      ! Calling arguments
      integer,intent(in) :: npl
      real(kind=8),dimension(npl),intent(in) :: coeff
      real(kind=8),dimension(2),intent(in) :: bounds
      real(kind=8),intent(in) :: h
      real(kind=8),external :: func
      real(kind=8),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: is, ie, i, ipl
      real(kind=8) :: x, xx, val_chebyshev, val_function, xxm1, xxm2, xxx, sigma, tau, error

      call f_routine(id='accuracy_of_chebyshev_expansion')

      sigma = 2.d0/(bounds(2)-bounds(1))
      tau = (bounds(1)+bounds(2))/2.d0

      is = nint(bounds(1)/h)
      ie = nint(bounds(2)/h)
      max_error = 0.d0
      mean_error = 0.d0
      do i=is,ie
          x = real(i,kind=8)*h
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
          !write(*,*) 'x, val_chebyshev, val_function', x, val_chebyshev, val_function
      end do
      mean_error = mean_error/real(ie-is+1,kind=8)

      call f_release_routine()

    end subroutine accuracy_of_chebyshev_expansion


    !!pure function x_power(x, power)
    !!  implicit none
    !!  real(kind=8),intent(in) :: x
    !!  real(kind=8),intent(in) :: power
    !!  real(kind=8) :: x_power
    !!  x_power = x**power
    !!end function x_power


end module foe_common
