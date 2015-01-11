!!****f* m_splines/splfit
!! NAME
!!  splfit
!!
!! FUNCTION
!!  Evaluate cubic spline fit to get function values on input set
!!  of ORDERED, UNFORMLY SPACED points.
!!  Optionally gives derivatives (first and second) at those points too.
!!  If point lies outside the range of arg, assign the extremal
!!  point values to these points, and zero derivative.
!! 
!! INPUTS
!!  arg(numarg)=equally spaced arguments (spacing delarg) for data 
!!   to which spline was fit.
!!  fun(numarg,2)=function values to which spline was fit and spline
!!   fit to second derivatives (from Numerical Recipes spline).
!!  ider=  see above
!!  newarg(numnew)=new values of arguments at which function is desired.
!!  numarg=number of arguments at which spline was fit.
!!  numnew=number of arguments at which function values are desired.
!!
!! OUTPUT
!!  derfun(numnew)=(optional) values of first or second derivative of function.
!!   This is only computed for ider=1 or 2; otherwise derfun not used.
!!  newfun(numnew)=values of function at newarg(numnew).
!!   This is only computed for ider=0 or 1.
!!
!! NOTES
!!        if ider=0, compute only the function (contained in fun)
!!        if ider=1, compute the function (contained in fun) and its first derivative (in derfun)
!!        if ider=2, compute only the second derivative of the function (in derfun)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine splfit(arg,derfun,fun,ider,newarg,newfun,numarg,numnew)

 use abi_defs_basis
 implicit none

 integer, intent(in) :: ider,numarg,numnew
 double precision, intent(in) :: arg(numarg),fun(numarg,2),newarg(numnew)
 double precision, intent(out) :: derfun(numnew),newfun(numnew)
 
 integer :: i,jspl
 double precision :: argmin,delarg,d,aa,bb,cc,dd
 character(len=500) :: msg

!argmin is smallest x value in spline fit; delarg is uniform spacing of spline argument
 argmin=arg(1)
 delarg=(arg(numarg)-argmin)/dble(numarg-1)

 jspl=-1

!Do one loop for no grads, other for grads:
 if (ider==0) then
! Spline index loop for no grads:
  do i=1,numnew
   if (newarg(i).ge.arg(numarg)) then
! MJV 6/4/2009 FIXME this message is never used
    write(msg,1000)char(10),i,newarg(i), &
&     jspl,char(10),numarg,arg(numarg),char(10),char(10),char(10)
1000 format(a1,' splfit: for arg number',i8,2x,'of value', &
&    1p,e12.4,1x,'jspl=',i8,a1,' is >= numarg=',i8,  &
&    '  (max arg(numarg)=',e12.4,')',a1,             &
&    ' Means function values are being requested outside',       &
&    ' range of data.',a1,' Function and slope will be set to',  &
&    ' values at upper end of data.',a1)

    newfun(i)=fun(numarg,1)

   else if (newarg(i).le.arg(1)) then
    newfun(i)=fun(1,1)

   else
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
    cc = aa*(aa**2-1.0d0)*(delarg**2/6.0d0)
    dd = bb*(bb**2-1.0d0)*(delarg**2/6.0d0)
    newfun(i)=aa*fun(jspl,1)+bb*fun(jspl+1,1)+cc*fun(jspl,2)+dd*fun(jspl+1,2)
   end if
  enddo

 else if(ider==1)then

! Spline index loop includes grads:
  do i=1,numnew

   if (newarg(i).ge.arg(numarg)) then
    newfun(i)=fun(numarg,1)
    derfun(i)=0.0d0

   else if (newarg(i).le.arg(1)) then
    newfun(i)=fun(1,1)
    derfun(i)=0.0d0

   else

!   cubic spline interpolation:
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
    cc = aa*(aa**2-1.0d0)*(delarg**2/6.0d0)
    dd = bb*(bb**2-1.0d0)*(delarg**2/6.0d0)
    newfun(i)=aa*fun(jspl,1)+bb*fun(jspl+1,1)+cc*fun(jspl,2)+dd*fun(jspl+1,2)
!   spline fit to first derivative:
!   note correction of Numerical Recipes sign error
    derfun(i) = (fun(jspl+1,1)-fun(jspl,1))/delarg +    &
&      (-(3.d0*aa**2-1.d0)*fun(jspl,2)+                 &
&        (3.d0*bb**2-1.d0)*fun(jspl+1,2)) * delarg/6.0d0

          end if
  enddo

 else if (ider==2) then

  do i=1,numnew

   if (newarg(i).ge.arg(numarg)) then
    derfun(i)=0.0d0

   else if (newarg(i).le.arg(1)) then
    derfun(i)=0.0d0

   else

!   cubic spline interpolation:
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
!   second derivative of spline (piecewise linear function)
    derfun(i) = aa*fun(jspl,2)+bb*fun(jspl+1,2)

   end if
  enddo

 end if

 return

end subroutine splfit
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline
!! NAME
!!  spline
!!
!! FUNCTION
!!  SPLINE (originally SPLINE_CUBIC_SET) computes the second derivatives
!!  of a cubic spline.
!!
!! INPUTS
!!    Input, integer N, the number of data points; N must be at least 2. 
!!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!!    spline will actually be linear. 
!!
!!    Input, double precision T(N), the knot values, that is, the points where data
!!    is specified.  The knot values should be distinct, and increasing.
!!
!!    Input, double precision Y(N), the data values to be interpolated.
!!
!!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!!
!! OUTPUT
!!    Output, double precision YPP(N), the second derivatives of the cubic spline.
!!    Work space, double precision DIAG(N) - should be removed ...
!!
!! SOURCE

subroutine spline( t, y, n, ybcbeg, ybcend, ypp )

!*******************************************************************************
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output, 
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to 
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) ) 
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1)) 
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      = 
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL) 
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)    
!
!    Boundary conditions must be applied at the first and last knots.  
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    07 February 1999
!    28 November 2004 XGonze : double precision
!                              make arguments similar to the Numeric Recipes routine
!                              also use algorithmics similar to the Numeric Recipes routine
!
!  Author:
!
!    John Burkardt
!    (XGonze got it from http://www.psc.edu/~burkardt/src/spline/spline.html)
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2. 
!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!    spline will actually be linear. 
!
!    Input, double precision T(N), the knot values, that is, the points where data
!    is specified.  The knot values should be distinct, and increasing.
!
!    Input, double precision Y(N), the data values to be interpolated.
!
!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!
!    Output, double precision YPP(N), the second derivatives of the cubic spline.
!
!    Work space, double precision DIAG(N) - should be removed ...
!
!
!    XG041127 : In the initial implementation, one had the control on
!     IBCBEG and IBCEND. Now, they are determined by the values
!     of YBCBEG, YBCEND. Option 2 has been disabled.
!
!    Input, integer IBCBEG, left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, integer IBCEND, right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.

  use abi_defs_basis
  use abi_interfaces_lowlevel
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: t(n)
  double precision, intent(in) :: y(n)
  double precision, intent(in) :: ybcbeg
  double precision, intent(in) :: ybcend

  double precision, intent(out) :: ypp(n)

  integer :: ibcbeg
  integer :: ibcend
  integer :: i,k
  double precision :: ratio,pinv
  double precision, allocatable :: tmp(:)
!
!  Check.
!
  if ( n <= 1 ) then
    write(std_out,* ) ' '
    write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
    write(std_out,* ) '  The number of knots must be at least 2.'
    write(std_out,* ) '  The input value of N = ', n
    call abi_wrtout(std_out,"Fatal error",'COLL')
    call abi_leave_new('COLL')
  end if

  allocate(tmp(n))

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write(std_out,* ) ' '
      write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
      write(std_out,* ) '  The knots must be strictly increasing, but'
      write(std_out,* ) '  T(',  i,') = ', t(i)
      write(std_out,* ) '  T(',i+1,') = ', t(i+1)
      call abi_wrtout(std_out,"Fatal error",'COLL')
      call abi_leave_new('COLL')
    end if
  end do
!
!  XG041127
  ibcbeg=1 ; ibcend=1
  if(ybcbeg>1.0d+30)ibcbeg=0
  if(ybcend>1.0d+30)ibcend=0
!
!  Set the first and last equations.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.d0
    tmp(1) = 0.d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = -0.5d0
    tmp(1) = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
  end if
  if ( ibcend == 0 ) then
    ypp(n) = 0.d0
    tmp(n) = 0.d0
  else if ( ibcend == 1 ) then
    ypp(n) = 0.5d0
    tmp(n) = (3.d0/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
  end if

!
!  Set the intermediate equations.
!
  do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0d0/(ratio*ypp(i-1) + 2.0d0)
   ypp(i) = (ratio-1.0d0)*pinv
   tmp(i)=(6.0d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&    /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0.d0)) tmp(i)=0.d0   !MT20050927
  enddo

! Solve the equations
  ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0d0)
  do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
  enddo

  deallocate(tmp)

  return
end subroutine spline
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_bicubic
!! NAME
!!  spline_bicubic
!!
!! FUNCTION
!!  Generates coefficients for bicubic spline interpolation.
!!
!! INPUTS
!!  n1 = length of first dimension
!!  n2 = length of second dimension
!!  x1 = positions on first dimension
!!  x2 = positions on second dimension
!!  y = function values on the (x1,x2) grid
!!  der1_x1 = first derivative of y wrt x1
!!  der1_x2 = first derivative of y wrt x2
!!  der2_x1x2 = second-order cross-derivative of y wrt x1x2
!!
!! OUTPUT
!!  spl_c = spline coefficients
!!
!! NOTES
!!  Adapted from Numerical Recipes and libbci.
!!
!! SOURCE

subroutine spline_bicubic(n1,n2,x1,x2,y,der1_x1,der1_x2,der2_x1x2,spl_c)

  use abi_defs_basis

  implicit none

  integer,intent(in)  :: n1,n2
  real(dp),intent(in) :: x1(n1),x2(n2),y(n1,n2)
  real(dp),intent(in) :: der1_x1(n1,n2),der1_x2(n1,n2),der2_x1x2(n1,n2)
  real(dp),intent(out):: spl_c(4,4,n1,n2)

  integer :: i1,i2
  real(dp) :: dx1,dx2,wt(16,16),z(16)

  data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4, &
&          8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6, &
&          2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3, &
&          2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0, &
&          -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0, &
&          -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2, &
&          -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

  ! Set coefficients for i1<n1 and i2<n2
  do i2 = 1,n2-1
    do i1 = 1,n1-1
      dx1 = x1(i1+1) - x1(i1)
      dx2 = x2(i2+1) - x2(i2)
      z(1)  = y(i1,i2)
      z(2)  = y(i1+1,i2)
      z(3)  = y(i1+1,i2+1)
      z(4)  = y(i1,i2+1)
      z(5)  = der1_x1(i1,i2) * dx1
      z(6)  = der1_x1(i1+1,i2) * dx1
      z(7)  = der1_x1(i1+1,i2+1) * dx1
      z(8)  = der1_x1(i1,i2+1) * dx1
      z(9)  = der1_x2(i1,i2) * dx2
      z(10) = der1_x2(i1+1,i2) * dx2
      z(11) = der1_x2(i1+1,i2+1) * dx2
      z(12) = der1_x2(i1,i2+1) * dx2
      z(13) = der2_x1x2(i1,i2) * dx1 * dx2
      z(14) = der2_x1x2(i1+1,i2) * dx1 * dx2
      z(15) = der2_x1x2(i1+1,i2+1) * dx1 * dx2
      z(16) = der2_x1x2(i1,i2+1) * dx1 * dx2
      z = matmul(wt,z)
      spl_c(:,:,i1,i2) = reshape(z,(/4,4/),order=(/2,1/))
    end do
  end do

! Set coefficients for i1=n1 and i2=n2 (valid only at the border)
  spl_c(:,:,n1,:) = 0
  spl_c(:,:,:,n2) = 0
  spl_c(1,1,n1,:) = y(n1,:)
  spl_c(1,1,:,n2) = y(:,n2)

end subroutine spline_bicubic
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_c
!! NAME
!!  spline_c
!!
!! FUNCTION
!!  Computes the spline of a complex function.
!!
!! INPUTS
!!  nomega_lo   = number of point in the non regular grid (e.g.  !logarithmic)
!!  nomega_li   = number of point in the regular grid on which the  spline is computed
!!  omega_lo    = value of freq on the 1st grid
!!  omega_li    = value of freq on the 2nd grid
!!  tospline_lo = function on the 1st grid
!!
!! OUTPUT
!!  splined_lo  = spline  (on the 2nd grid)
!!
!! SOURCE

subroutine spline_c( nomega_lo, nomega_li, omega_lo, omega_li, splined_li, tospline_lo)

 use abi_defs_basis
 implicit none

!Arguments --------------------------------------------
!scalars
 integer, intent(in) :: nomega_lo, nomega_li
 real(dp), intent(in) :: omega_lo(nomega_lo)
 real(dp), intent(in) :: omega_li(nomega_li)
 complex(dpc), intent(in) :: tospline_lo(nomega_lo)
 complex(dpc), intent(out) :: splined_li(nomega_li)

!Local variables---------------------------------------
!scalars
 complex(dpc) :: ybcbeg, ybcend
 complex(dpc), allocatable :: ysplin2_lo(:)

 ybcbeg=czero
 ybcend=czero

 allocate(ysplin2_lo(nomega_lo))
 call spline_complex(omega_lo, tospline_lo, nomega_lo, ybcbeg, ybcend, ysplin2_lo)
 call splint_complex( nomega_lo, omega_lo, tospline_lo,ysplin2_lo, nomega_li, omega_li, splined_li)
 deallocate(ysplin2_lo)

end subroutine spline_c
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_complex
!! NAME
!!  spline_complex
!!
!! FUNCTION
!!  spline_complex interfaces the usual spline routine in a case of a 
!!  complex function
!!
!! INPUTS
!!    Input, integer N, the number of data points; N must be at least 2. 
!!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!!    spline will actually be linear. 
!!
!!    Input, double precision T(N), the knot values, that is, the points where data
!!    is specified.  The knot values should be distinct, and increasing.
!!
!!    Input, complex Y(N), the data values to be interpolated.
!!
!!    Input, complex YBCBEG, YBCEND, the values to be used in the boundary
!!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!!
!! OUTPUT
!!    Output, complex YPP(N), the second derivatives of the cubic spline.
!!
!! SOURCE

subroutine spline_complex( t, y, n, ybcbeg, ybcend, ypp )

 use abi_defs_basis
 implicit none

 integer, intent(in) :: n
 real(dp), intent(in) :: t(n)
 complex(dpc), intent(in) :: y(n)
 complex(dpc), intent(in) :: ybcbeg
 complex(dpc), intent(in) :: ybcend
 complex(dpc), intent(out) :: ypp(n)

 real(dp), allocatable :: y_r(:)
 real(dp) :: ybcbeg_r
 real(dp) :: ybcend_r
 real(dp), allocatable :: ypp_r(:)
 real(dp), allocatable :: y_i(:)
 real(dp) :: ybcbeg_i
 real(dp) :: ybcend_i
 real(dp), allocatable :: ypp_i(:)

 allocate(y_r(n),ypp_r(n),y_i(n),ypp_i(n))
 y_r=real(y)
 y_i=imag(y)
 ybcbeg_r=real(ybcbeg)
 ybcbeg_i=imag(ybcbeg)
 ybcend_r=real(ybcend)
 ybcend_i=imag(ybcend)
 call spline( t, y_r, n, ybcbeg_r, ybcend_r, ypp_r )
 call spline( t, y_i, n, ybcbeg_i, ybcend_i, ypp_i )
 ypp=cmplx(ypp_r,ypp_i)
 deallocate(y_r,ypp_r,y_i,ypp_i)

 return

end subroutine spline_complex
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splint
!! NAME
!!  splint
!!
!! FUNCTION
!!  Compute spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): function on output mesh
!!  [ierr]=A non-zero value is used to signal that some points in xfit exceed xspline(nspline).
!!    The input value is incremented by the number of such points.
!!
!! SOURCE

subroutine splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit,ierr)


 use abi_defs_basis
 use abi_interfaces_lowlevel
 implicit none

 integer, intent(in) :: nfit, nspline
 integer, intent(out) :: ierr
 double precision, intent(in) :: xspline(nspline)
 double precision, intent(in) :: yspline(nspline)
 double precision, intent(in) :: ysplin2(nspline)
 double precision, intent(in) :: xfit(nfit)

 double precision, intent(out) :: yfit(nfit)


!local
 integer :: left,i,k,right,my_err
 double precision :: delarg,invdelarg,aa,bb

!source

 my_err=0

 left = 1
 do i=1, nfit
   yfit(i)=0.d0  ! Initialize for the unlikely event that rmax exceed r(mesh)
   !
   do k=left+1, nspline
     if(xspline(k) >= xfit(i)) then
       if(xspline(k-1) <= xfit(i)) then
         right = k
         left = k-1
       else
         if (k-1.eq.1 .and. i.eq.1) then
           call abi_wrtout(std_out,'xfit(1) < xspline(1)','COLL')
           call abi_leave_new('COLL')
           !my_err=my_err+1
           !exit
         else
           call abi_wrtout(std_out,'xfit not properly ordered','COLL')
           call abi_leave_new('COLL')
         end if
       end if
       delarg= xspline(right) - xspline(left)
       invdelarg= 1.0d0/delarg
       aa= (xspline(right)-xfit(i))*invdelarg
       bb= (xfit(i)-xspline(left))*invdelarg

       yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&               +( (aa*aa*aa-aa)*ysplin2(left) +         &
&                  (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0d0
       exit
     end if
   end do ! k
   !
   if (k==nspline+1) my_err=my_err+1 ! xfit not found 
 end do ! i

 ierr=my_err

end subroutine splint
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splint_complex
!! NAME
!!  splint_complex
!!
!! FUNCTION
!!  Interface to the usual splint to compute *complex* spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): complex function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): complex function on output mesh
!!
!! TODO
!! change double precision by real(dp) ( the same in splint.F90)
!!
!! SOURCE

subroutine splint_complex (nspline,xspline,yspline,ysplin2,nfit,xfit,yfit)

 use abi_defs_basis
 implicit none

 integer, intent(in) :: nfit, nspline
 real(dp), intent(in) :: xspline(nspline)
 complex(dpc), intent(in) :: yspline(nspline)
 complex(dpc), intent(in) :: ysplin2(nspline)
 real(dp), intent(in) :: xfit(nfit)
 complex(dpc), intent(out) :: yfit(nfit)

 integer :: ierr
 real(dp), allocatable :: ysplin2_r(:)
 real(dp), allocatable :: ysplin2_i(:)
 real(dp), allocatable :: yspline_r(:)
 real(dp), allocatable :: yspline_i(:)
 real(dp), allocatable :: yfit_r(:)
 real(dp), allocatable :: yfit_i(:)

 allocate(yspline_r(nspline))
 allocate(yspline_i(nspline))
 allocate(ysplin2_r(nspline))
 allocate(ysplin2_i(nspline))
 allocate(yfit_r(nfit))
 allocate(yfit_i(nfit))

!local

!source
 yspline_r=real(yspline)
 yspline_i=imag(yspline)
 ysplin2_r=real(ysplin2)
 ysplin2_i=imag(ysplin2)
 call splint (nspline,xspline,yspline_r,ysplin2_r,nfit,xfit,yfit_r,ierr)
 call splint (nspline,xspline,yspline_i,ysplin2_i,nfit,xfit,yfit_i,ierr)
 yfit=cmplx(yfit_r,yfit_i)
 deallocate(yspline_r,yspline_i,ysplin2_r,ysplin2_i,yfit_r,yfit_i)
 return

end subroutine splint_complex
!!***

