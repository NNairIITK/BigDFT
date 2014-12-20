!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_numeric
!! NAME
!!  m_paw_numeric
!!
!! FUNCTION
!!  Wrappers for various numeric operations (spline, sort, ...)
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (MT,TR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#include "libpaw.h"

module m_paw_numeric

 use defs_basis
 use m_errors
 USE_MEMORY_PROFILING

 implicit none

 private

!public procedures
 public:: paw_spline
 public:: paw_splint
 public:: paw_smooth
 public:: paw_sort_dp
 public:: jbessel
 public:: solvbes
 public:: jbessel_4spline
!!***

CONTAINS
!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/paw_spline
!! NAME
!!  paw_spline
!!
!! FUNCTION
!!  Computes the second derivatives of a cubic spline
!!
!! INPUTS
!!  * Input, integer N, the number of data points; N must be at least 2. 
!!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!!    spline will actually be linear. 
!!  * Input, real(dp) T(N), the knot values, that is, the points where data
!!    is specified.  The knot values should be distinct, and increasing.
!!  * Input, real(dp) Y(N), the data values to be interpolated.
!!  * Input, real(dp) YBCBEG, YBCEND, the values to be used in the boundary
!!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!!
!! OUTPUT
!!    Output, real(dp) YPP(N), the second derivatives of the cubic spline.
!!    Work space, real(dp) DIAG(N) - should be removed ...
!!
!! PARENTS
!!      eltfrxc3,m_atompaw,m_gaussfit,m_paw_pwaves_lmn,m_paw_pwij,m_pawpsp
!!      m_pawxmlps
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine paw_spline(t,y,n,ybcbeg,ybcend,ypp)

!*******************************************************************************
!  Discussion:
!    For data interpolation, the user must call SPLINE_CUBIC_SET to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output, 
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to 
!    evaluate the spline at any point.
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) ) 
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!    Since the first derivative of the spline is
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1)) 
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      = 
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!    or
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL) 
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!    - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)    
!    Boundary conditions must be applied at the first and last knots.  
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Author:
!    John Burkardt, modified by Xavier Gonze

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_spline'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: ybcbeg,ybcend
!arrays
 real(dp),intent(in) :: t(n),y(n)
 real(dp),intent(out) :: ypp(n)

!Local variables-------------------------------
!scalars
 integer :: ibcbeg,ibcend,i,k
 real(dp) :: ratio,pinv
!arrays
 real(dp),allocatable :: tmp(:)

! *************************************************************************

! DBG_ENTER("COLL")

!Check
 if (n<=1) then
   write(std_out,*) ' '
   write(std_out,*) 'SPLINE_CUBIC_SET - Fatal error!'
   write(std_out,*) '  The number of knots must be at least 2.'
   write(std_out,*) '  The input value of N = ', n
   MSG_ERROR("Fatal error")
 end if

 LIBPAW_ALLOCATE(tmp,(n))

 do i=1,n-1
   if (t(i)>=t(i+1)) then
     write(std_out,*) ' '
     write(std_out,*) 'SPLINE_CUBIC_SET - Fatal error!'
     write(std_out,*) '  The knots must be strictly increasing, but'
     write(std_out,*) '  T(',  i,') = ', t(i)
     write(std_out,*) '  T(',i+1,') = ', t(i+1)
     MSG_ERROR("Fatal error")
   end if
 end do

 ibcbeg=1;if(ybcbeg>1.0d+30)ibcbeg=0
 ibcend=1;if(ybcend>1.0d+30)ibcend=0

!Set the first and last equations
 if (ibcbeg==0) then
   ypp(1) = 0._dp
   tmp(1) = 0._dp
 else if ( ibcbeg == 1 ) then
   ypp(1) = -0.5_dp
   tmp(1) = (3._dp/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
 end if
 if (ibcend==0) then
   ypp(n) = 0._dp
   tmp(n) = 0._dp
 else if ( ibcend == 1 ) then
   ypp(n) = 0.5_dp
   tmp(n) = (3._dp/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
 end if

!Set the intermediate equations
 do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0_dp/(ratio*ypp(i-1) + 2.0_dp)
   ypp(i) = (ratio-1.0_dp)*pinv
   tmp(i)=(6.0_dp*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&        /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0._dp)) tmp(i)=0._dp
 end do

!Solve the equations
 ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0_dp)
 do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
 end do

 LIBPAW_DEALLOCATE(tmp)

! DBG_EXIT("COLL")

end subroutine paw_spline
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/paw_splint
!! NAME
!!  paw_splint
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
!! PARENTS
!!      m_atompaw,m_gaussfit,m_paw_finegrid,m_paw_pwaves_lmn,m_pawpsp
!!      m_pawxmlps
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine paw_splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_splint'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfit, nspline
 integer,optional,intent(out) :: ierr
!arrays
 real(dp),intent(in) :: xspline(nspline),yspline(nspline)
 real(dp),intent(in) :: ysplin2(nspline),xfit(nfit)
 real(dp),intent(out) :: yfit(nfit)

!Local variables-------------------------------
!scalars
 integer :: left,i,k,right,my_err
 real(dp) :: delarg,invdelarg,aa,bb
!arrays

! *************************************************************************

 DBG_ENTER("COLL")

 my_err=0
 left=1
 do i=1,nfit
   yfit(i)=0._dp  ! Initialize for the unlikely event that rmax exceed r(mesh)
   do k=left+1, nspline
     if(xspline(k) >= xfit(i)) then
       if(xspline(k-1) <= xfit(i)) then
         right = k
         left = k-1
       else
         if (k-1.eq.1 .and. i.eq.1) then
           MSG_ERROR('xfit(1) < xspline(1)')
         else
           MSG_ERROR('xfit not properly ordered')
         end if
       end if
       delarg= xspline(right) - xspline(left)
       invdelarg= 1.0_dp/delarg
       aa= (xspline(right)-xfit(i))*invdelarg
       bb= (xfit(i)-xspline(left))*invdelarg
       yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&               +( (aa*aa*aa-aa)*ysplin2(left) +         &
&                  (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0_dp
       exit
     end if
   end do ! k
   if (k==nspline+1) my_err=my_err+1 ! xfit not found 
 end do ! i
 if (present(ierr)) ierr=my_err

 DBG_EXIT("COLL")

end subroutine paw_splint
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/paw_smooth
!! NAME
!! paw_smooth
!!
!! FUNCTION
!! Smooth an array of given ordinates (y's) that are in order of
!! increasing abscissas (x's), but without using the abscissas themselves
!! supposed to be equally spaced.
!!
!! INPUTS
!!  it=number of abscissas to treat
!!  mesh=size of the array (number of abscissas)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  a(mesh)=array to be smoothed
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine paw_smooth(a,mesh,it)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_smooth'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: it,mesh
!arrays
 real(dp), intent(inout) :: a(mesh)

!Local variables-------------------------------
!scalars
 integer :: i,k
!arrays
 real(dp) :: asm(mesh)

! *************************************************************************

 DBG_ENTER("COLL")

 asm(1:4) = zero ! ?? Correct me ...
 do k=1,it
   asm(5)=0.2_dp*(a(3)+a(4)+a(5)+a(6)+a(7))
   asm(mesh-4)=0.2_dp*(a(mesh-2)+a(mesh-3)+a(mesh-4)+&
&                     a(mesh-5)+a(mesh-6))
   asm(mesh-3)=0.2_dp*(a(mesh-1)+a(mesh-2)+a(mesh-3)+&
&                     a(mesh-4)+a(mesh-5))
   asm(mesh-2)=0.2_dp*(a(mesh)+a(mesh-1)+a(mesh-2)+&
&                     a(mesh-3)+a(mesh-4))
   asm(mesh-1)=0.25_dp*(a(mesh)+a(mesh-1)+a(mesh-2)+a(mesh-3))
   asm(mesh)=1.0_dp/3.0_dp*(a(mesh)+a(mesh-1)+a(mesh-2))
   do i=6,mesh-5
     asm(i)=0.1_dp *a(i)+0.1_dp*(a(i+1)+a(i-1))+&
&           0.1_dp *(a(i+2)+a(i-2))+&
&           0.1_dp *(a(i+3)+a(i-3))+&
&           0.1_dp *(a(i+4)+a(i-4))+&
&           0.05_dp*(a(i+5)+a(i-5))
   end do
   do i=1,mesh
     a(i)=asm(i)
   end do
 end do

 DBG_EXIT("COLL")

end subroutine paw_smooth
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/paw_sort_dp
!! NAME
!!  paw_sort_dp
!!
!! FUNCTION 
!!  Sort real(dp) array list(n) into ascending numerical order using Heapsort
!!  algorithm, while making corresponding rearrangement of the integer
!!  array iperm. Consider that two real(dp) numbers
!!  within tolerance tol are equal.
!!
!! INPUTS
!!  n        intent(in)    dimension of the list
!!  tol      intent(in)    numbers within tolerance are equal
!!  list(n)  intent(inout) list of real(dp) numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation given the right ascending order
!!
!! PARENTS
!!      m_paw_finegrid
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine paw_sort_dp(n,list,iperm,tol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_sort_dp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: n
 real(dp), intent(in) :: tol
!arrays
 integer, intent(inout) :: iperm(n)
 real(dp), intent(inout) :: list(n)

!Local variables-------------------------------
!scalars
 integer :: l,ir,iap,i,j
 real(dp) :: ap
!arrays

! *************************************************************************

 DBG_ENTER("COLL")

!Accomodate case of array of length 1: already sorted!
 if (n==1) return

!Should not call with n<1
 if (n<1) then
   write(std_out,'(a,i12,/,a)') &
&   'paw_sort_dp has been called with array length n=',n,&
&   ' having a value less than 1.  This is not allowed.'
   MSG_ERROR('fatal error')
 end if

!Conduct the usual sort
 l=n/2+1 ; ir=n
 do ! Infinite do-loop
   if (l>1) then
     l=l-1
     ap=list(l)
     iap=iperm(l)
   else ! l<=1
     ap=list(ir)
     iap=iperm(ir)
     list(ir)=list(1)
     iperm(ir)=iperm(1)
     ir=ir-1
     if (ir==1) then
       list(1)=ap
       iperm(1)=iap
       exit   ! This is the end of this algorithm
     end if
   end if ! l>1
   i=l
   j=l+l
   do while (j<=ir) 
     if (j<ir) then
       if ( list(j)<list(j+1)-tol .or.  &
&          (list(j)<list(j+1)+tol.and.iperm(j)<iperm(j+1))) j=j+1
     endif
     if (ap<list(j)-tol.or.(ap<list(j)+tol.and.iap<iperm(j))) then
       list(i)=list(j)
       iperm(i)=iperm(j)
       i=j
       j=j+j
     else
       j=ir+1
     end if
   end do
   list(i)=ap
   iperm(i)=iap
 end do ! End infinite do-loop

 DBG_EXIT("COLL")

end subroutine paw_sort_dp
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/jbessel
!! NAME
!! jbessel
!!
!! FUNCTION
!! Compute spherical Bessel function j_l(x) and derivative(s)
!!
!! INPUTS
!!  ll=l-order of the Bessel function
!!  order=1 if first derivative is requested
!!        2 if first and second derivatives are requested
!!  xx=where to compute j_l
!!
!! OUTPUT
!!  bes= Bessel function j_l at xx
!!  besp= first derivative of j_l at xx (only if order>=1)
!!  bespp= second derivative of j_l at xx (only if order=2)
!!
!! PARENTS
!!      m_atompaw,m_paw_finegrid,m_paw_numeric,m_vcoul,pawtwdij_1,pawtwdij_2b
!!      pawtwdij_2e,psp11nl
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine jbessel(bes,besp,bespp,ll,order,xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jbessel'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,order
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: bes,besp,bespp

!Local variables ---------------------------------------
!scalars
 integer,parameter :: imax=40
 integer :: ii,il
 real(dp),parameter :: prec=1.d-15
 real(dp) :: besp1,fact,factp,factpp,jn,jnp,jnpp,jr,xx2,xxinv

! *********************************************************************

 DBG_ENTER("COLL")

 if (order>2) then 
   MSG_BUG("Wrong order in jbessel!")
 end if

 if (abs(xx)<prec) then
   bes=zero;if (ll==0) bes=one
   if (order>=1) then
     besp=zero;if (ll==1) besp=third
   end if
   if (order==2) then
     bespp=zero
     if (ll==0) bespp=-third
     if (ll==2) bespp=2._dp/15._dp
   end if
   return
 end if

 xxinv=one/xx
 if (order==0) then
   factp=zero
   factpp=zero
   jnp=zero
   jnpp=zero
 end if 

 if (xx<one) then
   xx2=0.5_dp*xx*xx
   fact=one
   do il=1,ll
     fact=fact*xx/dble(2*il+1)
   end do
   jn=one;jr=one;ii=0
   do while(abs(jr)>=prec.and.ii<imax)
     ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+1))
     jn=jn+jr
   end do
   bes=jn*fact
   if (abs(jr)>prec) then 
     MSG_ERROR('Bessel function did not converge!')
   end if
   if (order>=1) then
     factp=fact*xx/dble(2*ll+3)
     jnp=one;jr=one;ii=0
     do while(abs(jr)>=prec.AND.ii<imax)
       ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+3))
       jnp=jnp+jr
     end do
     besp=-jnp*factp+jn*fact*xxinv*dble(ll)
     if (abs(jr)>prec) then 
       MSG_ERROR('1st der. of Bessel function did not converge!')
     end if
   end if
   if (order==2) then
     factpp=factp*xx/dble(2*ll+5)
     jnpp=one;jr=one;ii=0
     do while(abs(jr)>=prec.AND.ii<imax)
       ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+5))
       jnpp=jnpp+jr
     end do
     besp1=-jnpp*factpp+jnp*factp*xxinv*dble(ll+1)
     if (abs(jr)>prec) then 
       MSG_ERROR('2nd der. of Bessel function did not converge !')
     end if
   end if
 else
   jn =sin(xx)*xxinv
   jnp=(-cos(xx)+jn)*xxinv
   do il=2,ll+1
     jr=-jn+dble(2*il-1)*jnp*xxinv
     jn=jnp;jnp=jr
   end do
   bes=jn
   if (order>=1) besp =-jnp+jn *xxinv*dble(ll)
   if (order==2) besp1= jn -jnp*xxinv*dble(ll+2)
 end if

 if (order==2) bespp=-besp1+besp*ll*xxinv-bes*ll*xxinv*xxinv

 DBG_EXIT("COLL")

end subroutine jbessel
!!***

!----------------------------------------------------------------------

!!****f* m_paw_numeric/solvbes
!! NAME
!! solvbes
!!
!! FUNCTION
!!    Find nq first roots of instrinsic equation:
!!               alpha.jl(Q) + beta.Q.djl/dr(Q) = 0
!!
!! INPUTS
!!  alpha,beta= factors in intrinsic equation
!!  ll= l quantum number
!!  nq= number of roots to find
!!
!! OUTPUT
!!  root(nq)= roots of instrinsic equation
!!
!! PARENTS
!!      m_atompaw
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

 subroutine solvbes(root,alpha,beta,ll,nq)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'solvbes'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: ll,nq
 real(dp) :: alpha,beta
!arrays
 real(dp) :: root(nq)

!Local variables-------------------------------
!scalars
 integer :: nroot
 real(dp),parameter :: dh=0.1_dp,tol=tol14
 real(dp) :: dum,hh,jbes,jbesp,qq,qx,y1,y2

! *************************************************************************

 DBG_ENTER("COLL")

 qq=dh;nroot=0

 do while (nroot<nq)
   call jbessel(jbes,jbesp,dum,ll,1,qq)
   y1=alpha*jbes+beta*qq*jbesp
   qq=qq+dh
   call jbessel(jbes,jbesp,dum,ll,1,qq)
   y2=alpha*jbes+beta*qq*jbesp

   do while (y1*y2>=zero)
     qq=qq+dh
     call jbessel(jbes,jbesp,dum,ll,1,qq)
     y2=alpha*jbes+beta*qq*jbesp
   end do

   hh=dh;qx=qq
   do while (hh>tol)
     hh=half*hh
     if (y1*y2<zero) then
       qx=qx-hh
     else
       qx=qx+hh
     end if
     call jbessel(jbes,jbesp,dum,ll,1,qx)
     y2=alpha*jbes+beta*qx*jbesp
   end do
   nroot=nroot+1
   root(nroot)=qx

 end do

 DBG_EXIT("COLL")

end subroutine solvbes
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/jbessel_4spline
!! NAME
!!  jbessel_4spline
!!
!! FUNCTION
!!  Compute spherical Bessel functions and derivatives. 
!!  A polynomial approximation is employed for q-->0.
!!  
!! INPUTS
!!  ll=l-order of the Bessel function
!!  tol=tolerance below which a Polynomial approximation is employed
!!   both for jl and its derivative (if required)
!!  order=1 if only first derivative is requested
!!        2 if first and second derivatives are requested
!!  xx=where to compute j_l
!!
!! OUTPUT
!!  bes=Spherical Bessel function j_l at xx
!!  besp= first derivative of j_l at xx (only if order>=1)
!!
!! TODO 
!! Remove inline definitions, they are obsolete in F2003
!!
!! PARENTS
!!      m_paw_pwij,m_pawpsp
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine jbessel_4spline(bes,besp,ll,order,xx,tol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jbessel_4spline'
!End of the abilint section

 implicit none


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jbessel_4spline'
!End of the abilint section

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,order
 real(dp),intent(in) :: xx,tol
 real(dp),intent(out) :: bes,besp

!Local variables ---------------------------------------
!scalars
 real(dp) :: bespp
 real(dp) :: arg,bes0a,bes0ap,bes0b,bes0bp,bes1a,bes1ap,bes1b,bes1bp
 real(dp) :: bes2a,bes2ap,bes2b,bes2bp,bes3a,bes3ap,bes3b,bes3bp
 character(len=500) :: msg
! *********************************************************************

! === l=0,1,2 and 3 spherical Bessel functions (and derivatives) ===
 bes0a(arg)=1.0_dp-arg**2/6.0_dp*(1.0_dp-arg**2/20.0_dp)
 bes0b(arg)=sin(arg)/arg
 bes1a(arg)=(10.0_dp-arg*arg)*arg/30.0_dp
 bes1b(arg)=(sin(arg)-arg*cos(arg))/arg**2
 bes2a(arg)=arg*arg/15.0_dp-arg**4/210.0_dp
 bes2b(arg)=((3.0_dp-arg**2)*sin(arg)-3.0_dp*arg*cos(arg))/arg**3
 bes3a(arg)=arg*arg*arg/105.0_dp-arg**5/1890.0_dp+arg**7/83160.0_dp
 bes3b(arg)=(15.0_dp*sin(arg)-15.0_dp*arg*cos(arg)-6.0_dp*arg**2*sin(arg)+arg**3*cos(arg))/arg**4
 bes0ap(arg)=(-10.0_dp+arg*arg)*arg/30.0_dp
 bes0bp(arg)=-(sin(arg)-arg*cos(arg))/arg**2
 bes1ap(arg)=(10.0_dp-3.0_dp*arg*arg)/30.0_dp
 bes1bp(arg)=((arg*arg-2.0_dp)*sin(arg)+2.0_dp*arg*cos(arg))/arg**3
 bes2ap(arg)=(1.0_dp-arg*arg/7.0_dp)*2.0_dp*arg/15.0_dp
 bes2bp(arg)=((4.0_dp*arg*arg-9.0_dp)*sin(arg)+(9.0_dp-arg*arg)*arg*cos(arg))/arg**4
 bes3ap(arg)=(1.0_dp/35-arg*arg/378.0_dp+arg**4/11880.0_dp)*arg*arg
 bes3bp(arg)=((-60.0_dp+27.0_dp*arg*arg-arg**4)*sin(arg)+(60.0_dp*arg-7.0_dp*arg**3)*cos(arg))/arg**5

 ! This is to test jbessel calculation without polynomial approximation for q-->0.
 ! call jbessel(bes,besp,bespp,ll,order,xx)
 ! RETURN

 if (order>2) then 
   MSG_ERROR("Wrong order in jbessel")
 end if

 select case (ll)
 case (0)
   if (xx<TOL) then
     bes=bes0a(xx)
     if (order>=1) besp=bes0ap(xx)
   else
     bes=bes0b(xx)
     if (order>=1) besp=bes0bp(xx)
   end if

 case (1)
  if (xx<TOL) then
    bes=bes1a(xx)
    if (order>=1) besp=bes1ap(xx)
  else
    bes=bes1b(xx)
    if (order>=1) besp=bes1bp(xx)
  end if

 case (2)
   if (xx<TOL) then
     bes=bes2a(xx)
     if (order>=1) besp=bes2ap(xx)
   else
     bes=bes2b(xx)
     if (order>=1) besp=bes2bp(xx)
   end if

 case (3)
   if (xx<TOL) then
     bes=bes3a(xx)
     if (order>=1) besp=bes3ap(xx)
   else
     bes=bes3b(xx)
     if (order>=1) besp=bes3bp(xx)
   end if

 case (4:)
   call jbessel(bes,besp,bespp,ll,order,xx)

 case default
   write(msg,'(a,i4)')' wrong value for ll = ',ll
   MSG_BUG(msg)
 end select

end subroutine jbessel_4spline
!!***

!----------------------------------------------------------------------

end module m_paw_numeric
!!***
