!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_paw_numeric
!! NAME
!!  m_paw_numeric
!!
!! FUNCTION
!!  Wrapper for splines operations
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2013 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"

module m_paw_numeric
    
 use defs_basis
 use m_pawrad, only: pawrad_type
 use m_pawtab, only: pawtab_type
 use m_errors
use interfaces_12_hide_mpi
use interfaces_14_hidewrite
use interfaces_16_hideleave
 use m_profiling

 implicit none

 private
 public:: paw_smooth
 public:: paw_spline
 public:: paw_splint
 public:: jbessel
 public:: solvbes
 public:: jbessel_4spline   ! Spherical Bessel functions and derivatives employing a polynomial approximation for q->0

!!***

CONTAINS
!===========================================================
!!***

!!****f* m_paw_numeric/paw_spline
!! NAME
!!  paw_spline
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
!! PARENTS
!!      atomden,calc_sig_cd_eet,calc_sigc_cd,calc_sigc_pole_cd,cc_derivatives
!!      denfgr,hirsh,init_bess_spl,init_occ_ent,integrho,m_atom,m_paw_pwij
!!      m_paw_slater,m_paw_toolbox,m_splines,optics_paw_core,pawdij0,pawinit
!!      pawkij,predict_string,psp10in,psp10nl,psp11nl,psp17in,psp1cc,psp1in
!!      psp1nl,psp2in,psp2nl,psp3in,psp3nl,psp4cc,psp5in,psp5nl,psp6cc,psp6in
!!      psp7calc,psp7cc,psp7cc_wvl,psp7in,psp7nl,psp7wvl1,psp8in,psp8lo,psp8nl
!!      psp9in,spline_paw_fncs,upf2abinit,vso_realspace_local
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_spline( t, y, n, ybcbeg, ybcend, ypp )

#if defined HAVE_LIBPAW_ABINIT
!use m_splines, only: spline
use interfaces_28_numeric_noabirule
#endif

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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_spline'
!End of the abilint section

  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: t(n)
  double precision, intent(in) :: y(n)
  double precision, intent(in) :: ybcbeg
  double precision, intent(in) :: ybcend

  double precision, intent(out) :: ypp(n)

#ifndef HAVE_LIBPAW_ABINIT
  integer :: ibcbeg
  integer :: ibcend
  integer :: i,k
  double precision :: ratio,pinv
  double precision, allocatable :: tmp(:)
#endif

#if defined HAVE_LIBPAW_ABINIT
 call spline( t, y, n, ybcbeg, ybcend, ypp )
#else

!
!  Check.
!
  if ( n <= 1 ) then
    write(std_out,* ) ' '
    write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
    write(std_out,* ) '  The number of knots must be at least 2.'
    write(std_out,* ) '  The input value of N = ', n
    call wrtout(std_out,"Fatal error",'COLL')
    call leave_new('COLL')
  end if

  ABI_ALLOCATE(tmp,(n))

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write(std_out,* ) ' '
      write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
      write(std_out,* ) '  The knots must be strictly increasing, but'
      write(std_out,* ) '  T(',  i,') = ', t(i)
      write(std_out,* ) '  T(',i+1,') = ', t(i+1)
      call wrtout(std_out,"Fatal error",'COLL')
      call leave_new('COLL')
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

  ABI_DEALLOCATE(tmp)

  return
#endif

end subroutine paw_spline
!!***

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
!!      atomden,calc_sig_cd_eet,calc_sigc_cd,calc_sigc_pole_cd,cc_derivatives
!!      denfgr,m_atom,m_paw_slater,m_paw_toolbox,m_splines,mkcore_inner
!!      mklocl_realspace,optics_paw_core,partial_dos_fractions,pawdij0,pawgylm
!!      pawkij,predict_string,psp17in,psp6cc,psp7calc,psp7cc,psp7in,psp9in
!!      spline_paw_fncs,vso_realspace_local,wffile,wvl_initro
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit,ierr)

#if defined HAVE_LIBPAW_ABINIT
! use m_splines, only: splint
use interfaces_28_numeric_noabirule
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_splint'
!End of the abilint section

 implicit none

 integer, intent(in) :: nfit, nspline
 integer,optional,intent(out) :: ierr
 double precision, intent(in) :: xspline(nspline)
 double precision, intent(in) :: yspline(nspline)
 double precision, intent(in) :: ysplin2(nspline)
 double precision, intent(in) :: xfit(nfit)

 double precision, intent(out) :: yfit(nfit)

#ifndef HAVE_LIBPAW_ABINIT
!local
 integer :: left,i,k,right,my_err
 double precision :: delarg,invdelarg,aa,bb
#endif

#if defined HAVE_LIBPAW_ABINIT
 if(present(ierr)) then
   call splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit,ierr)
 else
   call splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit)
 end if
#else

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
           MSG_ERROR('xfit(1) < xspline(1)')
!           call wrtout(std_out,'xfit(1,'COLL')
!           call leave_new('COLL')
           !my_err=my_err+1
           !exit
         else
           call wrtout(std_out,'xfit not properly ordered','COLL')
           call leave_new('COLL')
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

 if (PRESENT(ierr)) ierr=my_err
#endif

end subroutine paw_splint
!!***

subroutine paw_smooth(a,mesh,it)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_smooth'
!End of the abilint section

 implicit none
!Arguments
 integer, intent(in) :: it,mesh
 real(dp), intent(inout) :: a(mesh)
!Local variables
 real(dp) :: asm(mesh)
 integer :: i,k

 do k=1,it

 asm(5)=0.2d0*(a(3)+a(4)+a(5)+a(6)+a(7))
 asm(mesh-4)=0.2d0*(a(mesh-2)+a(mesh-3)+a(mesh-4)+&
&                 a(mesh-5)+a(mesh-6))
 asm(mesh-3)=0.2d0*(a(mesh-1)+a(mesh-2)+a(mesh-3)+&
&                 a(mesh-4)+a(mesh-5))
 asm(mesh-2)=0.2d0*(a(mesh)+a(mesh-1)+a(mesh-2)+&
&                 a(mesh-3)+a(mesh-4))
 asm(mesh-1)=0.25d0*(a(mesh)+a(mesh-1)+a(mesh-2)+a(mesh-3))
 asm(mesh)=1.0d0/3.0d0*(a(mesh)+a(mesh-1)+a(mesh-2))

 do i=6,mesh-5
 asm(i)=0.1d0*a(i)+0.1d0*(a(i+1)+a(i-1))+&
&         0.1d0*(a(i+2)+a(i-2))+&
&         0.1d0*(a(i+3)+a(i-3))+&
&         0.1d0*(a(i+4)+a(i-4))+&
&         0.05d0*(a(i+5)+a(i-5))
 enddo

 do i=1,mesh
 a(i)=asm(i)
 enddo

 enddo

end subroutine paw_smooth


!{\src2tex{textfont=tt}}
!!****f* m_paw_numeric/jbessel
!! NAME
!! jbessel
!!
!! FUNCTION
!! Compute spherical Bessel function j_l(x) and derivative(s)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2013 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      m_special_funcs,m_vcoul,pawgylm,pawshpfun,pawtwdij_1,pawtwdij_2b
!!      pawtwdij_2e,psp11nl,shapebes,solvbes
!!
!! CHILDREN
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

 if (order>2) then 
   call wrtout(std_out,"Wrong order in jbessel!",'COLL')
   call leave_new('COLL')
 end if

 if (abs(xx)<prec) then
   bes=zero;if (ll==0) bes=one
   if (order>=1) then
     besp=zero;if (ll==1) besp=third
   end if
   if (order==2) then
     bespp=zero
     if (ll==0) bespp=-third
     if (ll==2) bespp=2.d0/15.d0
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
   xx2=0.5d0*xx*xx
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
     call wrtout(std_out,'Bessel function did not converge!','COLL')
     call leave_new('COLL')
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
       call wrtout(std_out,'1st der. of Bessel function did not converge!','COLL')
       call leave_new('COLL')
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
       call wrtout(std_out,'2nd der. of Bessel function did not converge !','COLL')
       call leave_new('COLL')
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

end subroutine jbessel
!!***

!{\src2tex{textfont=tt}}
!!****f* m_paw_numeric/solvbes
!! NAME
!! solvbes
!!
!! FUNCTION
!!    Find nq first roots of instrinsic equation:
!!               alpha.jl(Q) + beta.Q.djl/dr(Q) = 0
!!
!! COPYRIGHT
!! Copyright (C) 1998-2013 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
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
!!      shapebes
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

end subroutine solvbes
!!***




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
!!      m_paw_pwij,pawpsp_nl
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

subroutine jbessel_4spline(bes,besp,ll,order,xx,tol)
!Arguments ---------------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jbessel_4spline'
!End of the abilint section

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
   call wrtout(std_out,"Wrong order in jbessel",'COLL')
   call leave_new('COLL')
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
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end select

end subroutine jbessel_4spline
!!***

!----------------------------------------------------------------------


end module m_paw_numeric
!!***
