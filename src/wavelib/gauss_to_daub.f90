!!****f* BigDFT/gauss_to_daub
!! FUNCTION
!!   Project gaussian functions in a mesh of Daubechies scaling functions
!!   Gives the expansion coefficients of :
!!     factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
!! INPUT
!!   hgrid    step size
!!   factor   normalisation factor
!!   gau_cen  center of gaussian function
!!   gau_a    parameter of gaussian
!!   n_gau    x**n_gau (polynomial degree)
!!   nmax     size of the grid
!!   nwork    size of the work array (ww) >= (nmax+1)*17
!!   periodic the flag for periodic boundary conditions
!!
!! OUTPUT
!!   n_left,n_right  interval where the gaussian is larger than the machine precision
!!   C(:,1)          array of scaling function coefficients:
!!   C(:,2)          array of wavelet coefficients:
!!   WW(:,1),WW(:,2) work arrays that have to be 17 times larger than C
!!   err_norm        normalisation error
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine gauss_to_daub(hgrid,factor,gau_cen,gau_a,n_gau,&!no err, errsuc
     nmax,n_left,n_right,c,err_norm,&                      !no err_wav. nmax instead of n_intvx
     ww,nwork,periodic)                         !added work arrays ww with dimension nwork
  use module_base
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: n_gau,nmax,nwork
  real(gp), intent(in) :: hgrid,factor,gau_cen,gau_a
  real(wp), dimension(0:nwork,2), intent(inout) :: ww 
  integer, intent(out) :: n_left,n_right
  real(gp), intent(out) :: err_norm
  real(wp), dimension(0:nmax,2), intent(out) :: c
  !local variables
  integer :: rightx,leftx,right_t,i0,i,k,length,j
  real(gp) :: a,z0,h,theor_norm2,x,r,coeff,r2,error,fac
  real(dp) :: cn2,tt
  real(wp) :: func
  integer, dimension(0:4) :: lefts,rights
  !include the convolutions filters
  include 'recs16.inc'! MAGIC FILTER  
  include 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
  include 'sym_16.inc'! WAVELET FILTERS

  !rescale the parameters so that hgrid goes to 1.d0  
  a=gau_a/hgrid

  i0=nint(gau_cen/hgrid) ! the array is centered at i0
  z0=gau_cen/hgrid-real(i0,gp)
  h=.125_gp*.5_gp

  !calculate the array sizes;
  !at level 0, positions shifted by i0 
  right_t= ceiling(15.d0*a)

  if (periodic) then
     !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
     !we limit however the folding to one cell on each side (it can be eliminated)
!!     lefts( 0)=max(i0-right_t,-nmax)
!!     rights(0)=min(i0+right_t,2*nmax)

     lefts( 0)=i0-right_t
     rights(0)=i0+right_t


     call gauss_to_scf()

     ! special for periodic case:
     call fold_tail
  else
     ! non-periodic: the Gaussian is bounded by the cell borders
     lefts( 0)=max(i0-right_t,   0)
     rights(0)=min(i0+right_t,nmax)

     call gauss_to_scf

     ! non-periodic: no tails to fold
     do i=0,length-1
        c(i+n_left,1)=ww(i       ,2) !n_left..n_right <->    0  ..  length-1
        c(i+n_left,2)=ww(i+length,2) !n_left..n_right <-> length..2*length-1
     end do
  endif

  !calculate the (relative) error
  cn2=0.0_dp
  do i=0,length*2-1
     tt=real(ww(i,2),dp)
     cn2=cn2+tt**2
  end do

  theor_norm2=valints(n_gau)*a**(2*n_gau+1)

  error=sqrt(abs(1.0_gp-real(cn2,gp)/theor_norm2))

  !write(*,*)'error, non scaled:',error
  !
  !RESCALE BACK THE COEFFICIENTS AND THE ERROR
  fac= hgrid**n_gau*sqrt(hgrid)*factor
  c=real(fac,wp)*c
  err_norm=error*fac

contains

  subroutine gauss_to_scf
    ! Once the bounds LEFTS(0) and RIGHTS(0) of the expansion coefficient array
    ! are fixed, we get the expansion coefficients in the usual way:
    ! get them on the finest grid by quadrature
    ! then forward transform to get the coeffs on the coarser grid.
    ! All this is done assuming nonperiodic boundary conditions
    ! but will also work in the periodic case if the tails are folded
    n_left=lefts(0)
    n_right=rights(0)
    length=n_right-n_left+1

    !print *,'nleft,nright',n_left,n_right

    do k=1,4
       rights(k)=2*rights(k-1)+m
       lefts( k)=2*lefts( k-1)-m
    enddo

    leftx = lefts(4)-n
    rightx=rights(4)+n  

    !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
  
    !corrected for avoiding 0**0 problem
    if (n_gau == 0) then
       do i=leftx,rightx
          x=real(i-i0*16,gp)*h
          r=x-z0
          r2=r/a
          r2=r2*r2
          r2=0.5_gp*r2
          func=real(dexp(-real(r2,kind=8)),wp)
          ww(i-leftx,1)=func
       enddo
    else
       do i=leftx,rightx
          x=real(i-i0*16,gp)*h
          r=x-z0
          coeff=r**n_gau
          r2=r/a
          r2=r2*r2
          r2=0.5_gp*r2
          func=real(dexp(-real(r2,kind=8)),wp)
          func=real(coeff,wp)*func
          ww(i-leftx,1)=func
       enddo
    end if

    call apply_w(ww(:,1),ww(:,2),&
         leftx   ,rightx   ,lefts(4),rights(4),h)

    call forward_c(ww(:,2),ww(:,1),&
         lefts(4),rights(4),lefts(3),rights(3)) 
    call forward_c(ww(:,1),ww(:,2),&
         lefts(3),rights(3),lefts(2),rights(2)) 
    call forward_c(ww(:,2),ww(:,1),&
         lefts(2),rights(2),lefts(1),rights(1)) 

    call forward(  ww(:,1),ww(:,2),&
         lefts(1),rights(1),lefts(0),rights(0)) 

    c=0.0_gp

  end subroutine gauss_to_scf

  subroutine fold_tail
    ! One of the tails of the Gaussian is folded periodically
    ! We assume that the situation when we need to fold both tails
    ! will never arise
    !implicit none


    !modification of the calculation.
    !at this stage the values of c are fixed to zero

    do i=n_left,n_right
       j=modulo(i,nmax+1)
       c(j,1)=c(j,1)+ww(i-n_left       ,2)
       c(j,2)=c(j,2)+ww(i-n_left+length,2)
    end do


!!    !write(*,*) 'I fold the tail'
!!    ! shift the resulting array and fold its periodic tails:
!!    if (n_left.ge.0) then
!!       if (n_right.le.nmax) then
!!          ! normal situation: the gaussian is inside the box
!!          do i=n_left,n_right
!!             c(i,1)=ww(i-n_left       ,2)
!!             c(i,2)=ww(i-n_left+length,2)
!!          enddo
!!       else
!!          ! the gaussian extends beyond the right border
!!
!!          ! the normal part:
!!          do i=n_left,nmax
!!             c(i,1)=ww(i-n_left       ,2)
!!             c(i,2)=ww(i-n_left+length,2)
!!          enddo
!!          ! the part of ww that goes beyond nmax 
!!          ! is shifted by nmax+1 to the left
!!          do i=nmax+1,n_right
!!             c(i-nmax-1,1)=ww(i-n_left       ,2)
!!             c(i-nmax-1,2)=ww(i-n_left+length,2)
!!          enddo
!!       endif
!!    else
!!       ! the gaussian extends beyond the left border
!!       ! the part of ww to the left of 0
!!       ! is shifted by nmax+1 to the right
!!       do i=n_left,-1
!!          c(i+nmax+1,1)=ww(i-n_left       ,2)
!!          c(i+nmax+1,2)=ww(i-n_left+length,2)
!!       enddo
!!       ! the normal part:
!!       do i=0,n_right
!!          c(i,1)=ww(i-n_left       ,2)
!!          c(i,2)=ww(i-n_left+length,2)
!!       enddo
!!    endif
  end subroutine fold_tail


end subroutine gauss_to_daub
!!***

!!****f* BigDFT/gauss_to_daub_k
!! FUNCTION
!!   Project gaussian functions in a mesh of Daubechies scaling functions
!!   Gives the expansion coefficients of :
!!     factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
!!   Multiply it for the k-point factor exp(Ikx)
!!   For this reason, a real (cos(kx)) and an imaginary (sin(kx)) part are provided 
!!   WARNING: in this version, we dephase the projector to wrt the center of the gaussian
!!            this should not have an impact on the results since the operator is unchanged
!! INPUT
!!   hgrid    step size
!!   factor   normalisation factor
!!   gau_cen  center of gaussian function
!!   gau_a    parameter of gaussian
!!   n_gau    x**n_gau (polynomial degree)
!!   nmax     size of the grid
!!   nwork    size of the work array (ww) >= (nmax+1)*17
!!   periodic the flag for periodic boundary conditions
!!   kval     value for the k-point
!!   ncplx    number of components in the complex direction (must be 2 if kval /=0)
!!
!! OUTPUT
!!   n_left,n_right  interval where the gaussian is larger than the machine precision
!!   C(:,1)          array of scaling function coefficients:
!!   C(:,2)          array of wavelet coefficients:
!!   WW(:,1),WW(:,2) work arrays that have to be 17 times larger than C
!!   err_norm        normalisation error
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine gauss_to_daub_k(hgrid,kval,ncplx,factor,gau_cen,gau_a,n_gau,&!no err, errsuc
     nmax,n_left,n_right,c,& 
     ww,nwork,periodic)      !added work arrays ww with dimension nwork
  use module_base
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: n_gau,nmax,nwork,ncplx
  real(gp), intent(in) :: hgrid,factor,gau_cen,gau_a,kval
  real(wp), dimension(0:nwork,2,ncplx), intent(inout) :: ww 
  integer, intent(out) :: n_left,n_right
  real(wp), dimension(ncplx,0:nmax,2), intent(out) :: c
  !local variables
  integer :: rightx,leftx,right_t,i0,i,k,length,j,icplx
  real(gp) :: a,z0,h,x,r,coeff,r2,fac,rk
  real(wp) :: func,cval,sval
  integer, dimension(0:4) :: lefts,rights
  !include the convolutions filters
  include 'recs16.inc'! MAGIC FILTER  
  include 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
  include 'sym_16.inc'! WAVELET FILTERS

  !rescale the parameters so that hgrid goes to 1.d0  
  a=gau_a/hgrid

  i0=nint(gau_cen/hgrid) ! the array is centered at i0
  z0=gau_cen/hgrid-real(i0,gp)
  h=.125_gp*.5_gp

  !calculate the array sizes;
  !at level 0, positions shifted by i0 
  right_t= ceiling(15.d0*a)

  !print *,'a,right_t',a,right_t,gau_a,hgrid

  !to rescale back the coefficients
  fac=hgrid**n_gau*sqrt(hgrid)*factor


  if (periodic) then
     !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
     !we limit however the folding to one cell on each side 
     !!(commented out)
     !!     lefts( 0)=max(i0-right_t,-nmax)
     !!     rights(0)=min(i0+right_t,2*nmax)

     lefts( 0)=i0-right_t
     rights(0)=i0+right_t
     
     
     call gauss_to_scf()
     
     ! special for periodic case:
     call fold_tail
  else
     ! non-periodic: the Gaussian is bounded by the cell borders
     lefts( 0)=max(i0-right_t,   0)
     rights(0)=min(i0+right_t,nmax)
     
     call gauss_to_scf
     
     !loop for each complex component
     do icplx=1,ncplx
        ! non-periodic: no tails to fold
        do i=0,length-1
           c(icplx,i+n_left,1)=fac*ww(i       ,2,icplx)
           c(icplx,i+n_left,2)=fac*ww(i+length,2,icplx) 
        end do
     end do
  endif


contains

  subroutine gauss_to_scf
    ! Once the bounds LEFTS(0) and RIGHTS(0) of the expansion coefficient array
    ! are fixed, we get the expansion coefficients in the usual way:
    ! get them on the finest grid by quadrature
    ! then forward transform to get the coeffs on the coarser grid.
    ! All this is done assuming nonperiodic boundary conditions
    ! but will also work in the periodic case if the tails are folded
    n_left=lefts(0)
    n_right=rights(0)
    length=n_right-n_left+1

    !print *,'nleft,nright',n_left,n_right

    do k=1,4
       rights(k)=2*rights(k-1)+m
       lefts( k)=2*lefts( k-1)-m
    enddo

    leftx = lefts(4)-n
    rightx=rights(4)+n  

    !loop for each complex component
    do icplx=1,ncplx

       !calculate the expansion coefficients at level 4, positions shifted by 16*i0 

       !corrected for avoiding 0**0 problem
       if (ncplx==1) then
          if (n_gau == 0) then
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                func=real(dexp(-real(r2,kind=8)),wp)
                ww(i-leftx,1,icplx)=func
             enddo
          else
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                coeff=r**n_gau
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                func=real(dexp(-real(r2,kind=8)),wp)
                func=real(coeff,wp)*func
                ww(i-leftx,1,icplx)=func
             enddo
          end if
       else if (icplx == 1) then
          if (n_gau == 0) then
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                rk=real(i,gp)*h
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                cval=real(cos(kval*rk),wp)
                func=real(dexp(-real(r2,kind=8)),wp)
                ww(i-leftx,1,icplx)=func*cval
             enddo
          else
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                rk=real(i,gp)*h
                coeff=r**n_gau
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                cval=real(cos(kval*rk),wp)
                func=real(dexp(-real(r2,kind=8)),wp)
                func=real(coeff,wp)*func
                ww(i-leftx,1,icplx)=func*cval
             enddo
          end if
       else if (icplx == 2) then
          if (n_gau == 0) then
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                rk=real(i,gp)*h
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                sval=real(sin(kval*rk),wp)
                func=real(dexp(-real(r2,kind=8)),wp)
                ww(i-leftx,1,icplx)=func*sval
             enddo
          else
             do i=leftx,rightx
                x=real(i-i0*16,gp)*h
                r=x-z0
                rk=real(i,gp)*h
                coeff=r**n_gau
                r2=r/a
                r2=r2*r2
                r2=0.5_gp*r2
                sval=real(sin(kval*rk),wp)
                func=real(dexp(-real(r2,kind=8)),wp)
                func=real(coeff,wp)*func
                ww(i-leftx,1,icplx)=func*sval
             enddo
          end if
       end if

       !print *,'here',gau_a,gau_cen,n_gau
       call apply_w(ww(:,1,icplx),ww(:,2,icplx),&
            leftx   ,rightx   ,lefts(4),rights(4),h)

       call forward_c(ww(:,2,icplx),ww(:,1,icplx),&
            lefts(4),rights(4),lefts(3),rights(3)) 
       call forward_c(ww(:,1,icplx),ww(:,2,icplx),&
            lefts(3),rights(3),lefts(2),rights(2)) 
       call forward_c(ww(:,2,icplx),ww(:,1,icplx),&
            lefts(2),rights(2),lefts(1),rights(1)) 

       call forward(  ww(:,1,icplx),ww(:,2,icplx),&
            lefts(1),rights(1),lefts(0),rights(0)) 

    end do


    c=0.0_gp

  end subroutine gauss_to_scf

  subroutine fold_tail
    ! One of the tails of the Gaussian is folded periodically
    ! We assume that the situation when we need to fold both tails
    ! will never arise
    !implicit none


    !modification of the calculation.
    !at this stage the values of c are fixed to zero
    do icplx=1,ncplx
       do i=n_left,n_right
          j=modulo(i,nmax+1)
          c(icplx,j,1)=c(icplx,j,1)+ww(i-n_left       ,2,icplx)
          c(icplx,j,2)=c(icplx,j,2)+ww(i-n_left+length,2,icplx)
       end do
    end do

    c=fac*c

  end subroutine fold_tail


end subroutine gauss_to_daub_k
!!***



!
!       APPLYING THE MAGIC FILTER ("SHRINK") 
!
subroutine apply_w(cx,c,leftx,rightx,left,right,h)
  use module_base
  implicit none
  integer, intent(in) :: leftx,rightx,left,right
  real(gp), intent(in) :: h
  real(wp), dimension(leftx:rightx), intent(in) :: cx
  real(wp), dimension(left:right), intent(out) :: c
  !local variables
  include 'recs16.inc'
  integer :: i,j
  real(wp) :: sqh,ci

  sqh=real(sqrt(h),wp)

  do i=left,right
     ci=0.0_wp
     do j=-n,n
        ci=ci+cx(i+j)*w(j)         
     enddo
     c(i)=ci*sqh
  enddo

end subroutine apply_w


!
!      FORWARD WAVELET TRANSFORM WITHOUT WAVELETS ("SHRINK")
!
subroutine forward_c(c,c_1,left,right,left_1,right_1)
  use module_base
  implicit none
  integer, intent(in) :: left,right,left_1,right_1
  real(wp), dimension(left:right), intent(in) :: c
  real(wp), dimension(left_1:right_1), intent(out) :: c_1
  !local variables
  integer :: i,i2,j
  real(wp) :: ci
  include 'sym_16.inc'

  !
  !      get the coarse scfunctions and wavelets
  !
  do i=left_1,right_1
     i2=2*i
     ci=0.0_wp
     do j=-m,m
        ci=ci+cht(j)*c(j+i2)
     enddo
     c_1(i)=ci
  enddo

end subroutine forward_c

!
!      CONVENTIONAL FORWARD WAVELET TRANSFORM ("SHRINK")
!
subroutine forward(c,cd_1,left,right,left_1,right_1)
  use module_base
  implicit none
  integer, intent(in) :: left,right,left_1,right_1
  real(wp), dimension(left:right), intent(in) :: c
  real(wp), dimension(left_1:right_1,2), intent(out) :: cd_1
  !local variables
  integer :: i,i2,j
  real(wp) :: ci,di
  include 'sym_16.inc'

  !
  !      get the coarse scfunctions and wavelets
  !
  do i=left_1,right_1
     i2=2*i
     ci=0.d0
     di=0.d0
     do j=-m,m
        ci=ci+cht(j)*c(j+i2)
        di=di+cgt(j)*c(j+i2)
     enddo
     cd_1(i,1)=ci
     cd_1(i,2)=di
  enddo

end subroutine forward
