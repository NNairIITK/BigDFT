!> @file
!!  Gaussian to Daubechies wavelets projection routines
!! @author
!!    Copyright (C) 2007-2014 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module gaussdaub
  use module_defs, only: wp,gp
  implicit none
  private
  
  !these parameters has to be modified to avoid naming confusions
  integer, parameter :: m=8,mm=m+2
  integer, parameter :: N=m
  integer :: iw !< index for initialization

  real(wp), dimension(-mm:mm), parameter :: CH = (/ &
       0.0_wp,0.0_wp,0.0_wp, &
       -0.0033824159510050025955_wp,-0.00054213233180001068935_wp, &
       0.031695087811525991431_wp,0.0076074873249766081919_wp, &
       -0.14329423835127266284_wp,-0.061273359067811077843_wp, &
       0.48135965125905339159_wp,0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp,0.0_wp,0.0_wp &
       /)
  ! Coefficients for wavelet transform (orthonormal wavelet)
  real(wp), dimension(-mm:mm), parameter :: cht=ch

  ! g coefficients from h coefficients
  real(wp), dimension(-mm:mm), parameter :: cg=[0.0_wp,((-1)**(iw+1)*cht(-iw),iw=-mm,mm-1)]
  real(wp), dimension(-mm:mm), parameter :: cgt=[0.0_wp,((-1)**(iw+1)*ch(-iw),iw=-mm,mm-1)]

  !magic filter coefficients for daubechies family
  real(wp), dimension(-n:n), parameter :: W = (/0.0_wp,&
       2.72734492911979659657715313017228e-6_wp,&
       -0.00005185986881173432922848639136911487_wp,&
       0.00049443227688689919192282259476750972_wp,&
       -0.00344128144493493857280881509686821861_wp,&
       0.01337263414854794752733423467013220997_wp,&
       -0.02103025160930381434955489412839065067_wp,&
       -0.0604895289196983516002834636_wp,0.9940415697834003993178616713_wp,&
       0.0612625895831207982195380597_wp,&
       0.02373821463724942397566389712597274535_wp,&
       -0.00942047030201080385922711540948195075_wp,&
       0.00174723713672993903449447812749852942_wp,&
       -0.00030158038132690463167163703826169879_wp,&
       0.00008762984476210559564689161894116397_wp,&
       -0.00001290557201342060969516786758559028_wp,&
       8.4334247333529341094733325815816e-7_wp&
       /)
  
  real(wp), dimension(0:6), parameter :: valints = (/& 
       1.7724538509055160272981674834_wp,0.8862269254527580136490837417_wp,&
       1.3293403881791370204736256126_wp,3.3233509704478425511840640314_wp,&
       11.63172839656744892914422411_wp,52.3427777845535201811490084949_wp,&
       287.8852778150443609963195467219_wp/)

  contains
    
    !convert a gaussian to one-dimensional functions
    subroutine gau_daub_1d(ncplx,ng,gau_cen,gau_a,n_gau,nmax,hgrid,factor,periodic,nres,&
         c,nwork,ww)
      implicit none
      logical, intent(in) :: periodic !< determine the bc
      integer, intent(in) :: nmax,nwork,nres
      !> 1 for real gaussians, 2 for complex ones 
      !! (to be generalized to k-points)
      integer, intent(in) :: ncplx
      integer, intent(in) :: ng !<number of different gaussians to convert
      !>principal quantum numbers of any of the gaussians
      integer, dimension(ng), intent(in) :: n_gau 
      !>multiplicative factors which have to be added to the different
      !!terms
      real(wp), dimension(ng), intent(in) :: factor
      !>standard deviations of the different gaussians (might be complex)
      real(wp), dimension(ng), intent(in) :: gau_a
      real(gp), intent(in) :: hgrid,gau_cen
      real(wp), dimension(ng,0:nmax,2), intent(inout) :: c
      real(wp), dimension(nwork,2), intent(inout) :: ww
      !local variables
      integer :: i0,right_t,inw,ig
      real(wp) :: x0
      integer, dimension(0:nres+1) :: lefts,rights !< use automatic arrays here, we should use parameters in the module
      real(gp), dimension(ng) :: a,fac,theor_norm2,error


      !first, determine the cutoffs where the 
      !bounds are to be calculated
      x0=gau_cen/hgrid
      i0=nint(x0) ! the array is centered at i0
      !here the are the quantities for any of the objects
      a=gau_a/hgrid
      right_t=ceiling(15.d0*maxval(a))
      !the multiplicative factors for any of the object
      do ig=1,ng
         fac(ig)=hgrid**n_gau(ig)*sqrt(hgrid)*factor(ig)
         !and the expected theoretical norm of the gaussian
         theor_norm2(ig)=valints(n_gau(ig))*a(ig)**(2*n_gau(ig)+1)
      end do
      
      !these are the limits of each of the convolutions
      call determine_bounds(nres,periodic,nmax,right_t,i0,&
         lefts,rights)

      !then check if the work array has a size which is big enough
      if (rights(nres+1) -lefts(nres+1) > nwork) then
         STOP 'gaustoisf'
      end if

      !fill the array of the high resolution values
      call gaus_highres(nres,ng,n_gau,a,x0,&
           lefts(nres+1),rights(nres+1),ww)

      !then come back to the original resolution level
      call magic_idwts(nres,ng,lefts,rights,nwork,ww,inw)

      !and retrieve the results and the error if desired
      call retrieve_results(periodic,ng,lefts(0),rights(0),nmax,ww(1,inw),&
           theor_norm2(1),fac(1),c,error(1))

    end subroutine gau_daub_1d

    !>   Project gaussian functions in a mesh of Daubechies scaling functions
    !!   Gives the expansion coefficients of :
    !!     factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
    !! INPUT
    !!   @param hgrid    step size
    !!   @param factor   normalisation factor
    !!   @param gau_cen  center of gaussian function
    !!   @param gau_a    parameter of gaussian
    !!   @param n_gau    x**n_gau (polynomial degree)
    !!   @param nmax     size of the grid
    !!   @param nwork    size of the work array (ww) >= (nmax+1)*17
    !!   @param periodic the flag for periodic boundary conditions
    !!
    !! OUTPUT
    !!   @param n_left,n_right  interval where the gaussian is larger than the machine precision
    !!   @param C(:,1)          array of scaling function coefficients:
    !!   @param C(:,2)          array of wavelet coefficients:
    !!   @param WW(:,1),WW(:,2) work arrays that have to be 17 times larger than C
    !!   @param err_norm        normalisation error
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
      real(gp) :: a,z0,h,theor_norm2,x,r,coff,r2,error,fac
      real(dp) :: cn2,tt
      real(wp) :: func
      integer, dimension(0:4) :: lefts,rights
      !include the convolutions filters
      !include 'recs16.inc' !< MAGIC FILTER  
      !include 'intots.inc' !< HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
      !include 'sym_16.inc' !< WAVELET FILTERS

      !rescale the parameters so that hgrid goes to 1.d0  
      a=gau_a/hgrid

      i0=nint(gau_cen/hgrid) ! the array is centered at i0
      z0=gau_cen/hgrid-real(i0,gp)
      h=.125_gp*.5_gp

      !calculate the array sizes;
      !at level 0, positions shifted by i0 
      right_t= ceiling(15.d0*a)

      ! initialise array
      c=0.0_gp


      if (periodic) then
         !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
         !we limit however the folding to one cell on each side (it can be eliminated)
         !!     lefts( 0)=max(i0-right_t,-nmax)
         !!     rights(0)=min(i0+right_t,2*nmax)

         lefts( 0)=i0-right_t
         rights(0)=i0+right_t


!!$         call gauss_to_scf()

         ! special for periodic case:
!!$         call fold_tail
      else
         ! non-periodic: the Gaussian is bounded by the cell borders
         lefts( 0)=max(i0-right_t,   0)
         rights(0)=min(i0+right_t,nmax)

!!$         call gauss_to_scf

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
      
!!$    contains
!!$
!!$      !> Once the bounds LEFTS(0) and RIGHTS(0) of the expansion coefficient array
!!$      !! are fixed, we get the expansion coefficients in the usual way:
!!$      !! get them on the finest grid by quadrature
!!$      !! then forward transform to get the coeffs on the coarser grid.
!!$      !! All this is done assuming nonperiodic boundary conditions
!!$      !! but will also work in the periodic case if the tails are folded
!!$      subroutine gauss_to_scf
!!$        n_left=lefts(0)
!!$        n_right=rights(0)
!!$        length=n_right-n_left+1
!!$
!!$        !print *,'nleft,nright',n_left,n_right
!!$
!!$        do k=1,4
!!$           rights(k)=2*rights(k-1)+m
!!$           lefts( k)=2*lefts( k-1)-m
!!$        enddo
!!$
!!$        leftx = lefts(4)-n
!!$        rightx=rights(4)+n  
!!$
!!$        !do not do anything if the gaussian is too extended
!!$        if (rightx-leftx > nwork) then
!!$           !STOP 'gaustodaub'
!!$           return
!!$        end if
!!$
!!$        !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
!!$
!!$        !corrected for avoiding 0**0 problem
!!$        if (n_gau == 0) then
!!$           do i=leftx,rightx
!!$              x=real(i-i0*16,gp)*h
!!$              r=x-z0
!!$              r2=r/a
!!$              r2=r2*r2
!!$              r2=0.5_gp*r2
!!$              func=safe_exp(-r2)
!!$              ww(i-leftx,1)=func
!!$           enddo
!!$        else
!!$           do i=leftx,rightx
!!$              x=real(i-i0*16,gp)*h
!!$              r=x-z0
!!$              coeff=r**n_gau
!!$              r2=r/a
!!$              r2=r2*r2
!!$              r2=0.5_gp*r2
!!$              func=safe_exp(-r2)
!!$              func=real(coeff,wp)*func
!!$              ww(i-leftx,1)=func
!!$           enddo
!!$        end if
!!$
!!$        call apply_w(1,ww(:,1),ww(:,2),&
!!$             leftx   ,rightx   ,lefts(4),rights(4),h)
!!$
!!$        call forward_c(1,ww(0,2),ww(0,1),&
!!$             lefts(4),rights(4),lefts(3),rights(3)) 
!!$        call forward_c(1,ww(0,1),ww(0,2),&
!!$             lefts(3),rights(3),lefts(2),rights(2)) 
!!$        call forward_c(1,ww(0,2),ww(0,1),&
!!$             lefts(2),rights(2),lefts(1),rights(1)) 
!!$
!!$        call forward(1,  ww(0,1),ww(0,2),&
!!$             lefts(1),rights(1),lefts(0),rights(0)) 
!!$
!!$
!!$      END SUBROUTINE gauss_to_scf
!!$
!!$
!!$      !> One of the tails of the Gaussian is folded periodically
!!$      !! We assume that the situation when we need to fold both tails
!!$      !! will never arise
!!$      subroutine fold_tail
!!$
!!$        !modification of the calculation.
!!$        !at this stage the values of c are fixed to zero
!!$
!!$        do i=n_left,n_right
!!$           j=modulo(i,nmax+1)
!!$           c(j,1)=c(j,1)+ww(i-n_left       ,2)
!!$           c(j,2)=c(j,2)+ww(i-n_left+length,2)
!!$        end do
!!$
!!$        !!
!!$        !!    !write(*,*) 'I fold the tail'
!!$        !!    ! shift the resulting array and fold its periodic tails:
!!$        !!    if (n_left.ge.0) then
!!$        !!       if (n_right.le.nmax) then
!!$        !!          ! normal situation: the gaussian is inside the box
!!$        !!          do i=n_left,n_right
!!$        !!             c(i,1)=ww(i-n_left       ,2)
!!$        !!             c(i,2)=ww(i-n_left+length,2)
!!$        !!          enddo
!!$        !!       else
!!$        !!          ! the gaussian extends beyond the right border
!!$        !!
!!$        !!          ! the normal part:
!!$        !!          do i=n_left,nmax
!!$        !!             c(i,1)=ww(i-n_left       ,2)
!!$        !!             c(i,2)=ww(i-n_left+length,2)
!!$        !!          enddo
!!$        !!          ! the part of ww that goes beyond nmax 
!!$        !!          ! is shifted by nmax+1 to the left
!!$        !!          do i=nmax+1,n_right
!!$        !!             c(i-nmax-1,1)=ww(i-n_left       ,2)
!!$        !!             c(i-nmax-1,2)=ww(i-n_left+length,2)
!!$        !!          enddo
!!$        !!       endif
!!$        !!    else
!!$        !!       ! the gaussian extends beyond the left border
!!$        !!       ! the part of ww to the left of 0
!!$        !!       ! is shifted by nmax+1 to the right
!!$        !!       do i=n_left,-1
!!$        !!          c(i+nmax+1,1)=ww(i-n_left       ,2)
!!$        !!          c(i+nmax+1,2)=ww(i-n_left+length,2)
!!$        !!       enddo
!!$        !!       ! the normal part:
!!$        !!       do i=0,n_right
!!$        !!          c(i,1)=ww(i-n_left       ,2)
!!$        !!          c(i,2)=ww(i-n_left+length,2)
!!$        !!       enddo
!!$        !!    endif
!!$      END SUBROUTINE fold_tail

    END SUBROUTINE gauss_to_daub

    !>   Project gaussian functions in a mesh of Daubechies scaling functions
    !!   Gives the expansion coefficients of :
    !!     factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
    !!   Multiply it for the k-point factor exp(Ikx)
    !!   For this reason, a real (cos(kx)) and an imaginary (sin(kx)) part are provided 
    !! INPUT
    !!   @param hgrid    step size
    !!   @param factor   normalisation factor
    !!   @param gau_cen  center of gaussian function
    !!   @param gau_a    parameter of gaussian
    !!   @param n_gau    x**n_gau (polynomial degree)
    !!   @param nmax     size of the grid
    !!   @param nwork    size of the work array (ww) >= (nmax+1)*17
    !!   @param periodic the flag for periodic boundary conditions
    !!   @param kval     value for the k-point
    !!   @param ncplx    number of components in the complex direction (must be 2 if kval /=0)
    !!
    !! OUTPUT
    !!   @param n_left,n_right  interval where the gaussian is larger than the machine precision
    !!   @param C(:,1)          array of scaling function coefficients:
    !!   @param C(:,2)          array of wavelet coefficients:
    !!   @param WW(:,1),WW(:,2) work arrays that have to be 17 times larger than C
    !!   @param err_norm        normalisation error
    !!@warning 
    !!  In this version, we dephase the projector to wrt the center of the gaussian
    !!  this should not have an impact on the results since the operator is unchanged
    subroutine gauss_to_daub_k(hgrid,kval,ncplx_w,ncplx_g,ncplx_k,&
         factor,gau_cen,gau_a,n_gau,&!no err, errsuc
         nstart,nmax,n_left,n_right,c,& 
         ww,nwork,periodic,gau_cut)      !added work arrays ww with dimension nwork
      use module_base
      !use gaussians, only: mp_exp
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: n_gau,nmax,nwork,nstart
      integer, intent(in) :: ncplx_w !size of the ww matrix
      integer, intent(in) :: ncplx_g !1 or 2 for simple or complex gaussians, respectively.
      integer, intent(in) :: ncplx_k !use 2 for k-points.
      real(gp), intent(in) :: hgrid,gau_cen,kval
      real(gp),dimension(ncplx_g),intent(in)::factor,gau_a
      real(wp), dimension(0:nwork,2,ncplx_w), intent(inout) :: ww 
      integer, intent(out) :: n_left,n_right
      real(wp), dimension(ncplx_w,0:nmax,2), intent(out) :: c
      real(gp), intent(in) :: gau_cut
      !local variables
      character(len=*), parameter :: subname='gauss_to_daub_k'
      integer :: rightx,leftx,right_t,i0,i,k,length,j,icplx
      real(gp) :: a1,a2,z0,h,x,r,coff,r2,rk,gcut
      real(gp) :: fac(ncplx_g)
      real(wp) :: func,cval,sval,cval2,sval2
      real(wp), dimension(:,:,:), allocatable :: cc
      integer, dimension(0:4) :: lefts,rights
      !include the convolutions filters
      !include 'recs16.inc'! MAGIC FILTER  
      !include 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
      !include 'sym_16.inc'! WAVELET FILTERS

      !rescale the parameters so that hgrid goes to 1.d0  
      !when calculating "r2" in gauss_to_scf 
      a1=gau_a(1)/hgrid
      if(ncplx_g==2) then
         a2=gau_a(2)*hgrid*hgrid
         cc = f_malloc((/ 1.to.ncplx_g, 0.to.nmax, 1.to.2 /),id='cc')
      end if
      i0=nint(gau_cen/hgrid) ! the array is centered at i0
      z0=gau_cen/hgrid-real(i0,gp)
      gcut=gau_cut/hgrid
      h=.125_gp*.5_gp

      !calculate the array sizes;
      !at level 0, positions shifted by i0 
      right_t= ceiling(15.d0*a1)

      !print *,'a,right_t',a1,right_t,gau_a,hgrid

      !to rescale back the coefficients
      fac(:)=hgrid**n_gau*sqrt(hgrid)*factor(:)

      !initialise array
      c=0.0_gp
      if(ncplx_g==2)cc=0.0_gp

      if (periodic) then
         !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
         !we limit however the folding to one cell on each side 
         !!(commented out)
         !     lefts( 0)=max(i0-right_t,-nmax)
         !     rights(0)=min(i0+right_t,2*nmax)

         lefts( 0)=i0-right_t
         rights(0)=i0+right_t

         call gauss_to_scf()

         ! special for periodic case:
         call fold_tail
      else
         ! non-periodic: the Gaussian is bounded by the cell borders
         lefts( 0)=max(i0-right_t,nstart)
         rights(0)=min(i0+right_t,nmax+nstart)

         call gauss_to_scf()

         n_left = n_left - nstart

         !loop for each complex component
         if(ncplx_g==1) then
            do icplx=1,ncplx_w
               ! non-periodic: no tails to fold
               do i=0,length-1
                  c(icplx,i+n_left,1)=ww(i       ,2,icplx)
                  c(icplx,i+n_left,2)=ww(i+length,2,icplx) 
               end do
            end do
         else !ncplx_g==2
            !use a temporary array cc instead
            do icplx=1,ncplx_w
               ! non-periodic: no tails to fold
               do i=0,length-1
                  cc(icplx,i+n_left,1)=ww(i       ,2,icplx)
                  cc(icplx,i+n_left,2)=ww(i+length,2,icplx) 
               end do
            end do
         end if
      endif

      ! Apply factor:
      if(ncplx_g==1) then
         c=fac(1)*c
      else
         c(1,:,:)=fac(1)*cc(1,:,:)-fac(2)*cc(2,:,:)
         c(2,:,:)=fac(1)*cc(2,:,:)+fac(2)*cc(1,:,:)

         call f_free(cc)
      end if

    contains

      !> Once the bounds LEFTS(0) and RIGHTS(0) of the expansion coefficient array
      !! are fixed, we get the expansion coefficients in the usual way:
      !! get them on the finest grid by quadrature
      !! then forward transform to get the coeffs on the coarser grid.
      !! All this is done assuming nonperiodic boundary conditions
      !! but will also work in the periodic case if the tails are folded

      subroutine gauss_to_scf
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

        !stop the code if the gaussian is too extended
        if (rightx-leftx > nwork) then
           !STOP 'gaustodaub'
           return
        end if

        if (ncplx_w==1) then
           !no kpts and real gaussians
           call gauss_to_scf_1()
        elseif(ncplx_k==2 .and. ncplx_g==1) then
           !kpts and real gaussians
           call gauss_to_scf_2()
        elseif(ncplx_k==1 .and. ncplx_g==2) then
           !no kpts and complex gaussians
           call gauss_to_scf_3()
        elseif(ncplx_k==2 .and. ncplx_g==2) then
           !kpts and complex gaussians
           call gauss_to_scf_4()
        endif

!!$        do icplx=1,ncplx_w
!!$           !print *,'here',gau_a,gau_cen,n_gau
!!$           call apply_w(ww(0,1,icplx),ww(0,2,icplx),&
!!$                leftx   ,rightx   ,lefts(4),rights(4),h)
!!$
!!$           call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
!!$                lefts(4),rights(4),lefts(3),rights(3)) 
!!$           call forward_c(ww(0,1,icplx),ww(0,2,icplx),&
!!$                lefts(3),rights(3),lefts(2),rights(2)) 
!!$           call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
!!$                lefts(2),rights(2),lefts(1),rights(1)) 
!!$
!!$           call forward(  ww(0,1,icplx),ww(0,2,icplx),&
!!$                lefts(1),rights(1),lefts(0),rights(0)) 
!!$
!!$        end do

      END SUBROUTINE gauss_to_scf

      ! Called when ncplx_w = 1
      subroutine gauss_to_scf_1

        !loop for each complex component
        !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
        !corrected for avoiding 0**0 problem
        icplx = 1
        if (n_gau == 0) then
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              r2=r/a1
              r2=r2*r2
              r2=0.5_gp*r2
              func=safe_exp(-r2)
              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
              ww(i-leftx,1,icplx)=func
           enddo
        else
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              coff=r**n_gau
              r2=r/a1
              r2=r2*r2
              r2=0.5_gp*r2
              func=safe_exp(-r2)
              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
              func=coff*func
              ww(i-leftx,1,icplx)=func
           enddo
        end if

      END SUBROUTINE gauss_to_scf_1

      ! Called when ncplx_k = 2 and ncplx_g = 1
      subroutine gauss_to_scf_2

        !loop for each complex component
        !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
        !corrected for avoiding 0**0 problem
        if (n_gau == 0) then
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              rk=real(i,gp)*h
              r2=r/a1
              r2=r2*r2
              r2=0.5_gp*r2
              cval=cos(kval*rk)
              func=safe_exp(-r2)
              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
              ww(i-leftx,1,1)=func*cval
              sval=sin(kval*rk)
              ww(i-leftx,1,2)=func*sval
           enddo
        else
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              rk=real(i,gp)*h
              coff=r**n_gau
              r2=r/a1
              r2=r2*r2
              r2=0.5_gp*r2
              cval=cos(kval*rk)
              func=safe_exp(-r2)
              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
              func=coff*func
              ww(i-leftx,1,1)=func*cval
              sval=sin(kval*rk)
              ww(i-leftx,1,2)=func*sval
           enddo
        end if

      END SUBROUTINE gauss_to_scf_2

      ! Called when ncplx_k = 1 and ncplx_g = 2
      ! no k-points + complex Gaussians
      subroutine gauss_to_scf_3

        if (n_gau == 0) then
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              if( abs(r)-gcut < 1e-8 ) then
                 r2=r*r
                 cval=cos(a2*r2)
                 sval=sin(a2*r2)
                 r2=0.5_gp*r2/(a1**2)
                 func=safe_exp(-r2)
                 ww(i-leftx,1,1)=func*cval
                 ww(i-leftx,1,2)=func*sval
              else
                 ww(i-leftx,1,1:2)=0.0_wp
              end if
           enddo
        else
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              if( abs(r)-gcut < 1e-8 ) then
                 r2=r*r
                 cval=cos(a2*r2)
                 sval=sin(a2*r2)
                 coff=r**n_gau
                 r2=0.5_gp*r2/(a1**2)
                 func=safe_exp(-r2)
                 func=coff*func
                 ww(i-leftx,1,1)=func*cval
                 ww(i-leftx,1,2)=func*sval
              else
                 ww(i-leftx,1,1:2)=0.0_wp
              end if
           enddo
        end if
      END SUBROUTINE gauss_to_scf_3

      ! Called when ncplx_k = 2 and ncplx_g = 2
      subroutine gauss_to_scf_4

        if (n_gau == 0) then
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              if( abs(r)-gcut < 1e-8 ) then
                 r2=r*r
                 cval=cos(a2*r2)
                 sval=sin(a2*r2)
                 rk=real(i,gp)*h
                 cval2=cos(kval*rk)
                 sval2=sin(kval*rk)
                 r2=0.5_gp*r2/(a1**2)
                 func=safe_exp(-r2)
                 ww(i-leftx,1,1)=func*(cval*cval2-sval*sval2)
                 ww(i-leftx,1,2)=func*(cval*sval2+sval*cval2)
              else
                 ww(i-leftx,1,1:2)=0.0_wp
              end if
           enddo
        else
           do i=leftx,rightx
              x=real(i-i0*16,gp)*h
              r=x-z0
              r2=r*r
              cval=cos(a2*r2)
              sval=sin(a2*r2)
              rk=real(i,gp)*h
              cval2=cos(kval*rk)
              sval2=sin(kval*rk)
              coff=r**n_gau
              r2=0.5_gp*r2/(a1**2)
              func=safe_exp(-r2)
              func=coff*func
              ww(i-leftx,1,1)=func*(cval*cval2-sval*sval2)
              ww(i-leftx,1,2)=func*(cval*sval2+sval*cval2)
           enddo
        end if
      END SUBROUTINE gauss_to_scf_4

      ! Original version
      !  subroutine gauss_to_scf
      !    n_left=lefts(0)
      !    n_right=rights(0)
      !    length=n_right-n_left+1
      !
      !    !print *,'nleft,nright',n_left,n_right
      !
      !    do k=1,4
      !       rights(k)=2*rights(k-1)+m
      !       lefts( k)=2*lefts( k-1)-m
      !    enddo
      !
      !    leftx = lefts(4)-n
      !    rightx=rights(4)+n  
      !
      !    !stop the code if the gaussian is too extended
      !    if (rightx-leftx > nwork) then
      !       !STOP 'gaustodaub'
      !       return
      !    end if
      !
      !    !loop for each complex component
      !    do icplx=1,ncplx
      !
      !       !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
      !
      !       !corrected for avoiding 0**0 problem
      !       if (ncplx==1) then
      !          if (n_gau == 0) then
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                ww(i-leftx,1,icplx)=func
      !             enddo
      !          else
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                coeff=r**n_gau
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                func=real(coeff,wp)*func
      !                ww(i-leftx,1,icplx)=func
      !             enddo
      !          end if
      !       else if (icplx == 1) then
      !          if (n_gau == 0) then
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                rk=real(i,gp)*h
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                cval=real(cos(kval*rk),wp)
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                ww(i-leftx,1,icplx)=func*cval
      !             enddo
      !          else
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                rk=real(i,gp)*h
      !                coeff=r**n_gau
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                cval=real(cos(kval*rk),wp)
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                func=real(coeff,wp)*func
      !                ww(i-leftx,1,icplx)=func*cval
      !             enddo
      !          end if
      !       else if (icplx == 2) then
      !          if (n_gau == 0) then
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                rk=real(i,gp)*h
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                sval=real(sin(kval*rk),wp)
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                ww(i-leftx,1,icplx)=func*sval
      !             enddo
      !          else
      !             do i=leftx,rightx
      !                x=real(i-i0*16,gp)*h
      !                r=x-z0
      !                rk=real(i,gp)*h
      !                coeff=r**n_gau
      !                r2=r/a
      !                r2=r2*r2
      !                r2=0.5_gp*r2
      !                sval=real(sin(kval*rk),wp)
      !                func=real(dexp(-real(r2,kind=8)),wp)
      !                func=real(coeff,wp)*func
      !                ww(i-leftx,1,icplx)=func*sval
      !             enddo
      !          end if
      !       end if
      !
      !       !print *,'here',gau_a,gau_cen,n_gau
      !       call apply_w(ww(0,1,icplx),ww(0,2,icplx),&
      !            leftx   ,rightx   ,lefts(4),rights(4),h)
      !
      !       call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
      !            lefts(4),rights(4),lefts(3),rights(3)) 
      !       call forward_c(ww(0,1,icplx),ww(0,2,icplx),&
      !            lefts(3),rights(3),lefts(2),rights(2)) 
      !       call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
      !            lefts(2),rights(2),lefts(1),rights(1)) 
      !
      !       call forward(  ww(0,1,icplx),ww(0,2,icplx),&
      !            lefts(1),rights(1),lefts(0),rights(0)) 
      !
      !    end do
      !
      !
      !  END SUBROUTINE gauss_to_scf


      !> One of the tails of the Gaussian is folded periodically
      !! We assume that the situation when we need to fold both tails
      !! will never arise
      subroutine fold_tail

        !modification of the calculation.
        !at this stage the values of c are fixed to zero
        !print *,'ncplx',ncplx,n_left,n_right,nwork,length
        do icplx=1,ncplx_w
           do i=n_left,n_right
              j=modulo(i,nmax+1)
              c(icplx,j,1)=c(icplx,j,1)+ww(i-n_left       ,2,icplx)
              c(icplx,j,2)=c(icplx,j,2)+ww(i-n_left+length,2,icplx)
           end do
        end do


      END SUBROUTINE fold_tail

    END SUBROUTINE gauss_to_daub_k


    subroutine gauss_c_to_daub_k(hgrid,kval,ncplx,gau_bf,ncs_s,factor , &
         gau_cen,gau_a, n_gau,&!no err, errsuc
         nmax,n_left,n_right,c,& 
         ww,nwork,periodic, hcutoff)      !added work arrays ww with dimension nwork
      use module_base
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: n_gau,nmax,nwork,ncs_s,ncplx
      real(gp), intent(in) :: hgrid,factor,gau_cen,gau_a,gau_bf
      real(wp), dimension(0:nwork,2,ncs_s, ncplx), intent(inout) :: ww 
      integer, intent(out) :: n_left,n_right
      real(wp), dimension(  ncs_s,ncplx,0:nmax,2), intent(out) :: c
      real(gp)  hcutoff

      !local variables
      real(gp), parameter :: pi=3.141592653589793_gp
      integer :: rightx,leftx,right_t,i0,i,k,length,j,ics, icplx
      real(gp) :: a,z0,h,x,r,coff,r2,fac
      real(wp) :: func,cval,sval
      integer, dimension(0:8) :: lefts,rights
      integer :: nrefinement, nforwards, ifwdtarget , ifwdsource, iswap
      real(gp) gau_kval, kval
      real(gp) cutoff, pishift

      !include the convolutions filters
      !include 'recs16.inc'! MAGIC FILTER  
      !include 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
      !include 'sym_16.inc'! WAVELET FILTERS

      !rescale the parameters so that hgrid goes to 1.d0  
      a=gau_a/hgrid
      gau_kval=gau_bf*hgrid*hgrid

      i0=nint(gau_cen/hgrid) ! the array is centered at i0

      z0=gau_cen/hgrid-real(i0,gp)
      cutoff= hcutoff /hgrid

      nrefinement=64
      nforwards=6

      h = (16 * .125_gp*.5_gp)/ nrefinement

      !calculate the array sizes;
      !at level 0, positions shifted by i0 
      right_t= ceiling(15.d0*a)

      !print *,'a,right_t',a,right_t,gau_a,hgrid

      !to rescale back the cofficients
      fac = hgrid**n_gau*sqrt(hgrid)*factor


      !initialise array
      c=0.0_gp

      if (periodic) then
         !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
         !we limit however the folding to one cell on each side 
         !!(commented out)
         !!     lefts( 0)=max(i0-right_t,-nmax)
         !!     rights(0)=min(i0+right_t,2*nmax)

         lefts( 0)=i0-right_t
         rights(0)=i0+right_t

         call gauss_c_to_scf()

         ! special for periodic case:
         call fold_tail
      else
         ! non-periodic: the Gaussian is bounded by the cell borders
         lefts( 0)=max(i0-right_t,   0)
         rights(0)=min(i0+right_t,nmax)

         call gauss_c_to_scf

         !loop for each complex component
         do icplx=1,ncplx
            do ics=1,ncs_s
               ! non-periodic: no tails to fold
               do i=0,length-1
                  c( ics,icplx,i+n_left,1)=fac*ww(i       ,2,ics, icplx)
                  c( ics,icplx,i+n_left,2)=fac*ww(i+length,2,ics, icplx) 
               end do
            end do
         end do
      endif


    contains

      subroutine gauss_c_to_scf
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

        do k=1,nforwards
           rights(k)=2*rights(k-1)+m
           lefts( k)=2*lefts( k-1)-m
        enddo

        leftx = lefts(nforwards)-n
        rightx=rights(nforwards)+n  

        !stop the code if the gaussian is too extended
        if (rightx-leftx > nwork) then
           STOP 'gaustodaub'
           return
        end if

        !loop for each complex component
        do icplx=1,ncplx
           pishift=-(icplx-1)*pi/2.0_gp
           do ics=1,ncs_s
              !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
              if( mod(nforwards,2)==0) then
                 ifwdtarget=1
                 ifwdsource=2
              else
                 ifwdtarget=2
                 ifwdsource=1
              endif


              if (ics == 1) then
                 if (n_gau == 0) then
                    do i=leftx,rightx
                       x=real(i-i0*nrefinement,gp)*h
                       sval=real(cos(kval*x+pishift))
                       r=x-z0
                       r2=r
                       r2=r2*r2
                       cval=real(cos(gau_kval*r2),wp)
                       r2=0.5_gp*r2/a/a
                       func=safe_exp(-r2)
                       if(abs(r)>cutoff) func=0
                       ww(i-leftx,ifwdtarget ,ics, icplx)=func*cval*sval
                    enddo
                 else
                    do i=leftx,rightx
                       x=real(i-i0*nrefinement,gp)*h
                       sval=real(cos(kval*x+pishift))
                       r=x-z0
                       coff=r**n_gau
                       r2=r
                       r2=r2*r2
                       cval=real(cos(gau_kval*r2),wp)
                       r2=0.5_gp*r2/a/a
                       func=safe_exp(-r2)
                       func=real(coff,wp)*func
                       if(abs(r)>cutoff) func=0
                       ww(i-leftx,ifwdtarget,ics, icplx)=func*cval*sval
                    enddo
                 end if
              else if (ics == 2) then
                 if (n_gau == 0) then
                    do i=leftx,rightx
                       x=real(i-i0*nrefinement,gp)*h
                       sval=real(cos(kval*x+pishift))
                       r=x-z0
                       r2=r
                       r2=r2*r2
                       cval=real(sin(gau_kval*r2),wp)
                       r2=0.5_gp*r2/a/a
                       func=safe_exp(-r2)
                       if(abs(r)>cutoff) func=0
                       ww(i-leftx,ifwdtarget,ics, icplx)=func*cval*sval
                    enddo
                 else
                    do i=leftx,rightx
                       x=real(i-i0*nrefinement,gp)*h
                       sval=real(cos(kval*x+pishift))
                       r=x-z0
                       coff=r**n_gau
                       r2=r
                       r2=r2*r2
                       cval=real(sin(gau_kval*r2),wp)
                       r2=0.5_gp*r2/a/a
                       func=safe_exp(-r2)
                       func=real(coff,wp)*func
                       if(abs(r)>cutoff) func=0
                       ww(i-leftx,ifwdtarget,ics, icplx)=func*cval*sval
                    enddo
                 end if
              end if

              !print *,'here',gau_a,gau_cen,n_gau

              iswap=ifwdsource
              ifwdsource=ifwdtarget
              ifwdtarget=iswap

!!$              call apply_w(ww(0,ifwdsource ,ics,icplx),ww(0,ifwdtarget ,ics,icplx),&
!!$                   leftx   ,rightx   ,lefts( nforwards),rights(nforwards  ),h)

              do i=nforwards,2,-1

                 iswap=ifwdsource
                 ifwdsource=ifwdtarget
                 ifwdtarget=iswap

!!$                 call forward_c(ww(0,ifwdsource ,ics, icplx),ww(0, ifwdtarget,ics, icplx),&
!!$                      lefts( i),rights( i),lefts(i-1),rights(i-1)) 

              enddo

              iswap=ifwdsource
              ifwdsource=ifwdtarget
              ifwdtarget=iswap

              if( ifwdsource .ne. 1) then
                 STOP ' ifwdsource .ne. 1  '
              endif

!!$              call forward(  ww(0,1,ics, icplx),ww(0,2,ics, icplx),&
!!$                   lefts(1),rights(1),lefts(0),rights(0)) 

           end do
        end do

      END SUBROUTINE gauss_c_to_scf

      subroutine fold_tail
        ! One of the tails of the Gaussian is folded periodically
        ! We assume that the situation when we need to fold both tails
        ! will never arise
        !implicit none


        !modification of the calculation.
        !at this stage the values of c are fixed to zero
        !print *,'ncs_s',ncs_s,n_left,n_right,nwork,length
        do icplx=1,ncplx
           do ics=1,ncs_s
              do i=n_left,n_right
                 j=modulo(i,nmax+1)
                 c(ics, icplx, j,1)=c(ics,icplx, j,1)+ww(i-n_left       ,2,ics, icplx)
                 c(ics, icplx, j,2)=c(ics,icplx,j,2)+ww(i-n_left+length,2,ics, icplx)
              end do
           end do
        end do


        do icplx=1,ncplx
           do ics=1,ncs_s
              do j=0,nmax
                 c(ics,icplx, j,1)=fac*c(ics, icplx , j, 1 )
                 c(ics,icplx, j,2)=fac*c(ics, icplx , j, 2 )
              enddo
           enddo
        end do

      END SUBROUTINE fold_tail

    END SUBROUTINE gauss_c_to_daub_k

    !> APPLYING THE MAGIC FILTER ("SHRINK") 
    subroutine apply_w(nc,cx,c,leftx,rightx,left,right,h)
      use module_base
      implicit none
      integer, intent(in) :: leftx,rightx,left,right,nc
      real(gp), intent(in) :: h
      real(wp), dimension(nc,leftx:rightx), intent(in) :: cx
      real(wp), dimension(nc,left:right), intent(out) :: c
      !local variables
      integer, parameter :: unroll=4
      integer :: i,j,ig,ngl,k,ngr
      real(wp) :: sqh
      real(wp), dimension(unroll) :: ci

      sqh=real(sqrt(h),wp)
      
      !limit for unrolling the convolution
      ngl=nc/unroll
      !remnant
      ngr=nc-ngl*unroll
      
      !!  !$omp parallel do default(shared) private(i,ci,j)
      do i=left,right
         do ig=0,ngl-1
            ci=0.0_wp
            do j=-n,n
               do k=1,unroll
                  ci(k)=ci(k)+cx(k+ig*unroll,i+j)*w(j)
               end do
            end do
            do k=1,unroll
               c(k+ig*unroll,i)=ci(k)*sqh
            end do
         end do
         do ig=1,ngr
            ci(ig)=0.0_wp
            do j=-n,n
               ci(ig)=ci(ig)+cx(ngl*unroll+ig,i+j)*w(j)
            end do
            c(ngl*unroll+ig,i)=ci(ig)*sqh
         end do
      end do
      !!  !$omp end parallel do

    END SUBROUTINE apply_w

    !> APPLYING THE INVERSE MAGIC FILTER ("GROW") 
    !!subroutine apply_inverse_w(cx,c,leftx,rightx,left,right,h)
    !!  use module_base
    !!  implicit none
    !!  integer, intent(in) :: leftx,rightx,left,right
    !!  real(gp), intent(in) :: h
    !!  real(wp), dimension(leftx:rightx), intent(in) :: cx
    !!  real(wp), dimension(left:right), intent(out) :: c
    !!  !local variables
    !!  include 'recs16.inc'
    !!  integer :: i,j
    !!  real(wp) :: sqh,ci
    !!
    !!  sqh=real(sqrt(h),wp)
    !!
    !!  do i=left,right
    !!     ci=0.0_wp
    !!     do j=-n,n
    !!        ci=ci+cx(i+j)*w(-j) !transposed MF         
    !!     enddo
    !!     c(i)=ci*sqh
    !!  enddo
    !!  
    !!END SUBROUTINE apply_inverse_w


    !> FORWARD WAVELET TRANSFORM WITHOUT WAVELETS ("SHRINK")
    subroutine forward_c(nc,c,c_1,left,right,left_1,right_1)
      implicit none
      integer, intent(in) :: left,right,left_1,right_1,nc
      real(wp), dimension(nc,left:right), intent(in) :: c
      real(wp), dimension(nc,left_1:right_1), intent(out) :: c_1
      !local variables
      integer, parameter :: unroll=4
      integer :: i,i2,j,ngl,ngr,k,ig
      real(wp), dimension(unroll) :: ci

      !limit for unrolling the convolution
      ngl=nc/unroll
      !remnant
      ngr=nc-ngl*unroll

      ! get the coarse scfunctions and wavelets
      !!  !$omp parallel do default(shared) private(i,i2,j,ci)
      do i=left_1,right_1
         i2=2*i
         do ig=0,ngl-1
            ci=0.0_wp
            do j=-m,m
               do k=1,unroll
                  ci(k)=ci(k)+cht(j)*c(k+ig*unroll,j+i2)
               end do
            enddo
            do k=1,unroll
               c_1(k+ig*unroll,i)=ci(k)
            end do
         end do
         do ig=1,ngr
            ci(ig)=0.0_wp
            do j=-m,m
               ci(ig)=ci(ig)+cht(j)*c(ig+ngl*unroll,j+i2)
            end do
            c_1(ig+ngl*unroll,i)=ci(ig)
         end do
      enddo
      !!  !$end parallel do

    END SUBROUTINE forward_c


    !>  CONVENTIONAL FORWARD WAVELET TRANSFORM ("SHRINK")
    subroutine forward(nc,c,cd_1,left,right,left_1,right_1)
      use module_base
      implicit none
      integer, intent(in) :: left,right,left_1,right_1,nc
      real(wp), dimension(nc,left:right), intent(in) :: c
      real(wp), dimension(nc,left_1:right_1,2), intent(out) :: cd_1
      !local variables
      integer, parameter :: unroll=4
      integer :: i,i2,j,ngr,ngl,ig,k
      real(wp), dimension(unroll) :: ci,di
      !include 'sym_16.inc'

      !limit for unrolling the convolution
      ngl=nc/unroll
      !remnant
      ngr=nc-ngl*unroll

      ! get the coarse scfunctions and wavelets
      do i=left_1,right_1
         i2=2*i
         do ig=0,ngl-1
            ci=0.d0
            di=0.d0
            do j=-m,m
               do k=1,unroll
                  ci(k)=ci(k)+cht(j)*c(k+ig*unroll,j+i2)
                  di(k)=di(k)+cgt(j)*c(k+ig*unroll,j+i2)
               end do
            end do
            do k=1,unroll
               cd_1(k+ig*unroll,i,1)=ci(k)
               cd_1(k+ig*unroll,i,2)=di(k)
            end do
         end do
         do ig=1,ngr
            ci(ig)=0.d0
            di(ig)=0.d0
            do j=-m,m
               ci(ig)=ci(ig)+cht(j)*c(ig+ngl*unroll,j+i2)
               di(ig)=di(ig)+cgt(j)*c(ig+ngl*unroll,j+i2)
            end do
            cd_1(ig+ngl*unroll,i,1)=ci(ig)
            cd_1(ig+ngl*unroll,i,2)=di(ig)
         end do
      enddo

    END SUBROUTINE forward


    subroutine determine_bounds(nres,periodic,nmax,right_t,i0,&
         lefts,rights)
      implicit none
      logical, intent(in) :: periodic !<boundary conditions
      integer, intent(in) :: nres !< number of level of resolution
      integer, intent(in) :: nmax !< size of the output array
      integer, intent(in) :: i0 !< center in the grid points
      integer, intent(in) :: right_t !< amplitude of the values
      integer, dimension(0:nres+1), intent(out) :: lefts !<left bounds
      integer, dimension(0:nres+1), intent(out) :: rights !<right bounds

      !local variables
      !>multiplicator after which the values of the gaussian become negligible
      !real(gp), parameter :: maxmult=15.0_gp 
      !!real(gp) :: right_t
      integer :: k
      !calculate the array sizes;
      !at level 0, positions shifted by i0 
      !right_t= ceiling(maxmult*a)

      !we expand the whole Gaussian in scfunctions and later fold one of its tails periodically
      !we limit however the folding to one cell on each side (it can be eliminated)
      if (periodic) then
         lefts( 0)=i0-right_t
         rights(0)=i0+right_t
      else
         ! non-periodic: the Gaussian is bounded by the cell borders
         lefts( 0)=max(i0-right_t,   0)
         rights(0)=min(i0+right_t,nmax)
      end if

      do k=1,nres
         rights(k)=2*rights(k-1)+m
         lefts( k)=2*lefts( k-1)-m
      enddo

      !last one, including magic filters
      lefts(nres+1) = lefts(nres)-n
      rights(nres+1)= rights(nres)+n  

    end subroutine determine_bounds

    !>inverse wavelet transforms to bring back the gaussian
    !!coefficients at the lowest resolution levels
    !! this set if idwt is combined with magic filters to 
    !! express correctly the wavefunction coefficients at highest resolution
    subroutine magic_idwts(nres,ng,lefts,rights,nwork,ww,inw)
      implicit none
      integer, intent(in) :: nres !< number of resolution levels
      !> size of the work array. Should be at least equal to 
      !! rights(nres+1)-lefts(nres+1). Otherwise the 
      !! transformation is simply skipped
      integer, intent(in) :: nwork,ng 
      integer, dimension(0:nres+1), intent(in) :: lefts
      integer, dimension(0:nres+1), intent(in) :: rights
      !> identify the index of the second dimension of
      !! the ww array where the final data have to be retrieved
      integer, intent(out) :: inw
      real(gp), dimension(ng,nwork,2), intent(inout) :: ww
      !local variables
      integer :: resolution,ires,inwo
      real(gp) :: h

      resolution=2**nres
      !the grid spacing is now given in terms of the level
      h=1.0_gp/(real(resolution,gp))

      !do not do anything if the gaussian is too extended
      !if (rightx-leftx > nwork) then
      !   !STOP 'gaustodaub'
      !   return
      !end if

      call apply_w(ng,ww(1,1,1),ww(1,1,2),&
           lefts(nres+1),rights(nres+1),lefts(nres),rights(nres),h)

      inw=2
      do ires=nres,2,-1
         inwo=3-inw
         call forward_c(ng,ww(1,1,inw),ww(1,1,inwo),&
              lefts(ires),rights(ires),lefts(ires-1),rights(ires-1)) 
         inw=3-inw
      end do
      call forward(ng,ww(1,1,inw),ww(1,1,3-inw),&
           lefts(1),rights(1),lefts(0),rights(0)) 

    end subroutine magic_idwts

    !> One of the tails of the Gaussian is folded periodically
    !! We assume that the situation when we need to fold both tails
    !! will never arise
    subroutine retrieve_results(periodic,ng,n_left,n_right,nmax,ww,theor_norm,fac,c,error)
      use module_base, only: f_memcpy,dot
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: n_left,n_right,nmax,ng
      !> the expected normalization of the gaussian. When different from zero,
      !! control the result and write the error in the error variable
      real(wp), intent(in) :: theor_norm 
      real(wp), intent(in) :: fac
      real(wp), dimension(n_left:n_right,2), intent(in) :: ww !< the values of the gaussian to be folded
      real(wp), intent(out) :: error !< the error in the expression of the function
      real(wp), dimension(0:nmax,2), intent(out) :: c !< final results
      !local variables
      integer :: i,j
      real(wp) :: cn2,tt

      if (periodic) then
         c=0.0_wp
         do i=n_left,n_right
            j=modulo(i,nmax+1)
            c(j,1)=c(j,1)+ww(i,1)
            c(j,2)=c(j,2)+ww(i,2)
         end do
      else
         !it might be that we have to transpose the final result
         call f_memcpy(src=ww,dest=c)
      end if

      if (theor_norm /= 0.0_wp) then
         !calculate the (relative) error
         cn2=dot(2*(n_right-n_left+1),ww(1,1),1,ww(1,1),1)
         error=sqrt(abs(1.0_gp-real(cn2,gp)/theor_norm))
      else
         error=0.0_gp
      end if

      !write(*,*)'error, non scaled:',error
      !
      !RESCALE BACK THE COEFFICIENTS AND THE ERROR
      c=real(fac,wp)*c
      error=error*fac

    end subroutine retrieve_results

    !> express a function defined as x**n_gau*exp(-(x-x0)**2/(2*a**2))
    pure subroutine gaus_highres(level,ng,n_gau,a,x0,leftx,rightx,ww)
      implicit none
      integer, intent(in) :: level !< resolution level of the evaluation
      integer, intent(in) :: ng !< number of gaussians
      integer, dimension(ng), intent(in) :: n_gau !< parameters of the gaussian
      integer, intent(in) :: leftx !<leftmost point of the grid
      integer, intent(in) :: rightx !<rightmost point of the grid
      real(gp), dimension(ng), intent(in) :: a !<standard deviations
      real(gp), intent(in) :: x0 !<gaussian center
      !>collcation values of the gaussian
      real(gp), dimension(ng,leftx:rightx), intent(out) :: ww 
      !local variables
      integer :: i,i0,resolution,ig
      real(gp) :: x,h,z0,r

      resolution=2**level
      !the grid spacing is now given in terms of the level
      h=1.0_gp/(real(resolution,gp))
      !the grid point closest to the center
      i0=nint(x0) ! the array is centered at i0
      z0=x0-real(i0,gp)
  
      !this part has to be extended to the different cases
      do i=leftx,rightx
         x=real(i-i0*resolution,gp)*h
         r=x-z0
         do ig=1,ng
            ww(ig,i)=collocate_gaussian(n_gau(ig),a(ig),r)
         end do
      end do
    end subroutine gaus_highres

    pure function collocate_gaussian(n_gau,a,r) result(func)
      use module_defs, only: safe_exp
      implicit none
      integer, intent(in) :: n_gau !<principal quantum number
      real(gp), intent(in) :: a !<standard deviation
      real(gp), intent(in) :: r !< function argument
      real(gp) :: func
      !local variables
      real(gp) :: r2

      r2=r/a
      r2=r2*r2
      r2=0.5_gp*r2
      func=safe_exp(-r2)
      func=coeff(n_gau,r)*func

    end function collocate_gaussian
    
    pure function coeff(n_gau,r)
      implicit none
      integer, intent(in) :: n_gau !<principal quantum number
      real(gp), intent(in) :: r !< function argument
      real(gp) :: coeff

      if (n_gau==0) then
         coeff=1.0_gp
      else
         coeff=r**n_gau
      end if
    end function coeff
end module gaussdaub
