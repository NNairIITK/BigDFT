!> @file
!!  Gaussian to Daubechies wavelets projection routines
!! @author
!!    Copyright (C) 2007-2014 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module to project gaussian functions to Daubechies wavelets
module gaussdaub
  use module_defs, only: wp,gp
  implicit none
  private
  
  !!@todo These parameters has to be modified to avoid naming confusions
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
  !> Coefficients for wavelet transform (orthonormal wavelet)
  real(wp), dimension(-mm:mm), parameter :: cht=ch

  !> g coefficients from h coefficients
  real(wp), dimension(-mm:mm), parameter :: cg=[0.0_wp,((-1)**(iw+1)*cht(-iw),iw=-mm,mm-1)]
  real(wp), dimension(-mm:mm), parameter :: cgt=[0.0_wp,((-1)**(iw+1)*ch(-iw),iw=-mm,mm-1)]

  !> Magic filter coefficients for daubechies family
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

  public :: gau_daub_1d,workspace_query_gau_daub_1d

  contains

    !> convert a shell in contracted Gaussian basis in wavelet
    !! given in a localization region
    subroutine gaussian_shell_to_daub(Gbasis,iter,factor,idir)
      use gaussians
      use psp_projectors_base
      implicit none
      integer, intent(in) :: idir !<derivative direction
      type(gaussian_basis_new), intent(in) :: Gbasis !<reference gaussian basis 
      real(gp), dimension(Gbasis%ncplx), intent(in) :: factor !<multipliciative factor
      !>iterator on the gaussians of the shell
      type(gaussian_basis_iter), intent(inout) :: iter
      !local variables
      integer, parameter :: nterm_max=20 !if GTH nterm_max=4
      !for the moment the number of gaussians per shell is 4 maximum
      integer, parameter :: ngau_max=4
      integer :: l,m,n,iex,iey,iez,idir2,nterm,iterm
      real(gp) :: fg_tmp
      real(gp), dimension(Gbasis%ncplx) :: coeff, expo,gau_ctmp
      !here the dimension of the n_gau array has to be inquired
      integer, dimension(2*iter%l-1) :: ng
      integer, dimension(ngau_max*nterm_max,3,2*iter%l-1) :: n_gau
      real(gp), dimension(Gbasis%ncplx,ngau_max*nterm_max,2*iter%l-1) :: gau_a
      real(gp), dimension(Gbasis%ncplx,ngau_max*nterm_max,2*iter%l-1) :: factors
      integer, dimension(3) :: nterm_arr
      integer, dimension(3,nterm_max,3) :: lxyz_arr
      real(gp), dimension(nterm_max,3) :: fac_arr

      
      !identify the quantum numbers
      l=iter%l
      n=iter%n

      ng=0
      iex=0
      iey=0
      iez=0
      !seq : 11 22 33 12 23 13
      if (idir == 4 .or. idir == 9) iex=1
      if (idir == 5 .or. idir == 7) iey=1
      if (idir == 6 .or. idir == 8) iez=1

      do while (gaussian_iter_next_gaussian(Gbasis, iter, coeff, expo))

         do m=1,2*l-1

            !set the exponents
            gau_ctmp(1:Gbasis%ncplx)=gau_c(Gbasis%ncplx,expo)
            fg_tmp=fgamma(Gbasis%ncplx,gau_ctmp(1),l,n)

            if (idir==0) then !normal projector calculation case
               idir2=1
               call calc_coeff_proj(l,n,m,nterm_max,nterm,&
                    n_gau(ng(m)+1,1,m),n_gau(ng(m)+1,2,m),&
                    n_gau(ng(m)+1,3,m),&
                    fac_arr)
               do iterm=1,nterm
                  factors(:,ng(m)+iterm,m)=factor(:)*fg_tmp*fac_arr(iterm,idir2)
               end do
            else !calculation of projector derivative
               idir2=mod(idir-1,3)+1
               call calc_coeff_derproj(l,n,m,nterm_max,gau_ctmp(1),&
                    nterm_arr,lxyz_arr,fac_arr)
               nterm=nterm_arr(idir2)
               do iterm=1,nterm
                  factors(:,ng(m)+iterm,m)=&
                       factor(:)*fg_tmp*fac_arr(iterm,idir2)
                  n_gau(ng(m)+iterm,1,m)=lxyz_arr(1,iterm,idir2)+iex
                  n_gau(ng(m)+iterm,2,m)=lxyz_arr(2,iterm,idir2)+iey
                  n_gau(ng(m)+iterm,3,m)=lxyz_arr(3,iterm,idir2)+iez
               end do
            end if
            do iterm=1,nterm
               gau_a(:,ng(m)+iterm,m)=gau_ctmp(:)
            end do
            ng(m)=ng(m)+nterm
         end do
      end do

!!$      !now fill the components for each of the shell numbers
!!$      do m=1,2*l-1
!!$         !here we can branch from the creation with tensor product decomposition to the one with actual 3d collocation
!!$         call gauss_to_daub_3d(periodic,ng(m),gau_cen,n_gau(1,1,m),kval,&
!!$              ncplx_a,gau_a(1,1,m),ncplx_f,factors(1,1,m),&
!!$              hgrid,nres,ncplx,nwork,ww,c,lr,phi(istart_c))
!!$         istart_c=istart_c+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx
!!$      end do

    end subroutine gaussian_shell_to_daub

    pure function gau_c(ncplx,gau_a)
      implicit none
      integer, intent(in) :: ncplx
      real(gp), dimension(ncplx), intent(in) :: gau_a
      real(gp), dimension(ncplx) :: gau_c
      !local variables
      real(gp) :: fpi,fgamma
      !this value can also be inserted as a parameter
      if (ncplx == 1) then
         gau_c(1) = 1._gp / sqrt(2._gp * gau_a(1))
      else
         gau_c(1) = 1._gp / sqrt(2._gp * gau_a(1))
         gau_c(2) = gau_a(2)
      end if
    end function gau_c

    pure function fgamma(ncplx,rloc,l,i)
      implicit none
      integer, intent(in) :: ncplx,l,i
      real(gp), intent(in) :: rloc
      real(gp) :: fgamma
      !local variables
      real(gp) :: fpi

      select case(ncplx)
      case(1)
         !fpi=pi^-1/4 pi^-1/2, pi^-1/4 comes from sqrt(gamma(x)) and pi^-1/2 from Ylm.
         !fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)
         fpi=0.42377720812375763_gp
         fgamma=sqrt(2.0_gp)*fpi/(sqrt(rloc)**(2*(l-1)+4*i-1))
      case(2)
         fpi=0.56418958354775628_gp
         select case(l)
         case(1)
            fgamma= 0.70710678118654757_gp !1.0/sqrt(2.0)
         case(2)
            fgamma= 0.8660254037844386_gp !sqrt(3)/2.0
         case(3)
            fgamma= 1.3693063937629153_gp  !sqrt(3*5)/(2.0*sqrt(2))
         case(4)
            fgamma= 2.5617376914898995_gp  !sqrt(7*5*3)/(4.0) 
         end select
         fgamma = fgamma * fpi
      end select
    end function fgamma

    !>convert a gaussian basis in daubechies basis set
    !! given by a localization region descriptor
    !! this routine should be prone to generalizations to non-orthorhombic
    !! cells
    subroutine gauss_to_daub_3d(periodic,ng,gau_cen,n_gau,kval,ncplx_a,gau_a,ncplx_f,factor,hgrid,nres,ncplx,nwork,ww,c,lr,phi)
      use yaml_strings, only: yaml_toa
      use dictionaries, only: f_err_throw
      use f_utils, only: f_zero
      use locregs, only: locreg_descriptors
      implicit none
      logical, dimension(3), intent(in) :: periodic !< determine the bc
      integer, intent(in) :: nwork,nres
      !> 1 for real gaussians, 2 for complex ones 
      !! (to be generalized to k-points)
      integer, intent(in) :: ncplx,ncplx_a,ncplx_f
      integer, intent(in) :: ng !<number of different gaussians to convert
      !>principal quantum numbers of any of the gaussians
      integer, dimension(ng,3), intent(in) :: n_gau 
      !> localization region descriptors of the resulting array
      type(locreg_descriptors), intent(in) :: lr
      !>multiplicative factors which have to be added to the different
      !!terms
      real(wp), dimension(ncplx_f,ng), intent(in) :: factor
      !>standard deviations of the different gaussians (might be complex)
      real(wp), dimension(ncplx_a,ng), intent(in) :: gau_a
      real(gp), dimension(3), intent(in) :: hgrid !< grid spacing in wavelets
      real(gp), dimension(3), intent(in) :: gau_cen !< center of the gaussian
      real(gp), dimension(3), intent(in) :: kval !< value of k-point
      !> work array, should be created by a workspace query
      real(wp), dimension(nwork,2), intent(inout) :: ww
      !> second work array for building the Gaussian
      real(wp), dimension(ncplx,ng,2,0:max(lr%d%n1,lr%d%n2,lr%d%n3),3) :: c
      !> buffer to put the wavelet components
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(inout) :: phi
      !local variables
      integer :: iseg,idir,n1p1,np,jj,j0,ii,i3,i2,i0,i1,i0jj,i,ig
      integer :: j1,ind_c,ind_f,mvctr_c
      real(wp) :: tt_c,tt_f1,tt_f2,tt_f3,tt_f4,tt_f5,tt_f6,tt_f7
      integer, dimension(3) :: nmax
      integer, dimension(3) :: nstart

      nmax(1)=lr%d%n1
      nmax(2)=lr%d%n2
      nmax(3)=lr%d%n3

      nstart(1)=lr%ns1
      nstart(2)=lr%ns2
      nstart(3)=lr%ns3

      do idir=1,3
         call gau_daub_1d(periodic(idir),ng,gau_cen(idir),n_gau(1,idir),&
              kval(idir),ncplx_a,gau_a,ncplx_f,factor,hgrid(idir),nres,nstart(idir),&
              nmax(idir),ncplx,c(1,1,1,0,idir),nwork,ww)
      end do

      n1p1=lr%d%n1+1
      np=n1p1*(lr%d%n2+1)
      mvctr_c=lr%wfd%nvctr_c
      !now the arrays for the gaussians are ready to be compressed in wavelet form
      do iseg=1,lr%wfd%nseg_c
         jj=lr%wfd%keyvglob(iseg)
         j0=lr%wfd%keyglob(1,iseg)
         j1=lr%wfd%keyglob(2,iseg)
         ii=j0-1
         i3=ii/(np)
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         i0jj=jj-i0
         do i=i0,i1
            tt_c=0.0_wp
            !real case only, the complex case will follow
            do ig=1,ng
               tt_c=tt_c+c(1,ig,1,i,1)*c(1,ig,1,i2,2)*c(1,ig,1,i3,3)
            end do
            ind_c=i+i0jj
            phi(ind_c)=tt_c
         end do
      end do
      
      ! Other terms: fine projector components
      !$omp do
      do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
         jj=lr%wfd%keyvglob(iseg)
         j0=lr%wfd%keyglob(1,iseg)
         j1=lr%wfd%keyglob(2,iseg)
         ii=j0-1
         i3=ii/(np)
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         i0jj=7*(jj-i0-1)+mvctr_c
         do i=i0,i1
            tt_f1=0.0_wp
            tt_f2=0.0_wp
            tt_f3=0.0_wp
            tt_f4=0.0_wp
            tt_f5=0.0_wp
            tt_f6=0.0_wp
            tt_f7=0.0_wp
            do ig=1,ng
               tt_f1=tt_f1+c(1,ig,2,i,1)*c(1,ig,1,i2,2)*c(1,ig,1,i3,3)
               tt_f2=tt_f2+c(1,ig,1,i,1)*c(1,ig,2,i2,2)*c(1,ig,1,i3,3)
               tt_f3=tt_f3+c(1,ig,2,i,1)*c(1,ig,2,i2,2)*c(1,ig,1,i3,3)
               tt_f4=tt_f4+c(1,ig,1,i,1)*c(1,ig,1,i2,2)*c(1,ig,2,i3,3)
               tt_f5=tt_f5+c(1,ig,2,i,1)*c(1,ig,1,i2,2)*c(1,ig,2,i3,3)
               tt_f6=tt_f6+c(1,ig,1,i,1)*c(1,ig,2,i2,2)*c(1,ig,2,i3,3)
               tt_f7=tt_f7+c(1,ig,2,i,1)*c(1,ig,2,i2,2)*c(1,ig,2,i3,3)
            end do
            ind_f=7*i+i0jj
            phi(ind_f+1)=tt_f1
            phi(ind_f+2)=tt_f2
            phi(ind_f+3)=tt_f3
            phi(ind_f+4)=tt_f4
            phi(ind_f+5)=tt_f5
            phi(ind_f+6)=tt_f6
            phi(ind_f+7)=tt_f7
         end do
      end do
      
    end subroutine gauss_to_daub_3d

    subroutine workspace_query_gau_daub_1d(&
         periodic,ng,gau_cen,kval,ncplx_a,gau_a,hgrid,nres,nstart,nmax,&
         nwork)
      use yaml_strings, only: yaml_toa
      use dictionaries, only: f_err_throw
      use f_utils, only: f_zero
      implicit none
      logical, intent(in) :: periodic !< determine the bc
      integer, intent(in) :: nmax,nres,nstart
      !> 1 for real gaussians, 2 for complex ones 
      !! (to be generalized to k-points)
      integer, intent(in) :: ncplx_a
      integer, intent(in) :: ng !<number of different gaussians to convert
      real(gp), intent(in) :: gau_cen !< center of the gaussian
      !>standard deviations of the different gaussians (might be complex)
      real(wp), dimension(ncplx_a,ng), intent(in) :: gau_a
      real(gp), intent(in) :: hgrid !< grid spacing in wavelets
      real(gp), intent(in) :: kval !< value of k-point
      integer, intent(out) :: nwork
      !local variables
      integer :: i0,ncplx_w
      real(wp) :: x0
      integer, dimension(0:nres+1) :: lefts,rights !< use automatic arrays here, we should use parameters in the module
      real(gp), dimension(ncplx_a) :: aval
      real(gp), dimension(ncplx_a,ng) :: a

      !first, determine the cutoffs where the 
      !bounds are to be calculated
      x0=gau_cen/hgrid
      i0=nint(x0) ! the array is centered at i0
      !here the are the quantities for any of the objects
      a=gau_a/hgrid
      aval=maxval(a,dim=2)

      !these are the limits of each of the convolutions
      !! here maxval does not have sense for a complex numbers
      call determine_bounds(nres,periodic,nstart,nmax,right_t(aval(1)),&
           i0,lefts,rights)

      if ( kval /= 0.0_gp .or. ncplx_a==2) then
         ncplx_w=2
      else
         ncplx_w=1
      end if

      !the total dimension of the work array should be
      nwork=ng*(rights(nres+1) -lefts(nres+1))*ncplx_w

    end subroutine workspace_query_gau_daub_1d

    !> provides the extreme from which the Gaussian function is
    !! assumed to be zero
    pure function right_t(aval)
      real(gp), intent(in) :: aval
      integer :: right_t

      right_t=max(ceiling(15.d0*aval),m+2)
    end function right_t
    
    !> Convert a gaussian to one-dimensional functions
    !! Gives the expansion coefficients of :
    !!   factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
    !! Multiply it for the k-point factor exp(Ikx)
    !! For this reason, a real (cos(kx)) and an imaginary (sin(kx)) part are provided 
    !! also, the value of factor and  gau_a might be complex.
    !! therefore the result has to be complex accordingly
    subroutine gau_daub_1d(periodic,ng,gau_cen,n_gau,kval,ncplx_a,gau_a,ncplx_f,factor,hgrid,nres,nstart,&
         nmax,ncplx,c,nwork,ww)
      use yaml_strings, only: yaml_toa
      use dictionaries, only: f_err_throw
      use f_utils, only: f_zero
      implicit none
      logical, intent(in) :: periodic !< determine the bc
      integer, intent(in) :: nmax,nwork,nres,nstart
      !> 1 for real gaussians, 2 for complex ones 
      !! (to be generalized to k-points)
      integer, intent(in) :: ncplx,ncplx_a,ncplx_f
      integer, intent(in) :: ng !<number of different gaussians to convert
      !>principal quantum numbers of any of the gaussians
      integer, dimension(ng), intent(in) :: n_gau 
      !>multiplicative factors which have to be added to the different
      !!terms
      real(wp), dimension(ncplx_f,ng), intent(in) :: factor
      !>standard deviations of the different gaussians (might be complex)
      real(wp), dimension(ncplx_a,ng), intent(in) :: gau_a
      real(gp), intent(in) :: hgrid !< grid spacing in wavelets
      real(gp), intent(in) :: gau_cen !< center of the gaussian
      real(gp), intent(in) :: kval !< value of k-point
      real(wp), dimension(ncplx,ng,2,0:nmax), intent(inout) :: c
      real(wp), dimension(nwork,2), intent(inout) :: ww
      !local variables
      integer :: i0,inw,ig,ncplx_w,nwork_tot,icplx,nst
      real(wp) :: x0
      integer, dimension(0:nres+1) :: lefts,rights !< use automatic arrays here, we should use parameters in the module
      real(gp), dimension(ncplx_a) :: aval
      real(gp), dimension(ncplx_a,ng) :: a
      real(gp), dimension(ncplx_f,ng) :: fac
      real(gp), dimension(ng) :: theor_norm2,error

      !should raise an error if the value ncplx is not consistent

      !first, determine the cutoffs where the 
      !bounds are to be calculated
      x0=gau_cen/hgrid
      i0=nint(x0) ! the array is centered at i0
      !here the are the quantities for any of the objects
      a=gau_a/hgrid
      aval=maxval(a,dim=2)

      !right_t=max(ceiling(15.d0*maxval(a)),m+2)
      !the multiplicative factors for any of the object
      do ig=1,ng
         do icplx=1,ncplx_f
            fac(icplx,ig)=hgrid**n_gau(ig)*sqrt(hgrid)*factor(icplx,ig)
         end do
         !and the expected theoretical norm of the gaussian
         theor_norm2(ig)=valints(n_gau(ig))*a(1,ig)**(2*n_gau(ig)+1)
      end do
      
      !these are the limits of each of the convolutions
      !! here maxval does not have sense for a complex numbers
      call determine_bounds(nres,periodic,nstart,nmax,right_t(aval(1)),&
           i0,&
         lefts,rights)

      if ( kval /= 0.0_gp .or. ncplx_a==2) then
         ncplx_w=2
      else
         ncplx_w=1
      end if

      !the total dimension of the work array should be
      nwork_tot=ng*(rights(nres+1) -lefts(nres+1))*ncplx_w
      
      !then check if the work array has a size which is big enough
      if (nwork_tot > nwork ) then
         call f_err_throw('The size of the work array ('//trim(yaml_toa(nwork))//&
              ') is too little to convert the gaussian in wavelets, put at least the value '//&
              trim(yaml_toa(nwork_tot)),&
              err_name='BIGDFT_RUNTIME_ERROR')
         return
      end if
      
      !fill the array of the high resolution values
      call gaus_highres(nres,ng,n_gau,ncplx_a,a,x0,kval,ncplx_w,&
           lefts(nres+1),rights(nres+1),ww)

      !then come back to the original resolution level
      call magic_idwts(nres,ncplx_w*ng,lefts,rights,nwork,ww,inw)

      nst=nstart
      if (periodic) nst=0
      !and retrieve the results and the error if desired
      call retrieve_results(periodic,ng,lefts(0),rights(0),nst,nmax,ncplx_w,&
           ww(1,inw),theor_norm2,ncplx_f,fac,ncplx,c,error)

    end subroutine gau_daub_1d

    !> Project gaussian functions in a mesh of Daubechies scaling functions
    !! Gives the expansion coefficients of :
    !!    factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
    subroutine gauss_to_daub(hgrid,factor,gau_cen,gau_a,n_gau,&!no err, errsuc
         nmax,n_left,n_right,c,err_norm,&                      !no err_wav. nmax instead of n_intvx
         ww,nwork,periodic)                         !added work arrays ww with dimension nwork
      use module_base
      implicit none
      logical, intent(in) :: periodic !< the flag for periodic boundary conditions
      integer, intent(in) :: n_gau    !< x**n_gau (polynomial degree)
      integer, intent(in) :: nmax     !< size of the grid
      integer, intent(in) :: nwork    !< size of the work array (ww) >= (nmax+1)*17
      real(gp), intent(in) :: hgrid   !< step size
      real(gp), intent(in) :: factor  !< normalisation factor
      real(gp), intent(in) :: gau_cen !< center of gaussian function
      real(gp), intent(in) :: gau_a   !< parameter of gaussian      
      real(wp), dimension(0:nwork,2), intent(inout) :: ww !< work arrays that have to be 17 times larger than C
      integer, intent(out) :: n_left,n_right !< interval where the gaussian is larger than the machine precision
      real(gp), intent(out) :: err_norm      !< normalisation error
      real(wp), dimension(0:nmax,2), intent(out) :: c !< c(:,1) array of scaling function coefficients
                                                      !! c(:,2) array of wavelet coefficients:
      !local variables
      integer :: right_t,i0,i,length
      real(gp) :: a,z0,h,theor_norm2,error,fac
      real(dp) :: cn2,tt
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


    !> Project gaussian functions in a mesh of Daubechies scaling functions
    !! Gives the expansion coefficients of :
    !!   factor*x**n_gau*exp(-(1/2)*(x/gau_a)**2)
    !! Multiply it for the k-point factor exp(Ikx)
    !! For this reason, a real (cos(kx)) and an imaginary (sin(kx)) part are provided 
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
      integer, intent(in) :: ncplx_w !< size of the ww matrix
      integer, intent(in) :: ncplx_g !< 1 or 2 for simple or complex gaussians, respectively.
      integer, intent(in) :: ncplx_k !< use 2 for k-points.
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

!!$        if (ncplx_w==1) then
!!$           !no kpts and real gaussians
!!$           call gauss_to_scf_1()
!!$        elseif(ncplx_k==2 .and. ncplx_g==1) then
!!$           !kpts and real gaussians
!!$           call gauss_to_scf_2()
!!$        elseif(ncplx_k==1 .and. ncplx_g==2) then
!!$           !no kpts and complex gaussians
!!$           call gauss_to_scf_3()
!!$        elseif(ncplx_k==2 .and. ncplx_g==2) then
!!$           !kpts and complex gaussians
!!$           call gauss_to_scf_4()
!!$        endif

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

!!!#!
!!!#!!!!!IMPLEMENTED
!!!#!      ! Called when ncplx_w = 1
!!!#!      subroutine gauss_to_scf_1
!!!#!
!!!#!        !loop for each complex component
!!!#!        !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
!!!#!        !corrected for avoiding 0**0 problem
!!!#!        icplx = 1
!!!#!        if (n_gau == 0) then
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              r2=r/a1
!!!#!              r2=r2*r2
!!!#!              r2=0.5_gp*r2
!!!#!              func=safe_exp(-r2)
!!!#!              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
!!!#!              ww(i-leftx,1,icplx)=func
!!!#!           enddo
!!!#!        else
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              coff=r**n_gau
!!!#!              r2=r/a1
!!!#!              r2=r2*r2
!!!#!              r2=0.5_gp*r2
!!!#!              func=safe_exp(-r2)
!!!#!              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
!!!#!              func=coff*func
!!!#!              ww(i-leftx,1,icplx)=func
!!!#!           enddo
!!!#!        end if
!!!#!
!!!#!      END SUBROUTINE gauss_to_scf_1
!!!#!!!!!
!!!#!
!!!#!
!!!#!!!!!IMPLEMENTED
!!!#!      ! Called when ncplx_k = 2 and ncplx_g = 1
!!!#!      subroutine gauss_to_scf_2
!!!#!
!!!#!        !loop for each complex component
!!!#!        !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
!!!#!        !corrected for avoiding 0**0 problem
!!!#!        if (n_gau == 0) then
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              rk=real(i,gp)*h
!!!#!              r2=r/a1
!!!#!              r2=r2*r2
!!!#!              r2=0.5_gp*r2
!!!#!              cval=cos(kval*rk)
!!!#!              func=safe_exp(-r2)
!!!#!              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
!!!#!              ww(i-leftx,1,1)=func*cval
!!!#!              sval=sin(kval*rk)
!!!#!              ww(i-leftx,1,2)=func*sval
!!!#!           enddo
!!!#!        else
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              rk=real(i,gp)*h
!!!#!              coff=r**n_gau
!!!#!              r2=r/a1
!!!#!              r2=r2*r2
!!!#!              r2=0.5_gp*r2
!!!#!              cval=cos(kval*rk)
!!!#!              func=safe_exp(-r2)
!!!#!              !func=mp_exp(h,i0*16*h+z0,0.5_gp/(a1**2),i,0,.true.)
!!!#!              func=coff*func
!!!#!              ww(i-leftx,1,1)=func*cval
!!!#!              sval=sin(kval*rk)
!!!#!              ww(i-leftx,1,2)=func*sval
!!!#!           enddo
!!!#!        end if
!!!#!
!!!#!      END SUBROUTINE gauss_to_scf_2
!!!#!!!!
!!!#!
!!!#!      ! Called when ncplx_k = 1 and ncplx_g = 2
!!!#!      ! no k-points + complex Gaussians
!!!#!      subroutine gauss_to_scf_3
!!!#!
!!!#!        if (n_gau == 0) then
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              if( abs(r)-gcut < 1e-8 ) then
!!!#!                 r2=r*r
!!!#!                 cval=cos(a2*r2)
!!!#!                 sval=sin(a2*r2)
!!!#!                 r2=0.5_gp*r2/(a1**2)
!!!#!                 func=safe_exp(-r2)
!!!#!                 ww(i-leftx,1,1)=func*cval
!!!#!                 ww(i-leftx,1,2)=func*sval
!!!#!              else
!!!#!                 ww(i-leftx,1,1:2)=0.0_wp
!!!#!              end if
!!!#!           enddo
!!!#!        else
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              if( abs(r)-gcut < 1e-8 ) then
!!!#!                 r2=r*r
!!!#!                 cval=cos(a2*r2)
!!!#!                 sval=sin(a2*r2)
!!!#!                 coff=r**n_gau
!!!#!                 r2=0.5_gp*r2/(a1**2)
!!!#!                 func=safe_exp(-r2)
!!!#!                 func=coff*func
!!!#!                 ww(i-leftx,1,1)=func*cval
!!!#!                 ww(i-leftx,1,2)=func*sval
!!!#!              else
!!!#!                 ww(i-leftx,1,1:2)=0.0_wp
!!!#!              end if
!!!#!           enddo
!!!#!        end if
!!!#!      END SUBROUTINE gauss_to_scf_3
!!!#!
!!!#!      ! Called when ncplx_k = 2 and ncplx_g = 2
!!!#!      subroutine gauss_to_scf_4
!!!#!
!!!#!        if (n_gau == 0) then
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              if( abs(r)-gcut < 1e-8 ) then
!!!#!                 r2=r*r
!!!#!                 cval=cos(a2*r2)
!!!#!                 sval=sin(a2*r2)
!!!#!                 rk=real(i,gp)*h
!!!#!                 cval2=cos(kval*rk)
!!!#!                 sval2=sin(kval*rk)
!!!#!                 r2=0.5_gp*r2/(a1**2)
!!!#!                 func=safe_exp(-r2)
!!!#!                 ww(i-leftx,1,1)=func*(cval*cval2-sval*sval2)
!!!#!                 ww(i-leftx,1,2)=func*(cval*sval2+sval*cval2)
!!!#!              else
!!!#!                 ww(i-leftx,1,1:2)=0.0_wp
!!!#!              end if
!!!#!           enddo
!!!#!        else
!!!#!           do i=leftx,rightx
!!!#!              x=real(i-i0*16,gp)*h
!!!#!              r=x-z0
!!!#!              r2=r*r
!!!#!              cval=cos(a2*r2)
!!!#!              sval=sin(a2*r2)
!!!#!              rk=real(i,gp)*h
!!!#!              cval2=cos(kval*rk)
!!!#!              sval2=sin(kval*rk)
!!!#!              coff=r**n_gau
!!!#!              r2=0.5_gp*r2/(a1**2)
!!!#!              func=safe_exp(-r2)
!!!#!              func=coff*func
!!!#!              ww(i-leftx,1,1)=func*(cval*cval2-sval*sval2)
!!!#!              ww(i-leftx,1,2)=func*(cval*sval2+sval*cval2)
!!!#!           enddo
!!!#!        end if
!!!#!      END SUBROUTINE gauss_to_scf_4
!!!#!
!!!#!      ! Original version
!!!#!      !  subroutine gauss_to_scf
!!!#!      !    n_left=lefts(0)
!!!#!      !    n_right=rights(0)
!!!#!      !    length=n_right-n_left+1
!!!#!      !
!!!#!      !    !print *,'nleft,nright',n_left,n_right
!!!#!      !
!!!#!      !    do k=1,4
!!!#!      !       rights(k)=2*rights(k-1)+m
!!!#!      !       lefts( k)=2*lefts( k-1)-m
!!!#!      !    enddo
!!!#!      !
!!!#!      !    leftx = lefts(4)-n
!!!#!      !    rightx=rights(4)+n  
!!!#!      !
!!!#!      !    !stop the code if the gaussian is too extended
!!!#!      !    if (rightx-leftx > nwork) then
!!!#!      !       !STOP 'gaustodaub'
!!!#!      !       return
!!!#!      !    end if
!!!#!      !
!!!#!      !    !loop for each complex component
!!!#!      !    do icplx=1,ncplx
!!!#!      !
!!!#!      !       !calculate the expansion coefficients at level 4, positions shifted by 16*i0 
!!!#!      !
!!!#!      !       !corrected for avoiding 0**0 problem
!!!#!      !       if (ncplx==1) then
!!!#!      !          if (n_gau == 0) then
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                ww(i-leftx,1,icplx)=func
!!!#!      !             enddo
!!!#!      !          else
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                coeff=r**n_gau
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                func=real(coeff,wp)*func
!!!#!      !                ww(i-leftx,1,icplx)=func
!!!#!      !             enddo
!!!#!      !          end if
!!!#!      !       else if (icplx == 1) then
!!!#!      !          if (n_gau == 0) then
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                rk=real(i,gp)*h
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                cval=real(cos(kval*rk),wp)
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                ww(i-leftx,1,icplx)=func*cval
!!!#!      !             enddo
!!!#!      !          else
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                rk=real(i,gp)*h
!!!#!      !                coeff=r**n_gau
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                cval=real(cos(kval*rk),wp)
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                func=real(coeff,wp)*func
!!!#!      !                ww(i-leftx,1,icplx)=func*cval
!!!#!      !             enddo
!!!#!      !          end if
!!!#!      !       else if (icplx == 2) then
!!!#!      !          if (n_gau == 0) then
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                rk=real(i,gp)*h
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                sval=real(sin(kval*rk),wp)
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                ww(i-leftx,1,icplx)=func*sval
!!!#!      !             enddo
!!!#!      !          else
!!!#!      !             do i=leftx,rightx
!!!#!      !                x=real(i-i0*16,gp)*h
!!!#!      !                r=x-z0
!!!#!      !                rk=real(i,gp)*h
!!!#!      !                coeff=r**n_gau
!!!#!      !                r2=r/a
!!!#!      !                r2=r2*r2
!!!#!      !                r2=0.5_gp*r2
!!!#!      !                sval=real(sin(kval*rk),wp)
!!!#!      !                func=real(dexp(-real(r2,kind=8)),wp)
!!!#!      !                func=real(coeff,wp)*func
!!!#!      !                ww(i-leftx,1,icplx)=func*sval
!!!#!      !             enddo
!!!#!      !          end if
!!!#!      !       end if
!!!#!      !
!!!#!      !       !print *,'here',gau_a,gau_cen,n_gau
!!!#!      !       call apply_w(ww(0,1,icplx),ww(0,2,icplx),&
!!!#!      !            leftx   ,rightx   ,lefts(4),rights(4),h)
!!!#!      !
!!!#!      !       call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
!!!#!      !            lefts(4),rights(4),lefts(3),rights(3)) 
!!!#!      !       call forward_c(ww(0,1,icplx),ww(0,2,icplx),&
!!!#!      !            lefts(3),rights(3),lefts(2),rights(2)) 
!!!#!      !       call forward_c(ww(0,2,icplx),ww(0,1,icplx),&
!!!#!      !            lefts(2),rights(2),lefts(1),rights(1)) 
!!!#!      !
!!!#!      !       call forward(  ww(0,1,icplx),ww(0,2,icplx),&
!!!#!      !            lefts(1),rights(1),lefts(0),rights(0)) 
!!!#!      !
!!!#!      !    end do
!!!#!      !
!!!#!      !
!!!#!      !  END SUBROUTINE gauss_to_scf
!!!#!

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



    pure subroutine determine_bounds(nres,periodic,nstart,nmax,right_t,i0,&
         lefts,rights)
      implicit none
      logical, intent(in) :: periodic !<boundary conditions
      integer, intent(in) :: nres !< number of level of resolution
      integer, intent(in) :: nstart !< starting point the output array
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
         lefts( 0)=max(i0-right_t,  nstart)
         rights(0)=min(i0+right_t,nmax+nstart)
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
      real(gp), dimension(nwork,2), intent(inout) :: ww
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

      call apply_w(ng,ww(1,1),ww(1,2),&
           lefts(nres+1),rights(nres+1),lefts(nres),rights(nres),h)

      inw=2
      do ires=nres,2,-1
         inwo=3-inw
         call forward_c(ng,ww(1,inw),ww(1,inwo),&
              lefts(ires),rights(ires),lefts(ires-1),rights(ires-1)) 
         inw=3-inw
      end do
      call forward(ng,ww(1,inw),ww(1,3-inw),&
           lefts(1),rights(1),lefts(0),rights(0)) 
      inw=3-inw
    end subroutine magic_idwts

    !> One of the tails of the Gaussian is folded periodically
    !! We assume that the situation when we need to fold both tails
    !! will never arise
    subroutine retrieve_results(periodic,ng,n_left,n_right,nst,nmax,ncplx_w,ww,theor_norm,ncplx_f,fac,ncplx,c,error)
      use module_base, only: f_memcpy,dot,f_zero
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: n_left,n_right,nmax,ng,ncplx_w,ncplx_f,ncplx,nst
      !> the expected normalization of the gaussian. When different from zero,
      !! control the result and write the error in the error variable
      real(wp), dimension(ng), intent(in) :: theor_norm 
      real(wp), dimension(ncplx_f,ng), intent(in) :: fac
      real(wp), dimension(ncplx_w,ng,n_left:n_right,2), intent(in) :: ww !< the values of the gaussian to be folded
      real(wp), dimension(ng), intent(out) :: error !< the error in the expression of the function
      real(wp), dimension(ncplx,ng,2,nst:nmax+nst), intent(out) :: c !< final results
      !local variables
      integer :: i,j,inl,ig
      real(wp) :: cn2,cfr,cfi

      !put to zero the results first
      call f_zero(c)
      call f_zero(error)

      if (ncplx_w == 1) then
         !for a real result retrieve also the normalization
         if (periodic) then
            do i=n_left,n_right
               j=modulo(i,nmax+1)
               do ig=1,ng
                  c(1,ig,1,j)=c(1,ig,1,j)+ww(1,ig,i,1)
                  c(1,ig,2,j)=c(1,ig,2,j)+ww(1,ig,i,2)
                  error(ig)=error(ig)+ww(1,ig,i,1)**2+ww(1,ig,i,2)**2
               end do
            end do
         else
            do i=n_left,n_right
               do ig=1,ng
                  c(1,ig,1,i)=ww(1,ig,i,1)
                  c(1,ig,2,i)=ww(1,ig,i,2)
                  error(ig)=error(ig)+ww(1,ig,i,1)**2+ww(1,ig,i,2)**2
               end do
            end do
         end if
         do ig=1,ng
            if (theor_norm(ig) /= 0.0_wp) then
               !calculate the (relative) error
               !cn2=dot(2*(n_right-n_left+1),ww(1,1),1,ww(1,1),1)
               error(ig)=sqrt(abs(1.0_gp-error(ig)/theor_norm(ig)))
            end if
         end do
      else
         if (periodic) then
            do i=n_left,n_right
               j=modulo(i,nmax+1)
               do ig=1,ng
                  c(1,ig,1,j)=c(1,ig,1,j)+ww(1,ig,i,1)
                  c(2,ig,1,j)=c(2,ig,1,j)+ww(2,ig,i,1)
                  c(1,ig,2,j)=c(1,ig,2,j)+ww(1,ig,i,2)
                  c(2,ig,2,j)=c(2,ig,2,j)+ww(2,ig,i,2)
               end do
            end do
         else
            do i=n_left,n_right
               do ig=1,ng
                  c(1,ig,1,i)=ww(1,ig,i,1)
                  c(2,ig,1,i)=ww(2,ig,i,1)
                  c(1,ig,2,i)=ww(1,ig,i,2)
                  c(2,ig,2,i)=ww(2,ig,i,2)
               end do
            end do
         end if
      end if

      !write(*,*)'error, non scaled:',error
      !
      !RESCALE BACK THE COEFFICIENTS AND THE ERROR
      if (ncplx_f == 1) then
         do i=nst,nmax+nst
            do inl=1,2
               do ig=1,ng
                  cfr=c(1,ig,inl,i)
                  c(1,ig,inl,i)=fac(1,ig)*cfr
               end do
            end do
         end do
         error=error*fac(1,:)         
      else if (ncplx_w == 1) then
         do i=nst,nmax+nst
            do inl=1,2
               do ig=1,ng
                  cfr=c(1,ig,inl,i)
                  c(1,ig,inl,i)=fac(1,ig)*cfr
                  c(2,ig,inl,i)=fac(2,ig)*cfr
               end do
            end do
         end do
      else
         do i=nst,nmax+nst
            do inl=1,2
               do ig=1,ng
                  cfr=c(1,ig,inl,i)
                  cfi=c(2,ig,inl,i)
                  c(1,ig,inl,i)=fac(1,ig)*cfr-fac(2,ig)*cfi
                  c(2,ig,inl,i)=fac(2,ig)*cfr+fac(1,ig)*cfi
               end do
            end do
         end do
      end if

    end subroutine retrieve_results

    !> express a function defined as x**n_gau*exp(-(x-x0)**2/(2*a**2))
    !! Multiply it for the k-point factor exp(Ikx)
    !! For this reason, a real (cos(kx)) and an imaginary (sin(kx)) part are provided 
    !! also, the value of a might be complex.
    !! therefore the result has to be complex accordingly
   
    pure subroutine gaus_highres(level,ng,n_gau,ncplx_a,a,x0,kval,ncplx,&
         leftx,rightx,ww)
      implicit none
      integer, intent(in) :: level !< resolution level of the evaluation
      integer, intent(in) :: ng !< number of gaussians
      !> complex exponents and final results
      integer, intent(in) :: ncplx,ncplx_a 
      integer, dimension(ng), intent(in) :: n_gau !< parameters of the gaussian
      integer, intent(in) :: leftx !<leftmost point of the grid
      integer, intent(in) :: rightx !<rightmost point of the grid
      !>standard deviations (exponent)
      real(gp), dimension(ncplx_a,ng), intent(in) :: a 
      real(gp), intent(in) :: x0 !<gaussian center
      real(gp), intent(in) :: kval !<k-point value
      !>collcation values of the gaussian
      real(gp), dimension(ncplx,ng,leftx:rightx), intent(out) :: ww 
      !local variables
      integer :: i,i0,resolution,ig
      real(gp) :: x,h,z0,r,rk

      resolution=2**level
      !the grid spacing is now given in terms of the level
      h=1.0_gp/(real(resolution,gp))
      !the grid point closest to the center
      i0=nint(x0) ! the array is centered at i0
      z0=x0-real(i0,gp)

      !print *,'XXXXXXXXX',ncplx_a,kval,leftx,rightx
      !this part has to be extended to the different cases
      if (ncplx_a==1) then
         if (kval==0.0_gp) then
            do i=leftx,rightx
               x=real(i-i0*resolution,gp)*h
               r=x-z0
               do ig=1,ng
                  ww(1,ig,i)=collocate_gaussian_r0(n_gau(ig),a(1,ig),r)
               end do
            end do
         else
            do i=leftx,rightx
               x=real(i-i0*resolution,gp)*h
               r=x-z0
               rk=real(i,gp)*h
               do ig=1,ng
                  ww(1:2,ig,i)=collocate_gaussian_rk(n_gau(ig),a(1,ig),&
                       r,kval,rk)
               end do
            end do
         end if
      else
         if (kval==0.0_gp) then
            do i=leftx,rightx
               x=real(i-i0*resolution,gp)*h
               r=x-z0
               do ig=1,ng
                  ww(1:2,ig,i)=collocate_gaussian_c0(n_gau(ig),a(1:2,ig),r)
               end do
            end do
         else
            do i=leftx,rightx
               x=real(i-i0*resolution,gp)*h
               r=x-z0
               rk=real(i,gp)*h
               do ig=1,ng
                  ww(1:2,ig,i)=collocate_gaussian_ck(n_gau(ig),a(1:2,ig),&
                       r,kval,rk)
               end do
            end do

         end if
      end if
    end subroutine gaus_highres

    !> complex exponent, k-point
    !! @warning: the meaning of the exponent is ambiguous:
    !! the real part corresponds to the standard deviation 
    !! and the imaginary part is the exponent
    !! this will have to be fixed in the future
    pure function collocate_gaussian_ck(n_gau,a,r,k,rk) result(func)
      use numerics, only: safe_exp
      implicit none
      integer, intent(in) :: n_gau !<principal quantum number
      real(gp), dimension(2), intent(in) :: a !<standard deviation
      real(gp), intent(in) :: r !< function argument
      real(gp), intent(in) :: k !< k-point coordinate
      real(gp), intent(in) :: rk !< argument of k-point value
      real(gp), dimension(2) :: func
      !local variables
      real(gp) :: r2,fn,cval,sval,cv2,sv2

      r2=r*r
      cval=cos(a(2)*r2)
      sval=sin(a(2)*r2)
      r2=r/a(1)
      r2=r2*r2
      r2=0.5_gp*r2
      cv2=cos(k*rk)
      sv2=sin(k*rk)        
      fn=safe_exp(-r2)
      !fn=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.true.)
      fn=coeff(n_gau,r)*fn
      func(1)=fn*(cval*cv2-sval*sv2)
      func(2)=fn*(cval*sv2+sval*cv2)

    end function collocate_gaussian_ck

    !> complex exponent, gamma point
    !! @warning: the meaning of the exponent is ambiguous:
    !! the real part corresponds to the standard deviation 
    !! and the imaginary part is the exponent
    !! this will have to be fixed in the future
    pure function collocate_gaussian_c0(n_gau,a,r) result(func)
      use numerics, only: safe_exp
      implicit none
      integer, intent(in) :: n_gau !<principal quantum number
      real(gp), dimension(2), intent(in) :: a !<standard deviation
      real(gp), intent(in) :: r !< function argument
      real(gp), dimension(2) :: func
      !local variables
      real(gp) :: r2,fn,cval,sval

      r2=r*r
      cval=cos(a(2)*r2)
      sval=sin(a(2)*r2)
      r2=r/a(1)
      r2=r2*r2
      r2=0.5_gp*r2
!!$      cval=cos(k*rk)
!!$      sval=sin(k*rk)        
      fn=safe_exp(-r2)
      fn=coeff(n_gau,r)*fn
      func(1)=fn*cval
      func(2)=fn*sval

    end function collocate_gaussian_c0

    !> real exponent, k-point
    pure function collocate_gaussian_rk(n_gau,a,r,k,rk) result(func)
      use numerics, only: safe_exp
      implicit none
      integer, intent(in) :: n_gau !<principal quantum number
      real(gp), intent(in) :: a !<standard deviation
      real(gp), intent(in) :: r !< function argument
      real(gp), intent(in) :: k !< k-point coordinate
      real(gp), intent(in) :: rk !< argument of k-point value
      real(gp), dimension(2) :: func
      !local variables
      real(gp) :: r2,fn,cval,sval

      r2=r/a
      r2=r2*r2
      r2=0.5_gp*r2
      cval=cos(k*rk)
      sval=sin(k*rk)        
      fn=safe_exp(-r2)
      fn=coeff(n_gau,r)*fn
      func(1)=fn*cval
      func(2)=fn*sval

    end function collocate_gaussian_rk

    !> real exponent, gamma point
    pure function collocate_gaussian_r0(n_gau,a,r) result(func)
      use numerics, only: safe_exp
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

    end function collocate_gaussian_r0
    
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
