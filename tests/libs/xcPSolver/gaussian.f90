!> @file

!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2013-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Program testing new ideas like momentum-preserving gaussian integrals
program MP_gaussian
  use module_base
  use gaussians
  use yaml_output
  use gaussdaub
  implicit none
  integer, parameter :: iplot=14,iunit=16 !< File unit for the plots
  integer, parameter :: nmoms=16          !< Number of calculated moments
  integer, parameter :: nstep=1           !< Number of resolution to calculate the moments
  integer, parameter :: nsigma=30         !< Number of different gaussian functions
  integer, parameter :: npts=16           !< Arrays from -npts to npts
  real(gp), parameter :: hgrid = 1.0_gp   !< Step grid
  integer, parameter :: nw=65536
  real(gp), dimension(0:nw,2,2) :: work
  integer :: i,j,imoms,pow,istep,isigma,unit
  integer, parameter :: itype_scf = 16
  integer :: n_scf,untplot,ml1,mu1
  integer :: istart, iend, i0, npf
  real(gp) :: pgauss,p0gauss,x0,x00,reference,max_phi,max_lag,scalar,tt,hh
  character(len=128) :: filename
  double precision, dimension(17), parameter :: MFdat=[&
       8.4334247333529341094733325815816d-7,&
       -0.1290557201342060969516786758559028d-4,&
       0.8762984476210559564689161894116397d-4,&
       -0.30158038132690463167163703826169879d-3,&
       0.174723713672993903449447812749852942d-2,&
       -0.942047030201080385922711540948195075d-2,&
       0.2373821463724942397566389712597274535d-1,&
       0.612625895831207982195380597d-1,&
       0.9940415697834003993178616713d0,&
       -0.604895289196983516002834636d-1,&
       -0.2103025160930381434955489412839065067d-1,&
       0.1337263414854794752733423467013220997d-1,&
       -0.344128144493493857280881509686821861d-2,&
       0.49443227688689919192282259476750972d-3,&
       -0.5185986881173432922848639136911487d-4,&
       2.72734492911979659657715313017228d-6,0.d0]
  real(gp), dimension(0:nmoms,3) :: moments
  real(gp), dimension(3,3,0:nmoms) :: avgmaxmin
  real(gp), dimension(:), allocatable :: fj_phi,fj_coll,fj_lag, psi_phi,psi_coll
  real(gp), dimension(:,:), allocatable :: f_mu, psi_mu
  real(gp), dimension(:), allocatable :: scf_dat,lag_dat
  real(gp), dimension(:), allocatable :: x_scf
  call f_lib_initialize()

!  call bacasable_valence()
! 
!  call f_lib_finalize()
!  stop

  pow=0
  unit=iunit+1
  untplot=iplot

  npf=2*npts-1
  hh=0.5d0*hgrid

  !pgauss=0.5_gp/((0.1_gp*hgrid)**2)!8.0e-3_dp*1.25_dp**(6*(8-1))
  !array where we have to write the value of the discretization
  fj_phi=f_malloc(-npf .to. npf,id='fj_phi')
  fj_lag=f_malloc(-npf .to. npf,id='fj_lag')
  fj_coll=f_malloc(-npf .to. npf,id='fj_coll')
  f_mu=f_malloc( [ -npts .to. npts, 1 .to. 2 ], id='f_mu')
  call initialize_real_space_conversion(npoints=2**8,isf_m=itype_scf,nmoms=0) !initialize the work arrays needed to integrate with isf

  psi_phi=f_malloc(-npf .to. npf,id='psi_phi')
  psi_coll=f_malloc(-npf .to. npf,id='psi_coll')
  psi_mu=f_malloc( [ -npts .to. npts, 1 .to. 2 ], id='psi_mu')

  call f_open_file(untplot,'Families.dat')
  call polynomial_exactness(npts,16,0,itype_scf,16,16,untplot,0,moments)
  call f_close(untplot)
  !stop

  !replot the smoothed magic filtered Daubechies
  call f_open_file(untplot,'daub_smoothed.dat')
  call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,8,MFdat)
  call f_close(untplot)

!  stop

  call f_open_file(unit,file='MultipolesError.dat')

  !Initialize the work array needed to integrate with the Lagrange polynomial
  n_scf = 2*itype_scf*(2**10)

  !Allocations for lagrange polynomial data array
  x_scf = f_malloc(0.to.n_scf)
  lag_dat = f_malloc(0.to.n_scf,id='lag_dat')
  scf_dat = f_malloc(0.to.n_scf,id='scf_dat')

  !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
  !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_dat)
  !call ISF_family(itype_scf,0,n_scf,n_range,x_scf,scf_dat)

  !call yaml_set_stream(record_length=150)
  call yaml_map('itype_scf',itype_scf)
  !call yaml_map('range', n_range)
  call yaml_map('number of points',n_scf)
  call yaml_map('step',x_scf(1)-x_scf(0))
  call yaml_map('dx',real(2*itype_scf,gp)/real(n_scf,gp))

  !Step grid for the integration
  !dx = real(2*itype_scf,gp)/real(n_scf,gp)
  !starting point for the x coordinate for integration
  !x  = real(-itype_scf,gp)-dx
!  write(iplot,'(a)') '#Abscissae   Interpolating_scaling_function   Lagrange_polynomial' 
!  do i=0,n_scf
!     !x=x+dx
!     lag_dat(i) = lag_sym(x_scf(i),itype_scf,0)
!     write(iplot,'(f19.13,1x,3(1pe23.15,1x))') x_scf(i),scf_dat(i),lag_dat(i)
!  end do

  !Do scalar product
  scalar = 0.d0
  do i=0,n_scf
     scalar = scalar + scf_dat(i)*lag_dat(i)
  end do
  scalar = scalar*(x_scf(1)-x_scf(0))
  call yaml_map('<lag_dat|scf-dat>',scalar)

  call yaml_sequence_open('<lag|iscf>')
  istart=-(itype_scf/2-1)
  iend=itype_scf/2-1
  do i0=-itype_scf/2+1,itype_scf/2-1
     scalar = 0.d0
     do i=0,n_scf
        scalar = scalar + scf_dat(i)*lagrange(x_scf(i),istart,iend,i0)
        !scalar = scalar + scf_dat(i)*lag_sym(x_scf(i),itype_scf,i0)
        !scalar = scalar + scf_dat(i)*x_scf(i)**i0
     end do
     scalar = scalar*(x_scf(1)-x_scf(0))
     !call yaml_map(trim(yaml_toa(istart))//trim(yaml_toa(iend))//trim(yaml_toa(i0)), scalar)
     !call yaml_map(trim(yaml_toa(i0)), (/ scalar,lag_sym(real(i0,gp),itype_scf,i0) /))
     call yaml_map(trim(yaml_toa(i0)), (/ scalar,lagrange(real(i0,gp),istart,iend,i0) /))
  end do
  call yaml_sequence_close()

  !Orbital is a gaussian (g0)
  isigma=10
  istep=5
  pgauss=pg(isigma,hgrid)
  x0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
  call evaluate_moments(nmoms,npf,hh,pgauss,pow,x0,psi_phi,psi_coll,fj_lag,moments)
  !Daubechies | gaussian
  call gauss_to_daub_k(hgrid,0.d0,1,1,1,1.d0,x0,sqrt(0.5_gp/pgauss),0,-npts,2*npts,ml1,mu1,&
      psi_mu,work,nw,.False.,1.d0) 
  p0gauss=pgauss
  x00=x0

  ! Calculate for different nsigma sigma
  do isigma=1,nsigma
     pgauss=pg(isigma,hgrid)
     call yaml_map('sigma/h',sqrt(0.5_gp/pgauss)/hgrid)

     avgmaxmin=0.d0
     avgmaxmin(3,:,:)=1.d100
     max_phi=0.0_gp
     max_lag=0.0_gp

     do istep=1,nstep
        x0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
        call evaluate_moments(nmoms,npf,hh,pgauss,pow,x0,fj_phi,fj_coll,fj_lag,moments)

        !Daubechies | gaussian
        call gauss_to_daub_k(hgrid,0.d0,1,1,1,1.d0,x0,sqrt(0.5_gp/pgauss),0,-npts,2*npts,ml1,mu1,&
            f_mu,work,nw,.False.,1.d0) 
        !call gau_daub_1d(.False.,1, x0+real(npts,gp)*hgrid, [ 0 ], 0.d0, 1, [ sqrt(0.5_gp/pgauss) ], &
        !     & 1, [ 1.d0 ], hgrid,4,2*npts,1,f_mu,nw,work)
        tt = sqrt(pi_param/(pgauss+p0gauss))*exp(-(x00-x0)**2*(pgauss*p0gauss)/(pgauss+p0gauss))
        call yaml_map('Projector norms', (/ &
              &  1.d0-hh*scpr(2*npf+1,fj_phi,psi_phi)/tt, &
              &  1.d0-hh*scpr(2*npf+1,fj_coll,psi_phi)/tt, &
              &  1.d0-scpr(4*npts+2,f_mu,psi_mu)/tt, &
              &  sqrt(0.5_gp/pgauss)/hgrid /),fmt='(1pg17.10)')
        !call yaml_map('Projector norms', (/ nrm2(2*npts+1,fj_phi(-npts),1)/tt, nrm2(2*npts+1,fj_coll(-npts),1)/tt, &
        !      &  nrm2(4*npts+2,f_mu(-npts,1),1)/tt, sqrt(0.5_gp/pgauss)/hgrid /),fmt='(1pg12.5)')

        call filename_test(x0,pgauss,hgrid,filename)
        untplot=iplot+istep
        call f_open_file(untplot,trim(filename)//'daub.dat')
        write(untplot,'(a,2(1pe20.8))') '#<daub|gaussian>',pgauss,sqrt(0.5_gp/pgauss)
        write(untplot,'(a)') '#j,f_mu(j)'
        do j=-npts,npts
           write(untplot,*) j,f_mu(j,:)
        end do
        call f_close(untplot)

        call filename_test(x0,pgauss,hgrid,filename)
        call f_open_file(untplot,trim(filename)//'coll.dat')
        call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,npts,fj_coll)
        call f_close(untplot)

        call filename_test(x0,pgauss,hgrid,filename)
        call f_open_file(untplot,trim(filename)//'isf.dat')
        call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,npts,fj_phi)
        call f_close(untplot)


        max_phi=max(max_phi,maxval(abs(fj_coll-fj_phi)))
        max_lag=max(max_lag,maxval(abs(fj_coll-fj_lag)))

        !calculate the average, maximum and minimum of each moment in function of the reference
        !j=1 use the elemental property of the mp_exp function with fj_phi
        !j=2 Lagrange polynomial
        !j=3 collocation array with fj_coll
        do j=1,3
           do imoms=0,nmoms
              reference=gauint0(pgauss,imoms+pow)
              !print *,j,imoms,reference,moments(imoms,j),abs((moments(imoms,j)-reference))
              if (reference /= 0.0_gp) then
                 !x^even
                 moments(imoms,j) = abs((moments(imoms,j)-reference)/reference)
              else
                 !x^odd
                 moments(imoms,j) = abs(moments(imoms,j))
              end if
              avgmaxmin(1,j,imoms) = avgmaxmin(1,j,imoms)+moments(imoms,j)/real(nstep,gp)
              avgmaxmin(2,j,imoms) = max(moments(imoms,j),avgmaxmin(2,j,imoms))
              avgmaxmin(3,j,imoms) = min(moments(imoms,j),avgmaxmin(3,j,imoms))
           end do
        end do
        call yaml_mapping_open('Normalised moments')
        call yaml_map('x0/h',x0/hgrid)
        call yaml_map('C',moments(:,2),fmt='(1pe14.5)')
        call yaml_map('MP',moments(:,1),fmt='(1pe14.5)')
        call yaml_mapping_close()
     end do
     call yaml_map('Results',reshape(avgmaxmin,[6,nmoms+1]),fmt='(1pe14.5)')
     !Plot fort.(iunit+1)
     avgmaxmin(2,:,:)=avgmaxmin(2,:,:)-avgmaxmin(3,:,:)
     write(unit,'(104(1pe14.5))') sqrt(0.5_gp/pgauss)/hgrid,avgmaxmin(1:2,1:3:2,:)
     call yaml_map('maxdiff' // trim(yaml_toa(isigma)), [ sqrt(0.5_gp/pgauss)/hgrid, max_phi, max_lag ] )
  end do
  call yaml_map('Results (phi)',reshape(avgmaxmin(1:3,1:3:2,:),[6,nmoms+1]),fmt='(1pe14.5)')
  call yaml_map('Results (lag)',reshape(avgmaxmin(1:3,2:3,:),[6,nmoms+1]),fmt='(1pe14.5)')

  call f_close(unit)

  call finalize_real_space_conversion()

  call f_free(fj_phi,fj_coll,fj_lag)
  call f_free(psi_phi,psi_coll)
  call f_free(f_mu)
  call f_free(psi_mu)
  call f_free(scf_dat,lag_dat,x_scf)
  call f_lib_finalize()


contains

  !> Determine pgauss
  function pg(i,h)
    implicit none
    integer, intent(in) :: i
    real(gp), intent(in) :: h
    real(gp) :: pg
    !pg=0.5_gp/((0.1_gp+0.1_gp*real(i-1,gp)*hgrid)**2)
    pg=1.0_gp/((0.1_gp+0.1_gp*real(i-1,gp)*h)**2)
  end function pg

  !> Scalar products
  function scpr(n,a,b) result(s)
    implicit none
    integer, intent(in) :: n
    real(gp), dimension(n), intent(in) :: a,b
    real(gp) :: s
    s = (dot(n,a(1),1,b(1),1))
  end function scpr

  !> Classify the quality of a multipole extraction in both cases
  subroutine evaluate_moments(nmoms,npts,hgrid,pgauss,pow,x0,fj_phi,fj_coll,fj_lag,moments)
    use module_base, only: gp,safe_exp,f_open_file,f_close,yaml_toa
    use gaussians, only: mp_exp, scfdotf
    implicit none
    !Arguments
    integer, intent(in) :: npts,pow,nmoms
    real(gp), intent(in) :: hgrid,pgauss,x0
    real(gp), dimension(0:nmoms,3), intent(out) :: moments
    real(gp), dimension(-npts:npts), intent(out) :: fj_phi,fj_coll,fj_lag
    integer, parameter :: iunit = 100
    !local variables
    integer :: j,unit
    integer :: istep = 0
    character(len=128) :: filename

    call filename_test(x0,pgauss,hgrid,filename)
    unit=iunit+istep
    call f_open_file(unit,trim(filename))

    !use the elemental property of the mp_exp function
    fj_phi=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.true.)
    !  fj_phi=scfdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
    !scfdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
    call moments_1d(2*npts+1,fj_phi,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,1))

    !Lagrange polynomial array
    fj_lag=lagdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
    call moments_1d(2*npts+1,fj_lag,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,2))

    !collocation array
    fj_coll=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.false.)

    call moments_1d(2*npts+1,fj_coll,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,3))

    write(unit,'(a)') '#Projection of a gaussian with iscf, collocation and lagrange polynomials method'
    write(unit,'(a)') '#j,fj_phi(j),fj_coll(j),fj_lag(j)'
    do j=-npts,npts
       write(unit,*) j,fj_phi(j),fj_coll(j)
    end do
    istep = istep + 1

    call f_close(unit)

  end subroutine evaluate_moments

  !> Get the test filename
  subroutine filename_test(x0,pgauss,hgrid,filename)
    use yaml_strings
    implicit none
    double precision, intent(in) :: x0,pgauss,hgrid
    character(len=128), intent(out) :: filename

    call f_strcpy(dest=filename,src=&
         'gau'//trim(adjustl(yaml_toa(x0,fmt='(f5.2)')))//&
         'p'//trim(adjustl(yaml_toa(pgauss,fmt='(f5.2)')))//&
         'h'//trim(adjustl(yaml_toa(hgrid,fmt='(f5.2)'))))


!     call f_strcpy(dest=filename,src=&
!         'gau'+x0**'(f5.2)'+'p'+pgauss**'(f5.2)'+'h'+hgrid**'(f5.2)')


  end subroutine filename_test

  !> Calculate the moments of an array with respect to a reference point 
  subroutine moments_1d(n,array,x0,h,nmoms,moments)
    use module_base, only:gp
    implicit none
    !Arguments
    integer, intent(in) :: nmoms,n
    real(gp), intent(in) :: x0 !< Reference point
    real(gp), intent(in) :: h  !< Grid spacing
    real(gp), dimension(n), intent(in) :: array
    real(gp), dimension(0:nmoms), intent(out) :: moments
    !local variables
    integer :: j,kh,kl
    !iteger :: k,nc,nh
    real(gp) :: x

    moments(0)=h*sum(array)
    do j=1,nmoms
       moments(j)=0.0_gp
       !should start from the middle, but start from the exterior
       kl=1
       kh=n
       loop_j: do while(kl <= kh)
          x=real(kl,gp)*h-x0
          moments(j)=moments(j)+x**j*array(kl)
          if (kh /= kl) then
             x=real(kh,gp)*h-x0
             moments(j)=moments(j)+x**j*array(kh)
          end if
          kl=kl+1
          kh=kh-1
       end do loop_j
       
!!$       do k=1,n
!!$          x=real(k,gp)*h-x0
!!$          !if (j==0) then
!!$          !   moments(j)=moments(j)+array(k)
!!$          !else
!!$          moments(j)=moments(j)+x**j*array(k)
!!$          !end if
!!$       end do
       moments(j)=moments(j)*h
    end do

  end subroutine moments_1d


  !> This function calculates the scalar product between a Lagrange polynomial and a 
  !! input function, which is a gaussian times a power centered
  !! @f$g(x) = (x-x_0)^{pow} e^{-pgauss (x-x_0)}@f$
  !! here pure specifier is redundant
  !! we should add here the threshold from which the 
  !! normal function can be evaluated
  elemental pure function lagdotf(j,hgrid,pgauss,x0,pow) result(gint)
    implicit none
    !Arguments
    integer, intent(in) :: j !<value of the input result in the hgrid reference
    integer, intent(in) :: pow
    real(gp), intent(in) :: hgrid,pgauss,x0
    real(gp) :: gint
    !local variables
    integer :: i
    real(gp) :: x,absci,fabsci,dx
    gint=0.0_gp

    !Step grid for the integration
    dx = real(2*itype_scf,gp)/real(n_scf,gp)
    !starting point for the x coordinate for integration
    x  = real(j-itype_scf+1,gp)-dx

    !the loop can be unrolled to maximize performances
    if (pow /= 0) then
       do i=0,n_scf
          x=x+dx
          absci = x*hgrid - x0
          !here evaluate the function
          fabsci = absci**pow
          absci = -pgauss*absci*absci
          fabsci = fabsci*dexp(absci)
          !calculate the integral
          gint = gint + lag_dat(i)*fabsci
          !       print *,'test',i,scf_dat(i),fabsci,pgauss,pow,absci
       end do
    else
       do i=0,n_scf
          x=x+dx
          absci = x*hgrid - x0
          !here evaluate the function
          absci = -pgauss*absci*absci
          fabsci = safe_exp(absci)
          !calculate the integral
          gint = gint + lag_dat(i)*fabsci
          !       print *,'test',i,scf_dat(i),fabsci,pgauss,pow,absci
       end do
    end if
    gint = gint*dx

  end function lagdotf


  !> Calculate the Lagrange polynomial
  !! @f$ l_{i_0}(x) = \prod_{i=istart, i \neq i_0}^{iend} \frac{x-i}{i0-i}$
  elemental function lagrange(x,istart,iend,i0) result(y)
    implicit none
    !Arguments
    real(gp), intent(in) :: x
    integer, intent(in) :: istart !< First point
    integer, intent(in) :: iend   !< Last point
    integer, intent(in) :: i0     !< Center (removed from the product)
    real(gp) :: y
    !Local variables
    integer :: i
    !if (x >= real(istart,kind=gp) .and. x<= real(iend,kind=gp)) then
    y = 1.d0
    do i=istart,iend
       if (i /= i0) then
          y = y *(x-real(i,gp))/real(i0-i,gp)
       end if
    end do
    !else
    !  y = 0.d0
    !end if
  end function lagrange


  !> Calculate the Lagrange polynomial
  !! @f$ l_{i_0}(x) = \prod_{i=istart, i \neq i_0}^{iend} \frac{x-i}{i0-i}$
  elemental function lag_sym(x,itypescf,i0) result(y)
    implicit none
    !Arguments
    real(gp), intent(in) :: x
    integer, intent(in) :: itypescf !< order of scf (should be even)
    integer, intent(in) :: i0     !< Center 
    real(gp) :: y
    !Local variables
    integer :: i,ii
    !if (abs(x-real(i0,gp)) < 1d-4) then
    !  y = 1.d0/0.0009765625
    !else
    !  y = 0.d0
    !end if
    !return
    ii = itypescf/2-1
    if (abs(x) <= real(ii,gp)) then
       y = 1.d0
       do i=-ii,ii
          !if (i /= i0) then
          if (i /= 0) then
             y = y *(x-real(i0+i,gp))/real(-i,gp)
          end if
       end do
    else
       y = 0.d0
    end if
  end function lag_sym

end program


!> Verify the property x**p = sum_j x_j**p \phi_j(x) for different families
!! of interpolating functions
subroutine polynomial_exactness(npts,itype_scf,nmoms,p,itype_scf_dual,nmoms_dual,iplot,n_l,fl)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: npts,itype_scf,nmoms,p,iplot,itype_scf_dual,nmoms_dual,n_l
  double precision, dimension(-n_l:n_l), intent(in) :: fl !< coefficients of the interpolations
  !local variables
  integer :: n_scf,j,q,i,n_range,n_ext,n_ranget,n_scft,n_lag,n_lagt,n_lagvdm,imin,imax
  real(gp) :: phival,scalar
  real(gp), dimension(0:p) :: deviations
  real(gp), dimension(:), allocatable :: scf_dat,x_scf,x_scft,scft_dat!,f_l
  real(gp), dimension(:,:), allocatable :: pol_dat,pol_sp,Aiq
  n_lag=0
  if (itype_scf<0)  n_lag=-itype_scf
  !Initialize the work array needed to integrate with the Lagrange polynomial
  n_scf = 2*(abs(itype_scf)+nmoms)*npts
  n_lagvdm=abs(itype_scf)/2
  imin=-n_lagvdm+1
  imax=n_lagvdm
  if (n_l > 0) then
     do i=-n_l,n_l
        if (abs(fl(i))>1.d-12) then
           imin=min(imin,i)
           imax=max(imax,i)
        end if
     end do
     n_lagvdm=imax-imin+1
     if (modulo(n_lagvdm,2) /=0) then
        n_lagvdm=n_lagvdm+1
        imin=imin-1
     end if
  end if
  call yaml_map('Lagrange chosen',[imin,imax,n_lagvdm])
  !Allocations for lagrange polynomial data array
  x_scf = f_malloc(0.to.n_scf)
  scf_dat = f_malloc(0.to.n_scf,id='scf_dat')
  pol_dat = f_malloc0([0.to.n_scf,0.to.p],id='pol_dat')
     !f_l = f_malloc0(0.to.n_scf,id='f_l')
     Aiq= f_malloc([0.to.2*n_lagvdm-1,imin.to.imax],id='Aiq')
     do i=0,0!5
        call invert_vandermonde(imin+i,imax+i,Aiq)
     end do

  !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
  !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_dat)
  if (n_lag==0) then
     call ISF_family(itype_scf,nmoms,n_scf,n_range,x_scf,scf_dat)
  else
     call lagrange_family(n_lag,n_scf,n_range,x_scf,scf_dat)
  end if

  !override scf_dat with lagrange polynomials
  !n_range=4

  pol_sp = f_malloc0([-n_range.to.n_range,0.to.p],id='pol_sp')

  !call yaml_set_stream(record_length=150)
  call yaml_map('itype_scf',itype_scf)
  call yaml_map('nmoms',nmoms)
  call yaml_map('range', n_range)
  call yaml_map('number of points',n_scf)
  call yaml_map('step',x_scf(1)-x_scf(0))
  call yaml_map('dx',real(2*(abs(itype_scf)+nmoms),gp)/real(n_scf,gp))

  if (n_l > 0) then
     do i=0,n_scf
        write(iplot,'(f19.13,1x,30(1pe26.15e3,1x))') x_scf(i),phi(x_scf(i),0),&
             lagrange(x_scf(i),-n_range+1,n_range-1,0),&
             !sum([(sum([(fl(j)*Aiq(q,j)*x_scf(i)**q,q=0,2*n_lagvdm-1)]),j=imin,imax)]),&
             sum([(phi(x_scf(i),j)*fl(j),j=-n_l,n_l)])
             !sum([(Aiq(q,0)*sum([(phi(x_scf(i),j)*fl(j),j=-n_l,n_l)]),q=0,2*n_lagvdm-1)])!,&
             !sum([(Aiq(q,0)*sum([(phi(x_scf(i),j)*fl(j),j=-n_l,n_l)]),q=0,2*n_lagvdm-1)])
        !(sum([(Aiq(q,j)*x_scf(i)**q,q=0,2*n_lagvdm)]),j=-n_lagvdm,n_lagvdm)
     end do
  else
     n_ext=n_range/2
     write(iplot,'(a)') '#Abscissae   Interpolating_scaling_function   Lagrange_polynomial' 
     do q=0,p
        !evaluate the maxima
        do j=-n_range,n_range
           do i=0,n_scf
              phival=phi(x_scf(i),j)
              pol_dat(i,q)=pol_dat(i,q)+(real(j,kind=8)/real(n_ext,kind=8))**q*phival
              !change of variables here
              if (n_lag == 0) then
                 pol_sp(j,q)=pol_sp(j,q)+&
                      ((x_scf(i)+real(j,kind=8))/real(n_ext,kind=8))**q*scf_dat(i)
              else
                 pol_sp(j,q)=pol_sp(j,q)+&
                      ((x_scf(i))/real(n_ext,kind=8))**q*phival
              end if
           end do
           pol_sp(j,q)=pol_sp(j,q)*(x_scf(1)-x_scf(0))
           !print *,'j,q',j,q,pol_sp(j,q),(real(j,kind=8)/real(n_ext,kind=8))**q
           pol_sp(j,q)=abs(pol_sp(j,q)-(real(j,kind=8)/real(n_ext,kind=8))**q)
        end do
     end do
     call yaml_map('Polynomial scalprods <p|phi_j> - j**p',[(maxval(pol_sp(:,q)),q=0,p)])
     !evaluate the maxima
     deviations=0.d0
     do i=0,n_scf
        forall (q=0:p) deviations(q)=max(deviations(q),abs((x_scf(i)/real(n_ext,kind=8))**q-pol_dat(i,q)))
        !write(iplot,'(f19.13,1x,30(1pe23.15,1x))') x_scf(i),phi(x_scf(i),0),pol_dat(i,:),lagrange(x_scf(i),-n_range+1,n_range,0)
        write(iplot,'(f19.13,1x,30(1pe23.15,1x))') x_scf(i),phi(x_scf(i),0),&
             lagrange(x_scf(i),-n_range+1,n_range-1,0),&
             sum([(Aiq(q,0)*x_scf(i)**q,q=0,2*n_lagvdm-1)]),&
             sum([(Aiq(q,0)*sum([(phi(x_scf(i),j)*dble(j**q),j=-n_range/2,n_range/2+1)]),q=0,2*n_lagvdm-1)])
        !(sum([(Aiq(q,j)*x_scf(i)**q,q=0,2*n_lagvdm)]),j=-n_lagvdm,n_lagvdm)

     end do
     call yaml_map('Polynomial exactness |p> - sum_j j**p |phi_j>',deviations)

     !Do scalar product with the dual function
     !Initialize the work array needed to integrate with the Lagrange polynomial
     n_lagt=0
     if (itype_scf_dual<0)  n_lagt=-itype_scf_dual

     n_scft = 2*(abs(itype_scf_dual)+nmoms_dual)*npts
     !Allocations for lagrange polynomial data array
     x_scft = f_malloc(0.to.n_scft,id='x_scft')
     scft_dat = f_malloc(0.to.n_scft,id='scf_dat')
     if (n_lagt==0) then
        call ISF_family(itype_scf_dual,nmoms_dual,n_scft,n_ranget,x_scft,scft_dat)
     else
        call lagrange_family(n_lagt,n_scft,n_ranget,x_scft,scft_dat)
     end if

     call yaml_map('itype_scf',itype_scf_dual)
     call yaml_map('nmoms',nmoms_dual)
     call yaml_map('range', n_ranget)
     call yaml_map('number of points',n_scft)
     call yaml_map('step',x_scft(1)-x_scft(0))
     call yaml_map('dx',real(2*(abs(itype_scf_dual)+nmoms_dual),gp)/real(n_scft,gp))

     do j=-n_ranget,n_ranget
        scalar = 0.d0
        do i=0,n_scft
           phival=phi(x_scft(i),j)
           scalar = scalar + scft_dat(i)*phival
        end do
        scalar = scalar*(x_scft(1)-x_scft(0))
        call yaml_map('<phi|phi_{'//trim(adjustl(yaml_toa(j)))//'}>',scalar)
     end do
     call f_free(x_scft,scft_dat)
  end if
 
  call f_free(x_scf,scf_dat)!,f_l)
  call f_free(pol_dat)
  call f_free(pol_sp)
  call f_free(Aiq)


  contains

    pure function phi(x0,i0)
      implicit none
      integer, intent(in) :: i0
      double precision, intent(in) :: x0
      double precision :: phi
      !local variables
      integer :: ival
      double precision :: y
      phi=0.d0
      y=x0-dble(i0)
      ival=minloc(abs(x_scf-y),1)+lbound(x_scf,1)-1
      !print *,'ival',x0,ival,scf_dat(ival),x_scf(ival)
      if (n_lag==0) then
         if (x_scf(ival)-y == 0.d0) phi=scf_dat(ival)
      else
         if (x_scf(ival)-y == 0.d0) phi=lagrange(x0,-n_lag,n_lag,i0)
         !if (y >= real(-2*n_lag) .and. y <= real(2*n_lag)) phi=lagrange(x0,-n_lag,n_lag,i0)
      end if

    end function phi

    !> Calculate the Lagrange polynomial
    !! @f$ l_{i_0}(x) = \prod_{i=istart, i \neq i_0}^{iend} \frac{x-i}{i0-i}$
    elemental function lagrange(x,istart,iend,i0) result(y)
      implicit none
      !Arguments
      real(gp), intent(in) :: x
      integer, intent(in) :: istart !< First point
      integer, intent(in) :: iend   !< Last point
      integer, intent(in) :: i0     !< Center (removed from the product)
      real(gp) :: y
      !Local variables
      integer :: i
      if (i0 >= istart .and. i0<= iend) then
         y = 1.d0
         do i=istart,iend
            if (i /= i0) then
               y = y *(x-real(i,gp))/real(i0-i,gp)
            end if
         end do
      else
         y = 0.d0
      end if
    end function lagrange

    subroutine lagrange_family(m,npts,n_range,x_scf,lag_dat)
      implicit none
      integer, intent(in) :: m
      integer, intent(in) :: npts
      integer, intent(out) :: n_range
      double precision, dimension(0:npts), intent(out) :: x_scf
      double precision, dimension(0:npts), intent(out) :: lag_dat
      !local variables
      integer :: i,ni

      n_range=m

      ni=2*n_range
      do i=0,npts
         x_scf(i) = real(2*i*ni,kind=8)/real(npts,kind=8)-real(ni,kind=8)
         lag_dat(i)=lagrange(x_scf(i),-n_range,n_range,0) !
      end do

    end subroutine lagrange_family
 
end subroutine polynomial_exactness

!>calculate the inverse of the vandermonde matrix and compare it to the lagrange polynomial
!! calculated with the usual formule
subroutine invert_vandermonde(istart,iend,Aiq)
  use dictionaries
  use yaml_output
  implicit none
  integer, intent(in) :: istart,iend
  double precision, dimension(0:iend-istart,istart:iend), intent(out) :: Aiq
  !local variables
  integer :: i,j,q,p,m,info,ishift
  double precision :: delta,deltadev
  integer, dimension(iend-istart+1) :: ipiv
  double precision, dimension(iend-istart+1) :: work
  external :: dgetrf,dgetri

  m=iend-istart+1
  !fill Vdm
  forall (i=istart:iend,q=0:m-1) Aiq(i-istart,q+istart)=dble(i**q)

  !lapack inverse
  call dgetrf(m,m,Aiq,m,ipiv,info)
  if (info /=0) call f_err_throw('Error in dgetrf, info='//trim(yaml_toa(info)))
  !fill the unit element of the LU factorization
!!$  do q=0,m-1
!!$     Aiq(q+istart,q)=1.d0
!!$  end do
  call dgetri(m,Aiq,m,ipiv,work,m,info)
  if (info /=0) call f_err_throw('Error in dgetri, info='//trim(yaml_toa(info)))

  !verify that the inverse is as such (plain dgemm)
  do ishift=0,0!istart-m,iend+m
     deltadev=0.d0
     do i=istart,iend
        do j=istart,iend
           delta=0.d0
           do q=0,m-1
              delta=delta+Aiq(q,j)*dble((i+ishift)**q)
           end do
           !print *,'i,j,delta',i,j,delta
           if (i==j) then
              deltadev=max(deltadev,abs(1.d0-delta))
           else
              deltadev=max(deltadev,abs(delta))
           end if
        end do
     end do
     call yaml_map('Deltadeviation, ishift='//trim(yaml_toa(ishift)),deltadev)
  end do
  deltadev=0.d0
  !the transposed
  do p=0,m-1
     do q=0,m-1
        delta=0.d0
        do i=istart,iend
           delta=delta+Aiq(p,i)*dble(i**q)
        end do
        if (p==q) then
           deltadev=max(deltadev,abs(1.d0-delta))
        else
           deltadev=max(deltadev,abs(delta))
        end if
        !print *,'p,q,delta',p,q,delta
     end do
  end do
  call yaml_map('Deltadeviation',deltadev)

  call yaml_map('Extremes',[istart,iend])
  call yaml_sequence_open('Inverse Vandermonde')
  do i=istart,iend
     call yaml_sequence(yaml_toa(Aiq(0:5,i),fmt='(1pe12.2)'),advance='no')
     call yaml_comment(yaml_toa(i))
  end do
  call yaml_sequence_close()
end subroutine invert_vandermonde


!!$
!!$!> takes a Gaussian of exponent pgauss an center x0 and discretize it on a grid of size 1 in units of sqrt(0.5*[pgauss])
!!$!! f(x)=fac*exp(-pgauss*(x-x0)**2)
!!$subroutine discretize_gaussian(nrange,fac,pgauss,x0,hgrid,filename)
!!$  use module_base
!!$  use gaussians
!!$  !not yet use dynamic_memory !let's start
!!$  implicit none
!!$  integer, intent(in) :: nrange
!!$  real(gp), intent(in) :: fac,pgauss,x0,hgrid
!!$  character(len=*), intent(in) :: filename
!!$  !local variables
!!$  integer, parameter :: itype_scf=16,npoints=2**6,nmoms=6
!!$  integer :: k,n_scf,n_range,j
!!$  real(gp) :: dx,x
!!$  real(gp), dimension(0:nmoms,3) :: moments !<to check all methods
!!$  real(dp), dimension(0:2048) :: fISF
!!$  real(gp), dimension(:), allocatable :: f_i,f_phi_i,f_phi_i2,x_scf,y_scf
!!$
!!$  !allocate the three gaussians arrays
!!$  allocate(f_i(-nrange:nrange))
!!$  allocate(f_phi_i(-nrange:nrange))
!!$  allocate(f_phi_i2(-nrange:nrange))
!!$
!!$  !Number of integration points = 2*itype_scf*n_points
!!$  n_scf=2*itype_scf*npoints
!!$
!!$  !Other allocations
!!$  allocate(x_scf(0:n_scf))
!!$  allocate(y_scf(0:n_scf))
!!$
!!$  !Build the scaling function
!!$  call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
!!$
!!$  !Step grid for the integration
!!$  dx = real(nrange,gp)/real(n_scf,gp)
!!$
!!$  !first, collocate the gaussian
!!$  do k=-nrange,nrange
!!$     f_i(k)=fac*exp(-pgauss*(hgrid*real(k,gp)-x0)**2)
!!$  end do
!!$
!!$  f_phi_i2=0.0_gp
!!$  f_phi_i=0.0_gp
!!$  !then, calculate its values integrated with SCF
!!$  !use the i2 array as a workspace
!!$  call my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,nrange,n_scf,x_scf,y_scf,f_phi_i,f_phi_i2)
!!$  call dscal(2*nrange+1,fac,f_phi_i,1)
!!$  !Only itype=8,14,16,20,24,30,40,50,60,100
!!$  call four_isf_coeffs(itype_scf,fISF)
!!$  f_phi_i2=0.0_gp
!!$  !then calculate the values via the analytic method (for the moment only on the second half)
!!$  !this routine is not useful for performances of precision and speed
!!$  !call my_analytic_integral(hgrid*sqrt(pgauss),x0/hgrid,&
!!$  !     nrange,itype_scf,f_phi_i2,fISF,64)
!!$  !  do k=-nrange,-1
!!$  !     f_phi_i2(k)=f_phi_i2(-k)
!!$  !  end do
!!$  !call dscal(2*nrange+1,fac,f_phi_i2,1)
!!$
!!$  !use daubechies wavelets for expressing the function
!!$  !  call gauss_to_daub(hgrid,fac,x0,sqrt(0.5/pgauss),0,&!no err, errsuc
!!$  !       2*nrange+1,n_left,n_right,c,err_norm,&         !no err_wav. nmax instead of n_intvx
!!$  !       ww,nwork,.false.)                             !added work arrays ww with dimension nwork
!!$  !create the isf array in the gaussian module
!!$  call initialize_real_space_conversion()
!!$  do j=-nrange,nrange
!!$     f_phi_i2(j)=scfdotf(j,hgrid,pgauss,x0,0)
!!$     !print *,'j',j,f_phi_i2(j)
!!$  end do
!!$  !use the elemental property of the scfdotf function
!!$  f_phi_i2=scfdotf((/(j,j=-nrange,nrange)/),hgrid,pgauss,x0,0)
!!$  call finalize_real_space_conversion('discretize_gaussian')
!!$
!!$  !then print the results
!!$  moments=0.0_dp
!!$  !open(unit=200, file = trim(filename), status = 'unknown')
!!$  do j=0,nmoms
!!$     open(unit=200+j)
!!$     do k=-nrange,nrange
!!$        x=real(k,gp)*hgrid-x0
!!$        moments(j,1)=moments(j,1)+x**j*f_i(k)
!!$        moments(j,2)=moments(j,2)+x**j*f_phi_i(k)
!!$        moments(j,3)=moments(j,3)+x**j*f_phi_i2(k)
!!$        write(200+j,'(i4,3(1pe26.17e3))')k,x**j*f_i(k),x**j*f_phi_i(k),x**j*f_phi_i2(k)
!!$     end do
!!$     moments(j,:) = moments(j,:)*hgrid
!!$     print '(a,i3,8(1pe14.5))',trim(filename),j,&
!!$          pgauss,moments(j,:),gauint0(pgauss,j),moments(j,:)-gauint0(pgauss,j)
!!$     close(unit=200+j)
!!$  end do
!!$
!!$  deallocate(f_i,f_phi_i,f_phi_i2)
!!$
!!$end subroutine discretize_gaussian
!!$
!!$
!!$!> Write a routine which performs the integration of a function on a given grid
!!$subroutine my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: n_range,n_scf,itype_scf
!!$  real(dp), intent(in) :: pgauss,hgrid,dx,x0
!!$  real(dp), dimension(0:n_scf), intent(in) :: x_scf
!!$  real(dp), dimension(0:n_scf), intent(in) :: y_scf
!!$  real(dp), dimension(-n_range:n_range), intent(inout) :: work
!!$  real(dp), dimension(-n_range:n_range), intent(inout) :: kernel_scf
!!$  !local variables
!!$  real(dp), parameter :: p0_ref = 1.0_dp
!!$  integer :: n_iter,i_kern,i
!!$  real(dp) :: p0_cell,p0gauss,absci,kern
!!$
!!$  !Step grid for the integration
!!$  !dx = real(n_range,dp)/real(n_scf,dp)
!!$
!!$  !To have a correct integration
!!$  p0_cell = p0_ref/(hgrid*hgrid)
!!$
!!$  !write(*,*) 'p0_cell = ', p0_cell 
!!$
!!$  !We calculate the number of iterations to go from pgauss to p0_ref
!!$  n_iter = nint((log(pgauss) - log(p0_cell))/log(4.0_dp))
!!$
!!$  write(*,*) 'n_iter = ', n_iter
!!$
!!$  if (n_iter <= 0 .or. x0 /= 0.d0)then
!!$     n_iter = 0
!!$     p0gauss = pgauss
!!$  else
!!$     p0gauss = pgauss/4._dp**n_iter
!!$  end if
!!$
!!$  !  write(*,*) 'p0gauss =',  p0gauss
!!$
!!$  !Stupid integration
!!$  !Do the integration with the exponential centered in i_kern
!!$  kernel_scf(:) = 0.0_dp
!!$  do i_kern=-n_range,n_range
!!$     kern = 0.0_dp
!!$     do i=0,n_scf
!!$        absci = x_scf(i) - real(i_kern,dp)
!!$        absci = x0-absci*hgrid !sign has changed here
!!$        absci = absci*absci
!!$        kern = kern + y_scf(i)*dexp(-p0gauss*absci)
!!$     end do
!!$     kernel_scf(-i_kern) = kern*dx
!!$     !write(17,*)i_kern,kernel_scf(i_kern)
!!$     !kernel_scf(-i_kern) = kern*dx
!!$     !if (abs(kern) < 1.d-18) then
!!$     !   !Too small not useful to calculate
!!$     !   exit
!!$     !end if
!!$  end do
!!$
!!$  !  do i=0,n_scf
!!$  !     write(18,*)i,x_scf(i),y_scf(i)
!!$  !  end do
!!$
!!$  !Start the iteration to go from p0gauss to pgauss
!!$  call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,work)
!!$
!!$END SUBROUTINE my_gauss_conv_scf
!!$
!!$
!!$!> Here alpha corresponds to sqrt(alpha) in mathematica
!!$!! the final result is fwork(j+m)-fwork(j-m)
!!$subroutine my_analytic_integral(alpha,x0,ntot,m,fwork,fISF,argument_nf)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: ntot,m
!!$  real(dp), intent(in) :: alpha,x0 !<x0 is the deviation wrt the grid spacing
!!$  real(dp), dimension(-ntot:ntot), intent(inout) :: fwork
!!$  real(dp), dimension(0:2048), intent(in) :: fISF
!!$  integer, intent(in) :: argument_nf 
!!$
!!$  !local variables
!!$  integer :: nf
!!$  real(dp), parameter :: pi=3.1415926535897932384_dp
!!$  logical :: flag,flag1,flag2
!!$  integer :: j,q,jz
!!$  real(dp) :: if,r1,r2,res,ypm,ymm,erfcpm,erfcmm,factor,re,ro,factorend
!!$
!!$  !write(*,*) fISF(1000), fISF16(1000)
!!$
!!$  !  if(present(argument_nf)) then 
!!$  nf=argument_nf
!!$  !  else
!!$  !     nf = 64 ! "default value"
!!$  !  endif
!!$
!!$  flag=.false.
!!$  factor=pi/real(2*m,dp)/alpha
!!$  factorend=sqrt(pi)/alpha/real(4*m,dp)*exp((alpha*x0)**2)
!!$  !fill work array
!!$  !the calculation for j=0 can be separated from the rest
!!$  !since it only contains error functions
!!$  !the values of the loop are not symmetric in principle
!!$  loop_nonzero: do j=-ntot,ntot
!!$     ypm=alpha*(real(j+m,dp)+x0)
!!$     ymm=alpha*(real(j-m,dp)+x0)
!!$     call derfcf(erfcpm,ypm) ! erfcpm = erfc(ypm)
!!$     call derfcf(erfcmm,ymm)
!!$
!!$     !assume nf even
!!$     res=0._dp
!!$     !reso=0._dp
!!$     do q=nf,2,-2
!!$        !the sign of q only influences the imaginary part
!!$        !so we can multiply by a factor of two
!!$
       ! call wofz_mod(alpha,m,q,-j-m,r1,if,flag1)
       ! call wofz_mod(alpha,m,q,-j+m,r2,if,flag2)
!!$        call GcplxInt(alpha,m,q,-j-m,-x0,r1,if,flag1)
       ! if (q==12 .and. j==15) then
       !    print *,'here',alpha,m,q,j,x0,r1,if,flag1
       !    stop
       ! end if
!!$        call GcplxInt(alpha,m,q,-j+m,-x0,r2,if,flag2)
!!$        re=r1-r2
!!$        flag=flag1 .or. flag2 .or. flag
!!$        !if (flag) then
!!$        !   print *,'here',r1,if,q,j,x0
!!$        !   stop 
!!$        !end if
       ! call wofz_mod(alpha,m,q-1,-j-m,r1,if,flag1)
       ! call wofz_mod(alpha,m,q-1,-j+m,r2,if,flag2)
!!$        call GcplxInt(alpha,m,q-1,-j-m,-x0,r1,if,flag1)
!!$        call GcplxInt(alpha,m,q-1,-j+m,-x0,r2,if,flag2)
!!$        ro=r1-r2
!!$        flag=flag1 .or. flag2 .or. flag
!!$        !if (flag) then
!!$        !   print *,'there',xo,y
!!$        !   stop 
!!$        !end if
!!$        !write(16,'(2(i4),6(1pe15.7))')j,q,re,ro,erfcmm-erfcpm
!!$        re=re*fISF(q)
!!$        ro=ro*fISF(q-1)
!!$        res=res+re-ro
!!$     end do
!!$     !q=0 
!!$     !fwork(j)=derf(y)+rese-reso
!!$     fwork(j)=erfcmm-erfcpm+2.0_dp*res!e-reso
!!$     fwork(j)=factorend*fwork(j)
!!$     !exit form the loop if it is close to zero
!!$     !     if (abs(fwork(j)) < 1.e-25_dp) exit loop_nonzero
!!$     !write(17,'(i4,8(1pe15.7))')j,derf(y),erfcsgn,rese,reso,derf(y)+rese-reso,&
!!$     !     -erfcsgn+rese-reso,erfcsgn+rese-reso
!!$  end do loop_nonzero
!!$
!!$  !check flag
!!$  if (flag) then
!!$     write (*,*)'value of alpha',alpha
!!$     stop 'problem occurred in wofz'
!!$  end if
!!$
!!$  !put to zero the rest
!!$  do jz=j+1,ntot
!!$     fwork(jz)=0.0_dp
!!$  end do
!!$
!!$END SUBROUTINE my_analytic_integral

!> determine the valence of any of the atoms and the corresponding symbol
subroutine bacasable_valence()
  use ao_inguess
  use module_base
  use yaml_output
  implicit none
  integer :: izatom,ival,ierr,stderr,stdout
  integer :: nzatom, nelpsp, npspcode,ixc
  logical :: exists
  character(len=2) :: symbol
  character(len=256) :: msg
  type(dictionary), pointer :: semicores
  real(gp), dimension(0:4,0:6) :: psppar

  semicores=>list_new(.item. ["Ru","Rh","Pd","In","Ir","Pt","Au","Tl"])

  stderr=f_get_free_unit(17)
  !open error queue, leave it as default queue
  call yaml_set_stream(unit=stderr,filename='errors') 
  stdout=6
  call yaml_set_stream(unit=stdout,setdefault=.false.,tabbing=0)
  ixc=-101130
  izatom=1
  call f_err_open_try()
  do while(izatom <= 86)
     ival=0
     find_symbol: do while (ival <= 30)
        !see if an atom exists with this value
        call atomic_info(izatom,ival,symbol=symbol)
        if (f_err_check()) then
           ierr=f_err_pop(add_msg=msg)
           call yaml_map('Error for'//trim(yaml_toa([izatom,ival])),msg)
        else
           !call yaml_map(symbol,[izatom,ival],unit=stdout)
           exit find_symbol
        end if
        ival=ival+1
     end do find_symbol
     if (ival /= 31) then
        !the symbol has been found therefore we can inspect the number of electrons
        if (symbol .in. semicores) then
           call f_strcpy(src=trim(symbol)//'_sc',dest=msg)
        else
           call f_strcpy(src=trim(symbol),dest=msg)
        end if
        call psp_from_data(trim(msg), nzatom, nelpsp, npspcode, ixc, psppar, exists)
        if (nzatom /= 0) call yaml_map(trim(msg),[nzatom, nelpsp, npspcode],unit=stdout)
     end if
     izatom=izatom+1
  end do
  call f_err_close_try()

  call yaml_close_stream(unit=stderr)

  call dict_free(semicores)
  
end subroutine bacasable_valence
