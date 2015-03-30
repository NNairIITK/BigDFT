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
  use gaussians, only: initialize_real_space_conversion, finalize_real_space_conversion,&
      &     mp_exp, gauint0,scfdotf
  use yaml_output
  implicit none
  integer, parameter :: itype_scf = 8 
  integer, parameter :: nmoms=16           !< Number of calculated moments
  integer, parameter :: nstep=1            !< Number of resolution to calculate the moments
  integer, parameter :: nsigma=1           !< Number of different gaussian functions
  integer, parameter :: npts=32            !< Arrays from -npts to npts
  real(gp), parameter :: hgrid = 1.0_gp    !< Step grid
  integer, parameter :: iplot=14,iunit=16  !< File unit for the plots
  integer :: j,imoms,pow,istep,isigma,unit
  real(gp) :: pgauss,x0,reference,max_phi
  real(gp), dimension(0:nmoms,2) :: moments
  real(gp), dimension(3,2,0:nmoms) :: avgmaxmin
  real(gp), dimension(:), allocatable :: fj_phi,fj_coll
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
  call f_lib_initialize()

  call bacasable_valence()

  call f_lib_finalize()
  stop

  pow=0
  unit=iunit+1
  untplot=iplot

  !pgauss=0.5_gp/((0.1_gp*hgrid)**2)!8.0e-3_dp*1.25_dp**(6*(8-1))
  !array where we have to write the value of the discretization
  fj_phi=f_malloc(-npts .to. npts,id='fj_phi')
  fj_coll=f_malloc(-npts .to. npts,id='fj_coll')

  call f_open_file(untplot,'Families.dat')
  call polynomial_exactness(npts,16,0,itype_scf,16,16,untplot,0,moments)
  call f_close(untplot)
  !stop

  !replot the smoothed magic filtered Daubechies
  call f_open_file(untplot,'daub_smoothed.dat')
  call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,8,MFdat)
  call f_close(untplot)
  !stop


  !initialize the work arrays needed to integrate with isf
  call initialize_real_space_conversion(npoints=2**8,isf_m=itype_scf,nmoms=6)

  call f_open_file(unit,file='MultipolesError.dat')

  call yaml_map('itype_scf',itype_scf)
  call yaml_map('nmoms',nmoms)

  !Evaluate the moments of the isf
  do imoms=0,nmoms
    moments(imoms,1) = scfdotf(0,hgrid,0.d0,0.d0,imoms)
  end do
  call yaml_map('Isf moments',moments(:,1))

  ! Calculate for different nsigma sigma
  do isigma=1,nsigma
     pgauss=0.5_gp/((0.7_gp+0.1_gp*(isigma-1)*hgrid)**2)
     call yaml_map('sigma/h',sqrt(0.5_gp/pgauss)/hgrid)

     avgmaxmin=0.d0
     avgmaxmin(3,:,:)=1.d100
     max_phi=0.0_gp

     do istep=1,nstep
        x0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp))*hgrid
        call evaluate_moments(nmoms,npts,hgrid,pgauss,pow,x0,fj_phi,fj_coll,moments)

        call filename_test(x0,pgauss,hgrid,filename)
        call f_open_file(untplot,trim(filename)//'coll.dat')
        call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,npts,fj_coll)
        call f_close(untplot)

        call filename_test(x0,pgauss,hgrid,filename)
        call f_open_file(untplot,trim(filename)//'isf.dat')
        call polynomial_exactness(npts,16,0,itype_scf,-8,0,untplot,npts,fj_phi)
        call f_close(untplot)


        max_phi=max(max_phi,maxval(abs(fj_coll-fj_phi)))

        call yaml_mapping_open('Normalised moments')
        call yaml_map('x0/h',x0/hgrid)
        call yaml_map('C  moments',moments(:,2),fmt='(1pe14.5)')
        call yaml_map('MP moments',moments(:,1),fmt='(1pe14.5)')
        call yaml_mapping_close()
        do j=1,2
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
     end do
     call yaml_map('Results',reshape(avgmaxmin,[6,nmoms+1]),fmt='(1pe14.5)')
     !Plot fort.(iunit+1)
     avgmaxmin(2,:,:)=avgmaxmin(2,:,:)-avgmaxmin(3,:,:)
     write(unit,'(104(1pe14.5))') sqrt(0.5_gp/pgauss)/hgrid,avgmaxmin(1:2,:,:)
     call yaml_map('maxdiff' // trim(yaml_toa(isigma)), (/ sqrt(0.5_gp/pgauss)/hgrid, max_phi /) )
  end do

  call f_close(unit)

  call finalize_real_space_conversion()

  call f_free(fj_phi,fj_coll)
  call f_lib_finalize()

contains


  !> Classify the quality of a multipole extraction in both cases
  subroutine evaluate_moments(nmoms,npts,hgrid,pgauss,pow,x0,fj_phi,fj_coll,moments)
    use module_base, only: gp,safe_exp,f_open_file,f_close,yaml_toa
    use gaussians, only: mp_exp, scfdotf
    implicit none
    !Arguments
    integer, intent(in) :: npts,pow,nmoms
    real(gp), intent(in) :: hgrid,pgauss,x0
    real(gp), dimension(0:nmoms,2), intent(out) :: moments
    real(gp), dimension(-npts:npts), intent(out) :: fj_phi,fj_coll
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
    call moments_1d(2*npts+1,fj_phi,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,1))

    !collocation array
    fj_coll=mp_exp(hgrid,x0,pgauss,(/(j,j=-npts,npts)/),pow,.false.)

    call moments_1d(2*npts+1,fj_coll,x0+hgrid*(npts+1),hgrid,nmoms,moments(0,2))

    write(unit,'(a)') '#Projection of a gaussian with iscf, collocation and lagrange polynomials method'
    write(unit,'(a)') '#j,fj_phi(j),fj_coll(j),fj_lag(j)'
    do j=-npts,npts
       write(unit,*) j,fj_phi(j),fj_coll(j)
    end do
    istep = istep + 1

    call f_close(unit)

  end subroutine evaluate_moments

  !>get the test filename
  subroutine filename_test(x0,pgauss,hgrid,filename)
    use yaml_strings
    implicit none
    double precision, intent(in) :: x0,pgauss,hgrid
    character(len=128), intent(out) :: filename

    call f_strcpy(dest=filename,src=&
         'gau'//trim(adjustl(yaml_toa(x0,fmt='(f5.2)')))//&
         'p'//trim(adjustl(yaml_toa(pgauss,fmt='(f5.2)')))//&
         'h'//trim(adjustl(yaml_toa(hgrid,fmt='(f5.2)'))))
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
    integer :: j,k,nh,kh,nc,kl
    real(gp) :: x

    moments(0)=h*sum(array)
    do j=1,nmoms
       moments(j)=0.0_gp
       !should start from the middle, but start from the exterior
       kl=1
          moments(j)=moments(j)+x**j*array(k)
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


end program

!> interpolate a function strarting from its coefficient in a given basis set
subroutine interpolated_function(ncoeff,fi,itype_scf,nmoms)
  implicit none
  integer, intent(in) :: ncoeff, itype_scf, nmoms
  double precision, dimension(-ncoeff:ncoeff), intent(in) :: fi

end subroutine interpolated_function


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
      integer :: ival,nrange
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