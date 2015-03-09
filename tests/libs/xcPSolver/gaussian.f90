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
  call f_lib_initialize()

  pow=0
  unit=iunit+1

  !pgauss=0.5_gp/((0.1_gp*hgrid)**2)!8.0e-3_dp*1.25_dp**(6*(8-1))
  !array where we have to write the value of the discretization
  fj_phi=f_malloc(-npts .to. npts,id='fj_phi')
  fj_coll=f_malloc(-npts .to. npts,id='fj_coll')
  call polynomial_exactness(npts,8,0,itype_scf,8,8,iplot)
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
     call yaml_map('Results',reshape(avgmaxmin,(/6,nmoms+1/)),fmt='(1pe14.5)')
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

  unit=iunit+istep
  call f_open_file(unit,'gau'//trim(adjustl(yaml_toa(x0,fmt='(f5.2)')))//&
       'p'//trim(adjustl(yaml_toa(pgauss,fmt='(f5.2)')))//&
       'h'//trim(adjustl(yaml_toa(hgrid,fmt='(f5.2)'))))

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
    integer :: j,k
    real(gp) :: x

  moments(0)=h*sum(array)
  do j=1,nmoms
       moments(j)=0.0_gp
       do k=1,n
          x=real(k,gp)*h-x0
          moments(j)=moments(j)+x**j*array(k)
       end do
       moments(j)=moments(j)*h
    end do

  end subroutine moments_1d


end program MP_gaussian


!> Verify the property x**p = sum_j x_j**p \phi_j(x) for different families
!! of interpolating functions
subroutine polynomial_exactness(npts,itype_scf,nmoms,p,itype_scf_dual,nmoms_dual,iplot)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: npts,itype_scf,nmoms,p,iplot,itype_scf_dual,nmoms_dual
  !local variables
  integer :: n_scf,j,q,i,n_range,n_ext,n_ranget,n_scft
  real(gp) :: phival,scalar
  real(gp), dimension(0:p) :: deviations
  real(gp), dimension(:), allocatable :: scf_dat,x_scf,x_scft,scft_dat
  real(gp), dimension(:,:), allocatable :: pol_dat,pol_sp

  !Initialize the work array needed to integrate with the Lagrange polynomial
  n_scf = 2*(itype_scf+nmoms)*npts

  !Allocations for lagrange polynomial data array
  x_scf = f_malloc(0.to.n_scf)
  scf_dat = f_malloc(0.to.n_scf,id='scf_dat')
  pol_dat = f_malloc0([0.to.n_scf,0.to.p],id='pol_dat')
  

  !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
  !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_dat)
  call ISF_family(itype_scf,nmoms,n_scf,n_range,x_scf,scf_dat)

  pol_sp = f_malloc0([-n_range.to.n_range,0.to.p],id='pol_sp')
 
  !call yaml_set_stream(record_length=150)
  call yaml_map('itype_scf',itype_scf)
  call yaml_map('nmoms',nmoms)
  call yaml_map('range', n_range)
  call yaml_map('number of points',n_scf)
  call yaml_map('step',x_scf(1)-x_scf(0))
  call yaml_map('dx',real(2*(itype_scf+nmoms),gp)/real(n_scf,gp))

  n_ext=n_range/2
  !Step grid for the integration
  !dx = real(2*itype_scf,gp)/real(n_scf,gp)
  !starting point for the x coordinate for integration
  !x  = real(-itype_scf,gp)-dx
  write(iplot,'(a)') '#Abscissae   Interpolating_scaling_function   Lagrange_polynomial' 
  do q=0,p
     !evaluate the maxima
     do j=-n_range,n_range
        do i=0,n_scf
           phival=phi(x_scf(i)-real(j,kind=8))
           pol_dat(i,q)=pol_dat(i,q)+(real(j,kind=8)/real(n_ext,kind=8))**q*phival
           pol_sp(j,q)=pol_sp(j,q)+((x_scf(i)+real(j,kind=8))/real(n_ext,kind=8))**q*scf_dat(i)
        end do
        pol_sp(j,q)=pol_sp(j,q)*(x_scf(1)-x_scf(0))
        pol_sp(j,q)=abs(pol_sp(j,q)-(real(j,kind=8)/real(n_ext,kind=8))**q)
     end do
  end do
  call yaml_map('Polynomial scalprods',[(maxval(pol_sp(:,q)),q=0,p)])
  !evaluate the maxima
  deviations=0.d0
  do i=0,n_scf
     !x=x+dx
     !lag_dat(i) = lag_sym(x_scf(i),itype_scf,0)
     forall (q=0:p) deviations(q)=max(deviations(q),abs((x_scf(i)/real(n_ext,kind=8))**q-pol_dat(i,q)))
     write(iplot,'(f19.13,1x,30(1pe23.15,1x))') x_scf(i),scf_dat(i),pol_dat(i,:)
  end do
  call yaml_map('Polynomial exactness',deviations)

  !Do scalar product with the dual function
  !Initialize the work array needed to integrate with the Lagrange polynomial
  n_scft = 2*(itype_scf_dual+nmoms_dual)*npts
  !Allocations for lagrange polynomial data array
  x_scft = f_malloc(0.to.n_scft)
  scft_dat = f_malloc(0.to.n_scft,id='scf_dat')
  call ISF_family(itype_scf_dual,nmoms_dual,n_scft,n_ranget,x_scft,scft_dat)

  call yaml_map('itype_scf',itype_scf_dual)
  call yaml_map('nmoms',nmoms_dual)
  call yaml_map('range', n_ranget)
  call yaml_map('number of points',n_scft)
  call yaml_map('step',x_scft(1)-x_scft(0))
  call yaml_map('dx',real(2*(itype_scf_dual+nmoms_dual),gp)/real(n_scft,gp))


  do j=-n_ranget,n_ranget
     scalar = 0.d0
     do i=0,n_scft
        phival=phi(x_scft(i)-real(j,kind=8))
        scalar = scalar + scft_dat(i)*phival
     end do
     scalar = scalar*(x_scft(1)-x_scft(0))
     call yaml_map('<phi|phi_{'//trim(adjustl(yaml_toa(j)))//'}>',scalar)
  end do
!!$
!!$  call yaml_sequence_open('<lag|iscf>')
!!$  istart=-(itype_scf/2-1)
!!$  iend=itype_scf/2-1
!!$  do i0=-itype_scf/2+1,itype_scf/2-1
!!$     scalar = 0.d0
!!$     do i=0,n_scf
!!$        scalar = scalar + scf_dat(i)*lagrange(x_scf(i),istart,iend,i0)
!!$        !scalar = scalar + scf_dat(i)*lag_sym(x_scf(i),itype_scf,i0)
!!$        !scalar = scalar + scf_dat(i)*x_scf(i)**i0
!!$     end do
!!$     scalar = scalar*(x_scf(1)-x_scf(0))
!!$     !call yaml_map(trim(yaml_toa(istart))//trim(yaml_toa(iend))//trim(yaml_toa(i0)), scalar)
!!$     !call yaml_map(trim(yaml_toa(i0)), (/ scalar,lag_sym(real(i0,gp),itype_scf,i0) /))
!!$     call yaml_map(trim(yaml_toa(i0)), (/ scalar,lagrange(real(i0,gp),istart,iend,i0) /))
!!$  end do
!!$  call yaml_sequence_close()
 
  call f_free(x_scf,scf_dat)
  call f_free(pol_dat)
  call f_free(pol_sp)

  contains

    pure function phi(x0)
      implicit none
      double precision, intent(in) :: x0
      double precision :: phi
      !local variables
      integer :: ival
      phi=0.d0
      ival=minloc(abs(x_scf-x0),1)+lbound(x_scf,1)-1
      !print *,'ival',x0,ival,scf_dat(ival),x_scf(ival)
      if (x_scf(ival) -x0 == 0.d0)  phi=scf_dat(ival)
      
    end function phi

end subroutine polynomial_exactness
