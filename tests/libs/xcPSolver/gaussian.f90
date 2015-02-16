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

    write(iunit+istep,'(a)') '#Projection of a gaussian with iscf, collocation method'
    write(iunit+istep,'(a)') '#j,fj_phi(j),fj_coll(j)'
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
