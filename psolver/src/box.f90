!> @file
!!    Modulefile for handling fundamental data structed and methods of the simulation box
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module box

  use PSbase

  private

  !>parameter for the definition of the bc
  integer, parameter :: FREE=0
  integer, parameter :: PERIODIC=1

  type, public :: cell
     logical :: orthorhombic !<true if the cell is orthorhombic
     integer, dimension(3) :: bc
     integer, dimension(3) :: ndims
     real(gp), dimension(3) :: hgrids
     real(gp), dimension(3) :: angrad !<angles between the dimensions in radiant
     !derived data
     real(gp) :: volume_element
     real(gp), dimension(3,3) :: habc !<primitive volume elements in the translation vectors direction
  end type cell

  public :: cell_r,cell_periodic_dims,minimum_distance,closest_r,square,cell_new

contains

  pure function cell_new(geocode,ndims,hgrids,angrad) result(mesh)
    use numerics, only: onehalf,pi
    use wrapper_linalg, only: det_3x3
    implicit none
    character(len=1), intent(in) :: geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    real(gp), dimension(3), intent(in), optional :: angrad
    type(cell) :: mesh
    !local variables
    real(gp) :: aa,cc,a2,cosang

    mesh%bc=FREE
    if (geocode /= 'F') mesh%bc(1)=PERIODIC
    if (geocode == 'P') mesh%bc(2)=PERIODIC
    if (geocode /= 'F') mesh%bc(3)=PERIODIC
    mesh%ndims=ndims
    mesh%hgrids=hgrids
    if (present(angrad)) then
       mesh%angrad=angrad
       !some consistency check on the angles should be performed
       !1) sum(angrad) < twopi
       !2) all(angrad) > 0.0_gp
       if (all(angrad==angrad(1)) .and. angrad(1) /= onehalf*pi) then
          !Treat the case of equal angles (except all right angles) :
          !generates trigonal symmetry wrt third axis
          cosang=cos(angrad(1))
          a2=2.0_gp/3.0_gp*(1.0_gp-cosang)
          aa=sqrt(a2)
          cc=sqrt(1.0_gp-a2)
          mesh%habc(1,1)=aa; mesh%habc(2,1)=0.0_gp; mesh%habc(3,1)=cc
          mesh%habc(1,2)=-0.5_gp*aa ; mesh%habc(2,2)=sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,2)=cc
          mesh%habc(1,3)=-0.5_gp*aa ; mesh%habc(2,3)=-sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,3)=cc
       else
          mesh%habc(:,:)=0.0_gp
          mesh%habc(1,1)=1.0_gp
          mesh%habc(1,2)=cos(angrad(3))
          mesh%habc(2,2)=sin(angrad(3))
          mesh%habc(1,3)=cos(angrad(2))
          mesh%habc(2,3)=(cos(angrad(1))-mesh%habc(1,2)*mesh%habc(1,3))/mesh%habc(2,2)
          mesh%habc(3,3)=sqrt(1.0_gp-mesh%habc(1,3)**2-mesh%habc(2,3)**2)
       end if
       !Rescale habc using hgrid
       mesh%habc(:,1)=hgrids*mesh%habc(:,1)
       mesh%habc(:,2)=hgrids*mesh%habc(:,2)
       mesh%habc(:,3)=hgrids*mesh%habc(:,3)
       !the volume element
       !Compute unit cell volume
       mesh%volume_element=det_3x3(mesh%habc)
    else
       mesh%angrad=onehalf*pi
       mesh%volume_element=product(mesh%hgrids)
    end if
    mesh%orthorhombic=all(mesh%angrad==onehalf*pi)
  end function cell_new

  !> returns a logical array of size 3 which is .true. for all the periodic dimensions
  pure function cell_periodic_dims(mesh) result(peri)
    implicit none
    type(cell), intent(in) :: mesh
    logical, dimension(3) :: peri
    !local variables

    peri= mesh%bc == PERIODIC

  end function cell_periodic_dims

  !>gives the value of the coordinate from the grid point
  elemental pure function cell_r(mesh,i,dim) result(t)
    implicit none
    integer, intent(in) :: i
    type(cell), intent(in) :: mesh
    integer, intent(in) :: dim
    real(gp) :: t

    t=mesh%hgrids(dim)*(i-1)
  end function cell_r

  function minimum_distance(mesh,v1,v2) result(d)
    use dictionaries, only: f_err_throw
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(cell), intent(in) :: mesh
    real(gp) :: d
    !local variables
    integer :: i
    real(gp) :: d2

    if (mesh%orthorhombic) then
       d2=0.0_gp
       do i=1,3
          d2=d2+min_dist(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               v1(i),v2(i))**2
       end do
       d=sqrt(d2)
    else
       call f_err_throw('Minimum distance not yet implemented for nonorthorhombic cells')
    end if

  end function minimum_distance

  !> Calculates the minimum difference between two coordinates
  pure function min_dist(bc,alat,r,r_old)
    implicit none
    integer, intent(in) :: bc
    real(gp), intent(in) :: r,r_old,alat
    real(gp) :: min_dist

    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
    min_dist=abs(r-r_old)
    if (bc==PERIODIC) then
       if (min_dist > 0.5_gp*alat) then
          if (r < 0.5_gp*alat) then
             min_dist=abs(r+alat-r_old)
          else
             min_dist=abs(r-alat-r_old)
          end if
       end if
    end if

  end function min_dist

  !> Calculates the minimum difference between two coordinates
  pure function r_wrap(bc,alat,r,c)
    implicit none
    integer, intent(in) :: bc
    real(gp), intent(in) :: r,c,alat
    real(gp) :: r_wrap

    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
    r_wrap=r-c
    if (bc==PERIODIC) then
       if (abs(r_wrap) > 0.5_gp*alat) then
          if (r < 0.5_gp*alat) then
             r_wrap=r+alat-c
          else
             r_wrap=r-alat-c
          end if
       end if
    end if

  end function r_wrap

  !>find the closest center according to the periodiciy of the
  !! box and provide the vector
  pure function closest_r(mesh,v,center) result(r)
    implicit none
    real(gp), dimension(3), intent(in) :: v,center
    type(cell), intent(in) :: mesh
    real(gp), dimension(3) :: r
    !local variables
    integer :: i

    if (mesh%orthorhombic) then
       do i=1,3
          r(i)=r_wrap(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               v(i),center(i))
       end do
    end if

  end function closest_r

  pure function square(mesh,v)
    implicit none
    real(gp), dimension(3), intent(in) :: v
    type(cell), intent(in) :: mesh
    real(gp) :: square

    if (mesh%orthorhombic) then
       square=v(1)**2+v(2)**2+v(3)**2
    end if

  end function square

end module box
