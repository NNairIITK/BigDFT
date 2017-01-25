!> @file
!!  Test of some functionalities of the numeric groups
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program numeric_check
  use futile
  use f_harmonics
  character(len=*), parameter :: input1=&
       "  {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the array for multipoles,"//&
       "  help_dict: {Allowed values: integer}}"
  character(len=*), parameter :: input2=&
       "  {name: boldify, shortname: b, default: None,"//&
       "  help_string: Boldify the string as a test,"//&
       "  help_dict: {Allowed values: string scalar}}"
  character(len=*), parameter :: input3=&
       "  {name: blinkify, shortname: l, default: None,"//&
       "  help_string: Make the string blinking,"//&
       "  help_dict: {Allowed values: string scalar}}"

  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2//f_cr//&
       '-'//input3


  integer :: n,i
  type(f_multipoles) :: mp
  type(dictionary), pointer :: options
  real(f_double), dimension(3) :: rxyz
  real(f_double), dimension(:), allocatable :: density
  
  call f_lib_initialize()
  call yaml_new_document()
  call yaml_argparse(options,inputs)
  n=options//'ndim'

  density=f_malloc(n,id='density')
  call f_random_number(density)

  rxyz=1.0_f_double

  !create random multipoles
  call f_multipoles_create(mp,2)

  do i=1,n
     call f_multipoles_accumulate(mp%Q,mp%lmax,rxyz,density(i))
  end do

  !here we may print the results of the multipole calculations
  call yaml_mapping_open('Multipoles of the array')
  call yaml_map('q0',sum(density))
  call yaml_mapping_close()

  call yaml_mapping_open('Calculated multipoles')
  call yaml_map('q0',mp%Q(0)%ptr)
  call yaml_map('q1',mp%Q(1)%ptr)
  call yaml_map('q2',mp%Q(2)%ptr)
  call yaml_mapping_close()

  call f_multipoles_free(mp)

  call f_free(density)
  call dict_free(options)

!!$  !test of the multipole preserving routine
!!$  !initialize the work arrays needed to integrate with isf
!!$  !names of the routines to be redefined
!!$  call initialize_real_space_conversion(isf_m=mp_isf_order)
!!$
!!$  boxit = box_iter(mesh,origin=rxyz,cutoff=cutoff)
!!$  call finalize_real_space_conversion()

  call f_lib_finalize()

end program numeric_check

!> Creates the charge density of a Gaussian function, to be used for the local part
!! of the pseudopotentials (gives the error function term when later processed by the Poisson solver).
subroutine gaussian_density(rxyz,rloc, zion, multipole_preserving, use_iterator,boxit,&
     mp_isf,nmpx, nmpy, nmpz, mpx, mpy, mpz, nrho, density)
  !use gaussians, only: mp_exp
  use box
  use f_precisions, only: dp=>f_double
  use multipole_preserving
  use dynamic_memory
  use numerics, only: pi
  use dictionaries, only: f_err_throw
  implicit none
  ! Calling arguments
  logical,intent(in) :: multipole_preserving, use_iterator
  integer,intent(in) :: nrho
  real(dp),intent(in) :: rloc
  integer,intent(in) :: zion !< ionic charge (integer!)
  integer,intent(in) :: mp_isf !< interpolating scaling function order for the multipole preserving
  integer,intent(in) :: nmpx, nmpy, nmpz !< sizes of the temporary arrays; if too small the code stops
  type(box_iterator), intent(inout) :: boxit
  real(dp), dimension(3) :: rxyz
  real(kind=8),dimension(0:nmpx),intent(inout) :: mpx !< temporary array for the exponetials in x direction
  real(kind=8),dimension(0:nmpy),intent(inout) :: mpy !< temporary array for the exponetials in y direction
  real(kind=8),dimension(0:nmpz),intent(inout) :: mpz !< temporary array for the exponetials in z direction
  real(dp),dimension(nrho),intent(inout) :: density
  ! Local variables
  real(dp),parameter :: mp_tiny = 1.e-30_dp
  logical :: perx, pery, perz,gox, goy, goz
  logical, dimension(3) :: peri
  integer :: i3s, n3pi,n1i,n2i,n3i
  real(dp) :: rlocinv2sq, charge, cutoff, xp, yp, zp, rx, ry, rz, hxh, hyh, hzh,fx,fy,fz
  integer :: i1, i2, i3, isx, iex, isy, iey, isz, iez, j1, j2, j3, ind

  call f_routine(id='gaussian_density')

  rx=boxit%oxyz(1)
  ry=boxit%oxyz(2)
  rz=boxit%oxyz(3)

  i3s=boxit%i3s
  n3pi=boxit%i3e-i3s+1

  hxh=boxit%mesh%hgrids(1)
  hyh=boxit%mesh%hgrids(2)
  hzh=boxit%mesh%hgrids(3)

  peri=cell_periodic_dims(boxit%mesh)
  perx=peri(1)
  pery=peri(2)
  perz=peri(3)

  n1i=boxit%mesh%ndims(1)
  n2i=boxit%mesh%ndims(2)
  n3i=boxit%mesh%ndims(3)

  rlocinv2sq=0.5_dp/rloc**2
  charge=real(zion,dp)/(2.0_dp*pi*sqrt(2.0_dp*pi)*rloc**3)

  !cutoff of the range
  cutoff=10.0_dp*rloc
  if (multipole_preserving) then
     !We want to have a good accuracy of the last point rloc*10
     !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=dp)
     cutoff=cutoff+max(hxh,hyh,hzh)*real(mp_isf,kind=dp)
  end if

  isx=boxit%nbox(1,1)
  isy=boxit%nbox(1,2)
  isz=boxit%nbox(1,3)
  iex=boxit%nbox(2,1)
  iey=boxit%nbox(2,2)
  iez=boxit%nbox(2,3)

  ! Check whether the temporary arrays are large enough
  if (iex-isx>nmpx) then
     call f_err_throw('Temporary array in x direction too small')
  end if
  if (iey-isy>nmpy) then
     call f_err_throw('Temporary array in y direction too small')
  end if
  if (iez-isz>nmpz) then
     call f_err_throw('Temporary array in z direction too small')
  end if

  do i1=isx,iex
     mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,multipole_preserving)
  end do
  do i2=isy,iey
     mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,multipole_preserving)
  end do
  do i3=isz,iez
     mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,multipole_preserving)
  end do

  if (use_iterator) then
     do while(box_next_z(boxit))
        fz=mpz(boxit%ibox(3)-boxit%nbox(1,3))
        do while(box_next_y(boxit))
           fy=mpy(boxit%ibox(2)-boxit%nbox(1,2))
           do while(box_next_x(boxit))
              fx=mpx(boxit%ibox(1)-boxit%nbox(1,1))
              xp=fx*fy*fz
              density(boxit%ind) = density(boxit%ind) - xp*charge
           end do
        end do
     end do
  else
     do i3=isz,iez
        zp = mpz(i3-isz)
        if (abs(zp) < mp_tiny) cycle
        call ind_positions_new(perz,i3,n3i,j3,goz) 
        do i2=isy,iey
           yp = zp*mpy(i2-isy)
           if (abs(yp) < mp_tiny) cycle
           call ind_positions_new(pery,i2,n2i,j2,goy)
           do i1=isx,iex
              xp = yp*mpx(i1-isx)
              if (abs(xp) < mp_tiny) cycle
              call ind_positions_new(perx,i1,n1i,j1,gox)
              if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1i+(j3-i3s)*n1i*n2i
                 density(ind)=density(ind)-xp*charge
              endif
           enddo
        enddo
     enddo
  end if

  call f_release_routine()

  contains

    !> Determine the index in which the potential must be inserted, following the BC
    !! Determine also whether the index is inside or outside the box for free BC
    pure subroutine ind_positions_new(periodic,i,ni,j,go)
      implicit none
      logical, intent(in) :: periodic
      integer, intent(in) :: i,ni
      logical, intent(out) :: go
      integer, intent(out) :: j

      if (periodic) then
         go=.true.
         j=modulo(i,ni)
      else
         j=i
         if (i >= -14 .and. i <= ni-15) then
            go=.true.
         else
            go=.false.
         end if
      end if

    END SUBROUTINE ind_positions_new


end subroutine gaussian_density
