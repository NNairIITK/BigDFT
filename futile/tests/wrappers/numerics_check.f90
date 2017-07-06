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


  !here some tests about the box usage
  call test_box_functions()

  call f_lib_finalize()

end program numeric_check


subroutine test_box_functions()
  use futile, gp=>f_double
  use box
  implicit none
  !local variables
  integer :: i,i1,i2,i3
  integer(f_long) :: t0,t1
  real(gp) :: totdot,res
  type(cell) :: mesh_ortho
  integer, dimension(3) :: ndims
  real(gp), dimension(:,:,:,:), allocatable :: v1,v2


  ndims=[100,100,100]

  mesh_ortho=cell_new('P',ndims,[1.0_gp,1.0_gp,1.0_gp])

  v1=f_malloc([3,ndims(1),ndims(2),ndims(3)],id='v1')
  v2=f_malloc([3,ndims(1),ndims(2),ndims(3)],id='v2')

  do i3=1,ndims(3)
     do i2=1,ndims(2)
        do i1=1,ndims(1)
           !the scalar product of these objects is 20.0
           v1(:,i1,i2,i3)=[1.0_gp,2.0_gp,3.0_gp]
           v2(:,i1,i2,i3)=[2.0_gp,3.0_gp,4.0_gp]
        end do
     end do
  end do

  totdot=0.0_gp
  t0=f_time()
  !$omp parallel do default(shared) &
  !$omp private(i,res)&
  !$omp reduction(+:totdot)
  do i3=1,ndims(3)
     do i2=1,ndims(2)
        do i1=1,ndims(1)
           res=dotp(mesh_ortho,v1(1,i1,i2,i3),v2(:,i1,i2,i3))
           res=res/20.0_gp
           totdot=totdot+res
           v2(:,i1,i2,i3)=res
        end do
     end do
  end do
  !$omp end parallel do
  t1=f_time()

  call yaml_map('TotDot',totdot)
  call yaml_map('TotSum',sum(v2))
  call yaml_map('Time spent in the loop (ns)',t1-t0)

  call f_free(v1)
  call f_free(v2)

end subroutine test_box_functions
