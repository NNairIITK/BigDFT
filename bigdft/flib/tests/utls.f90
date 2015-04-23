!> @file
!! Routine to tests f_utils module
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine f_utils_test()
  use f_precisions
  use f_utils
  use yaml_output
  use f_enums
  use dictionaries, only: f_loc
  implicit none
  !local variables
  type(f_enumerator) :: greetings=f_enumerator('Greetings',10) 
  type(f_enumerator) :: f1=f_enumerator('Ciao',1)              
  type(f_enumerator) :: f2=f_enumerator('Hello',2)             
  type(f_enumerator) :: f3=f_enumerator('Pizza',3)             
  integer :: unt,unt2,u
!  double precision :: t0
  integer(kind=8) :: i0,i1
  real(simple), dimension(3) :: r1
  real(double), dimension(3) :: r2
  real(quadruple), dimension(3) :: r4
  complex(simple), dimension(3) :: c1
  complex(double), dimension(3) :: c2
  complex(quadruple), dimension(3) :: c4
  integer(short), dimension(3) :: is
  integer(four), dimension(3) :: i4
  integer(long), dimension(3) :: il
  logical(byte), dimension(3) :: lb
  logical, dimension(3) :: l

  r4=real(10.0,quadruple)

!  call expq(r4(1),r4(2))

  r2=10.d0
!  call expq(r2(1),r2(2))

  !check the precisions
  call yaml_map('Kinds of reals',[f_loc(r1(2))-f_loc(r1(1)),f_loc(r2(2))-f_loc(r2(1)),f_loc(r4(2))-f_loc(r4(1))])
  call yaml_map('Kinds of complex',[f_loc(c1(2))-f_loc(c1(1)),f_loc(c2(2))-f_loc(c2(1)),f_loc(c4(2))-f_loc(c4(1))])
  call yaml_map('Kinds of integer',[f_loc(is(2))-f_loc(is(1)),f_loc(i4(2))-f_loc(i4(1)),f_loc(il(2))-f_loc(il(1))])
  call yaml_map('Kinds of logicals',[f_loc(lb(2))-f_loc(lb(1)),f_loc(l(2))-f_loc(l(1))])

!  print *,'test',r4(1),r4(2),'scale',scale(100.0_quadruple,2),scale(100.0_quadruple,4)!,exp(r4(1))
!  print *,'test',r2(1),r2(2),exp(r2(1)),exponent(r2(2)),scale(r2(2),exponent(r2(2)))!,scale(r2(1))

  !see if the attribute can live outside a given scope
  call f_enum_attr(f1,attr=greetings)
  call f_enum_attr(f2,attr=greetings)

  call yaml_map('Enum1 char',char(f1))
  call yaml_map('Enum1 int',int(f1))
  call yaml_map('Enum1 check',f1=='Ciao')

  call yaml_map('Greetings 1a',f1 .hasattr. 'Greetings') 
  call yaml_map('Greetings 1b',f1 .hasattr. 10)
  call yaml_map('Greetings 1c',f1 .hasattr. greetings) 

  call yaml_map('Greetings 2a',f2 .hasattr. 'Greetings') 
  call yaml_map('Greetings 2b',f2 .hasattr. 10)
  call yaml_map('Greetings 2c',f2 .hasattr. greetings) 

  call yaml_map('Greetings 3a',f3 .hasattr. 'Greetings') 
  call yaml_map('Greetings 3b',f3 .hasattr. 10)
  call yaml_map('Greetings 3c',f3 .hasattr. greetings) 

  !wait one second
  !t0=dble(f_time())*1.d-9
  i0=f_time()
  call yaml_map('Absolute time before pause (since epoch)',yaml_walltime_toa(i0))
  call f_pause(1)
  call yaml_map('Time spent after pause (s)',yaml_walltime_toa(f_time()-i0))!dble(f_time())*1.d-9-t0)
  !call yaml_map('Conversion of walltimes in standard form',yaml_walltime_toa(f_time()))
  !open files and get free_units
  unt=f_get_free_unit()
  call yaml_map('First unit which can be opened',unt)
  call f_open_file(unt,'Testopen')
  call yaml_map('Opened file in unit',unt)
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Second unit which can be opened',unt2)
  call f_open_file(unt2,'Testopen2')
  call yaml_map('Opened file in unit',unt2)
  call f_close(unt)
  call f_delete_file('Testopen')
  call f_delete_file('Testopen2')
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Third unit which can be opened',unt2)

  call yaml_mapping_open('Bastian test')
  u=f_get_free_unit()                                                   
  call yaml_map('First file',u)
  open(unit=u,file='test') 
  !call f_open_file(u,'test')
  call yaml_map('First file opened',u)
  unt=f_get_free_unit()                                                 
  call yaml_map('Second file',unt)
  open(unit=unt,file='test2')
  call yaml_map('Second file opened',unt)
  unt2=f_get_free_unit()                                      
  call yaml_map('Third file',unt2)
  open(unit=unt2,file='test3')     
  call yaml_map('Third file opened',unt2)
  call f_delete_file('test')
  call f_delete_file('test2')
  call f_delete_file('test3')
  call yaml_mapping_close()
  call yaml_map('If this value is 7 then all files have been correctly closed',f_get_free_unit())

end subroutine f_utils_test

!!$subroutine expq(x,res)
!!$  use f_precisions
!!$  implicit none
!!$  real(quadruple), intent(in) :: x
!!$  real(quadruple), intent(out) :: res
!!$  !local variables
!!$  real(quadruple), parameter :: log2=log(2.0_quadruple)
!!$  integer(long) :: i
!!$  real(quadruple) :: r,a
!!$
!!$  !easy case
!!$  if (x==real(0.0,quadruple)) then
!!$     res=real(1.0,quadruple)
!!$  else
!!$     !first determine the rest
!!$     r=x/log2
!!$     a=aint(r,quadruple)
!!$     r=r-a
!!$     r=2**r
!!$     i=int(a,long)
!!$     res=scale(r,i)
!!$  end if
!!$
!!$!  print *,'res',res,exp(x)
!!$  
!!$end subroutine expq
