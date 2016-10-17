!> @file
!!  Test of the buffer allocation informations examples of override
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program f_buffer_allocations
  use futile
  use wrapper_mpi
  implicit none
  real(f_double), dimension(:), pointer :: d1_ptr,d1_ptr_exotic

  call f_lib_initialize()
  call mpiinit()
  !allocate the pointer in the normal case
  d1_ptr=f_malloc_ptr(-1.to.34,id='d1_ptr')
  call output("Normal",d1_ptr)
  !now insert some extra information in the buffer
  !allocate the pointer in the normal case
  d1_ptr_exotic=f_malloc_ptr(-1.to.34,id='d1_ptr_exotic',&
       info='{Type: SHARED}')

  call output("Exotic",d1_ptr_exotic)

  call mpibarrier()

  call yaml_mapping_open('Status of the database',advance='no')
  call yaml_comment('Database status',hfill='~')
  call f_malloc_dump_status()
  call yaml_mapping_close(advance='no')
  call yaml_comment('',hfill='~')
  call f_free_ptr(d1_ptr)
  call f_free_ptr(d1_ptr_exotic)
  call mpifinalize()
  call f_lib_finalize()
  contains
  subroutine output(name, d1_ptr_in)
    character(len=*), intent(in):: name
    real(f_double), dimension(:), pointer, intent(in) :: d1_ptr_in
    call yaml_mapping_open(name//' pointer')
    call yaml_map('Shape',shape(d1_ptr_in))
    call yaml_map('Address',f_loc(d1_ptr_in))
    call yaml_map('Lbound',lbound(d1_ptr_in))
    call yaml_mapping_close()
  end subroutine
end program f_buffer_allocations
