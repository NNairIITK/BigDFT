!> @file
!! Manual on how to use plugins.
!! @author
!!    Copyright (C) 2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!
!! @defgroup FLIB_PLUGIN  Extension by plug-ins (flib)
!! @ingroup FLIB
!! @brief flib dynamic module load routines
!! @details
!! The library provides a way to dynamically load shared libraries at runtime.
!!
!! At any moment in the code, one can load a shared library with the following
!! statement. It takes as first argument the name of the library to load, without
!! the \c lib prefix and the \c .so suffix. The library should be in the current
!! working directory.
!! @snippet test_plugin.f90 load
!!
!! To be a valid loadable dynamic library, it should contains a symbol corresponding
!! to a routine called init().
!! @snippet pong.f90 init

!> Example of plugin usage to add capabilities.
program test
  use yaml_output
  use module_f_objects, except => f_object_signal_add_str, except2 => f_object_signal_connect
  use dictionaries

  implicit none
  
  integer :: ierr
  character(len = max_field_length) :: pong
  character(len = 256) :: mess

  call f_lib_initialize()

  ! Define here some signals.
  call f_object_add_signal("class", "ping", 1)

  call f_object_signal_connect("class", "ping", me, 1, ierr)

!!! [load]
  call plugin_load("pong", ierr)
!!! [load]
  call yaml_map("status", ierr)
  if (ierr /= 0) then
     call plugin_error(mess)
     call yaml_map("error", mess)
  end if

  if (f_object_signal_prepare("class", "ping")) then
     call f_object_signal_add_str("class", "ping", pong)
     call f_object_signal_emit("class", "ping")
     call yaml_map("ping", pong)
  else
     call yaml_map("ping", "no listeners")
  end if

  call f_lib_finalize()

contains

  subroutine me(p)
    use dictionaries
    character(len = max_field_length), intent(out) :: p

    write(p, "(A)") "no answer"
  end subroutine me

end program test
