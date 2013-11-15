!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_build_info
!! NAME
!!  m_build_info
!!
!! FUNCTION
!!  This module contains information about this particular version of ABINIT
!!  and its build parameters (useful for debugging).
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2013 ABINIT group (Yann Pouillon, Matteo Giantomassi)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

module m_build_info

  use defs_basis

  implicit none

! Try to prevent problems with the length of the argument 
! (should not exceed the line length of 132 character).
! So, always start the string to be replaced at the beginning of a line.

! Parameters set-up by Autoconf
  character(len=6),parameter :: abinit_version = &
&   "@VERSION@"
  character(len=*),parameter :: build_target   = &
&   "@ABINIT_TARGET@"

! More info on current version
  character(len=*),parameter :: version_major = &
&   "@ABINIT_VERSION_MAJOR@"
  character(len=*),parameter :: version_minor = &
&   "@ABINIT_VERSION_MINOR@"
  character(len=*),parameter :: version_micro = &
&   "@ABINIT_VERSION_MICRO@"
  character(len=*),parameter :: version_build = &
&   "@ABINIT_VERSION_BUILD@"

! Info on compilers. Try to prevent problems with the length of the argument 
! (should not exceed the line length of 132 character).
  character(len=*),parameter :: cc_info    = &
&   "@abi_cc_vendor@@abi_cc_version@"
  character(len=*),parameter :: cc_flags   = &
&   "@CFLAGS@"
  character(len=*),parameter :: cxx_info   = &
&   "@abi_cxx_vendor@@abi_cxx_version@"
  character(len=*),parameter :: cxx_flags  = &
&   "@CXXFLAGS@"
  character(len=*),parameter :: fc_info    = &
&   "@abi_fc_vendor@@abi_fc_version@"
  character(len=*),parameter :: fc_flags   = &
&   "@FCFLAGS@"
  character(len=*),parameter :: fc_ldflags = &
&   "@FC_LDFLAGS@"

! Info on optimizations
  character(len=*),parameter :: enable_debug    = &
&   "@enable_debug@"
  character(len=*),parameter :: enable_optim = &
&   "@enable_optim@"
  character(len=*),parameter :: cpu_info = &
&   "@abi_cpu_vendor@_@abi_cpu_model@"

! Info on MPI
  character(len=*),parameter :: enable_mpi       = &
&   "@enable_mpi@"
  character(len=*),parameter :: enable_mpi_io    = &
&   "@enable_mpi_io@"
  character(len=*),parameter :: enable_mpi_trace = &
&   "@enable_mpi_trace@"

! Info on GPU
  character(len=*),parameter :: enable_gpu    = &
&   "@enable_gpu@"

! Info on connectors / fallbacks
  character(len=*),parameter :: enable_connectors = &
&   "@enable_connectors@"
  character(len=*),parameter :: enable_fallbacks  = &
&   "@enable_fallbacks@"
  character(len=*),parameter :: dft_flavor        = &
&   "@lib_dft_flavor@"
  character(len=*),parameter :: fft_flavor        = &
&   "@lib_fft_flavor@"
  character(len=*),parameter :: linalg_flavor     = &
&   "@lib_linalg_flavor@"
  character(len=*),parameter :: math_flavor       = &
&   "@lib_math_flavor@"
  character(len=*),parameter :: timer_flavor      = &
&   "@lib_timer_flavor@"
  character(len=*),parameter :: trio_flavor       = &
&   "@lib_trio_flavor@"

! Info on experimental features
  character(len=*),parameter :: enable_bindings = &
&   "@enable_bindings@"
  character(len=*),parameter :: enable_exports = &
&   "@enable_exports@"
  character(len=*),parameter :: enable_gw_dpc = &
&   "@enable_gw_dpc@"
 
#if defined HAVE_BZR_BRANCH
! Info on Bazaar branch (if applies)
  character(len=*),parameter :: bzr_branch = &
&   "@bzr_branch@"
  character(len=*),parameter :: bzr_revno  = &
&   "@bzr_revno@"
  character(len=*),parameter :: bzr_clean  = &
&   "@bzr_clean@"
#endif

contains  !===========================================================
!!***

!!****f* ABINIT/m_build_info/dump_config
!! NAME
!!  dump_config
!!
!! FUNCTION
!!  Reports a printout of the information stored in m_build_infos,
!!  useful for error messages and debugging.
!!
!! INPUTS
!!  my_unit= Fortran unit number (optional, default is std_out)
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

subroutine dump_config(my_unit)

  use defs_basis

  implicit none

!Arguments ------------------------------------
  integer,optional,intent(in) :: my_unit

!Local variables-------------------------------
  integer :: unt

! *********************************************************************
 
! TODO: things that might be added through preprocessing options, e.g.
! date and time of compilation 

  unt = std_out
  if ( present(my_unit) ) unt = my_unit

  write(unt,*)
  write(unt,'(1x,a)') repeat('+',78)
  write(unt,*)
  write(unt,'(a)' )' === Build Information === '
  write(unt,'(2a)')'  Version       : ',trim(abinit_version)
  write(unt,'(2a)')'  Build target  : ',trim(build_target)
  write(unt,'(2a)')'  Build date    : ',trim(version_build)
  write(unt,*)
  write(unt,'(a)' )' === Compiler Suite === '
  write(unt,'(2a)')'  C compiler       : ',trim(cc_info)
  write(unt,'(2a)')'  CFLAGS           : ',trim(cc_flags)
  write(unt,'(2a)')'  C++ compiler     : ',trim(cxx_info)
  write(unt,'(2a)')'  CXXFLAGS         : ',trim(cxx_flags)
  write(unt,'(2a)')'  Fortran compiler : ',trim(fc_info)
  write(unt,'(2a)')'  FCFLAGS          : ',trim(fc_flags)
  write(unt,'(2a)')'  FC_LDFLAGS       : ',trim(fc_ldflags)
  write(unt,*)
  write(unt,'(a) ')' === Optimizations === '
  write(unt,'(2a)')'  Debug level        : ',trim(enable_debug)
  write(unt,'(2a)')'  Optimization level : ',trim(enable_optim)
  write(unt,'(2a)')'  Architecture       : ',trim(cpu_info)
  write(unt,*)
  write(unt,'(a) ')' === MPI === '
  write(unt,'(2a)')'  Parallel build : ',trim(enable_mpi)
  write(unt,'(2a)')'  Parallel I/O   : ',trim(enable_mpi_io)
  write(unt,'(2a)')'  Time tracing   : ',trim(enable_mpi_trace)
  write(unt,'(2a)')'  GPU support    : ',trim(enable_gpu)
  write(unt,*)
  write(unt,'(a) ')' === Connectors / Fallbacks === '
  write(unt,'(2a)')'  Connectors on : ',trim(enable_connectors)
  write(unt,'(2a)')'  Fallbacks on  : ',trim(enable_fallbacks)
  write(unt,'(2a)')'  DFT flavor    : ',trim(dft_flavor)
  write(unt,'(2a)')'  FFT flavor    : ',trim(fft_flavor)
  write(unt,'(2a)')'  LINALG flavor : ',trim(linalg_flavor)
  write(unt,'(2a)')'  MATH flavor   : ',trim(math_flavor)
  write(unt,'(2a)')'  TIMER flavor  : ',trim(timer_flavor)
  write(unt,'(2a)')'  TRIO flavor   : ',trim(trio_flavor)
  write(unt,*)
  write(unt,'(a)' )' === Experimental features === '
  write(unt,'(2a)')'  Bindings            : ',trim(enable_bindings)
  write(unt,'(2a)')'  Exports             : ',trim(enable_exports)
  write(unt,'(2a)')'  GW double-precision : ',trim(enable_gw_dpc)
  write(unt,*)
#if defined HAVE_BZR_BRANCH
  write(unt,'(a)' )' === Bazaar branch information === '
  write(unt,'(2a)')'  Branch ID : ',trim(bzr_branch)
  write(unt,'(2a)')'  Revision  : ',trim(bzr_revno)
  write(unt,'(2a)')'  Committed : ',trim(bzr_clean)
  write(unt,*)
#endif
  write(unt,'(1x,a)') repeat('+',78)
  write(unt,*)

end subroutine dump_config
!!***

end module m_build_info
!!***
