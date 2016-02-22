# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_PSOLVER],
[dnl Test for PSolver
  AC_ARG_WITH(psolver-libs, AS_HELP_STRING([--with-psolver-libs], [Give the linker flags for an external Psolver modules (default = None).]), ax_psolver_libdir=$withval, ax_psolver_libdir=)
  AC_ARG_WITH(psolver-incs, AS_HELP_STRING([--with-psolver-incs], [Give the compiler include flags for an external Psolver library (default = None).]), ax_psolver_incdir=$withval, ax_psolver_incdir=)
  
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  AC_REQUIRE([AX_FLIB])
  AC_REQUIRE([AX_LINALG])
  
  dnl Test the modules for compilation
  AC_MSG_CHECKING([for PSolver modules])
  FCFLAGS_SVG=$FCFLAGS
  if test -n "$ax_psolver_incdir" ; then
    FCFLAGS="$FCFLAGS $ax_psolver_incdir"
  elif test -n "$C_INCLUDE_PATH" ; then
    for path in ${C_INCLUDE_PATH//:/ }; do
      ax_psolver_incdir="$ax_psolver_incdir -I$path"
    done
    FCFLAGS="$FCFLAGS $ax_psolver_incdir"
  fi
  FCFLAGS="$FCFLAGS $LIB_FLIB_CFLAGS"
  AC_COMPILE_IFELSE([[program main
  use psbase
  use box
  use iobox
  use psbox
  use poisson_solver

  write(*,*) PS_getVersion()
end program]], withpsolvermod=yes, withpsolvermod=no)
  AC_MSG_RESULT($withpsolvermod)

  dnl Test the psolver library.
  AC_MSG_CHECKING([for PSolver library])
  LIBS_SVG=$LIBS
  if test -z "$ax_psolver_libdir" ; then
    ax_psolver_libdir="-lPSolver-1"
  fi
  LIBS="$ax_psolver_libdir $LINALG_LIBS $LIB_FLIB_LIBS $LIBS"
  AC_LINK_IFELSE(
    AC_LANG_PROGRAM([], [[
use Poisson_solver

type(coulomb_operator) :: kernel
real(dp), dimension(9) :: rhopot, potion
real(gp) :: eh

call H_Potential("G", kernel, rhopot, potion, eh, 0._dp, .false.)
]]),
    [ax_have_psolver=yes],
    [ax_have_psolver=no])
  AC_MSG_RESULT($ax_have_psolver)

  LIBS=$LIBS_SVG
  FCFLAGS=$FCFLAGS_SVG
  
  if test "$ax_have_psolver" = "yes" -a "$withpsolvermod" = "yes" ; then
    LIB_PSOLVER_CFLAGS=$ax_psolver_incdir
    LIB_PSOLVER_LIBS=$ax_psolver_libdir
    ax_have_psolver="yes"
  else
    ax_have_psolver="no"
  fi

  dnl LIB_XC_CFLAGS="-I/usr/include"
  dnl   PKG_CHECK_MODULES(LIB_XC, psolver >= 2.0, ax_have_psolver="yes", ax_have_psolver="no")
  
  AC_SUBST(LIB_PSOLVER_CFLAGS)
  AC_SUBST(LIB_PSOLVER_LIBS)

  AC_LANG_POP(Fortran)
])
