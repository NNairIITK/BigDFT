# -*- Autoconf -*-
#
# Copyright (c) 2015 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_FLIB],
[dnl Test for FLib
  AC_ARG_WITH(flib-libs, AS_HELP_STRING([--with-flib-libs], [Give the linker flags for an external FLib modules (default = None).]), ac_flib_libdir=$withval, ac_flib_libdir=)
  AC_ARG_WITH(flib-incs, AS_HELP_STRING([--with-flib-incs], [Give the compiler include flags for an external FLib library (default = None).]), ac_flib_incdir=$withval, ac_flib_incdir=)
  
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  
  dnl Test the modules for compilation
  AC_MSG_CHECKING([for FLib modules])
  FCFLAGS_SVG=$FCFLAGS
  if test -n "$ac_flib_incdir" ; then
    FCFLAGS="$FCFLAGS $ac_flib_incdir"
  elif test -n "$C_INCLUDE_PATH" ; then
    for path in ${C_INCLUDE_PATH//:/ }; do
      ac_flib_incdir="$ac_flib_incdir -I$path"
    done
    FCFLAGS="$FCFLAGS $ac_flib_incdir"
  fi
  AC_COMPILE_IFELSE([[program main
  use yaml_parse
  use yaml_output
  use f_utils
  use dynamic_memory
  use dictionaries

  call yaml_map("toto", "titi")
end program]], withflibmod=yes, withflibmod=no)
  AC_MSG_RESULT($withflibmod)

  dnl Test the library of flib.
  AC_MSG_CHECKING([for flib library])
  LIBS_SVG=$LIBS
  if test -z "$ac_flib_libdir" ; then
    ac_flib_libdir="-lflib-1"
  fi
  LIBS="$ac_flib_libdir $LIBS_SVG"
  AC_LINK_IFELSE(
    AC_LANG_PROGRAM([], [[
call f_lib_initialize()
]]),
    [ax_have_flib=yes],
    [ax_have_flib=no])
  if test $ax_have_flib != "yes" ; then
    dnl Static case, need to link with additional libs.
    ac_flib_libdir="$ac_flib_libdir -lyaml -lrt"
    LIBS="$ac_flib_libdir $LIBS_SVG"
    AC_LINK_IFELSE(
      AC_LANG_PROGRAM([], [[
call f_lib_initialize()
]]),
      [ax_have_flib=yes],
      [ax_have_flib=no])
  fi
  AC_MSG_RESULT($ax_have_flib)

  if test "$ax_have_flib" = "yes" -a "$withflibmod" = "yes" ; then
    LIB_FLIB_CFLAGS=$ac_flib_incdir
    LIB_FLIB_LIBS=$ac_flib_libdir
    ax_have_flib="yes"
  else
    ax_have_flib="no"
  fi

  dnl LIB_XC_CFLAGS="-I/usr/include"
  dnl   PKG_CHECK_MODULES(LIB_XC, flib >= 2.0, ax_have_flib="yes", ax_have_flib="no")
  
  AC_SUBST(LIB_FLIB_CFLAGS)
  AC_SUBST(LIB_FLIB_LIBS)

  dnl Try to find libflib-1.a for possible later inclusion.
  

  FCFLAGS=$FCFLAGS_SVG
  LIBS=$LIBS_SVG

  AC_LANG_POP(Fortran)
])