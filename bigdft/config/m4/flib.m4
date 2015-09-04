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
  FCFLAGS=$FCFLAGS_SVG

  dnl Test the library of flib.
  LIBS_SVG=$LIBS
  if test -n "$ac_flib_libdir" ; then
    LIBS="$LIBS $ac_flib_libdir"
  else
    ac_flib_libdir="-lflib-1"
  fi
  AC_CHECK_LIB(flib-1, f_lib_initialize, ac_use_flib=yes, ac_use_flib=no, -lyaml)
  LIBS=$LIBS_SVG
  
  if test "$ac_use_flib" = "yes" -a "$withflibmod" = "yes" ; then
    LIB_FLIB_CFLAGS=$ac_flib_incdir
    LIB_FLIB_LIBS=$ac_flib_libdir
    ac_use_flib="yes"
  else
    ac_use_flib="no"
  fi

  dnl LIB_XC_CFLAGS="-I/usr/include"
  dnl   PKG_CHECK_MODULES(LIB_XC, flib >= 2.0, ac_use_flib="yes", ac_use_flib="no")
  
  AC_SUBST(LIB_FLIB_CFLAGS)
  AC_SUBST(LIB_FLIB_LIBS)

  AC_LANG_POP(Fortran)
])