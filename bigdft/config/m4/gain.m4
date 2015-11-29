# -*- Autoconf -*-
#
# Copyright (c) 2015 BigDFT Group (Luigi Genovese)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_LIBGAIN],
[dnl Test for libGaIn
  AC_ARG_WITH(libgain-libs, AS_HELP_STRING([--with-libgain-libs], [Give the linker flags for an external libGaIn modules (default = None).]), ac_libgain_libdir=$withval, ac_libgain_libdir=)
  AC_ARG_WITH(libgain-incs, AS_HELP_STRING([--with-libgain-incs], [Give the compiler include flags for an external libGaIn library (default = None).]), ac_libgain_incdir=$withval, ac_libgain_incdir=)
  
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  
  dnl Test the modules for compilation
  AC_MSG_CHECKING([for libGaIn modules])
  FCFLAGS_SVG=$FCFLAGS
  if test -n "$ac_libgain_incdir" ; then
    FCFLAGS="$FCFLAGS $ac_libgain_incdir"
  elif test -n "$C_INCLUDE_PATH" ; then
    for path in ${C_INCLUDE_PATH//:/ }; do
      ac_libgain_incdir="$ac_libgain_incdir -I$path"
    done
    FCFLAGS="$FCFLAGS $ac_libgain_incdir"
  fi
  AC_COMPILE_IFELSE([[program main
  implicit none
  double precision, external :: cc_coulomb_cc
  write(*,*) cc_coulomb_cc(1.d0, 0.d0, 0,0,0, 1.d0, 0.d0, 0,0,0)
end program]], withlibgainmod=yes, withlibgainmod=no)
  AC_MSG_RESULT($withlibgainmod)
  FCFLAGS=$FCFLAGS_SVG

  dnl Test the library of libGaIn.
  LIBS_SVG=$LIBS
  if test -n "$ac_libgain_libdir" ; then
    LIBS="$LIBS $ac_libgain_libdir"
  else
    ac_libgain_libdir="-lxcf90 -lxc"
  fi
  AC_CHECK_LIB(GaIn, cc_coulomb_cc, ac_use_libgain=yes, ac_use_libgain=no, [])
  LIBS=$LIBS_SVG
  
  if test "$ac_use_libgain" = "yes" -a "$withlibgainmod" = "yes" ; then
    LIB_GAIN_CFLAGS=$ac_libgain_incdir
    LIB_GAIN_LIBS=$ac_libgain_libdir
    ac_use_libgain="yes"
  else
    ac_use_libgain="no"
  fi
  
  if test "$ac_use_libgain" = "no" ; then
    AC_MSG_ERROR([libGaIn is not available, install libGaIn and provide paths --with-libgain-libs --with-libgain-incs.])
  fi
  
  AC_SUBST(LIB_GAIN_CFLAGS)
  AC_SUBST(LIB_GAIN_LIBS)

  AC_LANG_POP(Fortran)
])
