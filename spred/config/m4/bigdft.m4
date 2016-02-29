# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_BIGDFT],
[dnl Test for Bigdft
  AC_ARG_WITH(bigdft-libs, AS_HELP_STRING([--with-bigdft-libs], [Give the linker flags for an external Bigdft modules (default = None).]), ax_bigdft_libs=$withval, ax_bigdft_libs=)
  AC_ARG_WITH(bigdft-incs, AS_HELP_STRING([--with-bigdft-incs], [Give the compiler include flags for an external Bigdft library (default = None).]), ax_bigdft_incdir=$withval, ax_bigdft_incdir=)
  
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])

  dnl try first with pkg-config
  PKG_CHECK_MODULES([BIGDFT],
                    [bigdft >= 1.8],
                    [ax_have_bigdft=yes],
                    [ax_have_bigdft=no])
  if test "$ax_have_bigdft" = "yes" ; then
    if test -z "${BIGDFT_CFLAGS// }" -a -n "$C_INCLUDE_PATH" ; then
      for path in ${C_INCLUDE_PATH//:/ }; do
        ax_bigdft_incdir="$ax_bigdft_incdir -I$path"
      done
      LIB_BIGDFT_CFLAGS=$ax_bigdft_incdir
    else
      LIB_BIGDFT_CFLAGS=$BIGDFT_CFLAGS
    fi
    LIB_BIGDFT_LIBS=$BIGDFT_LIBS
  fi

  dnl try by hand search if failed
  if test "$ax_have_bigdft" != "yes" ; then
    dnl Test the modules for compilation
    AC_MSG_CHECKING([for Bigdft modules])
    FCFLAGS_SVG=$FCFLAGS
    if test -n "$ax_bigdft_incdir" ; then
      FCFLAGS="$FCFLAGS $ax_bigdft_incdir"
    elif test -n "$C_INCLUDE_PATH" ; then
      for path in ${C_INCLUDE_PATH//:/ }; do
        ax_bigdft_incdir="$ax_bigdft_incdir -I$path"
      done
      FCFLAGS="$FCFLAGS $ax_bigdft_incdir"
    fi
    FCFLAGS="$FCFLAGS $LIB_FLIB_CFLAGS $LIB_PSOLVER_CFLAGS"
    AC_COMPILE_IFELSE([[program main
    use bigdft_run
  
    type(run_objects) :: run
    type(state_properties) :: outs
    integer :: info
    
    call bigdft_state(run, outs, info)
  end program]], withbigdftmod=yes, withbigdftmod=no)
    AC_MSG_RESULT($withbigdftmod)
  
    dnl Test the bigdft library.
    AC_MSG_CHECKING([for Bigdft library])
    LIBS_SVG=$LIBS
    if test -z "$ax_bigdft_libs" ; then
      ax_bigdft_libs="-lbigdft-1"
    fi
    LIBS="$ax_bigdft_libs $LIBS_SVG"
    AC_LINK_IFELSE(
      AC_LANG_PROGRAM([], [[
    use bigdft_run
  
    type(run_objects) :: run
    type(state_properties) :: outs
    integer :: info
    
    call bigdft_state(run, outs, info)
  ]]),
      [ax_have_bigdft=yes],
      [ax_have_bigdft=no])
    AC_MSG_RESULT($ax_have_bigdft)
  
    LIBS=$LIBS_SVG
    FCFLAGS=$FCFLAGS_SVG
    
    if test "$ax_have_bigdft" = "yes" -a "$withbigdftmod" = "yes" ; then
      LIB_BIGDFT_CFLAGS=$ax_bigdft_incdir
      LIB_BIGDFT_LIBS=$ax_bigdft_libs
      ax_have_bigdft="yes"
    else
      ax_have_bigdft="no"
    fi
  fi
  
  dnl LIB_XC_CFLAGS="-I/usr/include"
  dnl   PKG_CHECK_MODULES(LIB_XC, bigdft >= 2.0, ax_have_bigdft="yes", ax_have_bigdft="no")
  
  AC_SUBST(LIB_BIGDFT_CFLAGS)
  AC_SUBST(LIB_BIGDFT_LIBS)

  AC_LANG_POP(Fortran)
])
