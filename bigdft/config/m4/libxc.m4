# -*- Autoconf -*-
#
# Copyright (c) 2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_LIBXC],
[dnl Test for libXC
  AC_ARG_WITH(libxc-libs, AS_HELP_STRING([--with-libxc-libs], [Give the compiler include flags for an external libXC library (default = None).]), ac_libxc_libdir=$withval, ac_libxc_libdir=)
  AC_ARG_WITH(libxc-incs, AS_HELP_STRING([--with-libxc-incs], [Give the linker flags for an external libXC modules (default = None).]), ac_libxc_incdir=$withval, ac_libxc_incdir=)
  
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  
  dnl Run the detection.
  if test -n "$ac_libxc_libdir" -a -n "$ac_libxc_incdir" ; then
    dnl Test the modules for compilation
    AC_MSG_CHECKING([for libXC modules])
    FCFLAGS_SVG=$FCFLAGS
    FCFLAGS="$FCFLAGS $ac_libxc_incdir"
    AC_COMPILE_IFELSE([[program main
    use xc_f90_types_m
    use libxc_funcs_m
    use xc_f90_lib_m
  
    write(*,*) XC_FAMILY_GGA
  end program]], withlibxcmod=yes, withlibxcmod=no)
    AC_MSG_RESULT($withlibxcmod)
    FCFLAGS=$FCFLAGS_SVG
      
    dnl Test the library of libXC.
    LIBS_SVG=$LIBS
    LIBS="$LIBS $ac_libxc_libdir"
    AC_CHECK_LIB(xc, xc_f90_lda_vxc, ac_use_libxc=yes, ac_use_libxc=no)
    LIBS=$LIBS_SVG
    
    if test "$ac_use_libxc" = "yes" -a "$withlibxcmod" = "yes" ; then
      LIB_XC_CFLAGS=$ac_libxc_incdir
      LIB_XC_LIBS=$ac_libxc_libdir
      ac_use_libxc="yes"
    else
      ac_use_libxc="no"
    fi
  else
    LIB_XC_CFLAGS="-I/usr/include"
    PKG_CHECK_MODULES(LIB_XC, libxc >= 2.0, ac_use_libxc="yes", ac_use_libxc="no")
  fi
  
  if test "$ac_use_libxc" = "no" ; then
    AC_MSG_ERROR([libXC is not available, install libXC and provide paths --with-libxc-libs --with-libxc-incs.])
  fi
  
  AC_SUBST(LIB_XC_CFLAGS)
  AC_SUBST(LIB_XC_LIBS)

  AC_LANG_POP(Fortran)
])