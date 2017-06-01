# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_AMBERTOOLS],
[dnl test for the present of ambertools compiler (experimental)
  AC_PATH_PROG([NAB], [nab])
  AC_SUBST(NAB)
  AC_ARG_WITH(ambertools, AS_HELP_STRING([--with-ambertools],
              [Give the path of the ambertools libraries (example = /home/<username>/lib/). Do not use the -L before the path(es), just give the plain path.]),
              ac_ambertools_dir=$withval, ac_ambertools_dir=)
  ax_have_amber="no"
  if test -n "$ac_ambertools_dir" ; then
    LIB_AMBERTOOLS_LIBS="-L$ac_ambertools_dir -lsff -lnab -lpbsa -lcifparse -lrism -lfftw3 -larpack"
    if test x"$NAB" = "x" ; then
        AC_MSG_ERROR(["Path to ambertools library is given by user, but nab compiler was not found on the machine."])
    fi
    ax_have_amber="yes"
  fi
  AM_CONDITIONAL(HAVE_AMBERTOOLS, test $ax_have_amber = "yes")
])