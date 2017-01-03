# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_ALBORZ],
[dnl test for the presence of alborz
  AC_ARG_WITH(alborz, AS_HELP_STRING([--with-alborz],
              [Give the path of the alborz libraries (example = /home/<username>/alborz/). Do not use the -L before the path(es), just give the plain path.]),
              ac_alborz_dir=$withval, ac_alborz_dir=)
  ax_have_alborz="no"
  if test -n "$ac_alborz_dir" ; then
    LIB_ALBORZ_LIBS="$ac_alborz_dir"
    ax_have_alborz="yes"
  fi
  AM_CONDITIONAL(HAVE_ALBORZ, test $ax_have_alborz = "yes")
])
