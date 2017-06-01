# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_ARCHIVE],
[dnl Archive input file support
  AC_ARG_WITH([archives],
              AS_HELP_STRING([--with-archives], [add support of archives for output files.]),
              [ax_have_archive=$withval], [ax_have_archive="yes"])
  AC_ARG_WITH([archives-path],
              AS_HELP_STRING([--with-archives-path], [give a path to find libarchive.]),
              [ac_path_archive=$withval])
  if test x"$ax_have_archive" == x"yes" ; then
     if test x"$ac_path_archive" == x"" ; then
        PKG_CHECK_MODULES([LIB_ARCHIVE],
                          [libarchive >= 2.4],
                          [ax_have_archive=yes],
                          [ax_have_archive=warn])
     else
        LDFLAGS_SVG="$LDFLAGS"
        AC_LANG_PUSH(C)
        LDFLAGS="-L$ac_path_archive/lib"
        AC_CHECK_LIB([archive], [archive_read_data_block],
                     [ax_have_archive=yes], [ax_have_archive=warn])
        if test x"$ax_have_archive" = x"yes" ; then
           LIB_ARCHIVE_CFLAGS="-I$ac_path_archive/include"
           LIB_ARCHIVE_LIBS="-L$ac_path_archive/lib -larchive"
        fi
        AC_LANG_POP(C)
        LDFLAGS="$LDFLAGS_SVG"
     fi
     if test x"$ax_have_archive" = x"yes" ; then
        AC_DEFINE([HAVE_LIB_ARCHIVE], [1], [libarchive is linkable.])
     else
        ax_have_archive="warn"
        AC_MSG_WARN([libarchive is not available])
     fi
  fi
  AC_SUBST(LIB_ARCHIVE_CFLAGS)
  AC_SUBST(LIB_ARCHIVE_LIBS)
  AM_CONDITIONAL(HAVE_LIB_ARCHIVE, test x"$ax_have_archive" = x"yes")
])