# Define a macro to test flush() as intrinsict.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FC_FLUSH],
[
  AC_MSG_CHECKING([for flush(6) in Fortran.])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  cat > flushtest.f90 <<EOF
program test_flush

  implicit none

  write(*,"(A)") "yes"
  flush(6)

end program test_flush
EOF
  dnl Assume first that it should compile and run.
  ax_fc_flush="no"
  ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
    ax_fc_flush=`./flushtest.x 2> /dev/null`;
    if test "$?" != 0 ; then
      ax_fc_flush="no"
    fi
  fi
  dnl Assume second that it should compile and run with Intel option.
  FCFLAGS_SVG="$FCFLAGS"
  if test x"$ax_fc_flush" == x"no" ; then
    FCFLAGS="$FCFLAGS -assume noold_unit_star"
    ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_fc_flush=`./flushtest.x 2> /dev/null`;
      if test "$?" != 0 ; then
        ax_fc_flush="no"
      fi
    fi
  fi
  if test x"$ax_fc_flush" != x"yes" ; then
    FCFLAGS="$FCFLAGS_SVG"
  fi
  dnl Assume third that it should compile to have it.
  if test x"$ax_fc_flush" == x"no" ; then
    ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ax_fc_flush="yes"
    fi
  fi
  rm -f flushtest*
  if test x"$ax_fc_flush" == x"yes" ; then
    AC_DEFINE([HAVE_FC_FLUSH], [1], [Flush(6) can be used safely in fortran])
  fi
  AM_CONDITIONAL([HAVE_FC_FLUSH], [test x"$ax_fc_flush" == x"yes"])
  AC_LANG_POP(Fortran)

  AC_MSG_RESULT([$ax_fc_flush])
])
