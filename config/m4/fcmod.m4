# Define a macro to test module output of the fortran compiler.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FC_MOD],
[
  AC_MSG_CHECKING([for module output in Fortran.])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  ax_fc_mod_compile=no
  AC_COMPILE_IFELSE([
module modtest
  integer, public :: value
end module modtest
], [ax_fc_mod_compile=yes],
   [AC_MSG_FAILURE(Fortran compiler cannot compile modules.)])
  if test $ax_fc_mod_compile = "yes" ; then
    ax_fc_mod_name="unknown"
    if test -s modtest.mod ; then
      ax_fc_mod_ext="mod"
      ax_fc_mod_capitalize="no"
      ax_fc_mod_name="module"
      rm -f modtest.mod
    fi
    if test -s modtest.MOD ; then
      ax_fc_mod_ext="MOD"
      ax_fc_mod_capitalize="no"
      ax_fc_mod_name="module"
      rm -f modtest.MOD
    fi
    if test -s MODTEST.MOD ; then
      ax_fc_mod_ext="MOD"
      ax_fc_mod_capitalize="yes"
      ax_fc_mod_name="MODULE"
      rm -f MODTEST.MOD
    fi
    if test -s MODTEST.mod ; then
      ax_fc_mod_ext="mod"
      ax_fc_mod_capitalize="yes"
      ax_fc_mod_name="MODULE"
      rm -f MODTEST.mod
    fi
    if test $ax_fc_mod_name = "unknown" ; then
       AC_MSG_ERROR(Unknown module naming scheme for Fortran compiler.)
    fi  
  fi
  AC_MSG_RESULT([$ax_fc_mod_name.$ax_fc_mod_ext])
])
