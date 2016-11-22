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
[dnl Test for PSolver
AC_REQUIRE([AX_FLIB])
AC_REQUIRE([AX_LINALG])
AC_REQUIRE([AX_MPI])
AC_REQUIRE([AX_PSOLVER])
AX_PACKAGE([BIGDFT],[1.8],[-lbigdft-1],[$LIB_PSOLVER_LIBS $LIB_CHESS_LIBS $LIB_FUTILE_LIBS $LINALG_LIBS],[$LIB_FUTILE_CFLAGS $LIB_PSOLVER_CFLAGS],
             [program main
    use bigdft_run

    type(run_objects) :: run
    type(state_properties) :: outs
    integer :: info

    call bigdft_state(run, outs, info)
  end program
])
])
