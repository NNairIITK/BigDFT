# -*- Autoconf -*-
#
# Copyright (c) 2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_PYYAML],
[dnl Test for libXC
  dnl Python itself (set to ":" if not present).
  AC_REQUIRE([AM_PATH_PYTHON])

  AC_ARG_WITH([pyyaml-path],
              AS_HELP_STRING([--with-pyyaml-path], [give a path to find YAML python module.]),
              [AX_PYYAML_PATH=$withval], [AX_PYYAML_PATH=$pyexecdir])
  ax_have_pyyaml="no"
  if test "$PYTHON" != ":" ; then
    AC_MSG_CHECKING([for PyYAML and CLoader from])
    cat > pytest << EOF
import sys
sys.path.insert(0, "$AX_PYYAML_PATH")

try:
  import yaml

  try:
    l = yaml.CLoader
    print "0"
  except:
    print "2"
except:
  print "1"
EOF
    eval=$($PYTHON pytest)
    rm -f pytest
    case "$eval" in
      "0")
        AC_MSG_RESULT([$AX_PYYAML_PATH])
        ax_have_pyyaml="yes"
        ;;
      "2")
        AC_MSG_RESULT([no])
        AC_MSG_WARN([CLoader not available, check that the PyYAML python module was linked with the C-level YAML implementation.])
        ax_have_pyyaml="yes"
        ;;
      *)
        AC_MSG_RESULT([no])
        AC_MSG_WARN([PyYAML not available, use --with-pyyaml-path to specify the module location.])
        ax_have_pyyaml="no"
        ;;
    esac
  fi
  AM_CONDITIONAL(HAVE_PYYAML, test x"$ax_have_pyyaml" = x"yes")
  AC_SUBST(AX_PYYAML_PATH)
])