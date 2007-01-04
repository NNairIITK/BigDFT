# -*- Autoconf -*-
#
# Copyright (c) 2005-2006 The ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _ABI_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ABSOFT],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])

 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([ABSOFT_FC],1,[Define to 1 if you are using the ABSOFT Fortran compiler])
  fc_type="absoft"
  fc_version=`echo "${abi_result}" | sed -e 's/Pro Fortran //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ABSOFT



# _ABI_CHECK_FC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_COMPAQ],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])

 fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
 abi_result="${fc_info_string}"
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([COMPAQ_FC],1,[Define to 1 if you are using the COMPAQ Fortran compiler])
  fc_type="compaq"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V//;s/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_COMPAQ



# _ABI_CHECK_FC_FUJITSU(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_FUJITSU],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])

 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FUJITSU_FC],1,[Define to 1 if you are using the FUJITSU Fortran compiler])
  fc_type="fujitsu"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_FUJITSU



# _ABI_CHECK_FC_G95(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_G95],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the G95 Fortran compiler])

 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^G95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([GNU_FC],1,[Define to 1 if you are using the GNU Fortran compiler])
  AC_DEFINE([G95_FC],1,[Define to 1 if you are using the G95 Fortran compiler])
  fc_type="g95"
  fc_version=`echo ${abi_result} | sed -e 's/.*GCC //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_G95



# _ABI_CHECK_FC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GNU Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_GNU],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the GNU Fortran compiler])

 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^GNU Fortran 95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([GNU_FC],1,[Define to 1 if you are using the GNU Fortran compiler])
  fc_type="gnu"
  fc_version=`echo ${abi_result} | sed -e 's/.*GCC //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_GNU



# _ABI_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_IBM],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])

 fc_info_string=`$1 -qversion 2>&1`
 fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
 abi_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
  if test "${fc_garbage}" -gt 50; then
   AC_DEFINE([IBM_FC],1,[Define to 1 if you are using the IBM XL Fortran compiler])
   fc_type="ibm"
   fc_version="UNKNOWN"
   abi_result="yes"
  fi
 else
  AC_DEFINE([IBM_FC],1,[Define to 1 if you are using the IBM XL Fortran compiler])
  fc_type="ibm"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V\([[0-9\.]]*\) .*/\1/'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_IBM

# _ABI_CHECK_FC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_SUN],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the Sun WorkShop Fortran compiler])

 fc_info_string=`$1 -V 2>&1`
 abi_result=`echo "${fc_info_string}" | grep 'Sun Fortran 95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([SUN_FC],1,[Define to 1 if you are using the Sun WorkShop])
  fc_type="sun"
  fc_version="x"
  if test "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_SUN

# _ABI_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_INTEL],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the Intel Fortran compiler])

 fc_info_string=`$1 -v -V 2>&1 | sed -e '/^ifc: warning/d'`
 abi_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([INTEL_FC],1,[Define to 1 if you are using the Intel Fortran compiler])
  fc_type="intel"
  fc_version=`echo "${fc_info_string}" | grep '^Version' | sed -e 's/Version //;s/ .*//;s/ //g' | head -n 1`
  if test "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_INTEL



# _ABI_CHECK_FC_MIPSPRO(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_MIPSPRO],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])

 fc_info_string=`$1 -version 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([MIPSPRO_FC],1,[Define to 1 if you are using the MIPSpro Fortran compiler])
  fc_type="mipspro"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Version //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_MIPSPRO



# _ABI_CHECK_FC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PATHSCALE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])

 fc_info_string=`$1 -version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^PathScale'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([PATHSCALE_FC],1,[Define to 1 if you are using the PathScale Fortran compiler])
  fc_type="pathscale"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PATHSCALE



# _ABI_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PGI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])

 fc_info_string=`$1 -V 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^pgf90' | grep -v 'No files to process'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([PGI_FC],1,[Define to 1 if you are using the Portland Group Fortran compiler])
  fc_type="pgi"
  fc_version=`echo "${abi_result}" | sed -e 's/.* //; s/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PGI



 ##############################################################################



# ABI_CHECK_FORTRAN_EXIT()
# ------------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([ABI_CHECK_FORTRAN_EXIT],
[dnl Init
 fc_has_exit="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts exit()])

 dnl Try to compile a program calling exit()
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call exit(1)
  ]])], [fc_has_exit="yes"])
 AC_LANG_POP()

 if test "${fc_has_exit}" = "yes"; then
  AC_DEFINE([HAVE_FORTRAN_EXIT],1,
   [Define to 1 if your Fortran compiler supports exit()])
 fi

 AC_MSG_RESULT(${fc_has_exit})
]) # ABI_CHECK_FORTRAN_EXIT



# ABI_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ABI_PROG_FC],
[dnl Init
 if test "${fc_type}" = ""; then
  fc_type="UNKNOWN"
 fi
 if test "${fc_version}" = ""; then
  fc_version="UNKNOWN"
 fi
 fc_wrap="no"

 dnl Determine Fortran compiler type (the order is important)
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_G95(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_INTEL(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_GNU(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PATHSCALE(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PGI(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_COMPAQ(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_ABSOFT(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_MIPSPRO(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_FUJITSU(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_SUN(${FC})
 fi
 dnl Always keep that one at the end
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_IBM(${FC})
 fi

 dnl Schedule compiler info for substitution
 AC_SUBST(fc_type)
 AC_SUBST(fc_version)
 AC_SUBST(fc_wrap)
]) # ABI_PROG_FC



 ##############################################################################



# _ABI_TRICKS_LINALG()
# --------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([_ABI_TRICKS_LINALG],
[
])
