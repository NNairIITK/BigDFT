AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Chebyshev Sparse Solvers],[1.0],[stephan.mohr@bsc.es],[libCheSS])

AC_CONFIG_MACRO_DIR([m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE([-Wno-portability])

dnl Versioning
CHESS_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
CHESS_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
CHESS_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$CHESS_MICRO_VERSION" ; then
  CHESS_MICRO_VERSION=0
fi
CHESS_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(CHESS_MAJOR_VERSION)
AC_SUBST(CHESS_MINOR_VERSION)
AC_SUBST(CHESS_MICRO_VERSION)
AC_SUBST(CHESS_NANO_VERSION)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AX_FC_DEFS()
AX_FC_OPENMP()

dnl Set if the module name are capitalized or not with the FC compiler. 
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

if test -z "$FCFLAGS" ; then
  AC_SUBST(FCFLAGS, "-O2")
fi

dnl Look for ranlib and canonicalize it with the $host variable
AC_CHECK_TOOL(RANLIB, ranlib, :)
if test "$ranlib" = ":" ; then
  AC_MSG_ERROR(["No 'ranlib' program found."])
fi
dnl Look for ar and canonicalize it with the $host variable
AC_CHECK_TOOL(AR, ar, :)
if test "$ar" = ":" ; then
  AC_MSG_ERROR(["No 'ar' program found."])
fi

dnl Enable or not MPI.
AX_MPI()

dnl Python YAML support.
AX_PYYAML()

dnl The utilities executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_utilities=$enableval, ac_build_bin_utilities="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_utilities" = "yes")
AC_MSG_CHECKING([for CheSS binaries build])
AC_MSG_RESULT([$ac_build_bin_utilities])


AX_FLIB()

dnl Produce dynamic excecutables.
AX_DYNAMIC_LIBRARIES([], [f_lib_initialize], [$LIB_FUTILE_LIBS])

dnl Test accelerators
AX_ACC_CUDA()

AX_LINALG()
AX_ETSF_IO()

dnl Agregate all the optional libs.
LIB_CHESS_DEPS="$LIBETSFIO_LIBS $LIB_FUTILE_LIBS"
dnl if the futile libraries are dynamic add the dependencies
if test $ax_FUTILE_static != "yes" ; then
  LIB_CHESS_DEPS="$LIB_CHESS_DEPS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
fi
LIB_CHESS_PKG=""
if test $ax_build_dynamic != "yes" ; then
  LIB_CHESS_PKG="$LIB_CHESS_DEPS"
fi
AC_SUBST(LIB_CHESS_PKG)
AC_SUBST(LIB_CHESS_DEPS)

dnl Test for the developers' documentation.
AX_DOC()

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
chess.pc
doc/Makefile
doc/Doxyfile
sparsematrix/Makefile
src/configure.inc
src/Makefile
tests/Makefile
tests/sparsematrix/Makefile
])

AC_OUTPUT

echo "
Chebyshev Sparse Solvers configure summary:
================================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIB_CHESS_DEPS $LIBS

Build:
  Dynamic library:           $ax_build_dynamic

Options:
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With ETSF_IO:              $ax_have_etsf_io
   | include dir.:           $LIBETSFIO_INCLUDE

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$ax_have_etsf_io" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$ax_have_etsf_io" = x"warn" ; then
    echo "  ETSF_IO:                   disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi
