AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Spred - DFT over wavelets],[1.8],[Damien.Caliste@cea.fr],[spred])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE([-Wno-portability])
dnl To keep the object files in the same structure as the source files (pb does not work with libABINIT)
dnl AM_INIT_AUTOMAKE([--warnings=no-portability, subdir-objects])

dnl Versioning
SPRED_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
SPRED_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
SPRED_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$SPRED_MICRO_VERSION" ; then
  SPRED_MICRO_VERSION=0
fi
SPRED_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(SPRED_MAJOR_VERSION)
AC_SUBST(SPRED_MINOR_VERSION)
AC_SUBST(SPRED_MICRO_VERSION)
AC_SUBST(SPRED_NANO_VERSION)

dnl This is a hack for usage inside ABINIT.
dnl if test -n "$FCFLAGS" ; then
dnl   FCFLAGS=`echo $FCFLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi
dnl if test -n "$F90FLAGS" ; then
dnl   F90FLAGS=`echo $F90FLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
#	       ABI_PROG_FC()

dnl Fortran 77 (to be removed)
AC_PROG_F77([], [Fortran 77])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_DEFS()
AX_FC_OPENMP()

dnl Set if the module name are capitalized or not with the FC compiler. 
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

dnl Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS_f90" ; then
  FCFLAGS="$FCFLAGS $FCFLAGS_f90"
fi
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

dnl Check for python
AM_PATH_PYTHON([2.3], , [:])

dnl Search for pkg-config
PKG_PROG_PKG_CONFIG()

AX_MPI()
AX_LINALG()

dnl Libraries
AX_FLIB()
if test "$ax_have_FUTILE" != "yes" ; then
  AC_MSG_ERROR(["Cannot find Futile library."])
fi


dnl to be removed when MPI goes to flib.
AX_PSOLVER() 
if test "$ax_have_PSOLVER" != "yes" ; then
  AC_MSG_ERROR(["Cannot find PSolver library."])
fi
AX_BIGDFT()
if test "$ax_have_bigdft" != "yes" ; then
  AC_MSG_ERROR(["Cannot find BigDFT library."])
fi

dnl Look for GLib and friends for bindings.
AX_INTROSPECTION()
AX_DYNAMIC_LIBRARIES($ax_have_introspection, [f_lib_initialize], [$LIB_FUTILE_LIBS])

dnl Test libraries required for benching / testing.
AX_PYYAML()

dnl Agregate all the optional libs.
LIB_SPRED_DEPS="$LIB_BIGDFT_LIBS $LIB_FUTILE_LIBS $LINALG_LIBS $GLIB_LIBS $MPI_LDFLAGS $MPI_LIBS $PYTHON_LIBS"
AC_SUBST(LIB_SPRED_DEPS)

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/spred_baselib/Makefile
src/Makefile
src/art/Makefile
src/mhgps_files/Makefile
src/global_files/Makefile
tests/Makefile
])

AC_OUTPUT

echo "
Spred configure summary:
=======================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS_EXT_DEPS $LIBS

Options:
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With GLib and co.:         $ax_have_glib
   | linked libraries:       $GLIB_LIBS

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
