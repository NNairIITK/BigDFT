AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Poisson solver],[1.8],[Damien.Caliste@cea.fr],[psolver])

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
PSOLVER_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
PSOLVER_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
PSOLVER_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$PSOLVER_MICRO_VERSION" ; then
  PSOLVER_MICRO_VERSION=0
fi
PSOLVER_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(PSOLVER_MAJOR_VERSION)
AC_SUBST(PSOLVER_MINOR_VERSION)
AC_SUBST(PSOLVER_MICRO_VERSION)
AC_SUBST(PSOLVER_NANO_VERSION)

AC_SUBST(FLIB_MAJOR_VERSION, [$PSOLVER_MAJOR_VERSION])
AC_SUBST(FLIB_MINOR_VERSION, [$PSOLVER_MINOR_VERSION])
AC_SUBST(FLIB_MICRO_VERSION, [$PSOLVER_MICRO_VERSION])
AC_SUBST(FLIB_NANO_VERSION, [$PSOLVER_NANO_VERSION])

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the bindings and CXX for FFTW.
AC_PROG_CC()
AC_PROG_CXX()
AM_PROG_CC_C_O()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

dnl Fortran 77 (to be removed)
AC_PROG_F77([], [Fortran 77])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_FLUSH()
AX_FC_RECL_KIND()
AC_SUBST(RECL_INT_KIND, $ax_fc_recl_kind)
AX_FC_DEFS()
AX_FC_F2003()
AX_FC_OPENMP()

dnl Set if the module name are capitalized or not with the FC compiler. 
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)

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

dnl Enable or not MPI.
AX_MPI()

dnl C-level YAML support.
AX_YAML()

dnl Python YAML support.
AX_PYYAML()

dnl Produce dynamic excecutables.
AX_DYNAMIC_LIBRARIES()

dnl Test accelerators
AX_ACC_CUDA()
AX_ACC_OCL()
AX_ACC_MIC()

dnl Use external FLib.
AC_ARG_WITH(internal-flib, AS_HELP_STRING([--with-internal-flib], [Compile internal FLib (default = yes).]), internalflib=$withval, internalflib="yes")
if test x"$internalflib" != x"yes" ; then
  AX_FLIB()
  if test "$ac_use_flib" = "no" ; then
    AC_MSG_ERROR([FLib is not available, install FLib and provide paths --with-flib-libs --with-flib-incs.])
  fi
else
  LIB_FLIB_CFLAGS='-I$(top_builddir)/includes'
  AC_SUBST(LIB_FLIB_CFLAGS)
  LIB_FLIB_LIBS=""
  AC_SUBST(LIB_FLIB_LIBS)
fi
AM_CONDITIONAL(BUILD_FLIB, test "$internalflib" = "yes")

dnl Test the availability of clock_gettime() and rt library.
dnl This will be added in libs, on top, therefore care should be taken for a linker which require ordering
AC_LANG_PUSH(C)
AC_CHECK_HEADERS([time.h])
AC_CHECK_LIB([rt], [clock_gettime])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([strndup])
AC_LANG_POP(C)

dnl The Poisson solver library
AC_ARG_ENABLE(libpsolver, AS_HELP_STRING([--disable-libpsolver], [build the PSolver library (default = yes).]), ac_build_lib_PSolver=$enableval, ac_build_lib_PSolver="yes")

AM_CONDITIONAL(BUILD_LIB_PSOLVER, test "$ac_build_lib_PSolver" = "yes")
AC_MSG_CHECKING([for PSolver library build])
AC_MSG_RESULT([$ac_build_lib_PSolver])

dnl Linalg libraries.
if test $ac_build_lib_PSolver = "yes" ; then
AX_LINALG()
if test $ax_have_linalg != "yes" ; then
  AC_MSG_ERROR([Linear algebra libraries are not available.])
fi
fi

dnl Use the dgemmsy
if test $ac_build_lib_PSolver = "yes" ; then
ax_have_dgemmsy="no"
DGEMMSY_CPPFLAGS=
AC_ARG_ENABLE(dgemmsy, AS_HELP_STRING([--enable-dgemmsy], [Use dgemmsy (disabled by default).]),
                         ax_have_dgemmsy=$enableval, ax_have_dgemmsy="no")
if test x"$ax_have_dgemmsy" = "xyes" ; then
  AC_LANG_PUSH(C)

  DGEMMSY_CPPFLAGS="-msse3"
  CPPFLAGS_SVG=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $DGEMMSY_CPPFLAGS"
  AC_CHECK_HEADER([pmmintrin.h],
                  [ax_have_dgemmsy="yes"],
                  [ax_have_dgemmsy="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  AC_LANG_POP(C)
  if test "$ax_have_dgemmsy" = "no"; then
    DGEMMSY_CPPFLAGS=
    AC_MSG_WARN(["No SSE header found, dgemmsy disabled."])
  fi
fi
AC_SUBST(DGEMMSY_CPPFLAGS)
fi
AM_CONDITIONAL(USE_DGEMMSY, test "$ax_have_dgemmsy" = "yes")

dnl Support for ETSF_IO file format
if test $ac_build_lib_PSolver = "yes" ; then
AX_ETSF_IO()
else
AM_CONDITIONAL(HAVE_ETSF_IO, test "$ax_have_etsf_io" = "yes")
fi

dnl default installation directories
AC_ARG_WITH(moduledir,
            AS_HELP_STRING([--with-moduledir],
                           [installation directory for module files [[INCLUDEDIR]]]),
              ac_moduledir=$withval, ac_moduledir="no")
if test x"$ac_moduledir" != x"no" ; then
  moduledir=$ac_moduledir
else
  moduledir=${includedir}
fi
AC_SUBST(moduledir)

dnl Agregate all the optional libs.
LIBS_EXT_DEPS="$LIBCUDA_LIBS $LIBOCL_LIBS $LIBETSFIO_LIBS $LINALG_LIBS $LIB_FLIB_LIBS $GLIB_LIBS $MPI_LDFLAGS $MPI_LIBS"
AC_SUBST(LIBS_EXT_DEPS, [$LIBS_EXT_DEPS])

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
flib/Makefile
flib/src/Makefile
flib/src/f_utils.inc
flib/tests/Makefile
wrappers/Makefile
src/configure.inc
src/Makefile
tests/Makefile
])

AC_OUTPUT

echo "
Poisson-solver configure summary:
================================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS_EXT_DEPS $LIBS
  Linked dependencies:       $LIBS_SHORT_DEPS

Build:
  Library PSolver:           $ac_build_lib_PSolver
  Dynamic library:           $ax_build_dynamic
  Bindings:                  $ac_build_bindings
  User documentation:        $ac_have_pdflatex
  Devel. documentation:      $ac_devel_doc (Doxygen: $ac_have_doxygen)

Options:
  Max. Memory per task (GB): $MEMORY_LIMIT
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With Cuda support:         $ax_have_acc_cuda
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
  With OpenCL support:       $ax_have_acc_ocl
  With Intel MIC support:    $ax_have_acc_mic
  With dgemmsy support:      $ax_have_dgemmsy
  With ETSF_IO:              $ax_have_etsf_io
   | include dir.:           $LIBETSFIO_INCLUDE
  With libyaml:              $ax_have_yaml
   | C flags:                $LIB_YAML_CFLAGS
   | linked libraries:       $LIB_YAML_LIBS
  With GLib and co.:         $have_glib
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

if test x"$ac_use_mpi" = x"warn" -o x"$ac_use_etsfio" = x"warn" -o x"$ac_build_bindings" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$ac_use_etsfio" = x"warn" ; then
    echo "  ETSF_IO:                   disabled by tests"
  fi
  if test x"$ac_build_bindings" = x"warn" ; then
    echo "  C bindings:                disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
