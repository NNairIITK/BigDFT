dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - cluster], [0.2], [caliste@pcpm.ucl.ac.be], [bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl If F77 environement variable is not already set, it
dnl looks for a Fortran77 compiler and set F77 variable.
dnl AC_PROG_F77()

dnl If FC environement variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC()

dnl Set the FCFLAGS and FFLAGS variable
if test -n "$FFLAGS"; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS"; then
  AC_SUBST(FCFLAGS, "-O2")
fi

dnl Using libtool to handle libraries
dnl AC_CHECK_TOOL(LIBTOOL, libtool, :)
dnl if test "$libtool" = ":" ; then
dnl   AC_MSG_ERROR(["No 'libtool' program found."])
dnl fi
dnl AC_PROG_LIBTOOL
dnl AC_DEFUN([LT_SUPPORTED_TAG], [])
dnl LT_SUPPORTED_TAG([C])
dnl LT_SUPPORTED_TAG([FC])
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

dnl Look if a binary must be compiled
ac_build_exec="yes"
AC_ARG_ENABLE(build-binary,
              AS_HELP_STRING([--disable-build-binary],
                             [Don't create a binary (enable by default).]),
              ac_build_exec=$enableval, ac_build_exec="yes")
AM_CONDITIONAL(BUILD_BINARY, test "$ac_build_exec" = "yes")

dnl Look if a binary must be compiled
ac_build_lib="no"
AC_ARG_ENABLE(build-library,
              AS_HELP_STRING([--enable-build-library],
                             [Don't create a library (disable by default).]),
              ac_build_lib=$enableval, ac_build_lib="no")
AM_CONDITIONAL(BUILD_LIBRARY, test "$ac_build_lib" = "yes")

if [ test "$ac_build_lib" = "no" && test "$ac_build_exec" = "no" ] ; then
  AC_MSG_ERROR(["Nothing to build."])
fi

dnl Enable or not MPI.
AC_MSG_CHECKING([for builtin MPI capabilities])
ac_use_mpi="yes"
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--disable-mpi], [Don't use MPI (enable by default).]),
              ac_use_mpi=$enableval, ac_use_mpi="yes")
if test "$ac_use_mpi" = "yes"; then
  AC_DEFINE(USE_MPI)
  AC_SUBST(PARALLEL_OPTION, ".true.")
else
  AC_SUBST(PARALLEL_OPTION, ".false.")
fi
if ! test x"$MPI_INCLUDE" = x ; then
  MPI_INCLUDE=-I$MPI_INCLUDE
fi
AC_SUBST(MPI_INCLUDE)
if ! test x"$MPI_LDFLAGS" = x ; then
  MPI_LDFLAGS=-L$MPI_LDFLAGS
fi
AC_SUBST(MPI_LDFLAGS)
AM_CONDITIONAL(USE_MPI, test "$ac_use_mpi" = "yes")
if test "$ac_use_mpi" = "yes" ; then
  AC_SUBST(MPI_SOURCES, "")
else
  AC_SUBST(MPI_SOURCES, "  include \"MPIfake.f90\"")
fi
AC_MSG_RESULT([$ac_use_mpi])

dnl Use the simple or the optimised convolution code.
AC_MSG_CHECKING([for convolution optimisation])
ac_use_opti_convolut="yes"
AC_ARG_ENABLE(optimised-convolution, AS_HELP_STRING([--disable-optimised-convolution], [Use the simple convolution code (use the optimised code by default).]),
                         ac_use_opti_convolut=$enableval, ac_use_opti_convolut="yes")
AM_CONDITIONAL(USE_OPTI_CONVOLUT, test "$ac_use_opti_convolut" = "yes")
AC_MSG_RESULT([$ac_use_opti_convolut])
if test "$ac_use_opti_convolut" = "yes" ; then
  AC_SUBST(CONVOLUT_SOURCES, "convolut_optim_omp.f90")
else
  AC_SUBST(CONVOLUT_SOURCES, "convolut_simple.f90")
fi

dnl Get the lapack library
AC_ARG_WITH(lapack, AS_HELP_STRING([--with-lapack], [Give the filename of the Lapack library (default = /usr/lib/liblapack.a).]),
            ac_lapack_dir=$withval, ac_lapack_dir="/usr/lib/liblapack.a")
AC_SUBST(LAPACK_LIB, $ac_lapack_dir)

dnl Get the Blas library
AC_ARG_WITH(blas, AS_HELP_STRING([--with-blas], [Give the filename of Blas library (default = /usr/lib/libblas.a).]),
            ac_blas_dir=$withval, ac_blas_dir="/usr/lib/libblas.a")
AC_SUBST(BLAS_LIB, $ac_blas_dir)

dnl variables for autoconf
modulelibdir="$libdir"
AC_SUBST(modulelibdir)
modulemoddir="$includedir"
AC_SUBST(modulemoddir)

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
src/ABINIT-common/Makefile
src/parameters.h
lib/Makefile
lib/module_bigdft.f90
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
])
AC_OUTPUT

echo "
Basics:
  Prefix:                 $prefix
  Fortran90 compiler:     $FC
  Fortran90 flags:        $FCFLAGS
  Linker flags:           $LDFLAGS
  Binary build:           $ac_build_exec
  Library build:          $ac_build_lib

Options:
  With MPI:               $ac_use_mpi
   | include dir.:        $MPI_INCLUDE
   | lib. dir.:           $MPI_LDFLAGS
  With optimised conv.:   $ac_use_opti_convolut
  Blas lib.:              $ac_blas_dir
  Lapack lib.:            $ac_lapack_dir
"
