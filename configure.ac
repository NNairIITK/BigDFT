dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - cluster], [0.9], [caliste@pcpm.ucl.ac.be], [bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)

dnl We custom the module naming scheme depending on the compiler
dnl and the platform.
dnl Not done yet
AC_SUBST(MODULE_EXT, "mod")

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl Defining fortran language for tests.
AC_LANG(Fortran)

dnl We defined CC for the ABINIT part.
AC_PROG_CC()

dnl If FC environement variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC()
ABI_PROG_FC()

ABI_PROG_FC()

dnl Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
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
AC_MSG_CHECKING([for builtin MPI capabilities])
ac_use_mpi="yes"
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--disable-mpi], [Don't use MPI (enable by default).]),
              ac_use_mpi=$enableval, ac_use_mpi="yes")
if test "$ac_use_mpi" = "yes"; then
  AC_DEFINE(USE_MPI)
  case "$FC" in
    mpi*)
      dnl The compiler will pass everything needed
      MPI_LIBS=
      MPI_INCLUDE=
      MPI_LDFLAGS=
      ;;
    *)
      dnl We set all flags by hand
      if test x"$MPI_NAME" = x ; then
        MPI_NAME=mpi
      fi
      MPI_LIBS="-l$MPI_NAME -lf$MPI_NAME"
      if ! test x"$MPI_INCLUDE" = x ; then
        MPI_INCLUDE=-I$MPI_INCLUDE
      fi
      if ! test x"$MPI_LDFLAGS" = x ; then
        MPI_LDFLAGS=-L$MPI_LDFLAGS
      fi
      ;;
  esac
  AC_SUBST(PARALLEL_OPTION, ".true.")
else
  AC_SUBST(PARALLEL_OPTION, ".false.")
fi
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_INCLUDE)
AC_SUBST(MPI_LDFLAGS)
AM_CONDITIONAL(USE_MPI, test "$ac_use_mpi" = "yes")
AC_MSG_RESULT([$ac_use_mpi])

dnl Use the simple or the optimised convolution code.
AC_MSG_CHECKING([for convolution optimisation])
ac_use_opti_convolut="yes"
AC_ARG_ENABLE(optimised-convolution, AS_HELP_STRING([--disable-optimised-convolution], [Use the simple convolution code (use the optimised code by default).]),
                         ac_use_opti_convolut=$enableval, ac_use_opti_convolut="yes")
AM_CONDITIONAL(USE_OPTI_CONVOLUT, test "$ac_use_opti_convolut" = "yes")
AC_MSG_RESULT([$ac_use_opti_convolut])



dnl If the library is built, we doesn't need to check for
dnl libraries to link to (Lapack and Blas).
AC_ARG_ENABLE(library-only, AS_HELP_STRING([--enable-library-only], [build the BigDFT library instead of the excecutable program (default = no).]), library=yes, library=no)
AM_CONDITIONAL(BUILD_LIBRARY, test "$library" = "yes")
AC_MSG_CHECKING([for library compilation only])
if test "$library" = "yes" ; then
  AC_MSG_RESULT([$library])
else
  AC_MSG_RESULT([$library])

dnl Get the lapack library
lapack=yes
AC_ARG_WITH(lapack, AS_HELP_STRING([--with-lapack], [Link with Lapack library (default = yes).]), lapack=$withval, lapack=no)
AC_ARG_WITH(lapack-path, AS_HELP_STRING([--with-lapack-path], [Give the path of the Lapack library (default = /usr/lib).]),
            ac_lapack_dir=$withval, ac_lapack_dir=)
if test -n "$ac_lapack_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_lapack_dir"
fi
withlapack=no
if test "$lapack" = "yes" ; then
  AC_CHECK_LIB(lapack, dsysv, withlapack=yes, withlapack=no)
fi

dnl Get the Blas library
blas=yes
AC_ARG_WITH(blas, AS_HELP_STRING([--with-blas], [Link with Blas library (default = yes).]), blas=$withval, blas=no)
AC_ARG_WITH(blas-path, AS_HELP_STRING([--with-blas-path], [Give the path of the Blas library (default = /usr/lib).]),
            ac_blas_dir=$withval, ac_blas_dir=)
if test -n "$ac_blas_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_blas_dir"
fi
withblas=no
if test "$blas" = "yes" ; then
  AC_CHECK_LIB(blas, dcopy, withblas=yes, withblas=no)
fi

dnl Substitute Lapack and Blas by another specific library
withlinalg=no
AC_ARG_WITH(ext-linalg, AS_HELP_STRING([--with-ext-linalg],
            [Give the name of the library replacing Blas and Lapack (default = none sepcified).]),
            ac_linalg=$withval, ac_linalg=)
AC_ARG_WITH(ext-linalg-path, AS_HELP_STRING([--with-ext-linalg-path],
            [Give the path of the Linalg library (default = /usr/lib).]),
            ac_linalg_dir=$withval, ac_linalg_dir=)
if test -n "$ac_linalg_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_linalg_dir"
fi
if test -n "$ac_linalg" ; then
  AC_CHECK_LIB($ac_linalg, dsysv, withlinalg=yes, withlinalg=no)
  LIBS="-l$ac_linalg $LIBS"
fi

dnl Fallbacks for linear algebra libraries
if [ test $withlinalg = "no" && test $withblas = "no" ] ; then
  AC_CHECK_LIB(blas, dcopy, withblas=yes, withblas=no)
  if test $withblas = "no" ; then
    AC_MSG_ERROR(["No blas library. Use --with-blas-path or check your Blas installation."])
  fi
  LIBS="-lblas $LIBS"
fi
if [ test $withlinalg = "no" && test $withlapack = "no" ] ; then
  AC_CHECK_LIB(lapack, dsysv, withlapack=yes, withlapack=no)
  if test $withlapack = "no" ; then
    AC_MSG_ERROR(["No lapack library. Use --with-lapack-path or check your Lapack installation."])
  fi
  LIBS="-llapack $LIBS"
fi

fi # End of the if library compilation


dnl default installation directories
moduledir=${includedir}/$fc_type
AC_SUBST(moduledir)

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
src/ABINIT-common/Makefile
src/PSolver/Makefile
src/parameters.h
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
test/Makefile
test/C/Makefile
test/H/Makefile
test/NaCl/Makefile
test/SiH4/Makefile
])
AC_OUTPUT

echo "
Basics:
  Fortran90 compiler:      $FC
  Fortran90 compiler name: $fc_type
  Fortran90 flags:         $FCFLAGS
  Linker flags:            $LDFLAGS
  Linked libraries:        $LIBS

Options:
  Library making only:     $library
  With MPI:                $ac_use_mpi
   | MPI implementation:   $MPI_NAME
   | include dir.:         $MPI_INCLUDE
   | Linker flags:         $MPI_LDFLAGS
   | linked libraries:     $MPI_LIBS
  With optimised conv.:    $ac_use_opti_convolut

Installation paths:
  Prefix:                  $prefix
  Exec prefix:             $exec_prefix
  Binaries:                $bindir
  Static libraries:        $libdir
  Fortran modules:         $moduledir
"
