dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - DFT over wavelets], [0.9], [Damien.Caliste@cea.fr], [bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions
AM_CONFIG_HEADER(config.h)

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl This is a hack for usage inside ABINIT.
if test -n "$FCFLAGS" ; then
  FCFLAGS=`echo $FCFLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
fi
if test -n "$F90FLAGS" ; then
  F90FLAGS=`echo $F90FLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
fi

dnl Defining fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

dnl We defined CC for the ABINIT part.
AC_PROG_CC()

dnl If FC environement variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC()

ABI_PROG_FC()

dnl Set if the module name are capitalized or not with the FC compiler.
AC_MSG_CHECKING([for module extension for compiler '$fc_type'])
case "$fc_type" in
  "pathscale")
     capitalize_module=yes
     capitalize_module_ext="mod" ;;
  "mipspro")
     capitalize_module=yes
     capitalize_module_ext="mod" ;;
  *)
     capitalize_module=no
     capitalize_module_ext="mod" ;;
esac
AC_MSG_RESULT([$capitalize_module_ext])

dnl We custom the module naming scheme depending on the compiler
dnl and the platform.
dnl Not done yet
AM_CONDITIONAL(CAPITALIZE, test "$capitalize_module" = "yes")
AC_SUBST(MODULE_EXT, $capitalize_module_ext)

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
AC_MSG_CHECKING([for builtin MPI capabilities])
ac_use_mpi="yes"
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--disable-mpi], [Don't use MPI (enable by default).]),
              ac_use_mpi=$enableval, ac_use_mpi="yes")
if test "$ac_use_mpi" = "yes"; then
  case "$FC" in
    *mpi*)
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
AM_CONDITIONAL(USE_MPI, test "$ac_use_mpi" = "yes")
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


dnl AM_CONDITIONAL for the different libraries, set default values.
dnl the common library for XC computation.
ac_build_lib_ABINIT_common="yes"
dnl the pure PSolver library.
ac_build_lib_PSolver="yes"
dnl the BigDFT library (without PSolver)
ac_build_lib_BigDFT="yes"
dnl the standalone executable (required linking tests)
ac_build_bin_cluster="yes"


dnl disable different build
dnl The ABINIT common XC library
AC_ARG_ENABLE(libxc, AS_HELP_STRING([--enable-libxc], [build the XC library from ABINIT code (default = yes).]), ac_build_lib_ABINIT_common=$enableval, ac_build_lib_ABINIT_common="yes")
AM_CONDITIONAL(BUILD_LIB_XC, test "$ac_build_lib_ABINIT_common" = "yes")
AC_MSG_CHECKING([for XC library internal build])
AC_MSG_RESULT([$ac_build_lib_ABINIT_common])
dnl The PSolver library
AC_ARG_ENABLE(libpoissonsolver, AS_HELP_STRING([--enable-libpoissonsolver], [build the Poisson's solver library (default = yes).]), ac_build_lib_PSolver=$enableval, ac_build_lib_PSolver="yes")
AM_CONDITIONAL(BUILD_LIB_PSOLVER, test "$ac_build_lib_PSolver" = "yes")
AC_MSG_CHECKING([for Poisson's solver internal build])
AC_MSG_RESULT([$ac_build_lib_PSolver])
dnl The BigDFT library
AC_ARG_ENABLE(libbigdft, AS_HELP_STRING([--enable-libbigdft], [build the BigDFT library (default = yes).]), ac_build_lib_BigDFT=$enableval, ac_build_lib_BigDFT="yes")
AM_CONDITIONAL(BUILD_LIB_BIGDFT, test "$ac_build_lib_BigDFT" = "yes")
AC_MSG_CHECKING([for BigDFT library build])
AC_MSG_RESULT([$ac_build_lib_BigDFT])
dnl The cluster executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--enable-binaries], [build the binaries (default = yes).]), ac_build_bin_cluster=$enableval, ac_build_bin_cluster="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_cluster" = "yes")
AC_MSG_CHECKING([for BigDFT binaries build])
AC_MSG_RESULT([$ac_build_bin_cluster])

dnl Required path for ABINIT module inclusion
xc_path="../ABINIT-commmon"
AC_ARG_WITH(xc-module-path, AS_HELP_STRING([--with-xc-module-path], [give the path to the module file defining fortran type for XC (see ABINIT) (default = ../ABINIT-common).]), xc_path=$withval, xc_path="../ABINIT-common")
if test "$ac_build_lib_ABINIT_common" = "yes" ; then
  xc_path="../ABINIT-common"
fi
AC_SUBST(XC_INCLUDE, "-I"$xc_path)

dnl Required path for PSolver module inclusion
dnl xc_path="../PSolver"
dnl AC_ARG_WITH(xc-module-path, AS_HELP_STRING([--with-xc-module-path], [give the path to the module file defining fortran type for XC (see ABINIT) (default = ../ABINIT-common).]), xc_path=$withval, xc_path="../ABINIT-commmon")
dnl AC_SUBST(XC_INCLUDE, $xc_path)



dnl ---------------------------------------
dnl Test the linking part only if required.
if test "$ac_build_bin_cluster" = "yes" ; then

dnl Substitute Lapack and Blas by another specific library
withlinalg=no
AC_ARG_WITH(ext-linalg, AS_HELP_STRING([--with-ext-linalg],
            [Give the name of the libraries replacing Blas and Lapack (default = none specified). Use the -l before the name(s).]),
            ac_linalg=$withval, ac_linalg=)
AC_ARG_WITH(ext-linalg-path, AS_HELP_STRING([--with-ext-linalg-path],
            [Give the path of the other linear algebra libraries (default = -L/usr/lib). Use the -L before the path(es).]),
            ac_linalg_dir=$withval, ac_linalg_dir=)
if test -n "$ac_linalg_dir" ; then
  LDFLAGS="$LDFLAGS $ac_linalg_dir"
fi
if test -n "$ac_linalg" ; then
#  AC_LANG_CONFTEST([[program main
#                     call dsysv
#                     end program main]])
  AC_MSG_CHECKING([for dsysv in $ac_linalg])
  LIBS_OLD=$LIBS
  LIBS="$ac_linalg $LIBS"
  AC_LINK_IFELSE([[program main
  call dsysv
end]], withlinalg=yes, withlinalg=no)
#  AC_CHECK_LIB($ac_linalg, dsysv, withlinalg=yes, withlinalg=no)
  if test "$withlinalg" = "no"; then
    LIBS=$LIBS_OLD
  fi
  AC_MSG_RESULT([$withlinalg])
fi

if test "$withlinalg" = "no" ; then
dnl Get the Blas library
blas=yes
AC_ARG_WITH(blas, AS_HELP_STRING([--with-blas], [Link with Blas library (default = yes).]), blas=$withval, blas=yes)
AC_ARG_WITH(blas-path, AS_HELP_STRING([--with-blas-path], [Give the path of the Blas library (default = /usr/lib).]),
            ac_blas_dir=$withval, ac_blas_dir=)
if test -n "$ac_blas_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_blas_dir"
fi
withblas=no
if test "$blas" = "yes" ; then
  AC_CHECK_LIB(blas, dcopy, withblas=yes, withblas=no)
  if test "$withblas" = "yes"; then
    LIBS="-lblas $LIBS"
  fi
fi

dnl Get the lapack library
lapack=yes
AC_ARG_WITH(lapack, AS_HELP_STRING([--with-lapack], [Link with Lapack library (default = yes).]), lapack=$withval, lapack=yes)
AC_ARG_WITH(lapack-path, AS_HELP_STRING([--with-lapack-path], [Give the path of the Lapack library (default = /usr/lib).]),
            ac_lapack_dir=$withval, ac_lapack_dir=)
if test -n "$ac_lapack_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_lapack_dir"
fi
withlapack=no
if test "$lapack" = "yes" ; then
  AC_CHECK_LIB(lapack, dsysv, withlapack=yes, withlapack=no)
  if test "$withlapack" = "yes"; then
    LIBS="-llapack $LIBS"
  fi
fi

fi #end $withlinalg = no

dnl Fallbacks for linear algebra libraries
if [ test $withlinalg = "no" && test $withblas = "no" ] ; then
  AC_MSG_ERROR(["No blas library. Use --with-blas-path or check your Blas installation."])
fi
if [ test $withlinalg = "no" && test $withlapack = "no" ] ; then
  AC_MSG_ERROR(["No lapack library. Use --with-lapack-path or check your Lapack installation."])
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
src/convolutions/Makefile
src/profiling/Makefile
src/parameters.h
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
tests/Makefile
tests/C/Makefile
tests/Ca2/Makefile
tests/H/Makefile
tests/H-spin/Makefile
tests/Li+/Makefile
tests/NaCl/Makefile
tests/SiH4/Makefile
])
AC_OUTPUT

echo "
Basics:
  Fortran90 compiler:      $FC
  Fortran90 compiler name: $fc_type
  Fortran90 flags:         $FCFLAGS
  Linker flags:            $LDFLAGS
  Linked libraries:        $LIBS

Build:
  Library XC (from ABINIT) $ac_build_lib_ABINIT_common
  Library PSolver          $ac_build_lib_PSolver
  Library BigDFT           $ac_build_lib_BigDFT
  Binaries (cluster...)    $ac_build_bin_cluster

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
