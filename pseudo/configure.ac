AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Pseudo generation],[2.5],[Damien.Caliste@cea.fr],[pseudo])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
#	       ABI_PROG_FC()

dnl Defining fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

dnl Fortran 77 (to be removed)
AC_PROG_F77([], [Fortran 77])

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
ac_use_mpi="yes"
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--disable-mpi], [Don't use MPI (enabled by default).]),
              ac_use_mpi=$enableval, ac_use_mpi="yes")
if test "$ac_use_mpi" = "yes"; then
  dnl We set all flags by hand
  AC_ARG_WITH(mpi-include, AS_HELP_STRING([--with-mpi-include],
              [specify INCLUDE for MPI in case of non-built-in library. Use the -I before the path(es).]),
              ac_mpi_include=$withval, ac_mpi_include=$MPI_INCLUDE)
  AC_ARG_WITH(mpi-ldflags, AS_HELP_STRING([--with-mpi-ldflags],
              [specify LDFLAGS for MPI in case of non-built-in library. Use the -L before the path(es).]),
              ac_mpi_ldflags=$withval, ac_mpi_ldflags=$MPI_LDFLAGS)
  AC_ARG_WITH(mpi-libs, AS_HELP_STRING([--with-mpi-libs],
              [specify the libraries for MPI in case of non-built-in library. Use the -l before the name(s).]),
              ac_mpi_libs=$withval, ac_mpi_libs=$MPI_LIBS)

  LIBS_OLD=$LIBS
  LIBS="$LIBS $ac_mpi_libs"
  LDFLAGS_OLD=$LDFLAGS
  LDFLAGS="$LDFLAGS $ac_mpi_ldflags"
  FCFLAGS_OLD=$FCFLAGS
  FCFLAGS="$FCFLAGS $ac_mpi_include"

  AC_MSG_CHECKING([for MPI capabilities])
  AC_LINK_IFELSE([[program main
  include "mpif.h"

  call MPI_INIT
end]], withmpi=yes, withmpi=no)
  AC_MSG_RESULT([$withmpi])

  if test x"$withmpi" = x"yes" ; then
    dnl The compiler will pass everything needed
    MPI_LIBS=$ac_mpi_libs
    MPI_INCLUDE=$ac_mpi_include
    MPI_LDFLAGS=$ac_mpi_ldflags

    dnl Test for MPI2.
    AC_MSG_CHECKING([for MPI 2 capabilities])
    AC_LINK_IFELSE([[
program main
implicit none

integer :: ierr
real :: rho(3)

include "mpif.h"

call MPI_ALLREDUCE(MPI_IN_PLACE,rho,3,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
end
]], withmpi2=yes, withmpi2=no)
    AC_MSG_RESULT([$withmpi2])

    if test x"$withmpi2" = x"yes" ; then
      ac_with_mpi2="yes"
      AC_DEFINE([HAVE_MPI2], [1], [use MPI2 capabilities.])
    else
      ac_with_mpi2="no"
    fi

  else
    AC_MSG_WARN([MPI is not working and has been disabled.])
    AC_MSG_WARN([  If you don't want to see this warning,])
    AC_MSG_WARN([  set --disable-mpi in the configure line.])
    ac_use_mpi="Warning"
    MPI_LIBS=
    MPI_INCLUDE='-I$(top_builddir)/src/modules'
    MPI_LDFLAGS=
  fi

  LIBS="$LIBS_OLD"
  LDFLAGS="$LDFLAGS_OLD"
  FCFLAGS="$FCFLAGS_OLD"
else
  dnl We put an include in case of no MPI to have access to mpi fake.
  MPI_INCLUDE='-I$(top_builddir)/src/modules'
fi
AM_CONDITIONAL(USE_MPI, test "$ac_use_mpi" = "yes")
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_INCLUDE)
AC_SUBST(MPI_LDFLAGS)


dnl atom and pseudo executables
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_pseudo=$enableval, ac_build_bin_pseudo="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_pseudo" = "yes")
AC_MSG_CHECKING([for pseudo binaries build])
AC_MSG_RESULT([$ac_build_bin_pseudo])


dnl ---------------------------------------
dnl Test the linking part only if required.
if test "$ac_build_bin_pseudo" = "yes" ; then

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

dnl Test the libXC library
AC_ARG_WITH(libxc-path, AS_HELP_STRING([--with-libxc-path], [Give the path for an external libXC library, lib and module files (default = None).]),
            ac_libxc_dir=$withval, ac_libxc_dir=)
ac_use_libxc="no"
if test -n "$ac_libxc_dir" ; then
  LDFLAGS="$LDFLAGS -L$ac_libxc_dir"
  dnl Test the given implementation of libXC.
  AC_CHECK_LIB(xc, xc_f90_lda_vxc, withlibxc=yes, withlibxc=no)
  AC_CHECK_FILE($ac_libxc_dir/libxc_funcs_m.$capitalize_module_ext, funcs="yes", funcs="no")
  AC_CHECK_FILE($ac_libxc_dir/xc_f90_lib_m.$capitalize_module_ext, lib="yes", lib="no")
  AC_CHECK_FILE($ac_libxc_dir/xc_f90_types_m.$capitalize_module_ext, types="yes", types="no")
  if test "$withlibxc" = "yes" -a "$funcs" = "yes" -a "$lib" = "yes" -a "$types" = "yes"; then
    ac_use_libxc="yes"
    AC_SUBST(LIBXC_INCLUDE, "-I"$ac_libxc_dir)
    AC_SUBST(LIBXC_LIBS,"-L"$ac_libxc_dir" -lxc")
  fi
fi
FCFLAGS="$FCFLAGS $LIBXC_INCLUDE"
LIBS="$LIBXC_LIBS $LIBS"
dnl flags for F77 (to be removed)
FFLAGS="$FFLAGS $LIBXC_INCLUDE"

dnl Fallbacks for libXC library
if [ test $ac_use_libxc = "no" && test $ac_use_libxc = "no" ] ; then
  AC_MSG_ERROR(["No libxc library. Use --with-libxc-path or check your libxc installation."])
fi
fi # End of the if library compilation



dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
])

AC_OUTPUT

echo "
pseudo configure summary:
========================

Basics:
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Fortran90 compiler:        $FC
  Fortran90 compiler name:   $fc_type
  Fortran90 flags:           $FCFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS

Build:
  atom and pseudo binaries:  $ac_build_bin_pseudo

Options:
  With MPI:                  $ac_use_mpi
   | Include dir.:           $MPI_INCLUDE
   | Linker flags:           $MPI_LDFLAGS
   | Linked libraries:       $MPI_LIBS
   | MPI2 support:           $ac_with_mpi2
  With libXC:                $ac_use_libxc
   | include dir.:           $LIBXC_INCLUDE
   | linked libraries:       $LIBXC_LIBS

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"Warning" ; then
  echo "Warnings:
  MPI:                       disabled by tests
"
fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
sed -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
