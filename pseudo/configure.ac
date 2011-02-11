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


dnl Use the GPU computing (only NVIDIA cards with CUDA libraries)
AC_MSG_CHECKING([for GPU computing with CUDA libraries])
ac_use_cuda_gpu="no"
AC_ARG_ENABLE(cuda-gpu, AS_HELP_STRING([--enable-cuda-gpu], [Use CUDA implementation for GPU convolutions (disabled by default).]),
                         ac_use_cuda_gpu=$enableval, ac_use_cuda_gpu="no")
AM_CONDITIONAL(USE_CUDA_GPU, test "$ac_use_cuda_gpu" = "yes")
AC_MSG_RESULT([$ac_use_cuda_gpu])
if test x"$ac_use_cuda_gpu" = "xyes" ; then
  dnl Local CUDA usage (beside S_GPU).
  dnl Optional Cuda installation dir.
  AC_ARG_WITH(cuda-path, AS_HELP_STRING([--with-cuda-path], [give the path to the NVidia Cuda tools (default is /usr/local/cuda).]),
                         [ac_with_cuda_path=$withval], [ac_with_cuda_path="/usr/local/cuda"])
  CUDA_PATH=$ac_with_cuda_path
  AC_SUBST(CUDA_PATH)
  dnl Test the existence of things.
  AC_PATH_PROG([NVCC], [nvcc], [$ac_with_nvcc], [$PATH:${ac_with_cuda_path}/bin])
  AC_SUBST(NVCC)

  dnl Test the link on libcudart and libcublas.
  AC_LANG_PUSH(C)
  LIBS_OLD=$LIBS
  ac_cuda_lib_path="$ac_with_cuda_path/lib"
  LIBS="-L$ac_cuda_lib_path $LIBS"
  AC_CHECK_LIB(cublas, cublasSsymv, withlibcublas=yes, withlibcublas=no, [-lcudart])
  if test "$withlibcublas" = "no"; then
    ac_cuda_lib_path="$ac_with_cuda_path/lib64"
    LIBS="-L$ac_cuda_lib_path $LIBS_OLD"
    AC_CHECK_LIB(cublas, cublasScopy, withlibcublas=yes, withlibcublas=no, [-lcudart])
  fi
  if test "$withlibcublas" = "no"; then
    AC_MSG_ERROR(["No 'cublas' library found, link will fail."])
  fi
  LIBS=$LIBS_OLD
  AC_LANG_POP(C)

  dnl Add the flags.
  AC_ARG_WITH(cuda-cflags, AS_HELP_STRING([--with-nvcc-flags], [specify the flags for the NVidia Cuda Compiler (default is '-O3 --compiler-options -fno-strict-aliasing').]),
                         [ac_with_nvcc_flags=$withval],
                         [ac_with_nvcc_flags="no"])
  if test x"$ac_with_nvcc_flags" != x"no" ; then
    NVCC_FLAGS="$ac_with_nvcc_flags"
  else
    if test -z "NVCC_FLAGS" ; then
      NVCC_FLAGS="-O3 --compiler-options -fno-strict-aliasing"
    fi
  fi
  AC_SUBST(NVCC_FLAGS)
  LIBCUDA_LIBS="-L$ac_cuda_lib_path -lcudart -lcublas"
fi


dnl atom and pseudo executables
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_pseudo=$enableval, ac_build_bin_pseudo="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_pseudo" = "yes")
AC_MSG_CHECKING([for pseudo binaries build])
AC_MSG_RESULT([$ac_build_bin_pseudo])


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
  With Cuda GPU:             $ac_use_cuda_gpu
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
  With libXC:                $ac_use_libxc
   | internal built:         $ac_build_libxc
   | include dir.:           $LIBXC_INCLUDE

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
