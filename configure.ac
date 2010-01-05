

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - DFT over wavelets],[1.3.0-dev],[Damien.Caliste@cea.fr],[bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl This is a hack for usage inside ABINIT.
if test -n "$FCFLAGS" ; then
  FCFLAGS=`echo $FCFLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
fi
if test -n "$F90FLAGS" ; then
  F90FLAGS=`echo $F90FLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
fi

dnl We defined CC for the ABINIT part.
AC_PROG_CC()
AC_PROG_CXX()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

#	       ABI_PROG_FC()

dnl Defining fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

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

dnl Enable ndebug variable (disabled by default)
ac_use_debug="no"
AC_MSG_CHECKING([for debugging options])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug], [Debugging options (disabled by default)]),
              ac_use_debug=$enableval, ac_use_debug="no")
if test x"$ac_use_debug" = x"yes" ; then
  AC_DEFINE([HAVE_DEBUG], [1], [compile the code with debugging options])
fi
AC_MSG_RESULT([$ac_use_debug])

dnl Setup memory limit value (disabled by default)
ac_with_memory_limit="no"
AC_MSG_CHECKING([for memory limit in GB])
AC_ARG_WITH(memory-limit, AS_HELP_STRING([--with-memory-limit], [Memory limit value in GB. It stops EVERYTHING if some process passes such limit (disabled by default).]),
              ac_with_memory_limit=$withval, ac_with_memory_limit="no")
if test "${ac_with_memory_limit//\.*/}" -ge "0" 2> /dev/null ; then
  ac_memory_limit=$ac_with_memory_limit
  AC_DEFINE_UNQUOTED([MEMORY_LIMIT], [$ac_memory_limit], [Value in GB to stop processing.])
else
  ac_memory_limit="undefined"
fi
AC_MSG_RESULT([$ac_memory_limit])


dnl Use the simple or the optimised convolution code.
AC_MSG_CHECKING([for convolution optimisation])
ac_use_opti_convolut="yes"
AC_ARG_ENABLE(optimised-convolution, AS_HELP_STRING([--disable-optimised-convolution], [Use the simple convolution code (use the optimised code by default).]),
                         ac_use_opti_convolut=$enableval, ac_use_opti_convolut="yes")
AM_CONDITIONAL(USE_OPTI_CONVOLUT, test "$ac_use_opti_convolut" = "yes")
AC_MSG_RESULT([$ac_use_opti_convolut])


dnl Use the minima hopping method (build the executable global)
AC_MSG_CHECKING([for minima hopping])
ac_build_bin_minima_hopping="no"
AC_ARG_ENABLE(minima-hopping, AS_HELP_STRING([--enable-minima-hopping], [Compile the minima hopping code (need the file src/global.f90, not compiled by default).]),
                         ac_build_bin_minima_hopping=$enableval, ac_build_bin_minima_hopping="no")
AC_MSG_RESULT([$ac_build_bin_minima_hopping])
if test x"$ac_build_bin_minima_hopping" = "xyes" ; then
  dnl Check if src/global.f90 does exist.
  AC_CHECK_FILE(${srcdir}/src/global.f90, ac_have_global="yes", ac_have_global="no")
  if test x"$ac_have_global" = x"no" ; then
     AC_MSG_WARN(["No src/global.f90. This file is not freely available, minima hopping is disabled."])
     ac_build_bin_minima_hopping="no"
  fi
fi
AM_CONDITIONAL(USE_MINIMA_HOPPING, test "$ac_build_bin_minima_hopping" = "yes")


dnl Use the GPU computing (only NVIDIA cards with CUDA libraries)
AC_MSG_CHECKING([for GPU computing with CUDA libraries])
ac_use_cuda_gpu="no"
AC_ARG_ENABLE(cuda-gpu, AS_HELP_STRING([--enable-cuda-gpu], [Use CUDA implementation for GPU convolutions (disabled by default).]),
                         ac_use_cuda_gpu=$enableval, ac_use_cuda_gpu="no")
AM_CONDITIONAL(USE_CUDA_GPU, test "$ac_use_cuda_gpu" = "yes")
AC_MSG_RESULT([$ac_use_cuda_gpu])
if test x"$ac_use_cuda_gpu" = "xyes" ; then
  dnl Add no-second-underscore flag.
  dnl if `echo $FCFLAGS | grep -qvs -- "-fno-second-underscore"` ; then
  dnl    FCFLAGS=$FCFLAGS" -fno-second-underscore"
  dnl fi
  dnl Instead of second underscoring, test the compiler behaviour.
  AC_FC_WRAPPERS()
  dnl Optional Cuda installation dir.
  AC_ARG_WITH(cuda-path, AS_HELP_STRING([--with-cuda-path], [give the path to the NVidia Cuda tools (default is /usr/local/cuda).]),
                         [ac_with_cuda_path=$withval], [ac_with_cuda_path="/usr/local/cuda"])
  CUDA_PATH=$ac_with_cuda_path
  AC_SUBST(CUDA_PATH)
  dnl Test the existence of things.
  AC_PATH_PROG([NVCC], [nvcc], [$ac_with_nvcc], [$PATH:${ac_with_cuda_path}/bin])
  AC_SUBST(NVCC)

  dnl Add the flags.
  AC_ARG_WITH(cuda-cflags, AS_HELP_STRING([--with-cuda-cflags], [specify the flags for the NVidia Cuda Compiler (default is '-O3 --compiler-options -fno-strict-aliasing').]),
                         [ac_with_cuda_cflags=$withval],
                         [ac_with_cuda_cflags="no"])
  if test x"$ac_with_cuda_cflags" != x"no" ; then
    CUDA_FLAGS="$ac_with_cuda_cflags"
  else
    if test -z "$CUDA_FLAGS" ; then
      CUDA_FLAGS="-O3 --compiler-options -fno-strict-aliasing"
    fi
  fi
  AC_SUBST(CUDA_FLAGS)

  dnl Test if Cuda utils is installed.
  AC_ARG_WITH(lib-cutils, AS_HELP_STRING([--with-lib-cutils], [specifiy the directory for libcutil library (default is '$(top_srcdir)/src/lib/CUDA').]),
                         [ac_with_lib_cutils=$withval],
                         [ac_with_lib_cutils="no"])
  if test x"$ac_with_lib_cutils" = x"no" ; then
    LIB_CUTIL_PATH='$(top_srcdir)/src/lib/CUDA'
  else
    LIB_CUTIL_PATH=$ac_with_lib_cutils
  fi  
  AC_SUBST(LIB_CUTIL_PATH)

  if test -n "$GCC41HACK" ; then
    INC_GCC41CUDA="--compiler-bindir="$GCC41HACK
  else
    INC_GCC41CUDA="--compiler-bindir=/usr/bin"
  fi
  AC_SUBST(INC_GCC41CUDA)
fi


dnl AM_CONDITIONAL for the different libraries, set default values.
dnl the common routines from ABINIT.
ac_build_lib_ABINIT="yes"
dnl the pure PSolver library.
ac_build_lib_PSolver="yes"
dnl the BigDFT library (without PSolver)
ac_build_lib_BigDFT="yes"
dnl the standalone executable (required linking tests)
ac_build_bin_cluster="yes"


dnl disable different build
dnl The ABINIT common XC library
AC_ARG_ENABLE(libabinit, AS_HELP_STRING([--disable-libabinit], [build the ABINIT library from ABINIT code (default = yes).]), ac_build_lib_ABINIT=$enableval, ac_build_lib_ABINIT="yes")
AM_CONDITIONAL(BUILD_LIB_ABINIT, test "$ac_build_lib_ABINIT" = "yes")
AC_MSG_CHECKING([for ABINIT library internal build])
AC_MSG_RESULT([$ac_build_lib_ABINIT])
if test x"$ac_build_lib_ABINIT" = x"yes" ; then
  xc_path='$(top_builddir)/src/PSolver/ABINIT-common'
else
  dnl Required path for ABINIT module inclusion
  AC_ARG_WITH(xc-module-path, AS_HELP_STRING([--with-xc-module-path],
              [give the path to the module file defining XC (default = $(top_builddir)/src/PSolver/ABINIT-common).]),
              xc_path=$withval, xc_path=)
fi
AC_SUBST(XC_INCLUDE, "-I"$xc_path)
dnl The PSolver library
AC_ARG_ENABLE(libpoissonsolver, AS_HELP_STRING([--disable-libpoissonsolver], [build the Poisson's solver library (default = yes).]), ac_build_lib_PSolver=$enableval, ac_build_lib_PSolver="yes")
AM_CONDITIONAL(BUILD_LIB_PSOLVER, test "$ac_build_lib_PSolver" = "yes")
AC_MSG_CHECKING([for Poisson's solver internal build])
AC_MSG_RESULT([$ac_build_lib_PSolver])
dnl The BigDFT library
AC_ARG_ENABLE(libbigdft, AS_HELP_STRING([--disable-libbigdft], [build the BigDFT library (default = yes).]), ac_build_lib_BigDFT=$enableval, ac_build_lib_BigDFT="yes")
AM_CONDITIONAL(BUILD_LIB_BIGDFT, test "$ac_build_lib_BigDFT" = "yes")
AC_MSG_CHECKING([for BigDFT library build])
AC_MSG_RESULT([$ac_build_lib_BigDFT])
dnl The cluster executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_cluster=$enableval, ac_build_bin_cluster="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_cluster" = "yes")
AC_MSG_CHECKING([for BigDFT binaries build])
AC_MSG_RESULT([$ac_build_bin_cluster])


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

dnl Test the libXC library
internallibxc="yes"
AC_ARG_ENABLE(internal-libxc, AS_HELP_STRING([--disable-internal-libxc], [Do not build and link with internal libXC library (default = yes).]), internallibxc=$enval, internallibxc="yes")
AC_ARG_WITH(libxc-path, AS_HELP_STRING([--with-libxc-path], [Give the path for an external libXC library, lib and module files (default = None).]),
            ac_libxc_dir=$withval, ac_libxc_dir=)
ac_use_libxc="no"
ac_build_libxc="no"
if test "$internallibxc" = "yes" ; then
  ac_build_libxc="yes"
  ac_use_libxc="yes"
  AC_SUBST(LIBXC_INCLUDE, "-I\$(top_builddir)"/libXC/src)
  LDFLAGS="$LDFLAGS -L\$(top_builddir)/libXC/src"
else
  if test -n "$ac_libxc_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_libxc_dir/lib"
    
    dnl Test the given implementation of libXC.
    AC_CHECK_LIB(xc, xc_f90_lda_vxc, withlibxc=yes, withlibxc=no)
    AC_CHECK_FILE($ac_libxc_dir/include/libxc_funcs_m.$capitalize_module_ext, funcs="yes", funcs="no")
    AC_CHECK_FILE($ac_libxc_dir/include/xc_f90_lib_m.$capitalize_module_ext, lib="yes", lib="no")
    AC_CHECK_FILE($ac_libxc_dir/include/xc_f90_types_m.$capitalize_module_ext, types="yes", types="no")
    if test "$withlibxc" = "yes" -a "$funcs" = "yes" -a "$lib" = "yes" -a "$types" = "yes"; then
      ac_use_libxc="yes"
      AC_SUBST(LIBXC_INCLUDE, "-I"$ac_libxc_dir/include)
    fi
  fi
fi
AM_CONDITIONAL(BUILD_LIBXC, test x"$ac_build_libxc" = x"yes")
if test x"$ac_use_libxc" = x"yes" ; then
  AC_DEFINE([HAVE_LIBXC], [1], [compile the libXC section in the ABINIT code.])
  LIBS="-lxc $LIBS"
fi

dnl Test for the documentation.
AC_CHECK_PROG(ac_have_pdflatex, [pdflatex], [yes])
AM_CONDITIONAL(BUILD_PDFDOC, test x"$ac_have_pdflatex" = x"yes")
dnl default installation directories
AC_ARG_WITH(moduledir,
            AS_HELP_STRING([--with-moduledir],
                           [installation directory for module files [[INCLUDEDIR/FC_TYPE]]]),
              ac_moduledir=$withval, ac_moduledir="no")
if test x"$ac_moduledir" != x"no" ; then
  moduledir=$ac_moduledir
else
  moduledir=${includedir}/${fc_type}
fi
AC_SUBST(moduledir)
docdir="${docdir}"
AC_SUBST(docdir)

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
src/lib/Makefile
src/lib/fft/Makefile
src/modules/Makefile
src/profiling/Makefile
src/PSolver/ABINIT-common/Makefile
src/PSolver/Makefile
src/convolutions/Makefile
src/CUDA/Makefile
src/CUDA/S_GPU/Makefile
src/art/Makefile
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
tests/Makefile
tests/C/Makefile
tests/Ca2/Makefile
tests/H/Makefile
tests/H-spin/Makefile
tests/Li+/Makefile
tests/Mn2/Makefile
tests/NaCl/Makefile
tests/SiH4/Makefile
tests/H2O-CP2K/Makefile
tests/SiBulk/Makefile
tests/H2O-slab/Makefile
tests/Davidson-SiH4/Makefile
tests/GPU/Makefile
tests/MINHOP/Makefile
tests/H2-freq/Makefile
tests/NEB/Makefile
tests/NEB/input
tests/NEB/NEB_include.sh
tests/XabsCheb/Makefile
])
if test x"$ac_build_libxc" = x"yes" ; then
  # Call libXC ./configure script recursively.
  ac_configure_args=$ac_configure_args" --enable-fortran"
  AC_CONFIG_SUBDIRS([libXC])
fi
AC_OUTPUT

echo "
Basics:
  Fortran90 compiler:        $FC
  Fortran90 compiler name:   $fc_type
  Fortran90 flags:           $FCFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS

Build:
  Library ABINIT:            $ac_build_lib_ABINIT
  Library PSolver:           $ac_build_lib_PSolver
  Library BigDFT:            $ac_build_lib_BigDFT
  Main binaries (cluster...):$ac_build_bin_cluster
  Minima hopping binary:     $ac_build_bin_minima_hopping

Options:
  Debug version:             $ac_use_debug		
  With MPI:                  $ac_use_mpi
   | Include dir.:           $MPI_INCLUDE
   | Linker flags:           $MPI_LDFLAGS
   | Linked libraries:       $MPI_LIBS
   | MPI2 support:           $ac_with_mpi2
  With optimised conv.:      $ac_use_opti_convolut
  With Cuda GPU conv.:       $ac_use_cuda_gpu
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
   | GCC41HACK (gcc path):   $GCC41HACK
   | ext. cutil library:     $ac_with_lib_cutils
   |      cutil path:        $LIB_CUTIL_PATH
  With libXC:                $ac_use_libxc
   | internal built:         $ac_build_libxc
   | include dir.:           $LIBXC_INCLUDE

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
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
