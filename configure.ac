

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - DFT over wavelets],[1.5-dev.5],[Damien.Caliste@cea.fr],[bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE

dnl This is a hack for usage inside ABINIT.
dnl if test -n "$FCFLAGS" ; then
dnl   FCFLAGS=`echo $FCFLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi
dnl if test -n "$F90FLAGS" ; then
dnl   F90FLAGS=`echo $F90FLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi

dnl We defined CC for the ABINIT part.
AC_PROG_CC()
AC_PROG_CXX()

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
  MEMORY_LIMIT=$ac_with_memory_limit
else
  ac_with_memory_limit="undefined"
  MEMORY_LIMIT=0
fi
AC_SUBST(MEMORY_LIMIT)
AC_MSG_RESULT([$ac_with_memory_limit])


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

dnl Use the dgemmsy
ac_use_dgemmsy="no"
AC_ARG_ENABLE(dgemmsy, AS_HELP_STRING([--enable-dgemmsy], [Use dgemmsy (disabled by default).]),
                         ac_use_dgemmsy=$enableval, ac_use_dgemmsy="no")
if test x"$ac_use_dgemmsy" = "xyes" ; then
  AC_LANG_PUSH(C)

  CPPFLAGS_SVG=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS -msse3"
  AC_CHECK_HEADER([pmmintrin.h],
                  [ac_use_dgemmsy="yes"],
                  [ac_use_dgemmsy="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  AC_LANG_POP(C)
  if test "$ac_use_dgemmsy" = "no"; then
    AC_MSG_WARN(["No SSE header found, dgemmsy disabled."])
  fi
fi
AM_CONDITIONAL(USE_DGEMMSY, test "$ac_use_dgemmsy" = "yes")


dnl build the pseudo generation (build the executables atom and pseudo)
AC_MSG_CHECKING([for pseudo generation])
ac_build_bin_pseudo="no"
AC_ARG_ENABLE(pseudo, AS_HELP_STRING([--enable-pseudo], [Compile the pseudo code (not compiled by default).]),
                         ac_build_bin_pseudo=$enableval, ac_build_bin_pseudo="no")
AM_CONDITIONAL(BUILD_BINARIES_PSEUDO, test "$ac_build_bin_pseudo" = "yes")
AC_MSG_RESULT([$ac_build_bin_pseudo])


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

dnl Use the OpenCL support
ac_use_ocl="no"
AC_ARG_ENABLE(opencl, AS_HELP_STRING([--enable-opencl], [Use OpenCL implementation for GPU convolutions (disabled by default).]),
                         ac_use_ocl=$enableval, ac_use_ocl="no")
AC_MSG_CHECKING([for OpenCL support])
AC_MSG_RESULT([$ac_use_ocl])
if test x"$ac_use_ocl" = "xyes" ; then
  AC_LANG_PUSH(C)
  dnl Optional OpenCL installation dir.
  AC_ARG_WITH(ocl-path, AS_HELP_STRING([--with-ocl-path], [give the path to the OpenCL library (default is /usr).]),
                         [ac_ocl_path=$withval], [ac_ocl_path="/usr"])
  dnl Test the link on libOpenCL.
  dnl look for the header.
  CPPFLAGS_SVG=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS -I$ac_ocl_path/include"
  AC_CHECK_HEADER([CL/cl.h],
                  [ac_ocl_header="yes"],
                  [ac_ocl_header="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  if test "$ac_ocl_header" = "no"; then
    AC_MSG_WARN(["No 'OpenCL' header found, OpenCL parts disabled."])
    ac_use_ocl="no"
  else
    OCL_INCLUDE_PATH=$ac_ocl_path"/include"
    AC_SUBST(OCL_INCLUDE_PATH)
    dnl Now test for the library.
    LIBS_OLD=$LIBS
    LIBS="-L$ac_ocl_path/lib $LIBS"
    AC_CHECK_LIB(OpenCL, clCreateContext, withlibocl=yes, withlibocl=no)
    if test "$withlibocl" = "no"; then
      AC_MSG_WARN(["No 'OpenCL' library found, OpenCL parts disabled."])
      ac_use_ocl="no"
    else
dnl      LIBOCL_LIBS="-L$ac_ocl_path/lib -lOpenCL -lm -lrt -lstdc++"
      LIBOCL_LIBS="-L$ac_ocl_path/lib -lOpenCL -lm -lrt"
    fi
    LIBS=$LIBS_OLD
  fi
  AC_LANG_POP(C)
fi
AM_CONDITIONAL(USE_OCL, test "$ac_use_ocl" = "yes")

if test "$ac_use_ocl" = "yes" -o "$ac_use_cuda_gpu" = "yes" -o "$ac_use_dgemmsy" = "yes"; then
  dnl Instead of second underscoring, test the compiler behaviour.
  AC_FC_WRAPPERS()
fi

dnl Use ETSF_IO
ac_use_etsfio=no
AC_ARG_WITH(etsf-io, AS_HELP_STRING([--with-etsf-io], [Link with ETSF_IO library (default = no).]), etsfio=$withval, etsfio=no)
AC_ARG_WITH(etsf-io-path, AS_HELP_STRING([--with-etsf-io-path], [Give the path of the ETSF_IO installation (default = /usr).]),
            ac_etsfio_dir=$withval, ac_etsfio_dir=)
AC_ARG_WITH(netcdf-path, AS_HELP_STRING([--with-netcdf-path], [Give the path to NetCDF (required by ETSF_IO) (default = /usr).]),
            ac_netcdf_dir=$withval, ac_netcdf_dir=)
AC_ARG_WITH(netcdf-libs, AS_HELP_STRING([--with-netcdf-libs], [Give the library to link with NetCDF (required by ETSF_IO) (default = -lnetcdff -lnetcdf).]),
            ac_netcdf_libs=$withval, ac_netcdf_libs="-lnetcdff -lnetcdf")
if test "$etsfio" = "yes" ; then
   LDFLAGS_SVG="$LDFLAGS"
   LIBS_SVG="$LIBS"
   FCFLAGS_SVG="$FCFLAGS"
   if test -n "$ac_etsfio_dir" ; then
      LDFLAGS="$LDFLAGS -L$ac_etsfio_dir/lib"
      FCFLAGS="$FCFLAGS -I$ac_etsfio_dir/include"
   fi
   if test -n "$ac_netcdf_dir" -a x"$ac_netcdf_dir" != x"$ac_etsfio_dir" ; then
      LDFLAGS="$LDFLAGS -L$ac_netcdf_dir/lib"
      FCFLAGS="$FCFLAGS -I$ac_netcdf_dir/include"
   fi
   LIBS="$LIBS -letsf_io $ac_netcdf_libs"
   AC_MSG_CHECKING([for ETSF_IO library])
   AC_LINK_IFELSE([[
program main
  use etsf_io
  
  type(etsf_groups_flags) :: groups
  type(etsf_dims) :: dims
  logical :: lstat
  type(etsf_io_low_error) :: error_data

  call etsf_io_data_init("test", groups, dims, "test", "", lstat, error_data)
end]], ac_use_etsfio=yes, ac_use_etsfio=no)
   AC_MSG_RESULT([$ac_use_etsfio])
   LIBS="$LIBS_SVG"
   FCFLAGS="$FCFLAGS_SVG"
   if test "$ac_use_etsfio" = "yes"; then
      ac_libetsfio_libs="-letsf_io -lnetcdff -lnetcdf"
      if test -n "$ac_etsfio_dir" ; then
         AC_SUBST(LIBETSFIO_INCLUDE, "-I$ac_etsfio_dir/include")
      fi
      LIBETSFIO_LIBS="-letsf_io -lnetcdff -lnetcdf"
   fi
fi
AM_CONDITIONAL(HAVE_ETSF_IO, test "$ac_use_etsfio" = "yes")

dnl Test libraries required for benching.
AC_LANG_PUSH(C++)
AC_CHECK_LIB(config++, config_init, withconfig=yes, withconfig=no)
AM_CONDITIONAL(HAVE_LIBCONFIG, test "$withconfig" = "yes")
AC_LANG_POP(C++)

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
  ac_libxc_dir="$top_builddir"/libXC/src
  AC_SUBST(LIBXC_DIR, "\$(top_builddir)"/libXC/src)
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
      AC_SUBST(LIBXC_DIR, $ac_libxc_dir/include)
      AC_SUBST(LIBXC_INCLUDE, "-I"$ac_libxc_dir/include)
    fi
  fi
fi
AM_CONDITIONAL(BUILD_LIBXC, test x"$ac_build_libxc" = x"yes")
if test x"$ac_use_libxc" = x"yes" ; then
  AC_DEFINE([HAVE_LIBXC], [1], [use libXC in BigDFT.])
  LIBXC_LIBS="-lxc"
fi

dnl Test the ABINIT library
internallibabinit="yes"
AC_ARG_ENABLE(internal-libabinit, AS_HELP_STRING([--disable-internal-libabinit], [Do not build and link with internal libABINIT library (default = yes).]), internallibabinit=$enval, internallibabinit="yes")
AC_ARG_WITH(libabinit-path, AS_HELP_STRING([--with-libabinit-path], [Give the path for an external libABINIT library, lib and module files (default = None).]),
            ac_libabinit_dir=$withval, ac_libabinit_dir=)
ac_use_libabinit="no"
ac_build_libabinit="no"
if test "$internallibabinit" = "yes" ; then
  ac_build_libabinit="yes"
  ac_use_libabinit="yes"
  AC_SUBST(LIBABINIT_INCLUDE, "-I\$(top_builddir)"/libABINIT/src)
  LDFLAGS="$LDFLAGS -L\$(top_builddir)/libABINIT/src"
else
  if test -n "$ac_libabinit_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_libabinit_dir/lib"
    
    dnl Test the given implementation of libabinit.
    AC_CHECK_LIB(abinit, symbrav, withlibabinit=yes, withlibabinit=no)
    AC_CHECK_FILE($ac_libabinit_dir/include/ab6_moldyn.$capitalize_module_ext, moldyn="yes", moldyn="no")
    AC_CHECK_FILE($ac_libabinit_dir/include/ab6_symmetry.$capitalize_module_ext, sym="yes", sym="no")
    AC_CHECK_FILE($ac_libabinit_dir/include/libxc_functionals.$capitalize_module_ext, libxc="yes", libxc="no")
    if test "$withlibabinit" = "yes" -a "$moldyn" = "yes" -a "$sym" = "yes" -a "$libxc" = "yes"; then
      ac_use_libabinit="yes"
      AC_SUBST(LIBABINIT_INCLUDE, "-I"$ac_libabinit_dir/include)
    fi
  fi
fi
AM_CONDITIONAL(BUILD_LIBABINIT, test x"$ac_build_libabinit" = x"yes")
if test x"$ac_use_libabinit" = x"yes" ; then
  AC_DEFINE([HAVE_LIBABINIT], [1], [use libABINIT in BigDFT.])
  LIBABINIT_LIBS="-labinit"
fi

ac_build_libsgpu="no"
ac_use_libsgpu="no"
if test x"$ac_use_cuda_gpu" = "xyes" ; then
  dnl Test S_GPU library.
  internalsgpu="yes"
  AC_ARG_ENABLE(internal-libsgpu, AS_HELP_STRING([--disable-internal-libsgpu], [Do not build and link with internal S_GPU library (default = yes).]), internalsgpu=$enval, internalsgpu="yes")
  AC_ARG_WITH(libsgpu-path, AS_HELP_STRING([--with-libsgpu-path], [Give the path for an external S_GPU library, lib and module files (default = None).]),
            ac_libsgpu_dir=$withval, ac_libsgpu_dir=)
  if test "$internalsgpu" = "yes" ; then
    ac_build_libsgpu="yes"
    ac_use_libsgpu="yes"
    AC_SUBST(LIBSGPU_INCLUDE, "-I\$(top_srcdir)"/S_GPU/src)
  else
    if test -n "$ac_libsgpu_dir" ; then
      LDFLAGS="$LDFLAGS -L$ac_libsgpu_dir/lib"
    
      dnl Test the given implementation of sgpu.
      AC_CHECK_LIB(s_gpu, sg_create_stream, withsgpu=yes, withsgpu=no)
      AC_CHECK_FILE($ac_libsgpu_dir/include/s_gpu.h, header="yes", header="no")
      if test "$withsgpu" = "yes" -a "$header" = "yes" ; then
        ac_use_libsgpu="yes"
        AC_SUBST(LIBSGPU_INCLUDE, "-I"$ac_libsgpu_dir/include)
      fi
    fi
  fi
  if test x"$ac_use_libsgpu" = x"yes" ; then
    AC_DEFINE([HAVE_LIBSGPU], [1], [use S_GPU in BigDFT.])
    if test "$internalsgpu" = "yes" ; then
      LIBSGPU_LIBS="\$(top_builddir)/S_GPU/src/.libs/libs_gpu.a -lstdc++"
    else
      LIBSGPU_LIBS="-ls_gpu -lstdc++"
    fi
  fi
fi
AM_CONDITIONAL(BUILD_LIBSGPU, test x"$ac_build_libsgpu" = x"yes")

dnl Agregate all the optional libs.
LIBS="$LIBABINIT_LIBS $LIBXC_LIBS $LIBSGPU_LIBS $LIBCUDA_LIBS $LIBOCL_LIBS $LIBETSFIO_LIBS $LIBS"

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
libABINIT/src/Makefile
src/Makefile
src/modules/Makefile
src/modules/configure.inc
src/profiling/Makefile
src/PSolver/Makefile
src/convolutions/Makefile
src/CUDA/Makefile
src/OpenCL/Makefile
src/dgemmsy/Makefile
src/convolutions-c/Makefile
src/art/Makefile
utils/Makefile
utils/PSPfiles/Makefile
pseudo/Makefile
pseudo/src/Makefile
pseudo/tests/Makefile
doc/Makefile
tests/Makefile
tests/fft/Makefile
tests/PSolver/Makefile
tests/H2-freq/Makefile
tests/MINHOP/Makefile
tests/NEB/Makefile
tests/NEB/input
tests/NEB/NEB_include.sh
tests/XabsCheb/Makefile
tests/SPLSAD/Makefile
tests/DFT/Makefile
tests/OpenCL/Makefile
tests/CUDA/Makefile
])

# Call libXC ./configure script recursively.
# Call S_GPU ./configure script recursively.
AC_CONFIG_SUBDIRS([libXC S_GPU])


AC_OUTPUT


echo "
BigDFT configure summary:
========================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 compiler name:   $fc_type
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS

Build:
  Library ABINIT:            $ac_build_lib_ABINIT
  Library PSolver:           $ac_build_lib_PSolver
  Library BigDFT:            $ac_build_lib_BigDFT
  Main binaries (cluster...):$ac_build_bin_cluster
  Minima hopping binary:     $ac_build_bin_minima_hopping
  atom and pseudo binaries:  $ac_build_bin_pseudo

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
  With OpenCL support:       $ac_use_ocl
  With dgemmsy support:      $ac_use_dgemmsy
  With libXC:                $ac_use_libxc
   | internal built:         $ac_build_libxc
   | include dir.:           $LIBXC_INCLUDE
  With libABINIT:            $ac_use_libabinit
   | internal built:         $ac_build_libabinit
   | include dir.:           $LIBABINIT_INCLUDE
  With libS_GPU:             $ac_use_libsgpu
   | internal built:         $ac_build_libsgpu
   | include dir.:           $LIBSGPU_INCLUDE
  With ETSF_IO:              $ac_use_etsfio
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

if test x"$ac_use_mpi" = x"Warning" ; then
  echo "Warnings:
  MPI:                       disabled by tests
"
fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
sed -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
