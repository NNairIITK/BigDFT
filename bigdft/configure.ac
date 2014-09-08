

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - DFT over wavelets],[1.7.5],[Damien.Caliste@cea.fr],[bigdft])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE(-Wno-portability)

dnl Versioning
BIGDFT_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
BIGDFT_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
BIGDFT_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
BIGDFT_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(BIGDFT_MAJOR_VERSION)
AC_SUBST(BIGDFT_MINOR_VERSION)
AC_SUBST(BIGDFT_MICRO_VERSION)
AC_SUBST(BIGDFT_NANO_VERSION)

dnl This is a hack for usage inside ABINIT.
dnl if test -n "$FCFLAGS" ; then
dnl   FCFLAGS=`echo $FCFLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi
dnl if test -n "$F90FLAGS" ; then
dnl   F90FLAGS=`echo $F90FLAGS | sed 's/-qsuffix=[[a-zA-Z:0-9=]]*//'`
dnl fi

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the ABINIT part.
AC_PROG_CC()
AC_PROG_CXX()
AM_PROG_CC_C_O()

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
AX_FC_FLUSH()
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

dnl Search for pkg-config
PKG_PROG_PKG_CONFIG()

dnl Produce dynamic excecutables.
AC_ARG_ENABLE(dynamic-libraries, AS_HELP_STRING([--enable-dynamic-libraries],
                                               [Build dynamical libraries (disabled by default).]),
              ac_build_dynamic=$enableval, ac_build_dynamic="no")

dnl Basics for minimal bindings.
AX_FC_POINTER()

dnl Look for GLib and friends for bindings.
AC_ARG_ENABLE(bindings, AS_HELP_STRING([--enable-bindings], [Build C bindings (disabled by default).]),
              ac_build_bindings=$enableval, ac_build_bindings="no")
if test x"$ac_build_bindings" = x"yes" -a $F90_1D_POINTER_SIZE = 0 ; then
  ac_build_bindings="warn"
fi
AC_ARG_WITH(gobject, AS_HELP_STRING([--with-gobject],
                     [Build GObject bindings (default is yes).]),
            ac_use_glib=$withval, ac_use_glib=$ac_build_bindings)
if test x"$ac_use_glib" = x"yes"; then
  PKG_CHECK_MODULES(GLIB, glib-2.0 gobject-2.0 gthread-2.0 gio-2.0 >= 2.22, [have_glib=yes], [have_glib=no])
  PKG_CHECK_MODULES(GDBUS, glib-2.0 gobject-2.0 gthread-2.0 gio-2.0 >= 2.22 gio-unix-2.0 >= 2.30, [have_gdbus=yes], [have_gdbus=no])
  dnl Merge GLIB_LIBS and GDBUS_LIBS if GDBUS is available.
  if test "$have_gdbus" = "yes" ; then
    GLIB_LIBS=$GDBUS_LIBS
  fi
  PKG_CHECK_MODULES([GOBJECT_INTROSPECTION],
                    [gobject-introspection-1.0 >= 0.9],
                    [enable_introspection=yes],
                    [enable_introspection=no])
  if test "x$enable_introspection" = "xyes" ; then
    dnl girdir=$($PKG_CONFIG --variable=girdir gobject-introspection-1.0)
    dnl typelibsdir=$($PKG_CONFIG --variable=typelibdir gobject-introspection-1.0)
    girdir="$datadir/gir-1.0"
    typelibsdir="$libdir/girepository-1.0"
    AC_SUBST(girdir)
    AC_SUBST(typelibsdir)
    AC_SUBST([G_IR_SCANNER], [$($PKG_CONFIG --variable=g_ir_scanner gobject-introspection-1.0)])
    AC_SUBST([G_IR_COMPILER], [$($PKG_CONFIG --variable=g_ir_compiler gobject-introspection-1.0)])
    dnl Force dynamic library building.
    ac_build_dynamic="yes"
  fi
else
  have_glib=no
  have_gdbus=no
  enable_introspection=no
fi
dnl Test for library building tools.
if test x"$ac_build_dynamic" = x"yes" ; then
  AX_FLAG_PIC()
  if test -z "$ax_flag_pic" ; then
    AC_MSG_WARN(["No available position-independent code flag, bindings disabled."])
    ac_build_dynamic=no
    ac_build_bindings=warn
    have_glib=no
    have_gdbus=no
    enable_introspection=no
  else
    case $CFLAGS   in *"$ax_flag_pic"*) ;; *) CFLAGS="$CFLAGS $ax_flag_pic";; esac
    case $CXXFLAGS in *"$ax_flag_pic"*) ;; *) CXXFLAGS="$CXXFLAGS $ax_flag_pic";; esac
    case $FCFLAGS  in *"$ax_flag_pic"*) ;; *) FCFLAGS="$FCFLAGS $ax_flag_pic";; esac

    eval "set x $ac_configure_args"
    shift
    ac_configure_args=
    for ac_arg ; do
      case $ac_arg in
        CFLAGS=* | CXXFLAGS=* | FCFLAGS=*) ;;
        *) ac_configure_args="$ac_configure_args '$ac_arg'"
      esac
    done
    ac_configure_args="$ac_configure_args 'CFLAGS=$CFLAGS' 'CXXFLAGS=$CXXFLAGS' 'FCFLAGS=$FCFLAGS'"
  fi
else
  ax_flag_pic=""
fi
AM_CONDITIONAL([BUILD_DYNAMIC_LIBS], [test "x$ac_build_dynamic" = "xyes"])
AM_CONDITIONAL([WITH_GOBJECT_INTROSPECTION], [test "x$enable_introspection" = "xyes"])
if test x"$have_glib" = x"yes" ; then
  AC_DEFINE([HAVE_GLIB], [], [If set, we can call glib.h])
  AC_SUBST([GLIB_TRUE], [""])
  AC_SUBST([GLIB_END_TRUE], [""])
  AC_SUBST([GLIB_FALSE], ["/*"])
  AC_SUBST([GLIB_END_FALSE], ["*/"])
else
  AC_SUBST([GLIB_TRUE], ["/*"])
  AC_SUBST([GLIB_END_TRUE], ["*/"])
  AC_SUBST([GLIB_FALSE], [""])
  AC_SUBST([GLIB_END_FALSE], [""])
fi
AM_CONDITIONAL(HAVE_GLIB, test "$have_glib" = "yes")
if test x"$have_gdbus" = x"yes" ; then
  AC_DEFINE([HAVE_GDBUS], [], [If set, we can use gdbus])
fi
AM_CONDITIONAL(HAVE_GDBUS, test "$have_gdbus" = "yes")
AM_CONDITIONAL(BUILD_BINDINGS, test "$ac_build_bindings" = "yes")
if test x"$ac_build_bindings" = x"yes" ; then
  AC_SUBST([BUILD_BINDINGS_START], [""])
  AC_SUBST([BUILD_BINDINGS_END], [""])
else
  AC_SUBST([BUILD_BINDINGS_START], ["/*"])
  AC_SUBST([BUILD_BINDINGS_END], ["*/"])
fi

dnl Enable or not MPI.
AX_MPI()

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
if test `echo "$ac_with_memory_limit" | $SED "s/[.].*//"` -ge "0" 2> /dev/null ; then
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


dnl Use the dgemmsy
ac_use_dgemmsy="no"
DGEMMSY_CPPFLAGS=
AC_ARG_ENABLE(dgemmsy, AS_HELP_STRING([--enable-dgemmsy], [Use dgemmsy (disabled by default).]),
                         ac_use_dgemmsy=$enableval, ac_use_dgemmsy="no")
if test x"$ac_use_dgemmsy" = "xyes" ; then
  AC_LANG_PUSH(C)

  DGEMMSY_CPPFLAGS="-msse3"
  CPPFLAGS_SVG=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $DGEMMSY_CPPFLAGS"
  AC_CHECK_HEADER([pmmintrin.h],
                  [ac_use_dgemmsy="yes"],
                  [ac_use_dgemmsy="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  AC_LANG_POP(C)
  if test "$ac_use_dgemmsy" = "no"; then
    DGEMMSY_CPPFLAGS=
    AC_MSG_WARN(["No SSE header found, dgemmsy disabled."])
  fi
fi
AC_SUBST(DGEMMSY_CPPFLAGS)
AM_CONDITIONAL(USE_DGEMMSY, test "$ac_use_dgemmsy" = "yes")


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
  if test x"$ac_with_cuda_path" = x"no" ; then
    ac_with_cuda_path="/usr"
  fi
  CUDA_PATH=$ac_with_cuda_path
  AC_SUBST(CUDA_PATH)
  dnl Test the existence of things.
  AC_PATH_PROG([NVCC], [nvcc], [$ac_with_nvcc], [$PATH:${ac_with_cuda_path}/bin])
  AC_SUBST(NVCC)

  dnl look for the cublas header
  AC_LANG_PUSH(C)
  CPPFLAGS_SVG=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS -I$ac_with_cuda_path/include"
  AC_CHECK_HEADER([cublas.h],
                  [ac_cublas_header="yes"],
                  [ac_cublas_header="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  if test "$ac_cublas_header" = "no"; then
    AC_MSG_WARN(["No 'CUDA' header found, CUDA parts disabled."])
    ac_use_cuda_gpu="no"
  else
    CUDA_INCLUDE_PATH=$ac_with_cuda_path"/include"
    AC_SUBST(CUDA_INCLUDE_PATH)
    CUDA_INCLUDE="-I"$ac_with_cuda_path"/include"
    AC_SUBST(CUDA_INCLUDE)
  dnl Test the link on libcudart and libcublas.
  LIBS_OLD=$LIBS
  fi
  dnl look for the cublas library
  ac_cuda_lib_path="$ac_with_cuda_path/lib"
  LIBS="-L$ac_cuda_lib_path $LIBS_OLD"
  AC_CHECK_LIB(cublas, cublasSsymv, withlibcublas=yes, withlibcublas=no, [-lcudart])
  if test "$withlibcublas" = "no"; then
    AC_CHECK_LIB(cublas, cublasSsymv_v2, withlibcublas=yes, withlibcublas=no, [-lcudart])
  fi
  if test "$withlibcublas" = "no"; then
    ac_cuda_lib_path="$ac_with_cuda_path/lib64"
    LIBS="-L$ac_cuda_lib_path $LIBS_OLD"
    AC_CHECK_LIB(cublas, cublasScopy, withlibcublas=yes, withlibcublas=no, [-lcudart])
    if test "$withlibcublas" = "no"; then
      AC_CHECK_LIB(cublas, cublasScopy_v2, withlibcublas=yes, withlibcublas=no, [-lcudart])
    fi
  fi
  if test "$withlibcublas" = "no"; then
    AC_MSG_ERROR(["No 'cublas' library found, link will fail."])
  fi
  LIBS=$LIBS_OLD
  AC_LANG_POP(C)

  dnl Add the flags.
  AC_MSG_CHECKING([for NVCC flags])
  if test -z "$NVCC_FLAGS" ; then
    NVCC_FLAGS="-O3 --compiler-options '-fno-strict-aliasing $ax_flag_pic'"
  fi
  AC_SUBST(NVCC_FLAGS)
  AC_MSG_RESULT([$NVCC_FLAGS])
  LIBCUDA_LIBS="-lcudart -lcublas -lcufft"
  LDFLAGS="$LDFLAGS -L$ac_cuda_lib_path"
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
                         [ac_ocl_path=$withval], [ac_ocl_path=""])
  dnl Test the link on libOpenCL.
  dnl look for the header.
  CPPFLAGS_SVG=$CPPFLAGS
  if test -n "$ac_ocl_path" ; then
     CPPFLAGS="$CPPFLAGS -I$ac_ocl_path/include"
  fi
  AC_CHECK_HEADER([CL/cl.h],
                  [ac_ocl_header="yes"],
                  [ac_ocl_header="no"])
  CPPFLAGS=$CPPFLAGS_SVG
  if test "$ac_ocl_header" = "no"; then
    AC_MSG_WARN(["No 'OpenCL' header found, OpenCL parts disabled."])
    ac_use_ocl="no"
  else
    if test -n "$ac_ocl_path" ; then
       OCL_INCLUDE_PATH=$ac_ocl_path"/include"
       OCL_INCLUDE="-I"$ac_ocl_path"/include"
    fi
    AC_SUBST(OCL_INCLUDE_PATH)
    AC_SUBST(OCL_INCLUDE)
    dnl Now test for the library.
    LIBS_OLD=$LIBS
    if test -n "$ac_ocl_path" ; then
       LIBS="-L$ac_ocl_path/lib $LIBS"
    fi
    AC_CHECK_LIB(OpenCL, clCreateContext, withlibocl=yes, withlibocl=no)
    if test "$withlibocl" = "no"; then
      if test -n "$ac_ocl_path" ; then
        LIBS="-L$ac_ocl_path/lib64 $LIBS"
      fi
      AC_CHECK_LIB(OpenCL, clCreateContext, withlibocl=yes, withlibocl=no)
    fi
    if test "$withlibocl" = "no"; then
      AC_MSG_WARN(["No 'OpenCL' library found, OpenCL parts disabled."])
      ac_use_ocl="no"
    else
      if test -n "$ac_ocl_path" ; then
         LDFLAGS="$LDFLAGS -L$ac_ocl_path/lib"
      fi
      LIBOCL_LIBS="-lOpenCL -lm -lstdc++"
    fi
    LIBS=$LIBS_OLD
  fi
  AC_LANG_POP(C)
fi
AM_CONDITIONAL(USE_OCL, test "$ac_use_ocl" = "yes")

dnl Use the OpenCL support
ac_use_mic="no"
AC_ARG_ENABLE(intel-mic, AS_HELP_STRING([--enable-intel-mic], [Use MIC implementation for GPU convolutions (disabled by default).]),
                         ac_use_mic=$enableval, ac_use_mic="no")
AC_MSG_CHECKING([for Intel MIC support])
AC_MSG_RESULT([$ac_use_mic])
AC_ARG_WITH(intel-mic-libs, AS_HELP_STRING([--with-intel-mic-libs], [Link MIC executable with the following addional libraries (default = empty).]), ac_mic_libs=$withval, ac_mic_libs=no)
if test x"$ac_use_mic" = "xyes" ; then
  if test x"$ac_mic_libs" != "xno" ; then
    LDFLAGS_MIC="-offload-option,mic,ld,\"$ac_mic_libs\""
  fi
fi
AM_CONDITIONAL(USE_MIC, test "$ac_use_mic" = "yes")

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
      ac_etsfio_incs=""
      if test -n "$ac_etsfio_dir" ; then
         ac_etsfio_incs=$ac_etsfio_incs" -I$ac_etsfio_dir/include"
      fi
      if test -n "$ac_netcdf_dir" -a x"$ac_netcdf_dir" != x"$ac_etsfio_dir" ; then
         ac_etsfio_incs=$ac_etsfio_incs" -I$ac_netcdf_dir/include"
      fi
      AC_SUBST(LIBETSFIO_INCLUDE, $ac_etsfio_incs)
      LIBETSFIO_LIBS="-letsf_io_utils -letsf_io $ac_netcdf_libs"
   else
      ac_use_etsfio="warn"
   fi
fi
AM_CONDITIONAL(HAVE_ETSF_IO, test "$ac_use_etsfio" = "yes")

dnl Test libraries required for benching.
AC_LANG_PUSH(C++)
AC_CHECK_LIB(config++, config_init, withconfig=yes, withconfig=no)
AM_CONDITIONAL(HAVE_LIBCONFIG, test "$withconfig" = "yes")
AC_LANG_POP(C++)

dnl AM_CONDITIONAL for the different libraries, set default values.
dnl the BigDFT library.
ac_build_lib_BigDFT="yes"
dnl the standalone executable (required linking tests)
ac_build_bin_cluster="yes"


dnl disable different build
dnl The Poisson solver library
AC_ARG_ENABLE(libpsolver, AS_HELP_STRING([--disable-libpsolver], [build the PSolver library (default = yes).]), ac_build_lib_PSolver=$enableval, ac_build_lib_PSolver="yes")
dnl The BigDFT library
AC_ARG_ENABLE(libbigdft, AS_HELP_STRING([--disable-libbigdft], [build the BigDFT library (default = yes).]), ac_build_lib_BigDFT=$enableval, ac_build_lib_BigDFT="yes")
if test x"$ac_build_lib_PSolver" != x"yes" -a x"$ac_build_lib_BigDFT" == x"yes" ; then
  ac_build_lib_BigDFT="disabled"
fi
dnl The cluster executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_cluster=$enableval, ac_build_bin_cluster="yes")
if test x"$ac_build_lib_BigDFT" != x"yes" -a x"$ac_build_bin_cluster" == x"yes" ; then
  ac_build_bin_cluster="disabled"
fi

AM_CONDITIONAL(BUILD_LIB_PSOLVER, test "$ac_build_lib_PSolver" = "yes")
AC_MSG_CHECKING([for PSolver library build])
AC_MSG_RESULT([$ac_build_lib_PSolver])
AM_CONDITIONAL(BUILD_LIB_BIGDFT, test "$ac_build_lib_BigDFT" = "yes")
AC_MSG_CHECKING([for BigDFT library build])
AC_MSG_RESULT([$ac_build_lib_BigDFT])
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_cluster" = "yes")
AC_MSG_CHECKING([for BigDFT binaries build])
AC_MSG_RESULT([$ac_build_bin_cluster])

dnl ---------------------------------------
AX_LINALG()

dnl Error on the linking part only if required.
if test "$ac_build_bin_cluster" = "yes" ; then
  if [ test $withlinalg = "no" && test $withblas = "no" ] ; then
    AC_MSG_ERROR(["No blas library. Use --with-blas-path or check your Blas installation."])
  fi
  if [ test $withlinalg = "no" && test $withlapack = "no" ] ; then
    AC_MSG_ERROR(["No lapack library. Use --with-lapack-path or check your Lapack installation. Linalg: $LINALG_LIBS"])
  fi
fi # End of the if library compilation

dnl Archive input file support
ac_use_archive="yes"
AC_ARG_WITH([archives],
            AS_HELP_STRING([--with-archives], [add support of archives for output files.]),
            [ac_use_archive=$withval], [ac_use_archive=$ac_build_lib_BigDFT])
AC_ARG_WITH([archives-path],
            AS_HELP_STRING([--with-archives-path], [give a path to find libarchive.]),
            [ac_path_archive=$withval])
if test x"$ac_use_archive" != x"no" ; then
   if test x"$ac_path_archive" == x"" ; then
      PKG_CHECK_MODULES([LIB_ARCHIVE],
                        [libarchive >= 2.4],
                        [ac_use_archive=yes],
                        [ac_use_archive=warn])
   else
      LDFLAGS_SVG="$LDFLAGS"
      AC_LANG_PUSH(C)
      LDFLAGS="-L$ac_path_archive/lib"
      AC_CHECK_LIB([archive], [archive_read_data_block],
                   [ac_use_archive=yes], [ac_use_archive=warn])
      if test x"$ac_use_archive" = x"yes" ; then
         LIB_ARCHIVE_CFLAGS="-I$ac_path_archive/include"
         LIB_ARCHIVE_LIBS="-L$ac_path_archive/lib -larchive"
      fi
      AC_LANG_POP(C)
      LDFLAGS="$LDFLAGS_SVG"
   fi
   if test x"$ac_use_archive" = x"yes" ; then
      AC_DEFINE([HAVE_LIB_ARCHIVE], [1], [libarchive is linkable.])
   else
      AC_MSG_WARN([libarchive is not available])
   fi
fi
AC_LANG_PUSH(C)
AC_CHECK_FUNCS([strndup])
AC_LANG_POP(C)
AM_CONDITIONAL(HAVE_LIB_ARCHIVE, test x"$ac_use_archive" = x"yes")
AC_SUBST(LIB_ARCHIVE_CFLAGS)
AC_SUBST(LIB_ARCHIVE_LIBS)

dnl C-level YAML support.
AC_ARG_WITH([yaml-path],
            AS_HELP_STRING([--with-yaml-path], [give a path to find libyaml.]),
            [ac_path_yaml=$withval])

if test x"$ac_path_yaml" == x"" ; then
   ac_path_yaml="/usr"
fi

LDFLAGS_SVG="$LDFLAGS"
AC_LANG_PUSH(C)
LDFLAGS="-L$ac_path_yaml/lib"
AC_CHECK_LIB([yaml], [yaml_parser_parse],
             [ac_use_libyaml=yes], [ac_use_libyaml=no])
if test x"$ac_use_libyaml" = x"yes"; then
   if test x"$ac_path_yaml" != x"/usr" ; then
      LIB_YAML_CFLAGS="-I$ac_path_yaml/include"
      LIB_YAML_LIBS="-L$ac_path_yaml/lib "
   fi
   LIB_YAML_LIBS=$LIB_YAML_LIBS"-lyaml"
else
   AC_MSG_ERROR([libyaml is not available, install YAML and provide path --with-yaml-path.])
fi
AC_LANG_POP(C)
LDFLAGS="$LDFLAGS_SVG"

AC_SUBST(LIB_YAML_CFLAGS)
AC_SUBST(LIB_YAML_LIBS)

dnl Python YAML support.
AX_PYYAML()

dnl Test the libXC library
AX_LIBXC()

dnl Test the ABINIT library
AC_ARG_WITH(libabinit-path, AS_HELP_STRING([--with-libabinit-path],
                                           [Give the path for an external libABINIT library, lib and module files (default = /etsf).]),
            ac_libabinit_dir=$withval, ac_libabinit_dir="/etsf")
ac_use_libabinit="no"
dnl Test the given implementation of libabinit.
LDFLAGS="$LDFLAGS -L$ac_libabinit_dir/lib"
AC_CHECK_LIB(abinit, symbrav, withlibabinit=yes, withlibabinit=no)
AC_CHECK_FILE($ac_libabinit_dir/include/ab6_moldyn.$ax_fc_mod_ext, moldyn="yes", moldyn="no")
AC_CHECK_FILE($ac_libabinit_dir/include/m_ab6_symmetry.$ax_fc_mod_ext, sym="yes", sym="no")
dnl AC_CHECK_FILE($ac_libabinit_dir/include/libxc_functionals.$ax_fc_mod_ext, libxc="yes", libxc="no")
if test "$withlibabinit" = "yes" -a "$moldyn" = "yes" -a "$sym" = "yes"; then
  ac_use_libabinit="yes"
  ac_libabinit_incs="-I$ac_libabinit_dir/include"
else
  AC_MSG_ERROR([Missing libABINIT headers, consider using --with-liabinit-path.])
fi

dnl Test the libPAW library.
ac_internal_libpaw="yes"
AC_ARG_WITH(libpaw-path, AS_HELP_STRING([--with-libpaw-path],
                                        [Give the path for an external libPAW library, lib and module files (default = inside libABINIT).]),
            ac_internal_libpaw=$withval, ac_internal_libpaw="yes")
if test x"$ac_internal_libpaw" == x"yes" -o x"$ac_internal_libpaw" == x"internal" ; then
   ac_internal_libpaw="yes"
else
   dnl Given path should be tested here.
   ac_libpaw_incs="-I$ac_internal_libpaw/include"
   LDFLAGS="$LDFLAGS -L$ac_internal_libpaw/lib"
   LIBABINIT_LIBS="-lpaw $LIBABINIT_LIBS"
fi

AC_DEFINE([HAVE_LIBABINIT], [1], [use libABINIT in BigDFT.])
AC_SUBST(LIBABINIT_INCLUDE, "$ac_libpaw_incs $ac_libabinit_incs")
if test x"$ac_internal_libpaw" == x"yes" ; then
  LIBABINIT_LIBS="-labinit"
else
  LIBABINIT_LIBS="-lpaw -labinit"
fi

dnl Test S_GPU library.
ac_build_libsgpu="no"
ac_use_libsgpu="no"
AC_ARG_ENABLE(internal-libsgpu, AS_HELP_STRING([--disable-internal-libsgpu], [Do not build and link with internal S_GPU library (default = yes).]), ac_internal_sgpu=$enableval, ac_internal_sgpu="no")
ac_build_libsgpu=$ac_internal_sgpu
if test x"$ac_use_cuda_gpu" = "xyes" ; then
  dnl Test S_GPU library.
  ac_internal_sgpu="yes"
  AC_ARG_WITH(libsgpu-path, AS_HELP_STRING([--with-libsgpu-path], [Give the path for an external S_GPU library, lib and module files (default = None).]),
            ac_libsgpu_dir=$withval, ac_libsgpu_dir=)
  if test "$ac_internal_sgpu" = "yes" ; then
    ac_build_libsgpu="yes"
    ac_use_libsgpu="yes"
    AC_SUBST(LIBSGPU_INCLUDE, "-DCUBLAS_USE_THUNKING -I\$(top_srcdir)"/S_GPU/src)
  else
    if test -n "$ac_libsgpu_dir" ; then
      LDFLAGS="$LDFLAGS -L$ac_libsgpu_dir/lib"
    
      dnl Test the given implementation of sgpu.
      AC_CHECK_LIB(s_gpu, sg_create_stream, withsgpu=yes, withsgpu=no)
      AC_CHECK_FILE($ac_libsgpu_dir/include/s_gpu.h, header="yes", header="no")
      if test "$withsgpu" = "yes" -a "$header" = "yes" ; then
        ac_use_libsgpu="yes"
        AC_SUBST(LIBSGPU_INCLUDE, "-DCUBLAS_USE_THUNKING -I"$ac_libsgpu_dir/include)
      fi
    fi
  fi
  if test x"$ac_use_libsgpu" = x"yes" ; then
    AC_DEFINE([HAVE_LIBSGPU], [1], [use S_GPU in BigDFT.])
    if test "$ac_internal_sgpu" = "yes" ; then
      LIBSGPU_LIBS="\$(top_builddir)/S_GPU/src/.libs/libs_gpu.a -lstdc++"
    else
      LIBSGPU_LIBS="-ls_gpu -lstdc++"
    fi
  fi
fi
AM_CONDITIONAL(BUILD_LIBSGPU, test x"$ac_build_libsgpu" = x"yes")

dnl Test the availability of clock_gettime() and rt library.
AC_LANG_PUSH(C)
AC_CHECK_HEADERS([time.h])
AC_CHECK_LIB([rt], [clock_gettime])
AC_CHECK_FUNCS([clock_gettime])
AC_LANG_POP(C)

dnl Agregate all flags.
LDFLAGS="$LDFLAGS $LDFLAGS_MIC"

dnl Agregate all the optional libs.
LIBS_EXT_DEPS="$LIBCUDA_LIBS $LIBOCL_LIBS $LIBETSFIO_LIBS $LINALG_LIBS $LIB_ARCHIVE_LIBS $LIB_YAML_LIBS $GLIB_LIBS $MPI_LDFLAGS $MPI_LIBS"
AC_SUBST(LIBS_EXT_DEPS, [$LIBS_EXT_DEPS])
LIBS_SHORT_DEPS="$LIBABINIT_LIBS $LIB_XC_LIBS $LIBSGPU_LIBS $LIBS_EXT_DEPS"
AC_SUBST(LIBS_SHORT_DEPS, [$LIBS_SHORT_DEPS])

dnl Test for the documentation.
AC_CHECK_PROG(ac_have_pdflatex, [pdflatex], [yes], [no])
AM_CONDITIONAL(BUILD_PDFDOC, test x"$ac_have_pdflatex" = x"yes")
dnl Test for the developers' documentation.
AC_CHECK_PROG(ac_have_doxygen, [doxygen], [yes], [no])
AC_ARG_ENABLE(devel-doc,
              AS_HELP_STRING([--enable-devel-doc],
                             [compile documentation for developers using Doxygen.]),
              ac_devel_doc=$enableval, ac_devel_doc="no")
AM_CONDITIONAL(BUILD_DEVDOC, test x"$ac_have_doxygen" = x"yes" -a x"$ac_devel_doc" = x"yes")
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

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
bigdft.pc
flib/Makefile
flib/src/Makefile
flib/tests/Makefile
wrappers/Makefile
PSolver/Makefile
PSolver/src/configure.inc
PSolver/src/Makefile
PSolver/tests/Makefile
config/Doxyfile
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
src/bigdft-tool
src/Makefile
src/modules/Makefile
src/modules/bigdft_cst.h
src/modules/configure.inc
src/modules/configure.c
src/convolutions/Makefile
src/xcPSolver/Makefile
src/CUDA/Makefile
src/OpenCL/Makefile
src/convolutions-c/Makefile
src/bindings/Makefile
src/art/Makefile
src/tools/Makefile
src/tools/bader/Makefile
tests/Makefile
tests/fldiff_yaml.py
tests/report.py
tests/libs/Makefile
tests/libs/fft/Makefile
tests/libs/OP2P/Makefile
tests/libs/xc/Makefile
tests/libs/xcPSolver/Makefile
tests/libs/OpenCL/Makefile
tests/libs/CUDA/Makefile
tests/libs/IO-wf/Makefile
tests/libs/dsygv/Makefile
tests/libs/ftypes/Makefile
tests/libs/wrappers/Makefile
tests/libs/psp/Makefile
tests/libs/overlapgeneral/Makefile
tests/libs/coordinatetransform/Makefile
tests/DFT/Makefile
tests/DFT/cubic/Makefile
tests/DFT/linear/Makefile
tests/DFT/postSCF/Makefile
tests/overDFT/Makefile
tests/tutorials/Makefile
])

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
  Linked libraries:          $LIBS_DEPENDENCIES $LIBS
  Linked dependencies:       $LIBS_SHORT_DEPS

Build:
  Library flib:              yes
  Library PSolver:           $ac_build_lib_PSolver
  Library BigDFT:            $ac_build_lib_BigDFT
  Main binaries (bigdft,...):$ac_build_bin_cluster
  Bindings:                  $ac_build_bindings
  User documentation:        $ac_have_pdflatex
  Devel. documentation:      $ac_devel_doc (Doxygen: $ac_have_doxygen)

Options:
  Debug version:             $ac_use_debug
  Max. Memory per task (GB): $MEMORY_LIMIT
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With optimised conv.:      $ac_use_opti_convolut
  With Cuda GPU conv.:       $ac_use_cuda_gpu
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
  With OpenCL support:       $ac_use_ocl
  With Intel MIC support:    $ac_use_mic
  With dgemmsy support:      $ac_use_dgemmsy
  With libXC:                $ac_use_libxc
   | internal built:         $ac_build_libxc
   | include dir.:           $LIBXC_INCLUDE
  With libS_GPU:             $ac_use_libsgpu
   | internal built:         $ac_build_libsgpu
   | include dir.:           $LIBSGPU_INCLUDE
  With ETSF_IO:              $ac_use_etsfio
   | include dir.:           $LIBETSFIO_INCLUDE
  With libarchive:           $ac_use_archive
   | C flags:                $LIB_ARCHIVE_CFLAGS
   | linked libraries:       $LIB_ARCHIVE_LIBS
  With libyaml:              $ac_use_libyaml
   | internal built:         $ac_build_libyaml
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

if test x"$ac_use_mpi" = x"warn" -o x"$ac_use_archive" = x"warn" -o x"$ac_use_etsfio" = x"warn" -o x"$ac_build_bindings" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$ac_use_archive" = x"warn" ; then
    echo "  libarchive:                disabled by tests"
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
