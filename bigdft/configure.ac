AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([BigDFT - DFT over wavelets],[1.8],[Damien.Caliste@cea.fr],[bigdft])

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
BIGDFT_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
BIGDFT_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
BIGDFT_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$BIGDFT_MICRO_VERSION" ; then
  BIGDFT_MICRO_VERSION=0
fi
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

AC_LANG_PUSH(C)
AC_CHECK_FUNCS([strndup])
AC_LANG_POP(C)

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

dnl Libraries
AX_MPI()
AX_LINALG()
AX_FLIB()
AX_CHESS()
AX_PSOLVER()
if test "$ax_have_PSOLVER" != "yes" ; then
  AC_MSG_ERROR(["Cannot find PSolver library."])
fi


dnl Look for GLib and friends for bindings.
AX_FC_POINTER()
AX_INTROSPECTION()
AX_DYNAMIC_LIBRARIES($ax_have_introspection, [f_lib_initialize], [$LIB_FUTILE_LIBS])

AC_ARG_ENABLE(bindings, AS_HELP_STRING([--enable-bindings], [Build C bindings (disabled by default).]),
              ac_build_bindings=$enableval, ac_build_bindings="no")
if test x"$ac_build_bindings" = x"yes" -a $F90_1D_POINTER_SIZE = 0 ; then
  ac_build_bindings="warn"
fi
AM_CONDITIONAL(BUILD_BINDINGS, test "$ac_build_bindings" = "yes")
if test x"$ac_build_bindings" = x"yes" ; then
  AC_SUBST([BUILD_BINDINGS_START], [""])
  AC_SUBST([BUILD_BINDINGS_END], [""])
else
  AC_SUBST([BUILD_BINDINGS_START], ["/*"])
  AC_SUBST([BUILD_BINDINGS_END], ["*/"])
fi

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

dnl Accelerator support.
AX_ACC_CUDA()
AX_ACC_OCL()
AX_ACC_MIC()

dnl Test libraries required for benching / testing.
AX_PYYAML()
AC_LANG_PUSH(C++)
AC_CHECK_LIB(config++, config_init, withconfig=yes, withconfig=no)
AM_CONDITIONAL(HAVE_LIBCONFIG, test "$withconfig" = "yes")
AC_LANG_POP(C++)

dnl The cluster executable
AC_ARG_ENABLE(binaries, AS_HELP_STRING([--disable-binaries], [build the binaries (default = yes).]), ac_build_bin_cluster=$enableval, ac_build_bin_cluster="yes")
AM_CONDITIONAL(BUILD_BINARIES, test "$ac_build_bin_cluster" = "yes")
AC_MSG_CHECKING([for BigDFT binaries build])
AC_MSG_RESULT([$ac_build_bin_cluster])

dnl Error on the linking part only if required.
if test "$ac_build_bin_cluster" = "yes" -a $ax_have_linalg != "yes" ; then
  AC_MSG_ERROR(["No blas/lapack implementation."])
fi # End of the if library compilation

dnl Test libraries
AX_ETSF_IO()
AX_LIBXC()
AX_LIBGAIN()
AX_ARCHIVE()
AX_AMBERTOOLS()
dnl PEXSI
AC_ARG_WITH(pexsi, AS_HELP_STRING([--with-pexsi],
            [Give the link-line for PEXSI.]),
            ac_pexsi_linkline=$withval, ac_pexsi_linkline=)
ax_have_pexsi="no"
if test -n "$ac_pexsi_linkline" ; then
  PEXSI_LINKLINE="$ac_pexsi_linkline"
  ax_have_pexsi="yes"
fi
AM_CONDITIONAL(HAVE_PEXSI, test $ax_have_pexsi = "yes")
dnl cp2k
AC_ARG_WITH(cp2k, AS_HELP_STRING([--with-cp2k],
            [Give the link-line for cp2k.]),
            ac_cp2k_linkline=$withval, ac_cp2k_linkline=)
ax_have_cp2k="no"
if test -n "$ac_cp2k_linkline" ; then
  CP2K_LINKLINE="$ac_cp2k_linkline"
  ax_have_cp2k="yes"
fi
AM_CONDITIONAL(HAVE_CP2K, test $ax_have_cp2k = "yes")

dnl Test the ABINIT library
AC_ARG_WITH(libabinit-incs, AS_HELP_STRING([--with-libabinit-incs],
                                           [Give the compiler include flags for an external libABINIT library (default = None).]),
            ac_libabinit_incs=$withval, ac_libabinit_incs="")
AC_ARG_WITH(libabinit-incs, AS_HELP_STRING([--with-libabinit-incs],
                                           [Give the linker flags for an external libABINIT modules (default = None).]),
            ac_libabinit_libs=$withval, ac_libabinit_libs="")

dnl Test the given implementation of libabinit.
ac_use_libabinit="no"
AC_MSG_CHECKING([for libABINIT modules])
FCFLAGS_SVG=$FCFLAGS
if test -n "$ac_libabinit_incs" ; then
  FCFLAGS="$FCFLAGS $ac_libabinit_incs"
elif test -n "$C_INCLUDE_PATH" ; then
  for path in ${C_INCLUDE_PATH//:/ }; do
    ac_libabinit_incs="$ac_libabinit_incs -I$path"
  done
  FCFLAGS="$FCFLAGS $ac_libabinit_incs"
fi
AC_COMPILE_IFELSE([[program main
  use m_ab6_moldyn
  use m_ab6_symmetry
  use m_ab6_kpoints

  write(*,*) AB6_MAX_SYMMETRIES
end program]], withabinitmod=yes, withabinitmod=no)
FCFLAGS=$FCFLAGS_SVG
AC_MSG_RESULT($withabinitmod)

LIBS_SVG=$LIBS
if test -n "$ac_libabinit_libs" ; then
  LIBS="$LIBS $ac_libabinit_libs"
else
  ac_libabinit_libs="-labinit"
fi
AC_CHECK_LIB(abinit, abi_symbrav, withlibabinit=yes, withlibabinit=no)
LIBS=$LIBS_SVG

if test "$withlibabinit" = "yes" -a "$withabinitmod" = "yes"; then
  ac_use_libabinit="yes"
else
  AC_MSG_ERROR([Missing libABINIT, consider using --with-libabinit-[incs,libs].])
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
   LIB_ABINIT_LIBS="-lpaw $LIB_ABINIT_LIBS"
fi

AC_SUBST(LIBABINIT_INCLUDE, "$ac_libpaw_incs $ac_libabinit_incs")
if test x"$ac_internal_libpaw" == x"yes" ; then
  LIB_ABINIT_LIBS=$ac_libabinit_libs
else
  LIB_ABINIT_LIBS="-lpaw $ac_libabinit_libs"
fi

dnl Agregate all the optional libs.
LIB_BIGDFT_DEPS="$LIB_ABINIT_LIBS $LIB_XC_LIBS $LIB_GAIN_LIBS $CP2K_LINKLINE $PEXSI_LINKLINE $LIB_AMBERTOOLS_LIBS $LIBOCL_LIBS $LIBETSFIO_LIBS $LIB_ARCHIVE_LIBS $LIB_CHESS_DYNAMIC_LIBS $LIB_PSOLVER_DYNAMIC_LIBS $LIB_FUTILE_LIBS"
if test $ax_FUTILE_static != "yes" ; then
   LIB_BIGDFT_DEPS="$LIB_BIGDFT_DEPS $LIBCUDA_LIBS $LINALG_LIBS $GLIB_LIBS $MPI_LDFLAGS $MPI_LIBS $PYTHON_LIBS"
fi
dnl LIB_BIGDFT_DEPS="$LIB_ABINIT_LIBS $LIB_XC_LIBS $LIB_GAIN_LIBS $CP2K_LINKLINE $PEXSI_LINKLINE $LIB_AMBERTOOLS_LIBS $LIBCUDA_LIBS $LIBOCL_LIBS $LIBETSFIO_LIBS $LINALG_LIBS $LIB_ARCHIVE_LIBS $GLIB_LIBS $MPI_LDFLAGS $MPI_LIBS $PYTHON_LIBS"
LIB_BIGDFT_PKG=""
if test $ax_build_dynamic != "yes" ; then
  LIB_BIGDFT_PKG="$LIB_BIGDFT_DEPS"
fi
AC_SUBST(LIB_BIGDFT_PKG)
dnl LIB_BIGDFT_DEPS="$LIB_BIGDFT_DEPS $LIB_PSOLVER_LIBS $LIB_FUTILE_LIBS"
AC_SUBST(LIB_BIGDFT_DEPS)

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
docdir="${docdir}"
AC_SUBST(docdir)

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
bigdft.pc
config/Doxyfile
utils/Makefile
utils/PSPfiles/Makefile
doc/Makefile
src/bigdft-tool
src/Makefile
src/modules/Makefile
src/modules/configure.c
src/orbitals/Makefile
src/orbitals/configure.inc
src/xcPSolver/Makefile
src/CUDA/Makefile
src/OpenCL/Makefile
src/convolutions-c/Makefile
src/bindings/Makefile
src/bindings/bigdft_cst.h
src/geopt/Makefile
src/tools/Makefile
src/tools/bader/Makefile
tests/Makefile
tests/libs/Makefile
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
tests/libs/sparsematrix/Makefile
tests/benchs/Makefile
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
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIB_BIGDFT_DEPS $LIBS

Build:
  Main binaries (bigdft,...):$ac_build_bin_cluster
  Bindings:                  $ac_build_bindings
  User documentation:        $ac_have_pdflatex
  Devel. documentation:      $ac_devel_doc (Doxygen: $ac_have_doxygen)

Options:
  Debug version:             $ac_use_debug
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With optimised conv.:      $ac_use_opti_convolut
  With Cuda GPU conv.:       $ax_have_acc_cuda
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
  With OpenCL support:       $ax_have_acc_ocl
  With Intel MIC support:    $ax_have_acc_mic
  With AMBERTOOLS compiler:  $ax_have_amber
   | nab Compiler:           $NAB
  With CP2K:                 $ax_have_cp2k
  With PEXSI:                $ax_have_pexsi
  With libXC:                $ac_use_libxc
   | internal built:         $ac_build_libxc
   | include dir.:           $LIBXC_INCLUDE
  With ETSF_IO:              $ax_have_etsf_io
   | include dir.:           $LIBETSFIO_INCLUDE
  With libarchive:           $ax_have_archive
   | C flags:                $LIB_ARCHIVE_CFLAGS
   | linked libraries:       $LIB_ARCHIVE_LIBS
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

if test x"$ac_use_mpi" = x"warn" -o x"$ax_have_archive" = x"warn" -o x"$ax_have_etsf_io" = x"warn" -o x"$ac_build_bindings" = x"warn" -o x"$ax_have_scalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$ax_have_archive" = x"warn" ; then
    echo "  libarchive:                disabled by tests"
  fi
  if test x"$ax_have_etsf_io" = x"warn" ; then
    echo "  ETSF_IO:                   disabled by tests"
  fi
  if test x"$ac_build_bindings" = x"warn" ; then
    echo "  C bindings:                disabled by tests"
  fi
  if test x"$ax_have_scalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
