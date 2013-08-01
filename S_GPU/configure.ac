#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([S_GPU],[1.0],[jerome.reybert@imag.fr])
AC_CONFIG_SRCDIR([src/sg_common_def.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL

# Dynamic library ==================================
# For libtoolize, check number of things for libtool
#LT_INIT
#LT_INIT
#`AC_PROG_RANLIB' is rendered obsolete by `LT_INIT'
#AC_PROG_RANLIB
AM_PROG_LIBTOOL
# Turn off shared libraries during beta-testing, since they
# make the build process take too long.

#DX_INIT_DOXYGEN([S_GPU], [doc/doxy.conf], [doc])
#DX_DOXYGEN_FEATURE(ON)
#DX_HTML_FEATURE(ON)
#DX_PDF_FEATURE(ON)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([netdb.h netinet/in.h stddef.h stdlib.h sys/socket.h unistd.h ])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRERROR_R
AC_CHECK_FUNCS([bzero gethostname socket])



# Check for Python
AX_PYTHON_DEFAULT
AX_PYTHON_WITH               # Optional, allow the user to specify python path
AX_PYTHON_PATH               # Check for python
AX_PYTHON_MODULE(os)
AX_PYTHON_MODULE(sys)

AM_INIT_AUTOMAKE

# -------------------#
# Setup CUDA         #
# -------------------#

##### Cuda path #####

dnl The user can define the CUDA root path, ie. the path including the bin, lib ... directories
AC_ARG_WITH(cuda-path, AS_HELP_STRING([--with-cuda-path], [define the root path to the NVidia Cuda tools (eg. /usr/local/cuda/).]),
[ac_with_cuda_path=$withval], [ac_with_cuda_path="/usr/local/cuda"])

if test "x$ac_with_cuda_path" = "xno" ; then
dnl If the CUDA root path has not been defined, check if nvcc program is in $PATH
  AC_CHECK_PROG([NVCC], [nvcc], [nvcc])
else
dnl Else, check if nvcc program is in CUDA_PATH/bin directory
  AC_PATH_PROG([NVCC], [nvcc], [], [${ac_with_cuda_path}/bin/])
  AC_LDFLAGS_CUDA="-L${ac_with_cuda_path}/lib"
fi

dnl User can define the binary directory path of CUDA, if this one is diferent than CUDA_PATH/bin
AC_ARG_WITH(cuda-path-bin, AS_HELP_STRING([--with-cuda-path-bin], [define the binary directory path to the NVidia Cuda tools (eg. /usr/local/cuda/bin).]),
[ac_with_cuda_path_bin=$withval], [ac_with_cuda_path_bin="no"])

dnl If CUDA_BIN_PATH is defined, check if nvcc command exists here
if test "x$ac_with_cuda_path_bin" != "xno" ; then
  AC_PATH_PROG([NVCC], [nvcc], [], [${ac_with_cuda_path_bin}])
fi

dnl After all, if nvcc command has not been found, configure stops with an error
if test -z "$NVCC"; then
  AC_MSG_WARN([CUDA compiler not found. Use --with-cuda-path or --with-cuda-path-bin or check your CUDA installation.])
fi

dnl Copy value into Makefiles: $NVCC = @NVCC@
AC_SUBST(NVCC)
AC_SUBST(AC_LDFLAGS_CUDA)

##### Cuda flags #####

dnl The default flags passed to nvcc compiler AND linker (maybe be cool to split them into 2 variables)
default_nvcc_flags="-g -O3 -DHAVE_CONFIG_H -I. -I.. -arch sm_13"

dnl User can add some flags to nvcc
AC_ARG_WITH(nvcc-flags, AS_HELP_STRING([--with-nvcc-flags], [add flags for the NVidia Cuda compiler and linker (default is "-g -O3 -DHAVE_CONFIG_H -I. -I.. -arch sm_13".)]),
[ac_with_nvcc_flags=$withval],
[ac_with_nvcc_flags="no"])

dnl If so, there are add to the default
if test "x$ac_with_nvcc_flags" != "xno" ; then
  NVCCFLAGS="$default_nvcc_flags $ac_with_nvcc_flags "
else
dnl If user defined himself NVCCFLAGS, default flags are not set
  if test -z "$NVCCFLAGS" ; then
    NVCCFLAGS="$default_nvcc_flags "
  fi
fi

dnl User can define cuda library directory . It can be useful in some cases, eg. if nvcc is not able to link against lib64 directory
AC_ARG_WITH(cuda-path-lib, AS_HELP_STRING([--with-cuda-path-lib], [define the library directory path to the NVidia Cuda tools (eg. /usr/local/cuda/lib).]),
[ac_with_cuda_path_lib=$withval], [ac_with_cuda_path_lib="no"])

if test "x$ac_with_cuda_path_lib" != "xno" ; then
  NVCCFLAGS+=" -L$ac_with_cuda_path_lib"
fi

dnl User can define the C compiler used with nvcc. Useful if gcc default compiler is gcc-4.4 (incompatible with nvcc 2.3)
AC_ARG_WITH(nvcc-c-compiler, AS_HELP_STRING([--with-nvcc-c-compiler], [specify the C compiler used with nvcc (default is gcc). You must specify the complete path]),
[ac_with_nvcc_c_compiler=$withval])

if test "$ac_with_nvcc_c_compiler" 
then
  NVCCFLAGS+=" --compiler-bindir $ac_with_nvcc_c_compiler "
fi

dnl User can turn on device emulation for the compilation
AC_ARG_ENABLE([emu], AS_HELP_STRING([--enable-emu], [Turn on device emulation for CUDA.]),
[ac_enable_emu=$enableval], [ac_enable_emu="no"])

if test "x$ac_enable_emu" = "xyes"
then
   NVCCFLAGS+=" -deviceemu"
fi

AC_SUBST(NVCCFLAGS)

##### Check for nvcc compilation (and execution?) #####

if test -n "$NVCC" ; then
  AC_MSG_CHECKING([for nvcc compilation])
ac_compile_nvcc=no
    cat>conftest.cu<<EOF
    __global__ static void test_cuda() {
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;
    __syncthreads();
    }
EOF
    if $NVCC -c $NVCCFLAGS conftest.cu &> /dev/null
    then
        ac_compile_nvcc=yes
    fi
#    rm -f conftest.c conftest.o conftest.what conftest.ever
    AC_MSG_RESULT([$ac_compile_nvcc])

  if test "x$ac_compile_nvcc" = "xno" 
  then
    AC_MSG_WARN([CUDA compiler seems to have problem. Maybe check the gcc used for nvcc.])
    NVCC=""
  fi
fi

dnl Add the variables at the configure of configure --help
AC_ARG_VAR( [NVCC], [Cuda compiler command] )
AC_ARG_VAR( [NVCCFLAGS], [Cuda compiler flags] )


#====================== END CUDA PART =========================================

#====================#
#   Doc building     #
#====================#

# Documentation building

AC_MSG_CHECKING([whether to build documentation for this library])
AC_ARG_ENABLE([build-doc], AC_HELP_STRING([--enable-build-doc],
	      [build documentation for this library (default=no)]),
	      [ac_enable_build_doc=$enableval], [ac_enable_build_doc="no"])
AC_MSG_RESULT([$ac_enable_build_doc])
AM_CONDITIONAL([BUILD_DOC], [test $ac_enable_build_doc = yes])

if test $ac_enable_build_doc = yes; then

  AC_ARG_VAR([DOXYGEN], [path to doxygen binary])
  AC_ARG_VAR([LATEX], [path to latex binary])
  AC_ARG_VAR([DVIPS], [path to dvips binary])

  # doxygen
  AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])
  if test -z "$DOXYGEN"; then
    AC_MSG_ERROR([--enable-build-doc specified but could not find doxygen binary])
  fi

  # LaTeX output
  AC_CHECK_PROG([LATEX], [latex], [latex])
  if test -z "$LATEX"; then
    AC_MSG_WARN([latex not found, not generating LaTeX output])
    GENERATE_LATEX=NO
  else
    AC_CHECK_PROG([MAKEINDEX], [makeindex], [makeindex])
    if test -z "$MAKEINDEX"; then
      AC_MSG_WARN([latex not found, not generating LaTeX output])
      GENERATE_LATEX=NO
    else
      AC_CHECK_PROG([DVIPS], [dvips], [dvips])
      if test -z "$DVIPS"; then
	AC_MSG_WARN([dvips not found, not generating .ps output])
      fi
      GENERATE_LATEX=YES
    fi
  fi

  # LaTeX/PDF output
  AC_CHECK_PROG([PS2PDF], [ps2pdf], [ps2pdf])
  if test -z "$PS2PDF"; then
    AC_MSG_WARN([ps2pdf not found, the output PDF file wont have hyperlinks])
    PDF_HYPERLINKS=NO
  else
    PDF_HYPERLINKS=YES
  fi

  AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex])
  if test -z "$PDFLATEX"; then
    AC_MSG_WARN([pdflatex not found, latex will be use to generate the output PDF file])
    USE_PDFLATEX=NO
  else
    USE_PDFLATEX=YES
  fi



fi

#================================ END DOC BUILDING =======================#

##=======================#
##   Examples building   #
##=======================#
#
ac_enable_build_examples="no"
AM_CONDITIONAL([BUILD_EXAMPLES], [test $ac_enable_build_examples = yes])
#AC_MSG_CHECKING([wether to build examples for this library])
#AC_ARG_ENABLE([build-examples], AS_HELP_STRING([--enable-build-examples], [Build the examples for this library (default no).]),
#[ac_enable_build_examples=$enableval], [ac_enable_build_examples="no"])
#AM_CONDITIONAL([BUILD_EXAMPLES], [test $ac_enable_build_examples = yes])
#
#if test $ac_enable_build_examples = yes; then
#  AC_LANG(C)
#  AX_MPI
#  AC_PROG_FC
#  AC_LANG(Fortran)
#  AX_MPI
#  AC_ARG_WITH(mpi-path, AS_HELP_STRING([--with-mpi-path], [define the MPI path.]),
#    [ac_with_mpi_path=$withval], [ac_with_mpi_path="no"])
#
#  AC_LANG(C)
#if test "x$ac_with_mpi_path" = "xno" ; then
#  AC_CHECK_HEADER("mpi.h")
#  AC_CHECK_HEADER([mpi.h],[],AC_MSG_ERROR([Cannot find mpi.h Please specify the prefix directory of MPI with with-mpi-path=/path/to/mpi]))
#  AC_CHECK_LIB([mpi],[MPI_Init],[], AC_MSG_ERROR([Cannot find libmpi.so Please specify the prefix directory of MPI with with-mpi-path=/path/to/mpi]))
#else
#  AC_CHECK_HEADER([$ac_with_mpi_path/include/mpi.h], [], AC_MSG_ERROR([Cannot find $ac_with_mpi_path/include/mpi.h Please specify the prefix directory of MPI with with-mpi-path=/path/to/mpi]))
#  TMP_LD="$LDFLAGS"
#  LDFLAGS="$LDFLAGS -L$ac_with_mpi_path/lib/"
#  AC_LANG([C])
#  AC_CHECK_LIB([mpi],[MPI_Init],[], AC_MSG_ERROR([Cannot find $ac_with_mpi_path/lib/libmpi.so Please specify the prefix directory of MPI with with-mpi-path=/path/to/mpi]))
#  LDFLAGS="$TMP_LD"
#  MPI_INCLUDE_PATH=$ac_with_mpi_path/include/
#  MPI_LIBRARY_PATH=$ac_with_mpi_path/lib/
#fi
#
#if test "x$ac_with_cuda_path_lib" != "xno" ; then
#  CUDA_LIB_PATH="$ac_with_cuda_path_lib"
#fi
#
#
#AC_SUBST(CUDA_LIB_PATH)
#AC_SUBST(MPI_INCLUDE_PATH)
#AC_SUBST(MPI_LIBRARY_PATH)
#
#fi

#================================ END EXAMPLES BUILDING =======================#


AC_MSG_CHECKING([if trace is enabled or not])
AC_ARG_ENABLE([trace], AS_HELP_STRING([--enable-trace], [Enable trace functions for S_GPU (default no).]),
[ac_enable_trace=$enableval], [ac_enable_trace="no"])

if test x$ac_enable_trace = xyes
then
AC_DEFINE([TRACE_S_GPU], [], [Enable trace for S_GPU])
fi
AC_MSG_RESULT([$ac_enable_trace])

AC_MSG_CHECKING([if verbose is enabled or not])
AC_ARG_ENABLE([verbose], AS_HELP_STRING([--enable-verbose], [Enable verbose mode (default no).]),
[ac_enable_verbose=$enableval], [ac_enable_verbose="no"])

if test x$ac_enable_verbose = xyes
then
AC_DEFINE([VERBOSE_S_GPU], [], [Enable verbose for S_GPU])
fi
AC_MSG_RESULT([$ac_enable_verbose])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_FILES([
Makefile
src/Makefile
src/CUDA/Makefile
examples/Makefile
examples/c/Makefile
examples/c/simple/Makefile
examples/c/add_array/Makefile
examples/fortran/Makefile
examples/fortran/add_array/Makefile
examples/c/simple/Makefile.standalone
examples/c/add_array/Makefile.standalone
examples/fortran/add_array/Makefile.standalone
])

AC_OUTPUT


echo "
Basics:
C   compiler:              $CC
C   flags:                 $CCFLAGS
CXX compiler:              $CXX
CXX flags:                 $CXXFLAGS
Linker flags:              $LDFLAGS
Linked libraries:          $LIBS

NVidia Cuda Compiler:      $NVCC
NVCC flags:                $NVCCFLAGS
NVCC C Compiler:           $ac_with_nvcc_c_compiler

Cuda emulation:            $ac_enable_emu
Trace functions:           $ac_enable_trace
Verbose:                   $ac_enable_verbose

Installation paths:
Source code location:      ${srcdir}
Prefix:                    $prefix
Exec prefix:               $exec_prefix
Binaries:                  $bindir
Static libraries:          $libdir
"


cat >> ./libtool <<DELIM

# ### BEGIN LIBTOOL TAG CONFIG: CUDA

# The linker used to build libraries.
LD="/usr/bin/ld"

# Commands used to build an old-style archive.
old_archive_cmds="\\\$AR \\\$AR_FLAGS \\\$oldlib\\\$oldobjs~\\\$RANLIB \\\$oldlib"

# A language specific compiler.
CC="$NVCC"

# Is the compiler the GNU compiler?
with_gcc=yes

# Compiler flag to turn off builtin functions.
no_builtin_flag=" -fno-builtin"

# How to pass a linker flag through the compiler.
wl="-Xlinker "

# Additional compiler flags for building library objects.
pic_flag=" -fPIC -DPIC"

# Compiler flag to prevent dynamic linking.
link_static_flag="-static"

# Does compiler simultaneously support -c and -o options?
compiler_c_o="yes"

# Whether or not to add -lc for building shared libraries.
build_libtool_need_lc=no

# Whether or not to disallow shared libs when runtime libs are static.
allow_libtool_libs_with_static_runtimes=no

# Compiler flag to allow reflexive dlopens.
export_dynamic_flag_spec="\\\${wl}--export-dynamic"

# Compiler flag to generate shared objects directly from archives.
whole_archive_flag_spec="\\\${wl}--whole-archive\\\$convenience \\\${wl}--no-whole-archive"

# Whether the compiler copes with passing no objects directly.
compiler_needs_object="no"

# Create an old-style archive from a shared archive.
old_archive_from_new_cmds=""

# Create a temporary old-style archive to link instead of a shared archive.
old_archive_from_expsyms_cmds=""

# Commands used to build a shared archive.
archive_cmds="\\\$CC -shared \\\$predep_objects \\\$libobjs \\\$deplibs \\\$postdep_objects \\\$compiler_flags \\\${wl}-soname \\\$wl\\\$soname -o \\\$lib"
archive_expsym_cmds="\\\$CC -shared -nostdlib \\\$predep_objects \\\$libobjs \\\$deplibs \\\$postdep_objects \\\$compiler_flags \\\${wl}-soname \\\$wl\\\$soname \\\${wl}-retain-symbols-file \\\$wl\\\$export_symbols -o \\\$lib"

# Commands used to build a loadable module if different from building
# a shared archive.
module_cmds=""
module_expsym_cmds=""

# Whether we are building with GNU ld or not.
with_gnu_ld="yes"

# Flag that allows shared libraries with undefined symbols to be built.
allow_undefined_flag=""

# Flag that enforces no undefined symbols.
no_undefined_flag=""

# Flag to hardcode $libdir into a binary during linking.
# This must work even if $libdir does not exist
hardcode_libdir_flag_spec="\\\${wl}-rpath \\\${wl}\\\$libdir"

# If ld is used when linking, flag to hardcode $libdir into a binary
# during linking.  This must work even if $libdir does not exist.
hardcode_libdir_flag_spec_ld=""

# Whether we need a single "-rpath" flag with a separated argument.
hardcode_libdir_separator=""

# Set to "yes" if using DIR/libNAME${shared_ext} during linking hardcodes
# DIR into the resulting binary.
hardcode_direct=no

# Set to "yes" if using DIR/libNAME${shared_ext} during linking hardcodes
# DIR into the resulting binary and the resulting library dependency is
# "absolute",i.e impossible to change by setting ${shlibpath_var} if the
# library is relocated.
hardcode_direct_absolute=no

# Set to "yes" if using the -LDIR flag during linking hardcodes DIR
# into the resulting binary.
hardcode_minus_L=no

# Set to "yes" if using SHLIBPATH_VAR=DIR during linking hardcodes DIR
# into the resulting binary.
hardcode_shlibpath_var=unsupported

# Set to "yes" if building a shared library automatically hardcodes DIR
# into the library and all subsequent libraries and executables linked
# against it.
hardcode_automatic=no

# Set to yes if linker adds runtime paths of dependent libraries
# to runtime path list.
inherit_rpath=no

# Whether libtool must link a program against all its dependency libraries.
link_all_deplibs=no

# Fix the shell variable $srcfile for the compiler.
fix_srcfile_path=""

# Set to "yes" if exported symbols are required.
always_export_symbols=no

# The commands to list exported symbols.
export_symbols_cmds="\\\$NM \\\$libobjs \\\$convenience | \\\$global_symbol_pipe | \\\$SED 's/.* //' | sort | uniq > \\\$export_symbols"

# Symbols that should not be listed in the preloaded symbols.
exclude_expsyms="_GLOBAL_OFFSET_TABLE_|_GLOBAL__F[ID]_.*"

# Symbols that must always be exported.
include_expsyms=""

# Commands necessary for linking programs (against libraries) with templates.
prelink_cmds=""

# Specify filename containing input files.
file_list_spec=""

# How to hardcode a shared library path into an executable.
hardcode_action=immediate

# The directories searched by this compiler when creating a shared library.
compiler_lib_search_dirs=""

# Dependencies to place before and after the objects being linked to
# create a shared library.
predep_objects=""
postdep_objects=""
predeps=""
postdeps=""

# The library search path used internally by the compiler when linking
# a shared library.
compiler_lib_search_path=""

# ### END LIBTOOL TAG CONFIG: CUDA

DELIM
