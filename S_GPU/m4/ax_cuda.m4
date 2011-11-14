# AX_CUDA_RUN( CUDA_PROGRAM )
# -----------------
# Run a Python Test Program saving its output
# in ax_python_output and its condition code
# in ax_python_cc.

AC_DEFUN([AX_CUDA_RUN],
[
    AC_ARG_VAR( [NVCC], [Python Executable Path] )
    if test -z "$NVCC"
    then
        AC_MSG_ERROR([nvcc not found])
    else
        cat >cudatest.cu <<_ACEOF
$1
_ACEOF
        ax_cuda_output=`$NVCC cudatest.cu`
        ax_cuda_cc=$?
        rm cudatest.cu
        rm -f conftest.o conftest.lo
    fi
])


