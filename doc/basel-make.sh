#!/bin/sh

# All in one for compilation in Basel.

echo "Creating the Makefile.in."
./autogen.sh

echo "Configuring the package for paralell version with ifort"
export FC=/opt/intel/fc/9.0/bin/ifort   # Choose your Fortran Compiler
export MPI_INCLUDE=/opt/scali/include   # Set the include directory for MPI
export MPI_LDFLAGS=/opt/scali/lib       # Set the directory where to find MPI library
export LDFLAGS="/opt/intel/mkl72/lib/32/libmkl_ia32.a -lguide -lpthread -D_REENTRANT"
./configure --enable-mpi --without-lapack --without-blas --with-ext-linalg-path=/opt/intel/mkl72/lib/32 --with-ext-linalg=mkl_lapack --prefix=$HOME/usr



echo "Compile the sources."
make
