#!/bin/sh

# All in one for compilation in Basel.

echo "Creating the Makefile.in."
./autogen.sh

echo "Configuring the package for paralell version with ifort"
export FC=/opt/intel/fc/9.0/bin/ifort   # Choose your Fortran Compiler
export MPI_INCLUDE=/opt/scali/include   # Set the include directory for MPI
export MPI_LDFLAGS=/opt/scali/lib       # Set the directory where to find MPI library
export LDFLAGS="-L/opt/intel/mkl72/lib/32/ -lguide -lpthread -D_REENTRANT"
./configure --enable-mpi --with-lapack=/opt/intel/mkl72/lib/32/libmkl_lapack.a --with-blas=/opt/intel/mkl72/lib/32/libmkl_ia32.a --prefix=/home/caliste/usr

echo "Compile the sources."
make
