


                    ##########################
                    # Makefile for pseudo2.5 #
                    ##########################

##############  CUDA  ###########################################################################

# example: PGI cudafortran with cublas and OpenMPI on maia
 
CC=gcc
cuFlags= -Mcuda 
#-fast -Mcuda 
cuLibs=- /opt/pgi/linux86-64/2010/cuda/3.1/lib64
cuLink=-l cublas
cublasInc= -I /opt/pgi/linux86-64/2010/cuda/3.1/include


################################################################################################

# specify the MPI fortran compiler to be used
# and optionally, distinguish f90-, f77-compiler and linker

FC		=mpif90 ##pgfortran
F90C 		= $(FC) 
Linker 		= $(FC) 

# fortran compiler and linker flags

# example: maia PGI openMPI
FFlags 		= -Bdynamic  -fastsse -O2 -tp nehalem-64 -R/opt/intel/mkl/10.2.3.029/lib/em64t
F90Flags 	= $(FFlags)
LFlags 		= -mp -lpthread -pgf90libs


#library flags for linking to Lapack Blas etc

#example: maia INTEL
#Libs = -L/opt/intel/mkl/10.2.3.029/lib/em64t  -lmkl_lapack95_lp64 -lmkl_intel_lp64  -lguide -lmkl_intel_thread -lmkl_core

#example: maia PGI OpenMPI
Libs 		= -llapack -lblas -lmkl_intel_lp64 -lmkl_pgi_thread -lmkl_core  -L/opt/intel/mkl/10.2.3.029/lib/em64t 



#specify the sources directory of your libXC install
#for example
libXC		=/home/kernph/wilale00/comp.libxc-1.0

libxcInc	= -I$(libXC)/include/
libxcLib	= -L$(libXC)/lib/
libxcLink	= -lm -lxc 

################################################################################################

#   List of sources (objects) grouped by dependencies


OBJf77	=crtvh.o \
	penalty.o \
	amoeba.o \
	gatom.o \
	wave.o \
	wave2.o \
	wave3.o \
	detnp.o \
	resid.o \
	etot.o \
	pj2test.o \
	xpown.o \
	gamma.o \
	ppack.o \
	radgrid.o \
	zero.o \
	zbrent.o

OBJmain	= pseudo.o

OBJf90	= errorhandler.o

OBJxc	= xcfunction.o \
	driveXC.o

OBJat	= atom.o \
	atom.splines.o 

##############  CUDA  ###########################################################################

OBJcuda =cublas.fortran.bindings.o\
	ekin_gauss_wvlt.cublas.o

################################################################################################

# targets

all:		pseudo atom


pseudo:		$(OBJf77) $(OBJf90) $(OBJxc) $(OBJmain) $(OBJcuda)
		@echo "Linking pseudo with cudafortran kernels ... "
		$(Linker) $(libxcInc) $(LFlags) $(cuFlags) $(cuLibs)  $(libxcLib)  -o pseudo\
		$(OBJf77) $(OBJf90) $(OBJxc) $(OBJmain) \
		$(OBJcuda)\
		$(libxcLink) $(cuLink) $(Libs)
		@echo "done"

atom:		$(OBJat) $(OBJxc)
		@echo "Linking atom... "
		$(Linker) $(libxcInc)  $(LFlags) $(libxcLib)  -o atom  \
		$(OBJat) $(OBJxc)   \
		$(libxcLink) $(Libs)  
		@echo "done"


################################################################################################

# dependencies


$(OBJf77): %.o: ../%.f
		@echo f77 object...
		$(FC)   $(FFlags)  -c $< -o $@

$(OBJf90): %.o: ../%.f90
		@echo f90 objects...
		$(F90C) $(F90Flags)  -c $< -o $@

$(OBJmain): %.o: ../%.f
		@echo pseudo.o, needs libxc ...
		$(FC)  $(libxcInc)  $(FFlags) $(libxcLib) -c $< -o $@


$(OBJxc): %.o: ../%.f90
		@echo objects that need libxc ...
		$(F90C)  $(libxcInc)  $(F90Flags) $(libxcLib) -c $< -o $@


$(OBJat): %.o: ../%.f
		@echo objects for atom ...
		$(FC)  $(libxcInc)  $(FFlags) $(libxcLib) -c $< -o $@


##############  CUDA  ###########################################################################

cublas.fortran.bindings.o:      cublas.fortran.bindings.c
				$(CC) -O2 $(cublasInc) -c cublas.fortran.bindings.c

ekin_gauss_wvlt.cublas.o:	ekin_gauss_wvlt.cublas.f90\
				convol_kernel.f90\
				crt_kernel.f90
				$(FC) $(cuFlags) ekin_gauss_wvlt.cublas.f90 -c 

################################################################################################

.PHONY: clean	bruteclean

clean:
		rm -f *.o

bruteclean:
		cp Makefile ../
		rm ./*
		mv ../Makefile ./

