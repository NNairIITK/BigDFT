## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  fft \
	  PSolver \
	  DFT \
	  H2-freq \
	  NEB \
	  XabsCheb


#This will run the automatic tests also for minima hopping
if USE_MINIMA_HOPPING
SUBDIRS += MINHOP
endif

if USE_CUDA_GPU
SUBDIRS += CUDA
endif

if USE_OCL
SUBDIRS += OpenCL
endif

if USE_CUDA_GPU
SUBDIRS += CUDA
endif

if USE_OCL
SUBDIRS += OpenCL
endif

if HAVE_ETSF_IO
SUBDIRS += ETSF_IO
endif


if USE_MPI
mpirun_message=mpirun
else
mpirun_message=
endif

head_message: all-recursive
	@echo "=============================================="
	@echo "Type 'make check' in the directory tests      "

foot_message: $(mpirun_message) head_message
	@echo "=============================================="

mpirun: head_message
	@echo "     Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='mpirun -np 2'  "


EXTRA_DIST = \
	README \
	fldiff.py \
	report.py

CLEANFILES =

all: head_message $(mpirun_message) foot_message

check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi
