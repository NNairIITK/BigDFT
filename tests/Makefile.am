## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  H \
	  H-spin \
	  Li+ \
	  SiH4 \
	  NaCl \
	  C \
	  Ca2 \
	  Mn2  \
	  H2O-CP2K \
	  SiBulk \
	  H2O-slab \
	  Davidson-SiH4 \
	  H2-freq \
	  NEB \
	  NH3 \
	  XabsCheb \
	  FeHyb


#This will run the automatic tests also for GPU
if USE_CUDA_GPU
SUBDIRS += GPU
endif

#This will run the automatic tests also for minima hopping
if USE_MINIMA_HOPPING
SUBDIRS += MINHOP
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
	fldiff.py \
	inguess.dat \
	report.py \
	Refs

CLEANFILES =

all: head_message $(mpirun_message) foot_message

check:
	@python $(top_srcdir)/tests/report.py
