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
	  NEB


#This will run the automatic tests also for GPU
if USE_CUDA_GPU
SUBDIRS += GPU
endif

#This will run the automatic tests also for minima hopping
if USE_MINIMA_HOPPING
SUBDIRS += MINHOP
endif


EXTRA_DIST = \
	fldiff.py \
	inguess.dat \
	report.py \
	Refs

CLEANFILES =  

all:
	@echo "=============================================="
	@echo "Type 'make check' in the directory tests      "
	@echo "     Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='ompirun -np 2'  "
	@echo "=============================================="

check:
	@python $(top_srcdir)/tests/report.py
