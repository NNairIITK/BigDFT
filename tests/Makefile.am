## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  DFT \
	  overDFT \
	  tutorials \
	  libs

if USE_OCL
oclrun_message = oclrun
else
oclrun_message =
endif

if USE_MPI
mpirun_message = mpirun
else
mpirun_message =
endif

head_message: all-recursive
	@echo "======================================================="
	@echo "Type 'make check' in the directory tests      "

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "======================================================="

mpirun: head_message
	@echo "     Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo "     Use the environment variable run_ocl"
	@echo "     ex: export run_ocl='on' to use OpenCL acceleration"

EXTRA_DIST = \
	README \
	fldiff.py \
	report.py \
	styles images

CLEANFILES =

all: head_message foot_message

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_srcdir)/tests/report.py
