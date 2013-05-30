## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	README \
	fldiff.py \
	fldiff_yaml.py \
	report.py \
	styles images \
	yaml_hl.cfg \
	yaml_hl.py \
	parse_log.py \
	tols-BigDFT.yaml

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

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

head_message: all-recursive
	@echo "======================================================="
	@echo "Type 'make check' in the directory tests      "

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "======================================================="

mpirun: head_message
	@echo "     Use the environment variable run_parallel"
	@echo "         ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo "     Use the environment variable run_ocl"
	@echo "         ex: export run_ocl='on' to use OpenCL acceleration"
	@echo "     Use run_ocl='CPU' or 'ACC' to force use of hardware"
	@echo "         different than GPU"
	@echo "     and the environment variable ocl_platform"
	@echo "         ex: export ocl_platform='NVIDIA'"
	@echo "     and the environment variable ocl_devices"
	@echo "         ex: export ocl_devices='K20'"

CLEANFILES =

all: head_message foot_message

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
