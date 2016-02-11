## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	README \
	fldiff.py \
	styles \
	images \
	yaml_hl.cfg \
	yaml_hl.py \
	parse_log.py \
	tols-BigDFT.yaml \
	yaml_hl_remarks.cfg \
	fldiff_yaml.py.in \
	report.py.in

SUBDIRS = \
	  DFT \
	  overDFT \
	  tutorials \
	  libs \
	  benchs


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
	@echo "========================================================"
	@echo "To run the tests, type 'make check' (might take a while)"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "========================================================"

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

all: head_message foot_message report.py fldiff_yaml.py

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report.py: report.py.in
	@sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@

fldiff_yaml.py: fldiff_yaml.py.in
	@sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@

report:
	python $(top_builddir)/tests/report.py
