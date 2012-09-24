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

if BUILD_LIBYAML
pyyaml_build = pyyaml_build
pyyaml_clean = pyyaml_clean
else
pyyaml_build =
pyyaml_clean =
endif

pyyaml_build:
	cd $(top_builddir)/PyYAML-3.10 ; \
	rm -rf build ; \
	echo "[build_ext]" > setup.cfg ; \
	echo "define=YAML_DECLARE_STATIC" >> setup.cfg ; \
	echo "include_dirs=$(abs_top_srcdir)/yaml-0.1.4/include" >> setup.cfg ; \
	echo "library_dirs=$(top_builddir)/yaml-0.1.4/src/.libs" >> setup.cfg ; \
	python setup.py build ; \
	cd -

pyyaml_clean:
	cd $(top_builddir)/PyYAML-3.10 ; \
	python setup.py clean ; \
	cd -

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
	@echo "     use run_ocl='CPU' or 'ACC' to force use of hardware"
	@echo "     different than GPU"
	@echo "     and the environment variable ocl_platform"
	@echo "     ex: export ocl_platform='NVIDIA'"

EXTRA_DIST = \
	README \
	fldiff.py \
	fldiff_yaml.py \
	report.py \
	styles images \
	yaml_hl.cfg yaml_hl.py tols-BigDFT.yaml

CLEANFILES =

all: $(pyyaml_build) head_message foot_message

clean: $(pyyaml_clean) clean-recursive

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
