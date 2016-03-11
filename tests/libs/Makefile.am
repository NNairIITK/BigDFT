SUBDIRS = fft \
	xcPSolver \
	OP2P \
	psp \
	IO-wf \
	wrappers \
	overlapgeneral \
	coordinatetransform \
	sparsematrix

if USE_BLACS
SUBDIRS += dsygv
endif 

if USE_CUDA_GPU
SUBDIRS += CUDA
endif

if USE_OCL
SUBDIRS += OpenCL
endif

if USE_LIBXC
SUBDIRS += xc
endif

if BUILD_BINDINGS
SUBDIRS += ftypes
endif

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
