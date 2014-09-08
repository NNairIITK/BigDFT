SUBDIRS = fft \
	xcPSolver \
	OP2P \
	psp \
	IO-wf \
	wrappers \
	overlapgeneral \
	coordinatetransform \
	xc

if USE_BLACS
SUBDIRS += dsygv
endif 

if USE_CUDA_GPU
SUBDIRS += CUDA
endif

if USE_OCL
SUBDIRS += OpenCL
endif

if BUILD_BINDINGS
SUBDIRS += ftypes
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report: 
	python $(top_builddir)/tests/report.py
