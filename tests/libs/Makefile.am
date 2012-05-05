SUBDIRS = fft \
	xc \
	PSolver \
	OP2P \
	psp \
	IO-wf

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

check, report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi
