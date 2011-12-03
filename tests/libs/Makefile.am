SUBDIRS = fft \
	xc \
	PSolver \
	OP2P

if HAVE_ETSF_IO
SUBDIRS += ETSF_IO
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

check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi
