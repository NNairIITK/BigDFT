## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = xc.out.ref

CLEANFILES = xc.out xc.tmp.* xc.report

if USE_DGEMMSY
libs_dgemmsy = $(top_builddir)/src/dgemmsy/libdgemmsy.a
else
libs_dgemmsy =
endif
if USE_CUDA_GPU
libs_cuda=$(top_builddir)/src/CUDA/libGPU.a
else
libs_cuda=
endif
if USE_OCL
libs_ocl=$(top_builddir)/src/OpenCL/libOCL.a -lstdc++
else
libs_ocl=
endif

check_PROGRAMS = test_xc

test_xc_SOURCES = test_xc.f90 
test_xc_LDADD = -lbigdft $(libs_dgemmsy) $(libs_cuda) $(libs_ocl) -labinit -lbigdft
test_xc_LDFLAGS = -L$(top_builddir)/src

AM_FCFLAGS = -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@ 

xc.out: test_xc
	./test_xc > xc.out

check: xc.out xc.out.ref
	@cut -b -70 xc.out > xc.tmp.out
	@cut -b -70 $(srcdir)/xc.out.ref > xc.tmp.ref
	python $(top_srcdir)/tests/fldiff.py --discrepancy=9e-6 xc.tmp.out xc.tmp.ref | tee xc.report
	@rm -f xc.tmp.*
