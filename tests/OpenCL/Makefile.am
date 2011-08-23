EXTRA_DIST = \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

if USE_OCL
check_PROGRAMS = conv_check
my_tests = OpenCL.report
else
check_PROGRAMS =
my_tests =
endif

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


CLEANFILES = fort.1 conv_check.out OpenCL.report

AM_FCFLAGS = -I../../src/modules @MPI_INCLUDE@ -I. -I$(srcdir) @LIBABINIT_INCLUDE@ 
AM_CPPFLAGS = -I@OCL_INCLUDE_PATH@
AM_CFLAGS = -Wall -Werror

conv_check_SOURCES = $(top_srcdir)/src/OpenCL/conv_check.f90 
conv_check_LDADD =  $(top_builddir)/src/libbigdft.a \
	$(top_builddir)/src/PSolver/libpoissonsolver.a \
	$(top_builddir)/src/OpenCL/libbasec.a \
	$(top_builddir)/src/OpenCL/libbasef.a $(libs_dgemmsy) $(libs_cuda)
conv_check_LDFLAGS = -lrt -lstdc++


conv_check.out: conv_check
	echo 1 128 128 16384 16384 2 > fort.1
	./conv_check > conv_check.out

OpenCL.report: conv_check.out
	awk -f $(srcdir)/conv_check.awk $< | tee OpenCL.report

check: $(my_tests)
