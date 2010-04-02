EXTRA_DIST = autogen.sh \
	configure

ACLOCAL_AMFLAGS = -I config/m4

if BUILD_LIBABINIT
libABINIT_dir = libABINIT/src
else
libABINIT_dir = 
endif

if BUILD_LIBXC
libXC_dir = libXC
else
libXC_dir = 
endif

if BUILD_LIBSGPU
libSGPU_dir = S_GPU
else
libSGPU_dir = 
endif

if BUILD_BINARIES
tests_dir = tests
else
tests_dir = 
endif

SUBDIRS = $(libXC_dir) \
	$(libSGPU_dir) \
	$(libABINIT_dir) \
	src \
	doc \
	utils \
	$(tests_dir)

check:
	@python $(top_srcdir)/tests/report.py
