EXTRA_DIST = autogen.sh \
	configure

if BUILD_BINARIES
tests_dir = tests
else
tests_dir = 
endif

ACLOCAL_AMFLAGS = -I config/m4

if BUILD_LIBXC
libXC_dir = libXC
else
libXC_dir = 
endif

SUBDIRS = $(libXC_dir) \
	src \
	utils \
	$(tests_dir) \
	doc

