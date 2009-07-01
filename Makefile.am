EXTRA_DIST = autogen.sh \
	configure

if BUILD_BINARIES
tests_dir = tests
else
tests_dir = 
endif

ACLOCAL_AMFLAGS = -I config/m4

SUBDIRS = src \
	utils \
	$(tests_dir) \
	doc

