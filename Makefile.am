EXTRA_DIST = autogen.sh \
	configure \
	utils \
	test

if BUILD_BINARY
subdir_source = src
else
subdir_source = 
endif

if BUILD_LIBRARY
subdir_lib = lib
else
subdir_lib = 
endif



SUBDIRS = $(subdir_source) \
	$(subdir_lib)
