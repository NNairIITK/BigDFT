EXTRA_DIST = \
	$(convolution_includes)

convolution_includes = 	v.f90 \
	v_long.f90 \
	tree_long_shrink.f90 \
	long_shrink.f90 \
	tree_unrolled_grow.f90 \
	long_grow.f90

if USE_OPTI_CONVOLUT
convolution_source = convolut_simple.f90 \
	conv_per_optim.f90 \
	combined_grow.f90\
	combined_shrink.f90

else
convolution_source = convolut_simple.f90 \
	conv_per_optim.f90 \
	combined_grow.f90 \
	combined_shrink.f90

endif

AM_FCFLAGS = -I. -I$(srcdir)

noinst_LIBRARIES = libconvolutions.a

libconvolutions_a_SOURCES = $(convolution_source)

#dependencies

combined_grow.o: long_grow.f90 \
	tree_unrolled_grow.f90 \
	v.f90 \
	v_long.f90

combined_shrink.o: long_shrink.f90 \
	tree_long_shrink.f90 \
	v.f90 \
	v_long.f90

