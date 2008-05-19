EXTRA_DIST = \
	$(convolution_includes)

convolution_includes = 	v.inc \
	v_long.inc \
	v_17.inc 
#	tree_long_shrink.f90 \
#	long_shrink.f90 \
#	tree_unrolled_grow.f90 \
#	long_grow.f90

#GPU_INTERFACE if statement to be added
GPU_INTERFACE = interface_cuda_fake.f90
#	convolut_common_interface_cuda.f90 


if USE_OPTI_CONVOLUT
convolution_source = combined_shrink_optim.f90 \
	combined_grow_optim.f90 \
	combined_shrink_optim.f90 \
	convolut_ib_optim.f90 \
	conv_per_optim.f90 \
	convolut_common.f90 \
	$(GPU_INTERFACE)

else
convolution_source = combined_shrink_simple.f90 \
	combined_grow_simple.f90 \
	combined_shrink_simple.f90 \
	convolut_ib_simple.f90 \
	conv_per_simple.f90 \
	convolut_common.f90  \
	$(GPU_INTERFACE)
endif

AM_FCFLAGS = -I. -I$(srcdir) -I../modules

noinst_LIBRARIES = libconvolutions.a

libconvolutions_a_SOURCES = $(convolution_source)

#dependencies
combined_grow_simple.o : v_17.inc

combined_grow_optim.o : v_17.inc

combined_shrink_simple.o : v.inc

combined_shrink_optim.o : v.inc v_long.inc

convolution_prev.o : v.inc v_long.inc

