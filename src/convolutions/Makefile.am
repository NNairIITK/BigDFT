# To be installed:
noinst_LIBRARIES = libconvolutions.a

# Source code.
if USE_CUDA_GPU
CUDA_CONV = convolut_common_interface_cuda.f90
else
CUDA_CONV = interface_cuda_fake.f90
endif
if USE_OCL
OCL_CONV = convolut_common_interface_ocl.f90
else
OCL_CONV = interface_ocl_fake.f90
endif
if USE_OPTI_CONVOLUT
CONV_SOURCES = combined_shrink_optim.f90 \
	combined_grow_optim.f90 \
	growshrink_hyb_optim.f90 \
	convolut_ib_optim.f90 \
	conv_per_optim.f90 \
	convolut_optim_per.f90 \
	convolut_optim_slab.f90
else
CONV_SOURCES = combined_shrink_simple.f90 \
	combined_grow_simple.f90 \
	growshrink_hyb_simple.f90 \
	convolut_ib_simple.f90 \
	conv_per_simple.f90 \
	convolut_simple_per.f90 \
	convolut_simple_slab.f90
endif
libconvolutions_a_SOURCES = $(CUDA_CONV) $(OCL_CONV) $(CONV_SOURCES) \
	growshrink.f90 \
	growshrink_hyb_common.f90 \
	conv_per_common.f90 \
	convolut_common_slab.f90 \
	convolut_common_per.f90 \
	convolut_new_per.f90 \
	convolut_simple_per_k.f90 \
	convolut_simple_slab_k.f90
# convolut_ib_sse.c \ #to be verified, not working yet
EXTRA_DIST = v.inc v_long.inc v_17.inc

AM_FCFLAGS = -I$(srcdir) @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes
