noinst_LTLIBRARIES   = libcuda.la
libcuda_la_SOURCES   = cudafct.cu  \
	cuda_memory_management.cu \
	cudafct.h \
	cuda_memory_management.h

SUFFIXES = .cu

.cu.$(OBJEXT):
	$(NVCC) $(NVCCFLAGS) $(AM_CPPFLAGS) -c $< -o $@
.cu.lo:
	@PYTHON@ $(top_srcdir)/cudalt.py $@ $(NVCC) $(NVCCFLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS)\" -c $<