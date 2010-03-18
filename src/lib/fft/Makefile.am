## Process this file with automake to produce Makefile.in

OTHERS = fft_check.f90 fft_cache.f90

EXTRA_DIST = fft_check.out.ref fft_cache.gnuplot

check_PROGRAMS = fft_check fft_cache

fft_check_SOURCES = fft_check.f90 fft3d.f90 $(mpi_source)
fft_check_LDFLAGS =
fft_check_LDADD =

fft_cache_SOURCES = fft_cache.f90 fft3d.f90 $(mpi_source)
fft_cache_LDFLAGS =
fft_cache_LDADD =

CLEANFILES = fft_check.out fft_cache.out fft_cache.dat fft_columns.dat fft_perf.dat \
	fft_cache \
	module_fft_sg.@MODULE_EXT@ MODULE_FFT_SG.@MODULE_EXT@

fft_cache.o: fft3d.o
fft_check.o: fft3d.o

check:
	$(run_parallel) ./fft_check > fft_check.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.e-8 fft_check.out $(top_srcdir)/src/lib/fft/fft_check.out.ref
	@rm -f fft.cache 
	$(run_parallel) ./fft_cache 0
	$(run_parallel) ./fft_cache 6
	$(run_parallel) ./fft_cache 12
	$(run_parallel) ./fft_cache 24
	$(run_parallel) ./fft_cache 50
	$(run_parallel) ./fft_cache 75
	$(run_parallel) ./fft_cache 100
	python $(top_srcdir)/src/lib/fft/fft_cache.py
	
