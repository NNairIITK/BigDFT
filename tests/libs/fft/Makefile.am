## Process this file with automake to produce Makefile.in

OTHERS = fft_check.f90 fft_cache.f90

EXTRA_DIST = fft_check.out.ref fft_cache.gnuplot fft_cache.py

check_PROGRAMS = fft_check fft_cache

fft_check_SOURCES = fft_check.f90 $(top_srcdir)/wrappers/fft/fft3d.f90 $(mpi_source)
fft_check_LDFLAGS =
fft_check_LDADD = @LINALG_LIBS@

fft_cache_SOURCES = fft_cache.f90 $(top_srcdir)/wrappers/fft/fft3d.f90 $(mpi_source)
fft_cache_LDFLAGS =
fft_cache_LDADD = @LINALG_LIBS@

CLEANFILES = fft_check.out fft_cache.out fft_cache.dat fft_columns.dat fft_perf.dat \
	fft_cache fldiff.report \
	module_fft_sg.@MODULE_EXT@ MODULE_FFT_SG.@MODULE_EXT@

fft_cache.o: fft3d.o
fft_check.o: fft3d.o

fft_check.out: fft_check
	./fft_check > fft_check.out

fft_cache.dat: fft_cache
	./fft_cache 0
	./fft_cache 6
	./fft_cache 12
	./fft_cache 24
	./fft_cache 50
	./fft_cache 75
	./fft_cache 100

check: fft_check.out fft_cache.dat
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.e-8 \
		fft_check.out $(srcdir)/fft_check.out.ref | tee fldiff.report
	python $(srcdir)/fft_cache.py

