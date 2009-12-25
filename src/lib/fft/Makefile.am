## Process this file with automake to produce Makefile.in

OTHERS = fft_check.f90 fft_cache.f90

EXTRA_DIST = fft_check.out.ref fft_cache.gnuplot

SUBDIRS = $(sub_dirs_xc) 

check_PROGRAMS = fft_check

AM_FCFLAGS = -I../../PSolver

fft_check_SOURCES = fft_check.f90 fft3d.f90 $(mpi_source)
fft_check_LDFLAGS =
fft_check_LDADD = 

CLEANFILES = fft_check.out fft_cache.out fft_cache.dat fft_columns.dat fft_perf.dat \
			 fft_cache-temp.f90 fft3d-temp.f90 fft_cache

check:
	$(run_parallel) ./fft_check > fft_check.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.e-8 fft_check.out $(top_srcdir)/src/lib/fft/fft_check.out.ref
	@cp $(top_srcdir)/src/lib/fft/fft3d.f90 fft3d-temp.f90
	@cp $(top_srcdir)/src/lib/fft/fft_cache.f90 fft_cache-temp.f90
	python $(top_srcdir)/src/lib/fft/fft_cache.py $(FC) $(FCFLAGS)
