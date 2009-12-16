## Process this file with automake to produce Makefile.in

OTHERS = fft_check.f90

EXTRA_DIST = fft_check.out.ref

SUBDIRS = $(sub_dirs_xc) 

check_PROGRAMS = fft_check

AM_FCFLAGS = -I../../PSolver

fft_check_SOURCES = fft_check.f90 fft3d.f90 $(mpi_source)
fft_check_LDFLAGS =
fft_check_LDADD = 

CLEANFILES = fft_check.out

check:
	$(run_parallel) ./fft_check > fft_check.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.e-8 fft_check.out $(top_srcdir)/src/lib/fft/fft_check.out.ref
