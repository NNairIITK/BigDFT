## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  H \
	  H-spin \
	  Li+ \
	  SiH4 \
	  NaCl \
	  C \
	  Ca2 \
	  H2O-CP2K

EXTRA_DIST = \
	fldiff.py \
	inguess.dat \
	report.py \
	Refs

CLEANFILES =  

all:
	@echo "=============================================="
	@echo "Type 'make check' in the directory tests      "
	@echo "     Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='ompirun -np 2'  "
	@echo "=============================================="

check:
	@python $(top_srcdir)/tests/report.py
