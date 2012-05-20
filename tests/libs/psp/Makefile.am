## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = psp.out.ref

check_PROGRAMS = check_psp

check_psp_SOURCES = check_psp.f90
check_psp_FCFLAGS = -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@
check_psp_LDADD = @LIBS_DEPENDENCIES@

CLEANFILES = psp.out psp.report

check: psp.report

psp.report: $(srcdir)/psp.out.ref psp.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

psp.out: check_psp
	echo "test LDA" > $@
	@rm -f lda.out
	for psp in $(top_srcdir)/utils/PSPfiles/Krach-LDA/*-q* ; do \
	  ./check_psp $$psp >> lda.out ; \
	done
	@lang=C sort -k5 -g lda.out >> $@
	@rm -f lda.out
	echo "test PBE" >> $@
	@rm -f pbe.out
	for psp in $(top_srcdir)/utils/PSPfiles/Krach-PBE/*-q* ; do \
	  ./check_psp $$psp >> pbe.out ; \
	done
	@lang=C sort -k5 -g pbe.out >> $@
	@rm -f pbe.out
