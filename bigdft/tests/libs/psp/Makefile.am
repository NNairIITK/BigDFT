## Process this file with automake to produce Makefile.in

SUBDIRS =

EXTRA_DIST = psp.ref.yaml
SHORT_CHECK=
LONG_CHECK =psp.report.yaml

include $(srcdir)/../check.mk

check_PROGRAMS = check_psp

AM_FCFLAGS = @MPI_INCLUDE@ @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

check_psp_SOURCES = check_psp.f90
check_psp_LDADD = $(bigdft_library)

CLEANFILES = psp.out.yaml psp.report.yaml

check: psp.report.yaml
psp.report.yaml: psp.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $(abs_srcdir)/psp.ref.yaml -d $^ -o $@ -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=psp

psp.out.yaml: check_psp
	echo "Test LDA:" > $@
	for psp in $(top_srcdir)/utils/PSPfiles/Krach-LDA/*-q* ; do \
	  ./check_psp $$psp >> $@ ; \
	done
	echo "Test PBE:" >> $@
	for psp in $(top_srcdir)/utils/PSPfiles/Krach-PBE/*-q* ; do \
	  ./check_psp $$psp >> $@ ; \
	done
