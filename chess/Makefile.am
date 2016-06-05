EXTRA_DIST= m4

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = chess.pc

SUBDIRS = sparsematrix src doc

if HAVE_PYYAML
SUBDIRS += tests
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py

clean-local:
	rm -rf includes
