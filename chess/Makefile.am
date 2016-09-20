EXTRA_DIST= m4

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = chess.pc

SUBDIRS = sparsematrix src

if HAVE_PYYAML
SUBDIRS += tests
endif

SUBDIRS += doc

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py

clean-local:
	rm -rf includes
