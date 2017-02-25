EXTRA_DIST= m4 chess_environment.sh.in chessvars.sh.in

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = chess.pc

SUBDIRS = sparsematrix chesspexsi src

dist_bin_SCRIPTS = chessvars.sh chess_environment.sh


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
