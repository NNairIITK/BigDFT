EXTRA_DIST = config/m4

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

#Poisson Solver and wrappers
if BUILD_LIB_PSOLVER
SUBDIRS += wrappers PSolver
endif

# BigDFT itself and structure prediction library at the end
if BUILD_LIB_BIGDFT
SUBDIRS += src spred
endif

# Tests subdir is included only if binaries are compiled.
if BUILD_BINARIES
if HAVE_PYYAML
SUBDIRS += tests
endif
endif

SUBDIRS += \
	doc \
	utils

check report:
	python $(top_builddir)/tests/report.py

clean-local:
	-rm -rf includes

#lcov-analysis:
#	cp -p $(srcdir)/flib/src/{*-inc.f90,*-inc-ptr.f90} flib/src; \
#	cp -p $(srcdir)/flib/src/error_handling.f90 flib/src/ ; \
#	cp -p $(srcdir)/flib/tests/yamlout.f90 flib/tests/ ; \
#	cp -p $(srcdir)/PSolver/src/Finite*.inc PSolver/tests/ ; \
#	cp -p $(srcdir)/PSolver/src/*.inc PSolver/src/ ; \
#	cp -p $(srcdir)/PSolver/src/{PSolver_Main,createKernel}.f90 PSolver/src/ ; \
#	cp -p $(srcdir)/wrappers/*-inc.f90 wrappers/ ; \
#	cp -p $(srcdir)/src/convolutions/{*-inc.f90,*.inc} src/convolutions/ ; \
#	cp -p $(srcdir)/src/modules/*-inc.f90 src/modules/ ; \
#	cp -p $(srcdir)/src/wavelib/*.inc src/ ; \
#	lcov --capture --directory . --output-file doc/lcov.info ; \
#	genhtml --legend doc/lcov.info --output-dir doc/lcov-html --ignore-errors source
#	@echo "Works with gcc and the options '-fprofile-arcs -ftest-coverage'"
#	@echo "Type <Browser> doc/lcov-html/index.html"
