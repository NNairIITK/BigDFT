EXTRA_DIST = autogen.sh config/m4 PyYAML-3.10.tar.gz

ACLOCAL_AMFLAGS = -I config/m4

DISTCLEANFILES = config.inc bigdft.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bigdft.pc

SUBDIRS =

# This is required for testing.
if BUILD_LIBYAML
pyyaml_build = pyyaml_build
pyyaml_clean = pyyaml_clean
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
SUBDIRS += yaml-0.1.4
else
pyyaml_build =
pyyaml_clean =
endif

# Low level libraries provided by BigDFT team.
SUBDIRS += \
	flib \
	wrappers \
	PSolver

if BUILD_LIB_BIGDFT
# All BigDFT dependencies.
if BUILD_LIBXC
SUBDIRS += libxc-2.0.x
endif
if BUILD_LIBABINIT
SUBDIRS += libABINIT/src
endif
if BUILD_LIBSGPU
SUBDIRS += S_GPU
endif

# BigDFT itself.
SUBDIRS += src

# Additional tools around BigDFT.
if BUILD_BINARIES_PSEUDO
SUBDIRS += pseudo
endif
endif

SUBDIRS += \
	doc \
	utils

if BUILD_LIB_BIGDFT
if BUILD_BINARIES
SUBDIRS += tests
endif
endif

pyyaml_build:
	cd $(top_builddir)/PyYAML-3.10 ; \
	rm -rf build ; \
	echo "[build_ext]" > setup.cfg ; \
	echo "define=YAML_DECLARE_STATIC" >> setup.cfg ; \
	echo "include_dirs=$(abs_top_srcdir)/yaml-0.1.4/include" >> setup.cfg ; \
	echo "library_dirs=$(abs_top_builddir)/yaml-0.1.4/src/.libs" >> setup.cfg ; \
	python setup.py build ; \
	cd -

pyyaml_clean:
	cd $(top_builddir)/PyYAML-3.10 ; \
	python setup.py clean ; \
	cd -

all: $(pyyaml_build)

clean: $(pyyaml_clean) clean-recursive

check report:
	export PYTHONPATH=${PYTHONPATH}; export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ;python $(top_srcdir)/tests/report.py

clean-local:
	-rm -rf includes
