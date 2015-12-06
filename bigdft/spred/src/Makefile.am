## Process this file with automake to produce Makefile.in

# Main targets.
# ------------
dynlibdir           = $(libdir)
bin_PROGRAMS        = $(binaries_sources)
noinst_PROGRAMS     = $(local_binaries)
bin_SCRIPTS         = $(scripts_sources)
lib_LIBRARIES       = $(lib_spred)
module_DATA         = $(mod_spred)
dynlib_DATA         = $(spred_dynamic_library)

# Sub-directories.
# ---------------
SUBDIRS = \
	art \
	mhgps_files \
	global_files 

EXTRA_DIST = \
	NEB_driver.sh

# Build targets, binaries and libraries.
# -------------
noinst_LIBRARIES = libmain.a
# put all subs
subs_static_library =
out_static_library =  \
	$(top_builddir)/PSolver/src/libPSolver-1.a \
	$(top_builddir)/wrappers/libwrappers.a

if BUILD_FLIB
out_static_library += $(top_builddir)/flib/src/libflib-1.a
endif

out_static_library +=  \
	$(top_builddir)/src/libbigdft-1.a


if BUILD_DYNAMIC_LIBS
spred_dynamic_library = libspred-1.so.@BIGDFT_MINOR_VERSION@.0.@BIGDFT_MICRO_VERSION@
spred_soname_library = libspred-1.so.@BIGDFT_MINOR_VERSION@
spred_library = $(spred_soname_library)
AM_LDFLAGS = -Wl,-rpath=$(DESTDIR)$(dynlibdir)
else
spred_dynamic_library =
spred_soname_library =
spred_library = libspred-1.a $(out_static_library) @LIBS_SHORT_DEPS@ $(amber_ld) $(cp2k_ld) $(pexsi_ld)
endif
if BUILD_LIB_BIGDFT
lib_spred = libspred-1.a
else
lib_spred =
endif
if BUILD_BINARIES
binaries_sources = splsad NEB global globaltool mhgps mhgpstool
scripts_sources = NEB_driver.sh
EXTRA_PROGRAMS = bart #NEB_images
else
binaries_sources =
local_binaries =
scripts_sources =
sub_dirs_bin =
check_PROGRAMS =
endif

if HAVE_AMBERTOOLS
amber_ld= @LIB_AMBERTOOLS@  -lsff -lnab -lpbsa -lcifparse -lrism -lfftw3 -larpack
else
amber_ld =
endif

if HAVE_CP2K
cp2k_ld= @CP2K_LINKLINE@ @LINALG_LIBS@
else
cp2k_ld =
endif

if HAVE_PEXSI
pexsi_ld= @PEXSI_LINKLINE@ @LINALG_LIBS@
else
pexsi_ld =
endif

AM_FCFLAGS = -I. -I$(top_builddir)/includes \
	@LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ \
	@LIBETSFIO_INCLUDE@ @LIBSGPU_INCLUDE@ \
	@MPI_INCLUDE@
AM_CPPFLAGS = @LIB_ARCHIVE_CFLAGS@ @LIBSGPU_INCLUDE@

AM_CFLAGS = -I$(top_srcdir)/src/bindings -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules @GLIB_CFLAGS@

#objects which depends from module_types
high_level = images.f90

libmain_a_SOURCES = $(high_level) wrapper_neb.c

libspred_1_a_SOURCES =

tmp-libspred/extract.stamp: $(subs_static_library) $(out_static_library)
	test -e "tmp-libspred" || $(INSTALL) -d -m 755 tmp-libspred
	cd tmp-libspred ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libspred-1.a: $(libmain_a_OBJECTS) tmp-libspred/extract.stamp
	$(AR) $(ARFLAGS) $@ $(libmain_a_OBJECTS) tmp-libspred/*.o
	$(RANLIB) $@

$(spred_dynamic_library): $(libmain_a_OBJECTS) tmp-libspred/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(spred_library) -o $@ $(libmain_a_OBJECTS) tmp-libspred/*.o @LIBS_SHORT_DEPS@ $(amber_ld) $(cp2k_ld) $(pexsi_ld)
	@chmod a+x $@

libspred-1.so.@BIGDFT_MINOR_VERSION@: $(spred_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libspred-1.so

install-data-hook:
	if test -n "$(spred_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(spred_dynamic_library) ; \
	  ln -fs $(spred_dynamic_library) $(spred_library) ; \
	  ln -fs $(spred_library) libspred-1.so ; \
	fi

#NEB_images_SOURCES = NEB_images.f90
#NEB_images_LDADD = $(spred_library) @LIB_FLIB_LIBS@

splsad_SOURCES = splinedsaddle.f90
splsad_LDADD = $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

NEB_SOURCES = NEB.f90
NEB_LDADD = $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

bart_SOURCES = art.f90
bart_FCFLAGS = -I. -Iart
bart_LDADD = art/libart.a $(spred_library) @LINALG_LIBS@

mhgps_SOURCES = mhgps.f90
mhgps_FCFLAGS = -I. -Imhgps_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@
mhgps_LDADD = mhgps_files/libmhgps.a $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

mhgpstool_SOURCES = mhgpstool.f90
mhgpstool_FCFLAGS = -I. -Imhgps_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@
mhgpstool_LDADD = mhgps_files/libmhgps.a $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

global_SOURCES = global.f90
global_LDADD = $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

globaltool_SOURCES = globaltool.f90
globaltool_FCFLAGS = -I. -Iglobal_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@
globaltool_LDADD = global_files/libglobalf.a $(spred_library) @LIB_FLIB_LIBS@ @LINALG_LIBS@

CLEANFILES = *.@MODULE_EXT@ bindings_dbus.* introspect libspred-1.* input_variables_definition-inc.h
clean-local:
	rm -rf tmp-libspred

all: all-recursive

install: install-recursive

check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

NEB.o: images.o

