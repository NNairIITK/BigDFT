## Process this file with automake to produce Makefile.in

# Main targets.
# ------------
dynlibdir           = $(libdir)
bin_PROGRAMS        = splsad NEB global globaltool mhgps mhgpstool
EXTRA_PROGRAMS      = bart #NEB_images
bin_SCRIPTS         = NEB_driver.sh
lib_LIBRARIES       = libspred-1.a
dynlib_DATA         = $(spred_dynamic_library)
include_HEADERS     = bigdft_spred.h

# Sub-directories.
# ---------------
SUBDIRS = \
	spred_baselib \
	art \
	mhgps_files \
	global_files 

# put all subs
subs_static_library = \
	art/libart.a \
	mhgps_files/libmhgps.a \
	spred_baselib/libspredbasef.a \
	global_files/libglobalf.a 

if BUILD_DYNAMIC_LIBS
spred_dynamic_library = libspred-1.so.@SPRED_MINOR_VERSION@.0.@SPRED_MICRO_VERSION@
spred_soname_library = libspred-1.so.@SPRED_MINOR_VERSION@
spred_library = $(spred_soname_library)
AM_LDFLAGS = -Wl,-rpath=$(DESTDIR)$(dynlibdir)
else
spred_dynamic_library =
spred_soname_library =
spred_library = libspred-1.a @LIB_SPRED_DEPS@
endif
binaries_sources = splsad NEB global globaltool mhgps mhgpstool lowest_pathway


AM_FCFLAGS = -I. @MPI_INCLUDE@ @LIB_BIGDFT_CFLAGS@
AM_CFLAGS = @LIB_BIGDFT_CFLAGS@ @GLIB_CFLAGS@

libspred_1_a_SOURCES = images.f90 wrapper_neb.c

tmp-libspred/extract.stamp: $(subs_static_library)
	test -e "tmp-libspred" || $(INSTALL) -d -m 755 tmp-libspred
	cd tmp-libspred ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libspred-1.a: $(libspred_1_a_OBJECTS) tmp-libspred/extract.stamp
	$(AR) $(ARFLAGS) $@ $(libspred_1_a_OBJECTS) tmp-libspred/*.o
	$(RANLIB) $@

$(spred_dynamic_library): $(libspred_1_a_OBJECTS) tmp-libspred/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(spred_soname_library) -Wl,-rpath=$(DESTDIR)$(dynlibdir) -o $@ $(libspred_1_a_OBJECTS) tmp-libspred/*.o @LIB_SPRED_DEPS@
	@chmod a+x $@

$(spred_soname_library): $(spred_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libspred-1.so

install-data-hook:
	if test -n "$(spred_dynamic_library)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(spred_dynamic_library) ; \
	  ln -fs $(spred_dynamic_library) $(spred_soname_library) ; \
	  ln -fs $(spred_soname_library) libspred-1.so ; \
	fi

#NEB_images_SOURCES = NEB_images.f90
#NEB_images_LDADD = $(spred_library) @LIB_FUTILE_LIBS@

splsad_SOURCES = splinedsaddle.f90
splsad_LDADD = $(spred_library) @LIB_BIGDFT_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@

NEB_SOURCES = NEB.f90
NEB_LDADD = $(spred_library) @LIB_BIGDFT_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@

bart_SOURCES = art.f90
bart_FCFLAGS = -I. -Iart
bart_LDADD = $(spred_library) @LINALG_LIBS@

mhgps_SOURCES = mhgps.f90
mhgps_FCFLAGS = -I. -Ispred_baselib -Imhgps_files @LIB_BIGDFT_CFLAGS@
mhgps_LDADD = $(spred_library) @LIB_BIGDFT_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@
###mhgps_FCFLAGS = -I. -Ispred_baselib -Imhgps_files -I$(top_builddir)/includes @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@
###mhgps_LDADD = spred_baselib/libspredbasef.a mhgps_files/libmhgps.a $(spred_library) @LINALG_LIBS@
###mhgps_LDADD = $(spred_library) mhgps_files/libmhgps.a spred_baselib/libspredbasef.a @LINALG_LIBS@

mhgpstool_SOURCES = mhgpstool.f90
mhgpstool_FCFLAGS = -I. -Ispred_baselib -Imhgps_files @LIB_BIGDFT_CFLAGS@
##mhgpstool_LDADD = spred_baselib/libspredbasef.a $(spred_library) mhgps_files/libmhgps.a @LINALG_LIBS@
mhgpstool_LDADD = $(spred_library) @LIB_BIGDFT_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@

lowest_pathway_SOURCES = lowest_pathway.f90
lowest_pathway_FCFLAGS = 
lowest_pathway_LDADD =

global_SOURCES = global.f90
global_FCFLAGS = -I. -Ispred_baselib -Iglobal_files @LIB_BIGDFT_CFLAGS@
###global_LDADD = spred_baselib/libspredbasef.a $(spred_library) @LINALG_LIBS@
global_LDADD = $(spred_library) @LIB_BIGDFT_LIBS@ @LIB_FUTILE_LIBS@ @LINALG_LIBS@

globaltool_SOURCES = globaltool.f90
globaltool_FCFLAGS = -I. -Ispred_baselib -Iglobal_files @LIB_BIGDFT_CFLAGS@
##globaltool_LDADD = spred_baselib/libspredbasef.a  global_files/libglobalf.a $(spred_library)  @LINALG_LIBS@
globaltool_LDADD = $(spred_library) @LIB_FUTILE_LIBS@ @LINALG_LIBS@

CLEANFILES = *.@MODULE_EXT@ libspred-1.*
clean-local:
	rm -rf tmp-libspred

all: all-recursive

install: install-recursive

check:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

NEB.o: images.o

