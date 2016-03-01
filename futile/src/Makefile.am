# To be installed:
dynlibdir        = $(libdir)
noinst_LIBRARIES = $(futile_noinst)
lib_LIBRARIES    = $(futile_lib)
module_DATA      = $(modules)
dynlib_DATA      = $(futile_dynlib)
noinst_DATA      = $(futile_linklib)

# Static or dynamic library
futile_dynamic_library = libfutile-1.so.@FUTILE_MINOR_VERSION@.0.@FUTILE_MICRO_VERSION@
futile_library_soname  = libfutile-1.so.@FUTILE_MINOR_VERSION@

# We install futile library only if PSolver is deactivated.
futile_lib    = libfutile-1.a
if BUILD_DYNAMIC_LIBS
futile_dynlib = $(futile_dynamic_library)
futile_linklib = $(futile_library_soname)
else
futile_dynlib =
futile_linklib =
endif
futile_noinst =

libfutile_1_a_SOURCES = f_lib_highlev.f90

AM_FCFLAGS = -I$(top_builddir)/flib @MPI_INCLUDE@

tmp-libfutile/extract.stamp: $(top_builddir)/flib/libflib.a $(top_builddir)/wrappers/libwrappers.a
	test -e "tmp-libfutile" || $(INSTALL) -d -m 755 tmp-libfutile
	cd tmp-libfutile ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libfutile-1.a: $(libfutile_1_a_OBJECTS) tmp-libfutile/extract.stamp
	$(AR) $(ARFLAGS) $@ $(libfutile_1_a_OBJECTS) tmp-libfutile/*.o
	$(RANLIB) $@

# Dynamic library building
$(futile_dynamic_library): $(libfutile_1_a_OBJECTS) tmp-libfutile/extract.stamp
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(futile_library_soname) -o $@ $(libfutile_1_a_OBJECTS) tmp-libfutile/*.o @LIBS_EXT_DEPS@ $(LIBS)
	@chmod a+x $@

$(futile_library_soname): $(futile_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libfutile-1.so