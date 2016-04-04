# To be installed:
dynlibdir        = $(libdir)
lib_LIBRARIES    = libfutile-1.a
dynlib_DATA      = $(futile_dynlib)
noinst_DATA      = $(futile_linklib)

# Static or dynamic library
futile_dynamic_library = libfutile-1.so.@FUTILE_MINOR_VERSION@.0.@FUTILE_MICRO_VERSION@
futile_library_soname  = libfutile-1.so.@FUTILE_MINOR_VERSION@

# We install futile library only if PSolver is deactivated.
if BUILD_DYNAMIC_LIBS
futile_dynlib = $(futile_dynamic_library)
futile_linklib = $(futile_library_soname)
else
futile_dynlib =
futile_linklib =
endif

libfutile_1_a_SOURCES = f_lib_highlev.f90

AM_FCFLAGS = -I$(top_builddir)/flib @MPI_INCLUDE@

tmp-libfutile/extract.stamp: $(top_builddir)/flib/libflib.a $(top_builddir)/c-bindings/libbindings.a $(top_builddir)/wrappers/libwrappers.a
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

install-data-hook:
	if test -n "$(futile_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(futile_dynamic_library) ; \
	  ln -fs $(futile_dynamic_library) $(futile_library_soname) ; \
	  ln -fs $(futile_library_soname) libfutile-1.so ; \
	fi

CLEANFILES = libfutile-1.so*
clean-local:
	rm -rf tmp-libfutile
