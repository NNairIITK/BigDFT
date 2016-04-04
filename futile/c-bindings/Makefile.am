if HAVE_PYTHON
pylib = libpyimp.a
pymodule = futile.so
endif

noinst_LIBRARIES = libbindings.a $(pylib)
include_HEADERS  = dict.h futile.h futile_cst.h
pyexec_DATA      = $(pymodule)

if HAVE_PYTHON
pysource = python.c
endif

libbindings_a_SOURCES = misc.c futile.h dict.c dict.h dictf.f90 $(pysource)

libpyimp_a_SOURCES = pyobjs.c

AM_FCFLAGS = -I$(top_builddir)/flib
AM_CPPFLAGS = -I. @GLIB_CFLAGS@ @PYTHON_INCLUDES@ -D'PYTHON_EXECDIR="$(pyexecdir)"'

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
flib_libraries = $(top_builddir)/src/libfutile-1.a
LIBS += @LIB_YAML_LIBS@
endif

# Dynamic library building
futile.so: $(libpyimp_a_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -shared $(LDFLAGS) $(AM_LDFLAGS) -Wl,-soname=$@ -o $@ $(libpyimp_a_OBJECTS) $(flib_libraries) $(LIBS) @PYTHON_LIBS@
	@chmod a+x $@

install-data-hook:
	if test -n "$(pymodule)" ; then \
	  cd $(pyexecdir) ; \
	  chmod a+x $(pymodule) ; \
	fi
