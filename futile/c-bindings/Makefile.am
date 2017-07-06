noinst_LIBRARIES = libbindings.a
include_HEADERS = tree.h futile.h futileC_objects.h
#nodist_include_HEADERS  = futile_cst.h

libbindings_a_SOURCES = tree.c tree.h misc.c futile.h python.c futileC_objects.c treef.f90

AM_FCFLAGS = -I$(top_builddir)/dicts -I$(top_builddir)/flib
AM_CPPFLAGS = -I. -I$(top_builddir) -I$(top_builddir)/dicts -I$(top_srcdir)/dicts @GLIB_CFLAGS@ @PYTHON_INCLUDES@ -D'PYTHON_EXECDIR="$(pyexecdir)"'
AM_CFLAGS = -I$(top_builddir)/dicts
