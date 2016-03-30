noinst_LIBRARIES = libbindings.a
include_HEADERS = dict.h futile.h futile_cst.h

libbindings_a_SOURCES = misc.c futile.h dict.c dict.h dictf.f90

AM_FCFLAGS = -I$(top_builddir)/flib
AM_CPPFLAGS = -I. @GLIB_CFLAGS@