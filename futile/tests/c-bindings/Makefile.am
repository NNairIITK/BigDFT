check_PROGRAMS = dict

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
flib_libraries = $(top_builddir)/src/libfutile-1.a
LIBS += @LIB_YAML_LIBS@
endif

dict_SOURCES = dict.c
dict_LDADD = $(flib_libraries) @GLIB_LIBS@

AM_CPPFLAGS = -I$(top_builddir)/c-bindings -I$(top_srcdir)/c-bindings @GLIB_CFLAGS@