check_PROGRAMS = dict err

#we have to reduce the dict directory such as to move this test only there
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
else
#flib_libraries = $(top_builddir)/src/libfutile-1.a @LIBS_EXT_DEPS@
flib_libraries = $(top_builddir)/dicts/libdicts.a @LIBS_EXT_DEPS@
endif


dict_SOURCES = dict.c
dict_LDADD =  $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

err_SOURCES = err.c
err_LDADD =  $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

AM_CPPFLAGS = -I$(top_srcdir)/c-bindings -I$(top_builddir)/dicts -I$(top_srcdir)/dicts @GLIB_CFLAGS@
