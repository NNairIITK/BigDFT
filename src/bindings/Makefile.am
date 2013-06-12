# To be installed:
noinst_LIBRARIES = libbindings.a

# Source code.
if BUILD_BINDINGS
bindings_files = bindingsf.f90 \
	bindings_api.h \
	bindings.c bindings.h \
	bindings_atoms.c \
	bindings_regions.c \
	bindings_wavefunctions.c \
	bindings_localfields.c \
	bindings_signals.c bindings_signals.h \
	bindings_signals_inet.c
else
bindings_files = bindingsf_fake.f90
endif
if HAVE_GDBUS
bindings_files += bindings_signals_dbus.c
nodist_bindings_files = bindings_dbus.c bindings_dbus.h 
endif
libbindings_a_SOURCES = $(bindings_files)
nodist_libbindings_a_SOURCES = $(nodist_bindings_files)

EXTRA_DIST = bindings_dbus.xml

# Dependencies
bindings_dbus.h: bindings_dbus.c
bindings_dbus.c: bindings_dbus.xml
	gdbus-codegen --interface-prefix eu.etsf.bigdft.dbus. \
                      --generate-c-code bindings_dbus \
                      --c-namespace BigdftDBus \
                      --c-generate-object-manager $(srcdir)/bindings_dbus.xml

AM_FCFLAGS = @LIBABINIT_INCLUDE@ -I$(top_builddir)/includes
AM_CPPFLAGS = -I. -I../modules -I$(top_srcdir)/src/modules -I$(top_builddir)
AM_CFLAGS = @GLIB_CFLAGS@ @GDBUS_CFLAGS@

CLEANFILES = $(nodist_bindings_files)
