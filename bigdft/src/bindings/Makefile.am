# To be installed:
noinst_LIBRARIES = libbindings.a

# Source code.
if BUILD_BINDINGS
bindings_files = bindings.c \
	bindings_regions.c \
	bindings_wavefunctions.c \
	bindings_localfields.c \
	bindings_signals.c bindings_signals.h \
	bindings_signals_inet.c
if HAVE_GDBUS
bindings_files += bindings_signals_dbus.c
nodist_bindings_files = bindings_dbus.c bindings_dbus.h 
endif
else
bindings_files = bindingsf_fake.f90
endif
libbindings_a_SOURCES = bindingsf.f90 bindings_api.h bindings.h input_keys.h \
	wrappers.c wrapper_inputs.c wrapper_run.c \
	wrapper_atoms.c $(bindings_files)
nodist_libbindings_a_SOURCES = $(nodist_bindings_files)

EXTRA_DIST = bindings_dbus.xml

# Dependencies
bindings_dbus.h: bindings_dbus.c
bindings_dbus.c: bindings_dbus.xml
	gdbus-codegen --interface-prefix eu.etsf.bigdft.dbus. \
                      --generate-c-code bindings_dbus \
                      --c-namespace BigdftDBus \
                      --c-generate-object-manager $(srcdir)/bindings_dbus.xml

AM_FCFLAGS = -I$(top_builddir)/includes @LIBABINIT_INCLUDE@
AM_CPPFLAGS = -I. -I../modules -I$(top_srcdir)/src/modules @PYTHON_INCLUDES@
AM_CFLAGS = @GLIB_CFLAGS@ @GDBUS_CFLAGS@

CLEANFILES = $(nodist_bindings_files)
