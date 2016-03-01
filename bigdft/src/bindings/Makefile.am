# To be installed:
overridesdir = $(pyexecdir)/gi/overrides
noinst_LIBRARIES = libbindings.a
include_HEADERS = bigdft.h bigdft_input_keys.h
nodist_include_HEADERS = bigdft_cst.h
overrides_DATA = BigDFT.py

EXTRA_DIST = bigdft_cst.h.in

# Source code.
if BUILD_BINDINGS
bindings_files = bindings.c \
	bindings_regions.c \
	bindings_wavefunctions.c \
	bindings_localfields.c \
	bindings_signals.c bindings_signals.h \
	bindings_signals_inet.c
else
bindings_files = bindingsf_fake.f90
endif
libbindings_a_SOURCES = bindingsf.f90 bindings_api.h bindings.h input_keys.h \
	wrappers.c wrapper_inputs.c wrapper_run.c \
	wrapper_atoms.c $(bindings_files)
nodist_libbindings_a_SOURCES = $(nodist_bindings_files)

AM_FCFLAGS = -I$(top_builddir)/includes
AM_CPPFLAGS = -I. -I../modules -I$(top_srcdir)/src/modules @PYTHON_INCLUDES@
AM_CFLAGS = @GLIB_CFLAGS@

CLEANFILES = $(nodist_bindings_files)
