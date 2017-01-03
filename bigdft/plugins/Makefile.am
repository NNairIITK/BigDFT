pluginsdir = $(libdir)/bigdft

noinst_LIBRARIES = libvol.a
plugins_DATA = $(noinst_LIBRARIES:.a=.so)

AM_FCFLAGS = -I$(top_builddir)/includes @FC_PIC_FLAG@ \
	@LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@ \
	@LIBETSFIO_INCLUDE@ @MPI_INCLUDE@ @LIB_FUTILE_CFLAGS@

libvol_a_SOURCES = vol.f90

SUFFIXES = .so

.a.so:
	$(FC) @FC_BUILD_SHARED@ -o $@ $($(<:.a=_a_OBJECTS))