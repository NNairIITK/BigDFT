modules_source = base.f90 \
	types.f90 \
	interfaces.f90

AM_FCFLAGS = -I. -I$(srcdir)

noinst_LIBRARIES = libmodules.a

libmodules_a_SOURCES = $(modules_source)

#dependencies
interfaces.o: types.o
types.o: base.o

