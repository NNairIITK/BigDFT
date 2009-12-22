EXTRA_DIST = README

art_source = bigdft_forces.f90 \
	defs.f90 \
	diis.f90 \
	find_saddle.f90 \
	initialize.f90 \
	initialize_potential.f90 \
	lanczos.f90 \
	min_converge.f90 \
	random.f90 \
	read_parameters.f90 \
	restart.f90 \
	saddle_converge.f90 \
	utils.f90 \
	write_refconfig.f90


AM_FCFLAGS = -I. -I$(srcdir) -I../modules -I../

noinst_LIBRARIES = libart.a

libart_a_SOURCES = $(art_source)

#dependencies
read_parameters.o: defs.o

bigdft_forces.o: defs.o

find_saddle.o: defs.o lanczos.o

initialize.o: defs.o

initialize_potential.o: defs.o bigdft_forces.o

lanczos.o: defs.o random.o  bigdft_forces.o

restart.o: defs.o

min_converge.o: defs.o bigdft_forces.o

saddle_converge.o: defs.o bigdft_forces.o

utils.o: defs.o

write_refconfig.o: defs.o

diis.o: defs.o bigdft_forces.o random.o
