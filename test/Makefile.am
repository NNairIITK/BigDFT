## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	inguess.dat \
	input.dat \
	posinp

PPFILES = \
	psppar.H_lda \
	psppar.Si_lda

CLEANFILES = wavefunction*

#dependencies
$(PPFILES):
	ln -s ../utils/PSPfiles/$<

mpif.h:
	touch mpif.h
	echo "real*8 :: MPI_double_precision" >> mpif.h
	echo "real*8 :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h

cluster.o: parameters.h
test_forces.o: parameters.h
