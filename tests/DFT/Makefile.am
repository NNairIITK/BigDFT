## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  H \
	  C \
	  Li+ \
	  H2-spin \
	  SiH4 \
	  NH3 \
	  NaCl \
	  Ca2 \
	  Mn2 \
	  H2O-CP2K \
	  SiBulk \
	  H2O-slab \
	  Davidson-SiH4 \
	  Graphene \
	  FeHyb 


#This will run the automatic tests also for GPU
if USE_CUDA_GPU
SUBDIRS += GPU
endif

if USE_OCL
SUBDIRS += OCL
endif



if USE_MPI
mpirun_message=mpirun
else
mpirun_message=
endif


EXTRA_DIST = \
	README \
	inguess.dat

CLEANFILES =


check:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi
