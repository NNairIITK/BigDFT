## Process this file with automake to produce Makefile.in

OTHERS = PS_Check.f90 \
	PS_Program.f90 \
	PS_Exercise.f90 \
	PS_Integral.f90 \
	bench.csh \
	perfdata.t \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	COPYING \
	SUBDIRS \
	TODO \
	$(EXERCISE)

EXERCISE = Exercise/plot.gnuplot \
	Exercise/PS_Exercise.tex \
	Exercise/accF.pdf \
	Exercise/accS.pdf \
	Exercise/fandg.pdf \
	Exercise/PvsF.pdf \
	Exercise/rhoSpotS.pdf \
	Exercise/accP.pdf \
	Exercise/emepot.pdf \
	Exercise/rhoFpotF.pdf

EXTRA_DIST = PS_Check.ref.yaml tols.yaml

check_PROGRAMS = PS_Check PS_Program PS_Exercise 

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I$(top_builddir)/src/modules -I$(top_builddir)/src/PSolver -I$(top_srcdir)/src/PSolver -I. @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@  @MPI_INCLUDE@

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@


PS_Program_SOURCES = PS_Program.f90
PS_Program_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

PS_Exercise_SOURCES = PS_Exercise.f90
PS_Exercise_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@

#PS_Integral_SOURCES = PS_Integral.f90
#PS_Integral_LDADD = $(top_builddir)/src/libbigdft.a @LIBS_DEPENDENCIES@


CLEANFILES = malloc.prc time.prc PS_Check.out fldiff.report \
	module_fft_sg.@MODULE_EXT@ MODULE_FFT_SG.@MODULE_EXT@ \
	POISSON_SOLVER.@MODULE_EXT@ poisson_solver.@MODULE_EXT@ \
	module_base.@MODULE_EXT@

PS_Check.out.yaml: PS_Check PS_Program PS_Exercise 
	$(run_parallel) ./PS_Check 57 48 63  0 F  > PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63  1 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63 11 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 57 48 63 13 F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64  0 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64  1 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64 11 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 64 64 64 13 P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48  0 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48  1 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48 11 S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check 32 64 48 13 S >> PS_Check.out.yaml

check: PS_Check.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/PS_Check.ref.yaml -d PS_Check.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml

#	python $(top_srcdir)/tests/fldiff.py --mode=psolver --discrepancy=2.e-8 \
#		PS_Check.out $(srcdir)/PS_Check.out.ref | tee fldiff.report
