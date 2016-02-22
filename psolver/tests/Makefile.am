## Process this file with automake to produce Makefile.in

OTHERS = bench.csh \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
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

EXTRA_DIST = PS_Check.ref.yaml tols.yaml $(OTHERS)

check_PROGRAMS = PS_Check PS_Program PS_Exercise Generalized_PSolver

AM_FCFLAGS = -I. -I$(top_builddir)/includes -I$(top_srcdir)/src \
			 @LIB_FLIB_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/PSolver/src

f_libraries = $(top_builddir)/wrappers/libwrappers.a
if BUILD_FLIB
f_libraries += $(top_builddir)/flib/src/libflib-1.a
else
f_libraries += @LIB_FLIB_CFLAGS@
endif

if BUILD_DYNAMIC_LIBS
psolver_libraries = -lPSolver-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
psolver_libraries = $(top_builddir)/src/libPSolver-1.a $(f_libraries)
LIBS += @LIBS_EXT_DEPS@
endif

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = $(psolver_libraries) @LIB_FLIB_LIBS@ @LINALG_LIBS@

PS_Program_SOURCES = PS_Program.f90
PS_Program_LDADD = $(psolver_libraries) @LIB_FLIB_LIBS@

PS_Exercise_SOURCES = PS_Exercise.f90
PS_Exercise_LDADD = $(psolver_libraries) @LIB_FLIB_LIBS@

Generalized_PSolver_SOURCES = Generalized_PSolver.f90
Generalized_PSolver_LDADD = $(psolver_libraries) @LIB_FLIB_LIBS@


CLEANFILES = malloc.prc time.prc PS_Check.out.yaml PSolver.report.yaml

PS_Check.out.yaml: PS_Check PS_Program PS_Exercise Generalized_PSolver
	-$(run_parallel) ./PS_Check -n [57,48,63] -g F > PS_Check.out.yaml
	-$(run_parallel) ./PS_Check -n [64,64,64] -g P >> PS_Check.out.yaml
	-$(run_parallel) ./PS_Check -n [32,64,48] -g S >> PS_Check.out.yaml
	-$(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PI -i '{environment: {input_guess: No}}'
	cat log.yaml >> PS_Check.out.yaml
	-$(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PCG -i '{environment: {input_guess: No}}'
	cat log.yaml >> PS_Check.out.yaml
	cp PS_Check.out.yaml tmp
	cat tmp | grep -v 'Unable to read mpd.hosts' > PS_Check.out.yaml
	rm tmp

check: PS_Check.out.yaml PS_Check.ref.yaml tols.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $(srcdir)/PS_Check.ref.yaml -d PS_Check.out.yaml -t $(srcdir)/tols.yaml \
	        --output PSolver.report.yaml

#	python $(top_srcdir)/tests/fldiff.py --mode=psolver --discrepancy=2.e-8 \
#		PS_Check.out $(srcdir)/PS_Check.out.ref | tee fldiff.report
