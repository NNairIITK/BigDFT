## Process this file with automake to produce Makefile.in

OTHERS = overlap_ccs.txt \
		 hamiltonian_ccs.txt \
		 density_kernel_ccs.txt \
		 density_kernel_matmul_ccs.txt \
		 sparsematrix_metadata.dat \
		 mat_114_8986.txt \
		 mat_600_33620.txt \
		 mat_1202_77462.txt \
	     $(REFERENCES)

REFERENCES= driver_random.ref.yaml driver_foe.ref.yaml smatmul.ref.yaml

#EXERCISE = Exercise/plot.gnuplot \
#	Exercise/PS_Exercise.tex \
#	Exercise/accF.pdf \
#	Exercise/accS.pdf \
#	Exercise/fandg.pdf \
#	Exercise/PvsF.pdf \
#	Exercise/rhoSpotS.pdf \
#	Exercise/accP.pdf \
#	Exercise/emepot.pdf \
#	Exercise/rhoFpotF.pdf

EXTRA_DIST =  tols.yaml $(OTHERS)

check_PROGRAMS = driver_ccs driver_foe driver_single driver_random driver_eigenvalues smatmul

AM_FCFLAGS = -I. -I$(top_srcdir)/src -I$(top_builddir)/src @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
chess_library = $(top_builddir)/src/libCheSS-1.so.@CHESS_MINOR_VERSION@
else
chess_library = $(top_builddir)/src/libCheSS-1.a @LIB_CHESS_DEPS@
endif

driver_ccs_SOURCES = utilities.f90 driver_ccs.f90
driver_ccs_LDADD = $(chess_library) @LINALG_LIBS@

driver_SOURCES = driver.f90

driver_foe_SOURCES = utilities.f90 driver_foe.f90
driver_foe_LDADD = $(chess_library) @LINALG_LIBS@

driver_single_SOURCES = utilities.f90 driver_single.f90
driver_single_LDADD = $(chess_library) @LIB_FUTILE_LIBS@ @LINALG_LIBS@

driver_random_SOURCES = utilities.f90 driver_random.f90
driver_random_LDADD = $(chess_library) @LINALG_LIBS@

driver_eigenvalues_SOURCES = utilities.f90 driver_eigenvalues.f90
driver_eigenvalues_LDADD = $(chess_library) @LINALG_LIBS@

smatmul_SOURCES = smatmul.f90
smatmul_LDADD = $(chess_library) @LINALG_LIBS@

if BUILD_DYNAMIC_LIBS
driver_ccs_LDADD += @LIB_FUTILE_LIBS@
driver_foe_LDADD += @LIB_FUTILE_LIBS@
driver_single_LDADD += @LIB_FUTILE_LIBS@
driver_random_LDADD += @LIB_FUTILE_LIBS@
driver_eigenvalues_LDADD += @LIB_FUTILE_LIBS@
smatmul_LDADD += @LIB_FUTILE_LIBS@
endif

CLEANFILES = *.out.yaml *.report.yaml report_remarks *.mod time.yaml report

check: regtest report #PS_Basics.out.yaml PS_Check.out.yaml PS_Check.ref.yaml tols.yaml report

#test-random.report.yaml: $(srcdir)/driver_random.ref.yaml driver_random.out.yaml
#	python $(pythondir)/fldiff_yaml.py \
#		-r $< -d driver_random.out.yaml -t $(scrdir)/tols.yaml --output $@ 

regtest: driver_random driver_foe smatmul
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py \
	       --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir) \
	       --exclusive=${CHECK_ONLY}

report: driver_random.report.yaml driver_foe.report.yaml
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi


define F_REGTEST_INSTRUCTIONS
- driver_random:
   binary: driver_random
   runs:
    - $(run_parallel) ./driver_random \
        --sparsegen_method=random \
        --matgen_method=random \
        --nfvctr=100 \
        --nvctr=10000 \
        --condition_number=50 \
        --expo=-1 \
        --do_consistency_checks=.true.
    - $(run_parallel) ./driver_random \
         --sparsegen_method=random \
         --matgen_method=random \
         --nfvctr=2000 \
         --nvctr=100000 \
         --condition_number=100 \
         --expo=-0.5 \
         --do_consistency_checks=.true.
- driver_foe:
   binary: driver_foe
   runs:   
    - $(run_parallel) ./driver_foe \
        --sparsity_format=ccs \
        --matrix_format=serial_text \
        --metadata_file=$(srcdir)/sparsematrix_metadata.dat \
        --overlap_file=$(srcdir)/overlap_ccs.txt \
        --hamiltonian_file=$(srcdir)/hamiltonian_ccs.txt \
        --kernel_file=$(srcdir)/density_kernel_ccs.txt \
        --kernel_matmul_file=$(srcdir)/density_kernel_matmul_ccs.txt \
        --check_spectrum=yes \
		--betax=-1000
- smatmul:
   binary: smatmul
   runs:
    - $(run_parallel) ./smatmul \
        --filename=$(srcdir)/mat_114_8986.txt \
        --nit=10
    - $(run_parallel) ./smatmul \
        --filename=$(srcdir)/mat_600_33620.txt \
        --nit=100 \
        --verbosity=1
    - $(run_parallel) ./smatmul \
        --filename=$(srcdir)/mat_1202_77462.txt \
        --nit=20 --verbosity=1
- driver_eigenvalues:
   binary: driver_eigenvalues
   runs:
    - $(run_parallel) ./driver_eigenvalues \
        --metadata_file=$(srcdir)/sparsematrix_metadata-evals.dat \
        --overlap_file=$(srcdir)/overlap_sparse-evals.txt \
        --hamiltonian_file=$(srcdir)/hamiltonian_sparse-evals.txt \
        --kernel_file=$(srcdir)/density_kernel_sparse-evals.txt \
        --kernel_matmul_file=$(srcdir)/density_kernel_sparse_matmul-evals.txt \
        --iev_min=151 \
        --iev_max=156 \
        --fscale=1.e-2
endef
export F_REGTEST_INSTRUCTIONS

