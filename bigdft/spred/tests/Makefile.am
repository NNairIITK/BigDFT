# Give here the list of existing tests
SHORT_TESTDIRS = MINHOP

LONG_TESTDIRS = $(SHORT_TESTDIRS) \
	NEB \
	SPLSAD MHGPS_LJ_CONNECT

EXTRA_DIST = MH_taskgroups


# Give here the pseudo-potentials used by each test.
NEB.psp: psppar.H
SPLSAD.psp: HGH/psppar.N
MINHOP.psp: psppar.Si psppar.H

# Additional MINHOP in;
MINHOP-test.post-in:
	rm -f posloc* posimed* data/global.mon

# Additional NEB in.
NEB-test.post-in:
	if ! test -f NEB_include.sh ; then sed "s;%%abs_top_builddir%%;$(abs_top_builddir);g" NEB_include.sh.in > NEB_include.sh ; fi
	if ! test -L NEB_driver.sh ; then ln -fs $(abs_top_srcdir)/spred/src/NEB_driver.sh ; fi
	rm -f neb.NEB.restart
NEB-test.post-clean:
	rm -f NEB_include.sh NEB_driver.sh

# Additional reference for geopt.mon test
global.mon.out : test.minhop.out
	cp -f data/global.mon $@

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref | *.splsad.ref | *.minhop.ref  ) mode="--mode=bigdft";; \
          *.NEB.ref) mode="--mode=neb";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  NEB*)             prec="5e-5" ;; \
	  SPLSAD*)          prec="2e-8" ;; \
	  MINHOP*)          prec="5e-9" ;; \
          *)                prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name2=`basename $< .xabs.ref.yaml` ; \
	if test "$$name2" = "$<" ; then \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	else \
	name=`basename $< .xabs.ref.yaml`  ; fi;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ; \
	python $(abs_top_builddir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@

##if USE_OCL
##oclrun_message = oclrun
##else
##oclrun_message =
##endif
##
##if USE_MPI
##mpirun_message = mpirun
##else
##mpirun_message =
##endif
##
##
##head_message:
##	@echo "========================================================"
##	@echo "To run the tests, type 'make check' (might take a while)"
##
##foot_message: $(mpirun_message) $(oclrun_message) head_message
##	@echo "========================================================"
##
##mpirun: head_message
##	@echo "     Use the environment variable run_parallel"
##	@echo "         ex: export run_parallel='mpirun -np 2'  "
##
##oclrun: head_message $(mpirun_message)
##	@echo "     Use the environment variable run_ocl"
##	@echo "         ex: export run_ocl='on' to use OpenCL acceleration"
##	@echo "     Use run_ocl='CPU' or 'ACC' to force use of hardware"
##	@echo "         different than GPU"
##	@echo "     and the environment variable ocl_platform"
##	@echo "         ex: export ocl_platform='NVIDIA'"
##	@echo "     and the environment variable ocl_devices"
##	@echo "         ex: export ocl_devices='K20'"
##
##all: head_message foot_message
##

include $(srcdir)/check.mk
