# Give here the list of existing tests
TESTDIRS = Davidson-SiH4 \
	  TDDFT

EXTRA_DIST = \
	README

# Give here the pseudo-potentials used by each test.
Davidson-SiH4.psp: psppar.H psppar.Si
TDDFT.psp: psppar.N

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  TDDFT.out*)           prec="2.e-9" ;; \
          *)                    prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	if test -n "$$name" ; then \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log-$$name.yaml -t $(abs_srcdir)/../cubic/tols.yaml -o $@ ; \
	else \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log.yaml -t $(abs_srcdir)/../cubic/tols.yaml -o $@ ; \
	fi 

include $(srcdir)/../../check.mk
