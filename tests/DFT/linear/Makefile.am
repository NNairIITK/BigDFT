# Give here the list of existing tests
TESTDIRS = Locreg \
	   base \
	   mix \
	   cleanup \
       onlyhigh \
	   restartdisk

EXTRA_DIST = README \
	$(TESTDIRS) 

# Give here the pseudo-potentials used by each test.
base.psp: HGH/psppar.C extra/psppar.H
mix.psp: HGH/psppar.C extra/psppar.H
Locreg.psp : HGH-K/psppar.H HGH-K/psppar.Si

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  base.out*)                  prec="1.e-9" ;; \
      large.out*)                 prec="5.e-10" ;; \
      directmin.out*)             prec="1.1e-10" ;; \
	  standardnew.out*)           prec="2.e-10" ;; \
	  mix_base.out*)              prec="2.e-10" ;; \
	  mix_new.out*)               prec="1.1e-9" ;; \
	  mix_potential.out*)         prec="1.e-9" ;; \
	  onlyhigh.out*)            prec="1.1e-9" ;; \
	  cleanup.out*)            prec="1.1e-9" ;; \
          *)                          prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

include $(srcdir)/../../check.mk
