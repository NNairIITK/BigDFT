# Give here the list of existing tests
TESTDIRS = Locreg \
	   base \
	   mixpotential \
	   cleanup \
       onlyhigh \
	   restartdisk \
       directmin \
	   foe \
	   geopt \
	   testforces \
	   fructose

EXTRA_TESTDIRS = water \
                 molecule \
                 deltascf \
                 reformat \
		         paracetamol
#paracetamol move temporary

EXTRA_DIST = 

# Give here the pseudo-potentials used by each test.
base.psp: HGH/psppar.C extra/psppar.H
mixpotential.psp: HGH/psppar.C extra/psppar.H
Locreg.psp : HGH-K/psppar.H HGH-K/psppar.Si

#extra dependencies between tests in the same directory
restartdisk_read.out.out : restartdisk_charge.out.out

restartdisk_charge.out.post-out: 
	rm -rf data-read
	mv -f data-charge data-read

reformat_read.out.out : reformat_write.out.out

reformat_write.out.post-out: 
	rm -rf data-read
	mv -f data-write data-read

deltascf_read_charge.out.out : deltascf_write_neut.out.out

deltascf_write_neut.out.post-out: 
	rm -rf data-read_charge
	mv -f data-write_neut data-read_charge

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  base.out*)                  prec="1.e-9" ;; \
      large.out*)                     prec="5.e-10" ;; \
      directmin.out*)                 prec="1.1e-9" ;; \
	  standardnew.out*)           prec="2.e-10" ;; \
	  mixpotential.out*)          prec="1.e-9" ;; \
	  onlyhigh.out*)              prec="1.1e-9" ;; \
	  cleanup.out*)               prec="1.1e-8" ;; \
	  fructose.out*)              prec="9.0e-7" ;; \
	  geopt_diag.out*)            prec="2.e-9" ;; \
	  geopt_foe.out*)             prec="8.6e-7" ;; \
	  reformat_read.out*)         prec="1.e-9" ;; \
	  deltascf_read_charge.out*)  prec="1.e-9" ;; \
	  water_dirmin.out*)          prec="1.e-8" ;; \
	  water_mixdens.out*)         prec="1.1e-8" ;; \
	  water_mixpot.out*)          prec="1.1e-8" ;; \
	  water_foe.out*)             prec="1.e-8" ;; \
	  molecule_dmin.out*)         prec="1.e-9" ;; \
	  molecule_dmix.out*)         prec="1.e-9" ;; \
	  molecule_pmix.out*)         prec="1.1e-9" ;; \
	  molecule_foe.out*)          prec="1.e-9" ;; \
	  reformat_write.out*)        prec="3.e-10" ;; \
	  paracetamol.out*)           prec="2.e-10" ;; \
      *)                              prec="1.1e-10" ;; \
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
