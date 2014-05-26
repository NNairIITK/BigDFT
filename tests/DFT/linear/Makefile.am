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
	   fructose \
       rotatesmall \
       cdft \
       hybrid-directmin \
       experimentalmode \
	   linearcubic \
	   experimentalfoe \
	   sparsity \
	   explicit_locregcenters \
	   DIIShistory \
	   PBE

EXTRA_TESTDIRS = reformatcomplete \
    water \
	molecule \
	deltascf \
	reformat \
	paracetamol \
	rotate
#\
#	noortho


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

rotatesmall_read1.out.out : rotatesmall_write.out.out
rotatesmall_read2.out.out : rotatesmall_write.out.out

rotatesmall_write.out.post-out: 
	rm -rf data-read1/data-sih3
	mkdir -p data-read1
	cp -rf data-write data-read1/data-sih3
	cp -f write.xyz data-read1/sih3.xyz
	rm -rf data-read2/data-sih3
	mkdir -p data-read2
	cp -rf data-write data-read2/data-sih3
	cp -f write.xyz data-read2/sih3.xyz

rotate_read1.out.out : rotate_write.out.out
rotate_read2.out.out : rotate_write.out.out

rotate_write.out.post-out: 
	rm -rf data-read1/data-benzene
	mkdir -p data-read1
	cp -rf data-write data-read1/data-benzene
	cp -f write.xyz data-read1/benzene.xyz
	rm -rf data-read2/data-benzene
	mkdir -p data-read2
	cp -rf data-write data-read2/data-benzene
	cp -f write.xyz data-read2/benzene.xyz

cdft_ohoh.out.out : cdft_oh-.out.out
cdft_oh-.out.out : cdft_oh-a.out.out

cdft_oh-a.out.post-out: 
	rm -rf data-ohoh/data-oh-a
	mkdir -p data-ohoh
	mv -f data-oh-a data-ohoh/data-oh-a
	cp -f oh-a.xyz data-ohoh/oh-a.xyz

cdft_oh-.out.post-out: 
	rm -rf data-ohoh/data-oh-
	mv -f data-oh- data-ohoh/data-oh-
	cp -f oh-.xyz data-ohoh/oh-.xyz

reformatcomplete_read36.out.out : reformatcomplete_read45.out.out
reformatcomplete_read45.out.out : reformatcomplete_read35.out.out
reformatcomplete_read35.out.out : reformatcomplete_write.out.out

reformatcomplete_read45.out.post-out: 
	rm -rf data-read36/data-h2o
	mkdir -p data-read36
	mv -f data-read45/data-h2o data-read36/data-h2o
	cp -f write.xyz data-read36/h2o.xyz
reformatcomplete_read35.out.post-out: 
	rm -rf data-read45/data-h2o
	mkdir -p data-read45
	mv -f data-read35/data-h2o data-read45/data-h2o
	cp -f write.xyz data-read45/h2o.xyz
reformatcomplete_write.out.post-out: 
	rm -rf data-read35/data-h2o
	mkdir -p data-read35
	mv -f data-write data-read35/data-h2o
	cp -f write.xyz data-read35/h2o.xyz

linearcubic_cubic.out.out : linearcubic_linear.out.out

linearcubic_linear.out.post-out:
	rm -rf data-cubic
	mv data-linear data-cubic


# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  base.out*)                  prec="1.e-9" ;; \
      large.out*)                 prec="5.e-10" ;; \
      directmin.out*)             prec="2.6e-9" ;; \
	  standardnew.out*)           prec="2.e-10" ;; \
	  mixpotential.out*)          prec="1.e-9" ;; \
	  onlyhigh.out*)              prec="1.1e-9" ;; \
	  cleanup.out*)               prec="1.1e-8" ;; \
	  fructose.out*)              prec="1.0e-6" ;; \
	  geopt_diag.out*)            prec="2.e-9" ;; \
	  geopt_foe.out*)             prec="9.2e-7" ;; \
	  reformat_read.out*)         prec="1.e-9" ;; \
	  deltascf_read_charge.out*)  prec="1.e-9" ;; \
	  water_dirmin.out*)          prec="1.e-8" ;; \
	  water_mixdens.out*)         prec="1.1e-8" ;; \
	  water_mixpot.out*)          prec="1.1e-8" ;; \
	  water_foe.out*)             prec="3.e-8" ;; \
	  molecule_dmin.out*)         prec="1.1e-8" ;; \
	  molecule_dmix.out*)         prec="1.1e-8" ;; \
	  molecule_pmix.out*)         prec="1.1e-8" ;; \
	  molecule_foe.out*)          prec="5.4e-7" ;; \
	  foe.out*)                   prec="5.6e-10" ;; \
	  reformat_write.out*)        prec="3.e-10" ;; \
	  paracetamol.out*)           prec="2.1e-10" ;; \
	  rotatesmall.out*)           prec="1.1e-8" ;; \
	  test.testforces.out*)       prec="1.1e-9" ;; \
	  hybrid-directmin.out*)      prec="1.3e-9" ;; \
	  experimentalmode.out*)      prec="1.6e-6" ;; \
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
