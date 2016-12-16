SHORT_TESTDIRS = vol
LONG_TESTDIRS = $(SHORT_TESTDIRS)

EXTRA_DIST = 

vol-test.post-in:
	@cp $(abs_top_builddir)/plugins/libvol.so .

# Define the precision for specific directories.
%.report: %.ref
	python $(pythondir)/fldiff.py --mode=bigdft --discrepancy=1.1e-10 $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ; \
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

include $(srcdir)/../check.mk
