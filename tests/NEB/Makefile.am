EXTRA_DIST = \
	NEB_include.sh \
	input

psppar.H:
	ln -s $(top_srcdir)/utils/PSPfiles/$@

input.dat:
	ln -s $(srcdir)/$@

NEB_driver.sh:
	ln -s $(top_srcdir)/src/$@

check: psppar.H input.dat NEB_driver.sh NEB_include.sh
	cp -f $(srcdir)/input.dat .
	rm -f triH.NEB.it*
	$(top_builddir)/src/NEB < input | tee NEB.out
	python $(top_srcdir)/tests/fldiff.py NEB.out $(top_srcdir)/tests/Refs/NEB.out.ref | tee fldiff.report
	rm -rf triH.NEB.0*

CLEANFILES = triH.NEB.* velocities_file NEB.out NEB_driver.sh psppar.H gen_output_file fldiff.report
