EXTRA_DIST = \
	NEB_include.sh \
	input \
	input.dft
	start.xyz \
	stop.xyz

psppar.H:
	ln -s $(top_srcdir)/utils/PSPfiles/$@

input.dft:
	ln -s $(srcdir)/$@

NEB_driver.sh:
	ln -s $(top_srcdir)/src/$@

check: psppar.H input.dft NEB_driver.sh NEB_include.sh
	if ! test -f input.dft ; then cp -f ./input.dft .; fi
	rm -f triH.NEB.it*
	$(top_builddir)/src/NEB < input | tee NEB.out
	python $(top_srcdir)/tests/fldiff.py NEB.out $(top_srcdir)/tests/Refs/NEB.out.ref | tee fldiff.report
	rm -rf triH.NEB.0*

CLEANFILES = triH.NEB.* velocities_file NEB.out NEB_driver.sh psppar.H gen_output_file fldiff.report
