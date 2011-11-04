EXTRA_DIST = \
	NEB_include.sh \
	input \
	input.dft \
	start.xyz \
	stop.xyz \
	NEB.out.ref


CLEANFILES = triH.NEB.* velocities_file NEB.out NEB_driver.sh psppar.H gen_output_file fldiff.report malloc.prc


psppar.H:
	ln -s $(top_srcdir)/utils/PSPfiles/$@

NEB_driver.sh:
	ln -s $(top_srcdir)/src/$@

NEB.out: $(top_builddir)/src/NEB input NEB_include.sh psppar.H
	if ! test -f input.dft ; then cp -f $(srcdir)/input.dft .; fi
	rm -f triH.NEB.it*
	$(top_builddir)/src/NEB < input | tee NEB.out
	rm -rf triH.NEB.0*
	if ! test x"$(srcdir)" = x"." ; then rm -f input.dft ; fi

check: psppar.H input.dft NEB_driver.sh NEB_include.sh NEB.out
	python $(top_srcdir)/tests/fldiff.py --mode=neb --discrepancy=5e-5 NEB.out \
		$(top_srcdir)/tests/NEB/NEB.out.ref | tee fldiff.report
