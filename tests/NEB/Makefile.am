EXTRA_DIST = \
	NEB_include.sh \
	input \
	input.dft

psppar.H:
	ln -s $(top_srcdir)/utils/PSPfiles/$@

input.dft:
	ln -s $(srcdir)/$@

NEB_driver.sh:
	ln -s $(top_srcdir)/src/$@

check: psppar.H input.dft NEB_driver.sh NEB_include.sh
	cp -f $(srcdir)/input.dft .
	$(top_builddir)/src/NEB < input | tee NEB.out

CLEANFILES = triH.NEB.restart velocities_file NEB.out triH.NEB.it*.tar.bz2 NEB_driver.sh psppar.H gen_output_file input.dft

#python $(top_srcdir)/tests/fldiff.py --bigdft NaCl.out $(top_srcdir)/tests/Refs/NaCl.out.ref | tee fldiff.report

