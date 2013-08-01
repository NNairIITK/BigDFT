EXTRA_DIST = \
		README \
		ref/atom.out.ref\
		ref/fit.fixed.NLCC.out.ref\
		ref/fit.fixed.NLCC.proc.01.ref\
		ref/fit.free.NLCC.out.ref\
		ref/fit.free.NLCC.proc.01.ref\
		ref/fit.short.NLCC.out.ref\
		ref/fit.short.NLCC.proc.01.ref\
		ref/fit.without.FITPAR.ref\
		ref/fit.without.NLCC.out.ref\
		ref/fit.without.NLCC.proc.01.ref\
		ref/GTH.closedshell.out.ref\
		ref/GTH.polarized.out.ref\
		ref/nlcc.ref\
		ref/nlcc.short.ref\
		ref/psppar.nlcc.short.ref\
		ref/psppar.ref

check: report

clean:
	if ! test $(top_srcdir)/pseudo/tutorials/oxygen.NLCC -ef . ; then rm -f atom.dat ; fi
	rm -f *.out *.report input.* psppar* ae.* *.ae

report: atom.report
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

atom.out: atom.in $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@

atom.in:
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/atom.dat .

atom.report: atom.out
	python $(abs_top_srcdir)/tests/fldiff.py $<  $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/$<.ref | tee $@
	python $(abs_top_srcdir)/tests/fldiff.py atom.00.ae $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/atom.00.ae | tee atom.00.report
	python $(abs_top_srcdir)/tests/fldiff.py atom.01.ae $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/atom.01.ae | tee atom.01.report


