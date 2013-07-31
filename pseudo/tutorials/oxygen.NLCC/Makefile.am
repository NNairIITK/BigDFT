EXTRA_DIST = \
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
	\rm -f *.out

report: atom.report
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

atom.out: $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@

%.report: %.out
	python $(abs_top_srcdir)/tests/fldiff.py $<  $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/$<.ref | tee $@


