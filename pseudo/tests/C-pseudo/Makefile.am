EXTRA_DIST = \
			 input.fitpar \
			 psppar.HGH


REPORTS = \
		  atom.out.report \
		  pseudo.out.report

check: report

report: $(REPORTS)
	@if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

%.report: %
	python $(abs_top_srcdir)/tests/fldiff.py --mode=pseudo --discrepancy=1.e-4 \
			$<  $(abs_top_srcdir)/pseudo/tests/C-pseudo/$<.ref | tee $@

atom.in:
	rm -f input.weights
	cp $(top_srcdir)/pseudo/tests/C-pseudo/atom.dat .

atom.out: atom.in $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@
	mv psppar psppar.atom

pseudo.in: atom.out atom.out.report
	cp $(top_srcdir)/pseudo/tests/C-pseudo/input.pseudo .
	cp $(top_srcdir)/pseudo/tests/C-pseudo/input.fitpar .
	cp $(top_srcdir)/pseudo/tests/C-pseudo/psppar.HGH psppar

pseudo.out: pseudo.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

clean:
	if ! test $(top_srcdir)/pseudo/tests/C-pseudo -ef . ; then \
			rm -f atom.dat input.pseudo input.fitpar ; fi
	rm -f psppar *.out *.report ae.* *.ae psppar.* input.weights \
			local.pot *.gplt ps.* s.kernel.pot vertex.dump hostnames.*

