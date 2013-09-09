EXTRA_DIST = \
			 input.fitpar \
			 input.weights


REPORTS = \
		  atom.out.report \
		  pseudo.out.report

check: report

report: $(REPORTS)
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.report: %
	python $(abs_top_srcdir)/tests/fldiff.py --mode=pseudo --discrepancy=1.e-4 \
			$<  $(abs_top_srcdir)/pseudo/tests/H-pseudo/$<.ref | tee $@

atom.in:
	cp $(top_srcdir)/pseudo/tests/H-pseudo/atom.dat .

atom.out: atom.in $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@
	mv psppar psppar.atom
	mv input.weights input.weights.atom

pseudo.in: atom.out atom.out.report
	cp $(top_srcdir)/pseudo/tests/H-pseudo/input.pseudo .
	cp $(top_srcdir)/pseudo/tests/H-pseudo/input.fitpar .
	cp $(top_srcdir)/pseudo/tests/H-pseudo/input.weights .
	cp $(top_srcdir)/pseudo/tests/H-pseudo/psppar.GTH psppar

pseudo.out: pseudo.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

clean:
	if ! test $(top_srcdir)/pseudo/tests/H-pseudo -ef . ; then rm -f atom.dat input.pseudo ; fi
	rm -f *.out *.report input.* psppar ae.* *.ae \
			input.* local.pot *.gplt ps.* s.kernel.pot vertex.dump hostnames.*

