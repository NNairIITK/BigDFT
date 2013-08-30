EXTRA_DIST = \
		README \
		ref/atom.00.ae.ref\
		ref/atom.01.ae.ref\
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
		ref/input.fitpar.ref\
		ref/input.pseudo.ref\
		ref/input.weights.ref\
		ref/nlcc.gplt.ref\
		ref/orbitals.gplt.ref\
		ref/pseudo.nofit.out.ref \
		ref/pseudo.parallelnofit.out.ref \
		ref/pseudo.shortfit.out.ref \
		ref/psppar.atom.ref \
		ref/psppar.nlcc.short.ref

ATOM_OUTPUTS = \
			 atom.00.ae \
			 atom.01.ae \
			 input.pseudo \
			 input.weights \
			 input.fitpar \
			 psppar.atom

PSEUDO_OUTPUTS = \
				 orbitals.gplt \
				 nlcc.gplt

REPORTS = \
		  atom.out.report $(ATOM_OUTPUTS:=.report) \
		  pseudo.nofit.out.report $(PSEUDO_OUTPUTS:=.report) \
		  pseudo.shortfit.out.report \
		  pseudo.parallelnofit.out.report \
		  pseudo.NLCCnofit.out.report \
		  pseudo.NLCCfit.out.report \
		  pseudo.NLCCfree.out.report \
		  pseudo.NLCCfree2.out.report

check: report

report: $(REPORTS)
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.report: %
	python $(abs_top_srcdir)/tests/fldiff.py --mode=pseudo --discrepancy=1.e-4 \
			$<  $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/$<.ref | tee $@

$(ATOM_OUTPUTS): atom.out
$(PSEUDO_OUTPUTS): pseudo.nofit.out
nlcc.gplt: pseudo.NLCCnofit.out

atom.out: atom.in $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@
	cp psppar psppar.atom

atom.in:
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/atom.dat .

pseudo.nofit.out: pseudo.nofit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.nofit.in: atom.out atom.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.GTH psppar.GTH
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.GTH psppar

pseudo.shortfit.in: pseudo.nofit.out pseudo.nofit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.shortfit.pseudo input.pseudo

pseudo.shortfit.out: pseudo.shortfit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.parallelnofit.in: pseudo.shortfit.out pseudo.shortfit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.parallelnofit.pseudo  input.pseudo
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.parallelnofit.weights input.weights

pseudo.parallelnofit.out: pseudo.parallelnofit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.NLCCnofit.in: pseudo.parallelnofit.out pseudo.parallelnofit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.NLCCnofit.pseudo input.pseudo
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.addNLCC psppar.addNLCC
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.addNLCC psppar

pseudo.NLCCnofit.out: pseudo.NLCCnofit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.NLCCfit.in: pseudo.NLCCnofit.out pseudo.NLCCnofit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.NLCCfit.pseudo input.pseudo

pseudo.NLCCfit.out: pseudo.NLCCfit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.NLCCfree.in: pseudo.NLCCfit.out pseudo.NLCCfit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.NLCCfree.fitpar input.fitpar

pseudo.NLCCfree.out: pseudo.NLCCfree.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.NLCCfree2.in: pseudo.NLCCfree.out pseudo.NLCCfree.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.NLCCfree2.pseudo input.pseudo

pseudo.NLCCfree2.out: pseudo.NLCCfree2.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

clean:
	if ! test $(top_srcdir)/pseudo/tutorials/oxygen.NLCC -ef . ; then rm -f atom.dat psppar.GTH ; fi
	rm -f *.out *.report input.* psppar ae.* *.ae \
			input.* local.pot *.gplt ps.* s.kernel.pot vertex.dump hostnames.*

