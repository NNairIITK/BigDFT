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

ATOM_OUTPUTS = \
			 atom.00.ae \
			 atom.01.ae \
			 input.pseudo

PSEUDO_OUTPUTS = \
				 orbitals.gplt

REPORTS = \
		  atom.out.report $(ATOM_OUTPUTS:=.report) \
		  pseudo.nofit.out.report $(PSEUDO_OUTPUTS:=.report) \
		  pseudo.shortfit.out.report

check: report

report: $(REPORTS)
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.report: %
	python $(abs_top_srcdir)/tests/fldiff.py --mode=pseudo $<  $(abs_top_srcdir)/pseudo/tutorials/oxygen.NLCC/ref/$<.ref | tee $@

$(ATOM_OUTPUTS): atom.out
$(PSEUDO_OUTPUTS): pseudo.nofit.out

atom.out: atom.in $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@

atom.in:
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/atom.dat .

pseudo.nofit.out: pseudo.nofit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.nofit.in: atom.out atom.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.GTH psppar.GTH
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/psppar.GTH psppar

pseudo.shortfit.out: pseudo.shortfit.in $(abs_top_builddir)/pseudo/src/pseudo
	$(run_parallel) $(abs_top_builddir)/pseudo/src/pseudo > $@

pseudo.shortfit.in: pseudo.nofit.out pseudo.nofit.out.report
	cp $(top_srcdir)/pseudo/tutorials/oxygen.NLCC/input.shortfit.pseudo input.pseudo

clean:
	if ! test $(top_srcdir)/pseudo/tutorials/oxygen.NLCC -ef . ; then rm -f atom.dat psppar.GTH ; fi
	rm -f *.out *.report input.* psppar ae.* *.ae \
			input.* local.pot *.gplt ps.* s.kernel.pot vertex.dump hostnames.*

