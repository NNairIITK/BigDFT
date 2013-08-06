# Give here the list of existing tests
TESTDIRS = \
Ag_q1/pade\
Ag_q11/pade\
Al/pade\
Ar/pade\
As/pade\
At/pade\
Au_q1/pade\
Au_q11/pade\
B/pade\
Ba_q10/pade\
Ba_q2/pade\
Be_q2/pade\
Be_q4/pade\
Bi/pade\
Br/pade\
C/pade\
Ca_q10/pade\
Ca_q2/pade\
Cd_q12/pade\
Cd_q2/pade\
Ce_q12/pade\
Cl/pade\
Co_q17/pade\
Co_q9/pade\
Cr_q14/pade\
Cr_q6/pade\
Cs_q1/pade\
Cs_q9/pade\
Cu_q1/pade\
Cu_q11/pade\
Dy_q20/pade\
Er_q22/pade\
Eu_q17/pade\
F/pade\
Fe_q16/pade\
Fe_q8/pade\
Ga_q13/pade\
Ga_q3/pade\
Gd_q18/pade\
Ge/pade\
H/pade\
He/pade\
Hg_q12/pade\
Hg_q2/pade\
Ho_q21/pade\
I/pade\
In_q13/pade\
In_q3/pade\
Ir_q17/pade\
Ir_q9/pade\
K_q1/pade\
K_q9/pade\
Kr/pade\
La_q11/pade\
Li_q1/pade\
Li_q3/pade\
Lu_q25/pade\
Mg_q10/pade\
Mg_q2/pade\
Mn_q15/pade\
Mn_q7/pade\
Mo_q14/pade\
Mo_q6/pade\
N/pade\
N/pade.nonrel\
N/pade.pol\
N/PBE\
N/PBE.nonrel\
N/PBE.pol\
N/VWN\
Na_q1/pade\
Na_q9/pade\
Nb_q13/pade\
Nb_q5/pade\
Nd_q14/pade\
Ne/pade\
Ni_q10/pade\
Ni_q18/pade\
O/pade\
Os_q16/pade\
Os_q8/pade\
P/pade\
Pb/pade\
Pd_q10/pade\
Pd_q18/pade\
Pm_q15/pade\
Po/pade\
Pr_q13/pade\
Pt_q10/pade\
Pt_q18/pade\
Rb_q1/pade\
Rb_q9/pade\
Re_q15/pade\
Re_q7/pade\
Rh_q17/pade\
Rh_q9/pade\
Rn/pade\
Ru_q16/pade\
Ru_q8/pade\
S/pade\
Sb/pade\
Sc_q11/pade\
Sc_q3/pade\
Se/pade\
Si/pade\
Sm_q16/pade\
Sn/pade\
Sr_q10/pade\
Sr_q2/pade\
Ta_q13/pade\
Ta_q5/pade\
Tb_q19/pade\
Te/pade\
Ti_q12/pade\
Ti_q4/pade\
Tl_q13/pade\
Tl_q3/pade\
Tm_q23/pade\
V_q13/pade\
V_q5/pade\
W_q14/pade\
W_q6/pade\
Xe/pade\
Yb_q24/pade\
Y_q11/pade\
Y_q3/pade\
Zn_q12/pade\
Zn_q2/pade\
Zr_q12/pade\
Zr_q4/pade


# Generic part, do not touch.
CHECKS = $(TESTDIRS:=/atom.check)
CLEANS = $(TESTDIRS:=/atom.clean)

check: $(CHECKS) check_and_report

clean: $(CLEANS)

check_and_report: $(CHECKS)
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

atom.out: $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@

%.report: %.out
	python $(abs_top_srcdir)/tests/fldiff.py --discrepancy=2.e-8 --mindigits=2 $< $<.ref | tee $@

%.clean:
	@dir=`dirname $@` ; \
        if test x"$(srcdir)" = x"." ; then \
          rm -f $$dir.* $$dir/*.out $$dir/*.report ; \
        else \
          rm -rf $$dir.* $$dir ; \
        fi ; \
        echo "Clean in "$$dir

%.check:
	@dir=`dirname $@` ; \
        if ! test x"$(srcdir)" = x"." ; then \
          if [ ! -d $$dir ] ; then mkdir -p $$dir ; fi ; \
          for i in $(srcdir)/$$dir/* ; do cp -f $$i $$dir; done ; \
        fi ; \
        cd $$dir && $(MAKE) -f ../../Makefile atom.report
	touch $@

