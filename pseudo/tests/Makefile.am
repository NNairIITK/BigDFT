# Give here the list of existing tests
TESTDIRS = \
Ag_q1/atom\
Ag_q11/atom\
Al/atom\
Ar/atom\
As/atom\
At/atom\
Au_q1/atom\
Au_q11/atom\
B/atom\
Ba_q10/atom\
Ba_q2/atom\
Be_q2/atom\
Be_q4/atom\
Bi/atom\
Br/atom\
C/atom\
Ca_q10/atom\
Ca_q2/atom\
Cd_q12/atom\
Cd_q2/atom\
Ce_q12/atom\
Cl/atom\
Co_q17/atom\
Co_q9/atom\
Cr_q14/atom\
Cr_q6/atom\
Cs_q1/atom\
Cs_q9/atom\
Cu_q1/atom\
Cu_q11/atom\
Dy_q20/atom\
Er_q22/atom\
Eu_q17/atom\
F/atom\
Fe_q16/atom\
Fe_q8/atom\
Ga_q13/atom\
Ga_q3/atom\
Gd_q18/atom\
Ge/atom\
H/atom\
He/atom\
Hg_q12/atom\
Hg_q2/atom\
Ho_q21/atom\
I/atom\
In_q13/atom\
In_q3/atom\
Ir_q17/atom\
Ir_q9/atom\
K_q1/atom\
K_q9/atom\
Kr/atom\
La_q11/atom\
Li_q1/atom\
Li_q3/atom\
Lu_q25/atom\
Mg_q10/atom\
Mg_q2/atom\
Mn_q15/atom\
Mn_q7/atom\
Mo_q14/atom\
Mo_q6/atom\
N/atom\
Na_q1/atom\
Na_q9/atom\
Nb_q13/atom\
Nb_q5/atom\
Nd_q14/atom\
Ne/atom\
Ni_q10/atom\
Ni_q18/atom\
O/atom\
Os_q16/atom\
Os_q8/atom\
P/atom\
Pb/atom\
Pd_q10/atom\
Pd_q18/atom\
Pm_q15/atom\
Po/atom\
Pr_q13/atom\
Pt_q10/atom\
Pt_q18/atom\
Rb_q1/atom\
Rb_q9/atom\
Re_q15/atom\
Re_q7/atom\
Rh_q17/atom\
Rh_q9/atom\
Rn/atom\
Ru_q16/atom\
Ru_q8/atom\
S/atom\
Sb/atom\
Sc_q11/atom\
Sc_q3/atom\
Se/atom\
Si/atom\
Sm_q16/atom\
Sn/atom\
Sr_q10/atom\
Sr_q2/atom\
Ta_q13/atom\
Ta_q5/atom\
Tb_q19/atom\
Te/atom\
Ti_q12/atom\
Ti_q4/atom\
Tl_q13/atom\
Tl_q3/atom\
Tm_q23/atom\
V_q13/atom\
V_q5/atom\
W_q14/atom\
W_q6/atom\
Xe/atom\
Yb_q24/atom\
Y_q11/atom\
Y_q3/atom\
Zn_q12/atom\
Zn_q2/atom\
Zr_q12/atom\
Zr_q4/atom


# Generic part, do not touch.
CHECKS = $(TESTDIRS:=/atom.check)
CLEANS = $(TESTDIRS:=/atom.clean)

check: $(CHECKS) report

clean: $(CLEANS)

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

atom.out: $(abs_top_builddir)/pseudo/src/atom
	$(run_parallel) $(abs_top_builddir)/pseudo/src/atom > $@

%.report: %.ref
	python $(abs_top_srcdir)/tests/fldiff.py $*".out" $< | tee $@

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
        cd $$dir && $(MAKE) -f ../../Makefile atom.out
	touch $@

if USE_MPI
  mpirun_message=mpirun
else
  mpirun_message=
endif

EXTRA_DIST =
