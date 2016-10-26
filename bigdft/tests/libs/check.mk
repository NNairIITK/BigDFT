#write here a portable way to run only few test in the Makefile.
#taken as the third solution of the interesting webpage http://gallium.inria.fr/blog/portable-conditionals-in-makefiles/
#thanks to this trick we can test only few test by typing the command
# make check checkonly_that=C O2-Spin etc. CHECK_MODE=custom
#if CHECK_MODE has not been defined promote it to long check
#the following variable is true if we are inside a bzr branch, false otherwise
#if CHECK_MODE has not been defined promote it to long check
#the following variable is true if we are inside a bzr branch, false otherwise
CHECK_MODE_INTERNAL = $(shell if test -f $(top_srcdir)/branchfile; then echo true; else echo false;fi)

#this is the policy to be given in the case of explicit CHECK_MODE
checkonlyfoo_short_true= short
checkonlyfoo_short_false= short
checkonlyfoo_long_true= long
checkonlyfoo_long_false= long
#this is what would happen if the CHECK_MODE variable is not defined
checkonlyfoo__true= long
checkonlyfoo__false= short

checkonly_short=$(SHORT_CHECK)
checkonly_long=$(LONG_CHECK)
#this fixes the default value, if the CHECK_MODE is badly set
checkonly_=$(SHORT_CHECK)

present_checks := ${checkonly_${checkonlyfoo_${CHECK_MODE}_${CHECK_MODE_INTERNAL}}}

check: ${present_checks}
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
