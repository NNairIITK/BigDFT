SUBDIRS = src \
	tests

EXTRA_DIST = COPYING \
	README

check report:
	export PYTHONPATH=${PYTHONPATH}; export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ;python $(top_srcdir)/tests/report.py

