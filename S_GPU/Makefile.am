SUBDIRS = \
src \
examples

#if BUILD_EXAMPLES
#SUBDIRS += examples
#endif

EXTRA_DIST = cudalt.py autogen.sh

#include doc/doxyinclude.am

if BUILD_DOC 
doxygen:
	$(DOXYGEN) $(DOXYFLAGS) doc/doxy.conf

.PHONY: doxygen

#html_DATA = doc/html/*.html doc/html/*.css doc/html/*.png doc/html/*.gif \
#	    doc/html/*.chm doc/html/*.qch
#pdf_DATA = doc/latex/refman.pdf
#ps_DATA = doc/latex/refman.ps
#dvi_DATA = doc/latex/refman.dvi
#EXTRA_DIST = $(html_DATA) $(pdf_DATA) $(ps_DATA)

endif # BUILD_DOC


ACLOCAL_AMFLAGS = -I m4
