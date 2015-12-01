EXTRA_DIST = \
	basel-make.sh \
	coding-rules \
	halFFT.tex \
	cluster.tex \
	compress.pdf \
	figrid.pdf \
	Manual/grid100.pdf \
	Manual/grid200.pdf \
	Manual/grid300.pdf \
	Manual/grid400.pdf \
	Manual/BigDFT-manual.tex \
	Manual/sym_16.pdf \
	html/logo_header.png \
	html/extra.css

CLEANFILES = *.aux *.toc *.log *.xml

if BUILD_PDFDOC
pdfs = cluster.pdf \
	halFFT.pdf \
	BigDFT-manual.pdf
else
pdfs =
endif

if BUILD_DEVDOC
devdocs = doxygen-doc
else
devdocs =
endif

doc_DATA = $(pdfs) $(devdocs)
noinst_DATA = tests.xml

doxygen-doc: Doxygen/html

Doxygen/html: Doxygen.stamp web.stamp

Doxygen.stamp: $(top_builddir)/config/Doxyfile
	\rm -rf Doxygen;\
	doxygen $(top_builddir)/config/Doxyfile; \
	touch Doxygen.stamp

web.stamp: Doxygen.stamp
	touch web.stamp

cluster.pdf: cluster.tex compress.pdf figrid.pdf
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/cluster.tex

halFFT.pdf: halFFT.tex
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/halFFT.tex

BigDFT-manual.pdf: Manual/BigDFT-manual.tex Manual/grid100.pdf Manual/grid200.pdf Manual/grid300.pdf Manual/grid400.pdf Manual/sym_16.pdf
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/BigDFT-manual.tex
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/BigDFT-manual.tex

tests.xml:
	@echo "<?xml version=\"1.0\" encoding=\"iso-8859-15\"?>" > tests.xml ; \
	echo "<bigdft-tests>" >> tests.xml ; \
	for dir in $(top_srcdir)/tests/* ; do \
	  if test -d $$dir -a \
                  $$(basename $$dir) != "styles" -a \
                  $$(basename $$dir) != "images" -a \
                  $$(basename $$dir) != "tutorials" ; then \
	    echo "<dir name=\"$$(basename $$dir)\">" ; \
	    if test ! -f $$dir/README ; then \
	      echo "not documented yet." ; \
	    else \
	      head -n1 $$dir/README ; \
	    fi ; \
	    for subdir in $$dir/* ; do \
	      if test -d $$subdir ; then \
	        echo "<test name=\"$$(basename $$subdir)\">" ; \
	        if test ! -f $$subdir/README ; then \
		  echo "not documented yet." ; \
		else \
		  head -n1 $$subdir/README ; \
		fi ; \
		for file in `find $$subdir -name input\*` ; do \
		  echo "  <input name=\"$$(basename $$file)\" />" ; \
		done ; \
		for file in `find $$subdir -name posinp\*` ; do \
		  echo "  <posinp name=\"$$(basename $$file)\" />" ; \
		done ; \
		echo "</test>" ; \
	      fi ; \
	    done ; \
	    echo "</dir>" ; \
	  fi ; \
	done >> tests.xml ; \
	echo "</bigdft-tests>" >> tests.xml


clean-local:
	rm -rf Doxygen lcov* *.aux *.log cluster.pdf halFFT.pdf BigDFT-manual.pdf *.toc Doxygen *.stamp doxygen.warnings

