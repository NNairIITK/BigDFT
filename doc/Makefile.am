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
	Manual/manual.tex \
	Manual/sym_16.pdf \
	html/commun.css \
	html/doxygen.css \
	html/footer.html \
	html/header.html \
	html/logo_header.png \
	html/menu_fond.png \
	html/menu_footer.png \
	html/menu_header.png \
	html/original.css \
	html/point.png \
	html/sample.png


if BUILD_PDFDOC
pdfs = cluster.pdf \
	halFFT.pdf \
	manual.pdf
else
pdfs =
endif

if BUILD_DEVDOC
devdocs = Doxygen/html
else
devdocs =
endif

doc_DATA = $(pdfs) $(devdocs)

Doxygen/html: Doxygen.stamp web.stamp

Doxygen.stamp: $(top_builddir)/config/Doxyfile $(srcdir)/html/header.html $(srcdir)/html/footer.html
	\rm -rf Doxygen;\
	doxygen $(top_builddir)/config/Doxyfile; \
	touch Doxygen.stamp

web.stamp: Doxygen.stamp
	if [ ! -d Doxygen/images ] ; then mkdir Doxygen/images ; fi
	\cp $(srcdir)/html/*.png Doxygen/images
	if [ ! -d Doxygen/styles ] ; then mkdir Doxygen/styles ; fi
	\cp $(srcdir)/html/*.css Doxygen/styles
	touch web.stamp

cluster.pdf: cluster.tex compress.pdf figrid.pdf
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/cluster.tex

halFFT.pdf: halFFT.tex
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/halFFT.tex

manual.pdf: Manual/manual.tex Manual/grid100.pdf Manual/grid200.pdf Manual/grid300.pdf Manual/grid400.pdf Manual/sym_16.pdf
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/manual.tex
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/manual.tex

clean:
	rm -rf *.aux *.log cluster.pdf halFFT.pdf manual.pdf *.toc Doxygen *.stamp

