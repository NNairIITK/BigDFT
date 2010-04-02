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
	Manual/sym_16.pdf

if BUILD_PDFDOC
pdfs = \
	cluster.pdf \
	halFFT.pdf \
	manual.pdf
else
pdfs =
endif

doc_DATA = $(pdfs)

robodoc:
	export OLDPWD=`pwd`; cd $(top_srcdir); \
	pwd; robodoc --rc config/robodoc/robodoc-html.rc\
			--css ./config/robodoc/robodoc.css\
	        --doc $${OLDPWD}/www

cluster.pdf: cluster.tex compress.pdf figrid.pdf
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/cluster.tex

halFFT.pdf: halFFT.tex
	TEXINPUTS=$(top_srcdir)/doc:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/halFFT.tex

manual.pdf: Manual/manual.tex Manual/grid100.pdf Manual/grid200.pdf Manual/grid300.pdf Manual/grid400.pdf Manual/sym_16.pdf
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/manual.tex
	TEXINPUTS=$(top_srcdir)/doc/Manual:$${TEXINPUTS} pdflatex $(top_srcdir)/doc/Manual/manual.tex

clean:
	rm -f *.aux *.log cluster.pdf halFFT.pdf manual.pdf *.toc

