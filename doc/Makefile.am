EXTRA_DIST = \
	basel-make.sh \
	coding-rules \
	halFFT.tex \
	cluster.tex \
	compress.ps \
	figrid.ps


pdf: cluster.pdf halFFT.pdf

# doc/robodoc-html will contain all the documentation files made by robodoc
robodoc:
	export OLDPWD=`pwd`; cd $(top_srcdir); \
	pwd; robodoc --rc config/robodoc/robodoc-html.rc\
			--css ./config/robodoc/robodoc.css\
	        --doc $${OLDPWD}/www

cluster.pdf: cluster.tex
	pdflatex $(top_srcdir)/doc/cluster.tex

halFFT.pdf: halFFT.tex
	pdflatex $(top_srcdir)/doc/halFFT.tex

clean:
	rm -f *.aux *.log

