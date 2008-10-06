EXTRA_DIST = \
	basel-make.sh \
	coding-rules \
	halFFT.tex \
	cluster.tex \
	compress.ps \
	figrid.ps


# doc/robodoc-html will contain all the documnentation files made by robodoc
robodoc:
	export OLDPWD=`pwd`; cd $(top_srcdir); \
	pwd; robodoc --rc config/robodoc/robodoc-html.rc\
			--css ./config/robodoc/robodoc.css\
	        --doc $${OLDPWD}/www

cluster.pdf: cluster.tex
	pdflatex $(top_srcdir)/doc/cluster.tex

HalFFT.pdf: HalFFT.tex
	pdflatex $(top_srcdir)/doc/HalFFT.tex

clean:
	\rm *.aux *.log

