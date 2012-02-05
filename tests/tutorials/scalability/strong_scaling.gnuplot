#!/usr/bin/gnuplot -persist
set style fill  pattern 1 border -1
set style rectangle back fc  lt -3 fillstyle  solid 1.00 border -1
set key title ""
set key outside right bottom vertical Left reverse enhanced autotitles columnhead nobox
set key invert samplen 4 spacing 1 width 0 height 0 
set style increment default
unset style line
unset style arrow
set style histogram rowstacked title  offset character 0, 0, 0
set pointsize 1
set style data histograms
set style function lines
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0
set ytics autofreq 
set y2tics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0
set y2tics autofreq 
set xlabel "No. of cores"
set ylabel "Percent" 
set y2label "Speedup" 
GNUTERM = "wxt"
	set y2range [1:]
#set size 1.4,1
ncores(mpi,omp)=mpi*((omp==0?1:omp))
set title "Run analysis, strong scaling"
plot [:] [:100] 'strong.dat' u (100*($4/$9)):xtic(1) ls 1 t 'Comms',\
	'' u (100*($6/$9)) ls 2 t 'LinAlg',\
	'' u (100*($5/$9)) ls 4 t 'Conv',\
	'' u (100*($8/$9)) ls 7 t 'Potential',\
	'' u (100*($7/$9)) ls 3 t 'Other',\
	'' u 0:($10/$9) w lp lt 2 linecolor 1 lw 3.5 pt 7 ps 0.9 t 'Speedup' axis x1y2 ,\
	'' u 0:(100*($10/$9/ncores($2,$3))*ncores($11,$12)) w lp lt 1 linecolor 3 lw 3.5 pt 7 ps 0.9 t 'Efficiency (%)'
#    EOF
