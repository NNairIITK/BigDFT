set term wxt 1
set logscale x
plot "< awk '$1 == 0*1024'    fft_cache.dat" using 2:3 title "1k"    w lp 1, \
     "< awk '$1 == 6*1024'    fft_cache.dat" using 2:3 title "6k"    w lp 2, \
     "< awk '$1 == 12*1024'   fft_cache.dat" using 2:3 title "12k"   w lp 3, \
     "< awk '$1 == 24*1024'   fft_cache.dat" using 2:3 title "24k"   w lp 4, \
     "< awk '$1 == 50*1024'   fft_cache.dat" using 2:3 title "50k"   w lp 5, \
     "< awk '$1 == 75*1024'   fft_cache.dat" using 2:3 title "75k"   w lp 6, \
     "< awk '$1 == 100*1024'  fft_cache.dat" using 2:3 title "100k"  w lp 7


set term wxt 2
set logscale x
plot "fft_columns.dat" using 1:($2/$9) title "0k"    w lp 1, \
     "fft_columns.dat" using 1:($3/$9) title "6k"    w lp 2, \
     "fft_columns.dat" using 1:($4/$9) title "12k"   w lp 3, \
     "fft_columns.dat" using 1:($5/$9) title "20k"   w lp 4, \
     "fft_columns.dat" using 1:($6/$9) title "50k"   w lp 5, \
     "fft_columns.dat" using 1:($7/$9) title "75k"   w lp 6, \
     "fft_columns.dat" using 1:($8/$9) title "100k"  w lp 7

#plot "fft_columns.dat" using 1:($2/$9) title "0k"    w lp 1, \
#     "fft_columns.dat" using 1:($8/$9) title "100k"  w lp 7

set term wxt 3
unset logscale x
plot "fft_perf.dat" using 1:2 w lp
