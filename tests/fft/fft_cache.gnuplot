set term wxt 1
set logscale x
set title "FFT Time versus the FFT size for different cache size"
set xlabel "FFT size array"
set ylabel "Time (s)"
plot "< awk '$1 == 0*1024'    fft_cache.dat" using 2:4 title "no cache"    w lp 1, \
     "< awk '$1 == 6*1024'    fft_cache.dat" using 2:4 title "6kB"    w lp 2, \
     "< awk '$1 == 12*1024'   fft_cache.dat" using 2:4 title "12kB"   w lp 3, \
     "< awk '$1 == 24*1024'   fft_cache.dat" using 2:4 title "24kB"   w lp 4, \
     "< awk '$1 == 50*1024'   fft_cache.dat" using 2:4 title "50kB"   w lp 5, \
     "< awk '$1 == 75*1024'   fft_cache.dat" using 2:4 title "75kB"   w lp 6, \
     "< awk '$1 == 100*1024'  fft_cache.dat" using 2:4 title "100kB"  w lp 7


set term wxt 2
set logscale x
set title "Performance/Best Performance versus the FFT size for different cache size"
set xlabel "FFT size array"
set ylabel "Perf/Best"
plot "fft_columns.dat" using 1:($2/$9) title "no cache"    w lp 1, \
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
set title "Total time for FFT from 50 to 500 versus the cache size"
set xlabel "Cache size (kB)
set ylabel "Total time (s)"
plot "fft_perf.dat" using ($1/1024):2 notitle w lp 
