set logscale x
plot "< awk '$1 == 1*1024'    fft_cache.dat" using 2:3 title "1k"    w lp 1, \
     "< awk '$1 == 6*1024'    fft_cache.dat" using 2:3 title "6k"    w lp 2, \
     "< awk '$1 == 10*1024'   fft_cache.dat" using 2:3 title "10k"   w lp 3, \
     "< awk '$1 == 20*1024'   fft_cache.dat" using 2:3 title "20k"   w lp 4, \
     "< awk '$1 == 50*1024'   fft_cache.dat" using 2:3 title "50k"   w lp 5, \
     "< awk '$1 == 100*1024'  fft_cache.dat" using 2:3 title "100k"  w lp 6, \
     "< awk '$1 == 500*1024'  fft_cache.dat" using 2:3 title "500k"  w lp 7, \
     "< awk '$1 == 1000*1024' fft_cache.dat" using 2:3 title "1000k" w lp 8

