set logscale x
plot "fft_columns.dat" using 1:($2/$8) title "1k"    w lp 1, \
     "fft_columns.dat" using 1:($3/$8) title "6k"    w lp 2, \
     "fft_columns.dat" using 1:($4/$8) title "10k"   w lp 3, \
     "fft_columns.dat" using 1:($5/$8) title "20k"   w lp 4, \
     "fft_columns.dat" using 1:($6/$8) title "50k"   w lp 5, \
     "fft_columns.dat" using 1:($7/$8) title "100k"  w lp 6, \
     "fft_columns.dat" using 1:($8/$8) title "500k"  w lp 7, \
     "fft_columns.dat" using 1:($9/$8) title "1000k" w lp 8

