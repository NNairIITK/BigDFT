set term png size 300,200 ; set ou 'avg_y.png'
#set term png "Bold" 18 ; set ou 'avg_y.png'
#set xl 'molecule axis (origined at Oxsygen) (bohr)'
set xl 'Water molecule axis '
set yl 'Averaged elec. density'
unset xtics; unset ytics
p './electronic_density_avg_y' u ($2-8.42):3 lw 2 w lp t '' 
