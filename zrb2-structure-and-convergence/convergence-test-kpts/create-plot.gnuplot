set terminal pngcairo size 900,600
set output 'plot-zrb2-energy-vs-kpts.png'

set xlabel 'k-points N N N 0 0 0'
set ylabel 'Total Energy (Ry)'
set y2label 'Total Energy (eV)'
set title 'QE convergence test: zrb2 k-points'

# Set up right y-axis (y2) with appropriate scaling
set y2tics
set ytics nomirror

# Plot: left axis (Ry), right axis (eV)
plot 'energy-vs-kpts.dat' using 1:2 notitle with linespoints lw 3 pt 7 ps 2 lc "#00000000" axes x1y1,      '' using 1:3 notitle with lines lt -1 lc "#FFFF0000" axes x1y2
