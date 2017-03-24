set xlabel "tau"
set ylabel "etta"
set yrange [0:100]
set xrange [0:100]
set palette rgbformulae 30,31,32
splot 'velocity.txt' matrix w pm3d ti ""