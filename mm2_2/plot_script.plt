reset
set xlabel "z"
set ylabel "t"
set yrange [0:10]
set xrange [0:50]
set palette rgbformulae 30,31,32
set hidden3d
set pm3d interpolate 2,2
splot 'answer.txt' matrix w pm3d ti ""
