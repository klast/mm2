set terminal gif animate delay 10 size 1024,768
set output 'profile.gif'
set xrange [0:41]
set yrange [0:1.1]
stats 'delta.txt' nooutput
do for [i=1:int(STATS_blocks)] {
    plot 'delta.txt' index (i-1) title '' with lines
}