set terminal gif animate delay 10 size 1024,768
set output 'profile.gif'
set xrange [0:2]
set yrange [-0.6:0.6]
stats 'animation_data.txt' nooutput
do for [i=1:int(STATS_blocks)] {
    plot 'animation_data.txt' index (i-1) title '' with lines
}