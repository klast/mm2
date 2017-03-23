set terminal gif animate delay 20 size 1024,768
set output 'velocity.gif'
stats 'velocity.txt' nooutput
do for [i=1:int(STATS_blocks)] {
    plot 'velocity.txt' using 1:2:3:4 index (i-1) with vectors head filled lt 2 title ''
}

set terminal gif animate delay 20 size 1024,768
set output 'velocity_x.gif'
stats 'velocity_x.txt' nooutput
do for [i=1:int(STATS_blocks)] {
    plot 'velocity_x.txt' using 1:2:3:4 index (i-1) with vectors head filled lt 2
}

set terminal gif animate delay 20 size 1024,768
set output 'velocity_y.gif'
stats 'velocity_y.txt' nooutput
do for [i=1:int(STATS_blocks)] {
    plot 'velocity_y.txt' using 1:2:3:4 index (i-1) with vectors head filled lt 2
}