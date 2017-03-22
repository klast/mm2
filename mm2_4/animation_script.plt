set terminal gif animate delay 100
set output 'velocity_animation.gif'
stats 'velocity_field.txt' nooutput

do for [i=1:int(STATS_blocks)] {
    plot 'velocity_field.txt' using 1:2:3:4 index (i-1) with vectors head filled lt 2
}

