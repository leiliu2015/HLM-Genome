reset
fx = ARG1
mlwf = 3
mlw = 2
mps = 1.2
mpss = 1.0
# load '/home/kapok/Downloads/gnuplotPalettes/ylgnbu.pal'
set size ratio 0.64
set size ratio 0.60


set border linewidth 1.0
#set key right top samplen 0.1 spacing 1.0 maxrows 2 font "Arial-Italic"
set key out top samplen 0.1 maxrows 12

# Axes
set style line 11 lc rgb '#808080' lt 1
#set border 3 back #ls 11
set tics nomirror in scale 0.75
# Grid
set style line 12 lc rgb 'gray' dt (5,5) lw 1
set grid back ls 12

set xlabel 'Chromosome Size [Mb]' font "Arial-Italic, 18" offset 0,0
set xrange [30:270]
set xtics 0,50,300
set mxtics 5

set ylabel 'Squa. Rad. Dist.' offset 1,0 font "Arial-Italic, 18"
set ylabel 'Squ. Rad. Dist.' offset 1,0 font "Arial-Italic, 18"
set yrange [25:67]
set ytics 0,20,80
set mytics 4
#set ytics format "10^{%T}"

set bmargin at screen 0.27
set lmargin at screen 0.27
set rmargin at screen 0.87

###set bar 0
###plot fx u ($12/1000000.0):8:7:9 w error ls 3 lw mlw pt 4 ps 0 notitle, \
###     fx u ($12/1000000.0):8 w p ls 6 lw mlwf pt 6 ps mps notitle

set bar 0
plot fx u ($12/1000000.0):8:7:9 w error ls 3 lw mlw pt 4 ps 0 notitle, \
     fx u ($12/1000000.0):8 w p ls 6 lw mlwf pt 6 ps mps notitle, \
     fx u ($2==5  ? $12/1000000.0 : NaN):8 w p lc rgb 'cyan' pt 7 ps mpss notitle, \
     fx u ($2==18 ? $12/1000000.0 : NaN):8 w p lc rgb 'magenta' pt 7 ps mpss notitle, \
     fx u ($2==19 ? $12/1000000.0 : NaN):8 w p lc rgb 'orange' pt 7 ps mpss notitle





