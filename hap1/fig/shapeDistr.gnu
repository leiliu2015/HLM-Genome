reset
fx = ARG1
mlwf = 3
mlw = 2
mps = 1.2
mpss = 1.0
load '../pals/ylgnbu.pal'
set size ratio 0.64
set size ratio 0.6


set border linewidth 1.0
#set key right top samplen 0.1 spacing 1.0 maxrows 2 font "Arial-Italic"
#set key out top samplen 0.1 maxrows 12
unset key
unset colorbox
set tics front

# Axes
set style line 11 lc rgb '#808080' lt 1
#set border 3 back #ls 11
set tics nomirror in scale 0.75
## Grid
#set style line 12 lc rgb 'gray' dt (5,5) lw 1
#set grid back ls 12

set xlabel 'Shape factor' font "Arial-Italic, 18" offset 0,0
set xrange [-0.1:0.25]
set xtics -0.1,0.1,0.3
set mxtics 10

set ylabel 'Asphericity' offset 1,0 font "Arial-Italic, 18"
set yrange [0:0.25]
set ytics 0,0.1,0.3
set mytics 10
#set ytics format "10^{%T}"

set bmargin at screen 0.27
set lmargin at screen 0.27
set rmargin at screen 0.87

#set bar 0
#plot fx u 15:8:14:16:7:9 w xyerrorbars ls 3 lw mlw pt 4 ps 0 notitle, \
#     fx u 15:8 w p ls 6 lw mlwf pt 6 ps mps notitle, \
#     fx u ($2==5  ? $15 : NaN):8 w p lc rgb 'cyan' pt 7 ps mpss notitle, \
#     fx u ($2==18 ? $15 : NaN):8 w p lc rgb 'magenta' pt 7 ps mpss notitle, \
#     fx u ($2==19 ? $15 : NaN):8 w p lc rgb 'orange' pt 7 ps mpss notitle

dx = 0.01
set arrow 1 from first 0, graph 0 to first 0, graph 1 nohead ls 1 lw mlwf lc rgb 'grey' dt 3 front
plot fx matrix u ($1*dx-0.25+0.5*dx):($2*dx+0.5*dx):(log10($3)) w image notitle




