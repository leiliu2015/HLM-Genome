fa = ARG1
fb = ARG2
N = 23

set xrange [-0.5:N-0.5]
set xtics 0,5,25
set mxtics 5
do for [c=0:4] {
  set xtics add(sprintf("%d",c*5+1) c*5)
}
set xtics offset 0,0.25
set xlabel 'chr.' font "Arial-Italic, 18" offset 0,0.2

set yrange [N-0.5:-0.5]
set ytics 0,5,25
set mytics 5
do for [c=0:4] {
  set ytics add(sprintf("%d",c*5+1) c*5)
}
set ytics offset 0.5,0
set ylabel 'chr.' font "Arial-Italic, 18" offset 1.5,0

unset key
set size square
set tics front out scale 0.75 nomirror
load '../pals/inferno-inv.pal'
unset colorbox

set bmargin at screen 0.27
set lmargin at screen 0.27
set rmargin at screen 0.57

plot fa matrix u 1:2:($1<$2 ? $3 : NaN) w image pixels notitle, \
     fb matrix u 1:2:($1>$2 ? $3 : NaN) w image pixels notitle


