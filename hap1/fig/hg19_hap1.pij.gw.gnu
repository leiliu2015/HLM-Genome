apx = ARG1
load apx

fb  = ARG2
lnq = ARG3

unset key
set size square
set tics front out scale 0.3
set palette defined(0'white', 1'red')
set xtics font ',6' rotate by 270
set ytics font ',6'
set yrange [N-0.5:-0.5]

if (lnq == 1) {
  set cbrange [-3.7:-1.2]
  set cbtics in
} else {
  set cbrange [0:0.02]
  set cbtics in
}
load '../pals/inferno-inv.pal'

if (lnq == 1) {
  plot fa matrix u 1:2:($3>0 & $1<$2 ? log10($3) : NaN) w image, \
       fb matrix u 1:2:($3>0 & $1>$2 ? log10($3) : NaN) w image
} else {
  plot fa matrix u 1:2:($3>0 & $1<$2 ? $3 : NaN) w image, \
       fb matrix u 1:2:($3>0 & $1>$2 ? $3 : NaN) w image
}
