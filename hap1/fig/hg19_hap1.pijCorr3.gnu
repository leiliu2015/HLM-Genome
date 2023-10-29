fx = ARG1

array logTics[10] = [0.00000,0.30103,0.47712,0.60206,0.69897,0.77815,0.84510,0.90309,0.95424,1.00000]

set tics front in scale 0.4
unset key
set size square

set xrange [0:500]
set yrange [0:500]
set xlabel 'p@^{Hi-C}_{ij}'
set ylabel 'p@^{Model}_{ij}'

load '../pals/ylgnbu.pal'
#load '/home/kapok/Downloads/gnuplotPalettes/jet.pal'
set cbtics in scale 1
set logscale cb

set xtics 0,100,500 nomirror
set ytics 0,100,500 nomirror
do for [c=0:5] {
  # major tics
  if (c<4) {
    set xtics add(sprintf("10^{%d}",c-5) c*100)
    set ytics add(sprintf("10^{%d}",c-5) c*100)
  } else {
    if (c==4) {
      set xtics add('0.1' c*100)
      set ytics add('0.1' c*100)
    } else {
      set xtics add('1' c*100)
      set ytics add('1' c*100)
    }
  }

  # minor tics
  if (c<5) {
    do for [t=2:9] {
      set xtics add('' c*100+logTics[t]*100)
      set ytics add('' c*100+logTics[t]*100)
    }
  }

  # grids
  if (c>0 & c<5) {
    set arrow from first c*100,graph 0 to first c*100,graph 1 nohead dt 3 front
    set arrow from graph 0,first c*100 to graph 1,first c*100 nohead dt 3 front
  }
}


set term wxt 3
plot fx matrix u 2:1:($3>0 ? $3 : NaN) w image


