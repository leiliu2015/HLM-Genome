set size ratio 0.2
set tics front
unset key

fx = ARG1
array chroN[23] = [ 250, 244, 199, 192, 181, 172, 160, 147, 162, 136, 136, 134, 116, 108, 103, 91, 82, 79, 60, 64, 49, 32, 156 ]


set yrange [-1:1]
set ytics -1,1,1; set mytics 1
#set grid xtics front
#set grid ytics front
do for [c=-1:1] {
  set ytics add(''c*1.)
}
###do for [c=0:6] {
###  set xtics add(''c*50)
###}
set xtics nomirror
set ytics nomirror
#set xlabel 'genomic position' offset 0,1
#set ylabel '' offset 2,0

mlw = 1
mps = 0.7
set arrow 1 from graph 0, first 0 to graph 1, first 0 nohead ls 8 lw 1 lc rgb 'black' front
set key left top
set border 3

###set xrange [0:250]
###set xtics 0,50,300; set mxtics 5
###set ylabel 'chr1' offset 3,0
###set label 1 'Lamin' at first 10, first 1 front
###set label 2 'radial' at first 10, first -1 front
###plot fx index 1 u 1:($3>0 ? $3/4.0 : NaN) w impulse ls 1 lw mlw notitle, \
###     fx index 1 u 1:($4>0 ? -2.*$4 : NaN) w impulse ls 2 lw mlw notitle

set term wxt 0 size 1200, 700
set multiplot layout 6,4
do for [c=0:22] {
  if (c==22) {
    set ylabel 'chrX' offset 3,0
  } else {
    set ylabel sprintf("chr%d",c+1) offset 3,0
  }

  xmx = chroN[c+1]
  set xrange [0:xmx]
  if (xmx < 100) {
    set xtics 0,20,100; set mxtics 2
  } else {
    set xtics 0,50,300; set mxtics 5
  }

  plot fx index c u 1:($3>0 ? $3/4.0 : NaN) w impulse ls 1 lw mlw notitle, \
       fx index c u 1:($4>0 ? -2.*$4 : NaN) w impulse ls 2 lw mlw notitle
}
unset multiplot
