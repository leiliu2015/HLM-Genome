fx = ARG1

set tics front in scale 0.4
set key right top samplen 0

array chroName[23] = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
xmx = 300
ymi = 3.0e-4
set xrange [1:xmx]
set yrange [ymi:0.3]
set logscale xy 10
do for [i=-3:-2] {
  set ytics add(sprintf("10^{%d}",i) 10**(i))
}
set grid xtics back
set grid ytics back
mlw = 2
mlwf= 3
mps = 1

set term wxt 1 size 1600,800
set multiplot layout 4, 6
do for [c=0:22] {
  plot fx index c u 1:2 w p lw mlw pt 7 ps mps notitle, \
       fx index c u 1:3 w l lw mlwf t sprintf("chr%s",chroName[c+1])
}
unset multiplot


