apx = ARG1
load apx

fb  = ARG2
idx = ARG3 + 0
jdx = ARG4 + 0
lnq = ARG5

N1 = chroN[idx]
Z1 = chroZ[idx]
N2 = chroN[jdx]
Z2 = chroZ[jdx]

set xrange [-0.5:N1-0.5]
set yrange [N2-0.5:-0.5]
set size ratio -1
set tics front out scale 0.3
#set xtics 0,40,320 nomirror
#set ytics 0,40,320 nomirror
#set mxtics 4
#set mytics 4
set xtics nomirror
set ytics nomirror
unset key
set palette defined(0'white', 1'red'); # load '/home/kapok/Downloads/gnuplotPalettes/jet.pal'

if (lnq == 1) {
  set cbrange [-3.7:-1.2]
  set cbtics in
} else {
  set cbrange [0:0.02]
  set cbtics in
}
load '../pals/inferno-inv.pal'

if (idx == jdx) {
  if (lnq == 1) {
    plot fa matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 & ($1-Z1)<($2-Z2) ? log10($3) : NaN) w image, \
         fb matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 & ($1-Z1)>($2-Z2) ? log10($3) : NaN) w image
  } else {
    plot fa matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 & ($1-Z1)<($2-Z2) ? $3 : NaN) w image, \
         fb matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 & ($1-Z1)>($2-Z2) ? $3 : NaN) w image
  }
} else {
  if (lnq == 1) {
    set term wxt 1; set title 'P_{obs}' offset 0,-0.5; plot fa matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 ? log10($3) : NaN) w image
    set term wxt 2; set title 'P_{fit}' offset 0,-0.5; plot fb matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 ? log10($3) : NaN) w image
  } else {
    set cbrange [1.0e-4:1.0e-3]; set cbtics 0,2.0e-4,1.0e-3
    set term wxt 1; set title 'P_{obs}' offset 0,-0.5; plot fa matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 ? $3 : NaN) w image pixels
    set term wxt 2; set title 'P_{fit}' offset 0,-0.5; plot fb matrix u ($1-Z1):($2-Z2):((Z1<$1<=(Z1+N1)) & (Z2<$2<=(Z2+N2)) & $3>0 ? $3 : NaN) w image pixels
  }
}



