import os
import sys
import numpy as np

# mean/median radial postion of each chromosome
# ---
if not len(sys.argv) == 2:
    print('usage:: python core_K2cfg_r2Chrx.py xxx.r2')
    sys.exit()
fx = str(sys.argv[1])
fy = fx + 'Chrx'
chrNames = ["%d"%(i) for i in range(1,23)] + ['X']

if not os.path.isfile(fx):
    print('Cannot find ' + fx)
    sys.exit()
else:
    r2 = []
    cx = []
    with open(fx) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                cx.append( int(lt[0]) )
                r2.append( float(lt[4]) ) # \sigma_{ii}
    r2 = np.array(r2)
    cx = np.array(cx, dtype=int)

    fw = open(fy, 'w')
    L = max(cx) + 1
    for c in range(0, L):
        ix = np.where(cx==c)[0]

        r2_c = r2[ix]
        r2_i = np.where( np.isfinite(r2_c) )[0]
        rs = r2_c[r2_i]

        lt = "chr: %2s " % (chrNames[c])
        lt += "<r^2_i>: %11.5e " % (np.mean(rs))
        lt += "quantile: "
        for q in [5, 25, 50, 75, 95]:
            lt += "%11.5e " % (np.quantile(rs, q/100.))
        #print(lt)
        fw.write(lt+'\n')
    fw.close()


