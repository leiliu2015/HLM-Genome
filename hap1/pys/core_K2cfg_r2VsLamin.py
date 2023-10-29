import os
import sys
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Compare r^{2}_{i} with LaminB2_{i}
if not len(sys.argv) == 3:
    print('usage:: python core_K2cfg_r2VsLamin.py xxx.r2 lamin.txt')
    sys.exit()
fx = str(sys.argv[1])
fb = str(sys.argv[2])
chrNames = ["chr%d"%(i+1) for i in range(22)] + ['chrX']
fy = fx + 'VsLamin'

# Read in Lamin DamID signal
if not os.path.isfile(fb):
    print('Cannot find '+fb)
    sys.exit()
else:
    cx = []
    lx = []
    with open(fb) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                cx.append( int(lt[0]) )
                lx.append(float(lt[2]))
    cx = np.array(cx, dtype=int)
    lx = np.array(lx)
    L  = int(max(cx)+1)
    Nc = np.array([sum(cx==c) for c in range(L)], dtype=int)
    Zg = np.array([sum(Nc[:c]) for c in range(L)],dtype=int)
    Ng = sum(Nc)

# Read in <r^2_i>
if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    r2 = np.zeros(Ng) + np.nan
    with open(fx) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                r2[int(lt[2])] = float(lt[4])

    nr = np.zeros((Ng, 2)) # Normalize <r^2_i> genome-wide or per chain
    mx = np.nanmean(r2)
    nr[:,0] = np.log2( r2/mx )
    for c in range(0, L):
        ia = Zg[c]
        ib = Zg[c] + Nc[c]
        mx = np.nanmean(r2[ia:ib])
        nr[ia:ib,1] = np.log2( r2[ia:ib]/mx )

# <r^2_i> vs. Lamin


if True:
    mask = np.where(np.isfinite(lx) * np.isfinite(r2) == True)[0] # both finite

    notes = ["log2(r^2_i/<r^2>_{gw})", "log2(r^2_i/<r^2>_{pc})"]
    ct = ''
    for k in range(2):
        pc, pp =  pearsonr(lx[mask], nr[mask,k])
        sc, sp = spearmanr(lx[mask], nr[mask,k])
        ct += "#    gw pearson: %+12.5e %11.5e spearman: %+12.5e %11.5e " % (pc, pp, sc, sp) + ";%s\n"%(notes[k])
    #print(ct)
    fw = open(fy, 'w'); fw.write(ct)

    for c in range(0, L):
        ia = Zg[c]
        ib = Zg[c] + Nc[c]
        mask = np.where(np.isfinite(lx[ia:ib]) * np.isfinite(r2[ia:ib]) == True)[0] + Zg[c]

        pc, pp =  pearsonr(lx[mask], nr[mask,0])
        sc, sp = spearmanr(lx[mask], nr[mask,0])
        ct = "# %5s pearson: %+12.5e %11.5e spearman: %+12.5e %11.5e " % (chrNames[c], pc, pp, sc, sp)
        #print(ct)

        fw.write(ct+'\n'); fw.write("#i_pc i_gw laminB2 log2(r^2_i/<r^2>_{gw})\n")
        for i in range(ia, ib):
            lt  = "%5d %5d " % (i-Zg[c], i)
            lt += "%+12.5e "%(lx[i])   if np.isfinite(lx[i])   else "%12s "%('NaN')
            lt += "%+12.5e "%(nr[i,0]) if np.isfinite(nr[i,0]) else "%12s "%('NaN')
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()


