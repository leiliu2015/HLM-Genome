import os
import sys
from numpy import *

def saveMxTx(fn, mx, ct):
    N, M = shape(mx)
    fw= open(fn, 'w')
    ct+= "#shape: %d %d " % (N, M)
    ct+= "min: %11.5e " % (nanmin(mx))
    ct+= "max: %11.5e " % (nanmax(mx))
    fw.write(ct+'\n')
    for i in range(0, N):
        lt = ''
        for j in range(0, M):
            lt += "%11s " % ('NaN') if isnan(mx[i,j]) else "%11.5e "%(mx[i,j])
        fw.write(lt+'\n')
    fw.close()
    return

# Calculate mean inter-chromosome contact probability 
# <p_{ij}>_{i \in chrA, j \in chrB} LL@July.28.2023
#
if not len(sys.argv) == 3:
    print('usage:: python hg19_hap1.pij.interChain.py xxx_wm.gnu-bc.apx xxx.P_fit.sbs')
    sys.exit()
fz = str(sys.argv[1])
fx = str(sys.argv[2])
fy = fx + '.interChain'

if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    cm = []
    with open(fx) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                cm.append( array(list(map(float, lt))) )
    cm = array(cm)

if not os.path.isfile(fz):
    print('Cannot find '+fz)
    sys.exit()
else:
    with open(fz) as f:
        for line in f:
            if 'chroN' in line:
                lt = line.strip().split()[4:-1]
                Nc = array([int(s.strip(',')) for s in lt])
                break
    L = len(Nc)
    Ng = sum(Nc)
    Zg = [sum(Nc[:c]) for c in range(L)]

    if not Ng == len(cm):
        print('Check matrix shape')
        sys.exit()

    P_interX = diag(nan*ones(L))
    for u in range(0, L):
        for v in range(0, L):
            if not u == v:
                ia = Zg[u]
                ib = Zg[u] + Nc[u]
                ja = Zg[v]
                jb = Zg[v] + Nc[v]
                nx = 0.0
                for i in range(ia, ib):
                    for j in range(ja, jb):
                        if isfinite(cm[i,j]) and isfinite(cm[j,i]):
                            P_interX[u,v] += cm[i,j]
                            nx += 1
                P_interX[u,v] /= float(nx)
    saveMxTx(fy, P_interX, '')

