import os
import sys
import numpy as np

def readMx(fx):
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        mx = []
        with open(fx) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    mx.append( list(map(float, lt)) )
        mx = np.array(mx)
    return mx

def xyzSample(K):
    d = np.sum(K, axis=0) + np.diag(K)
    L = np.diag(d) - K
    lam, Q = np.linalg.eigh(L)

    X = np.zeros((N, 3))
    for k in range(3):
        X[1:,k] = np.random.randn(N-1) / np.sqrt(lam[1:]) # sqrt(1.0/3/lam[1:]) ?

    R = np.dot(Q, X)
    R = R - np.mean(R, axis=0) # Necessary when diag(K) != 0
    return R

# Onami's method. Note that here GAM = 1/(2*\gamma) = <r^2_{ij}>,
# which is different from \gamma used in Hyeon's paper in 2019.
def K2M(K):
    d = np.sum(K, axis=0) + np.diag(K)
    L = np.diag(d) - K
    b = L[0:N-1,N-1]

    M = np.zeros((N, N))
    Q = np.linalg.inv( L[0:N-1,0:N-1]-b-np.reshape(b,(-1,1))+L[N-1,N-1] )
    M[0:N-1,0:N-1] = 0.5*(Q + np.transpose(Q))

    p = np.sum(M[0:N-1], axis=1)
    M[N-1,0:N-1] = -p
    M[0:N-1,N-1] = -p
    M[N-1,N-1] = np.sum(p)
    return M

# Calculate r^2_i based on explicit configurations without ad hoc confinement
if not len(sys.argv) == 4:
    print('usage:: python core_K2cfg_r2.py K_fit xxx.zi Ncfgs[>=0]')
    sys.exit()
fk = str(sys.argv[1])
fz = str(sys.argv[2])
ns = int(sys.argv[3])
np.random.seed(1274)

if not os.path.isfile(fz):
    print('Cannot find '+fz)
    sys.exit()
else:
    ic = [] # particle chain index
    iw = [] # particle index with missing
    with open(fz) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                ic.append( int(lt[0]) )
                iw.append( int(lt[2]) )

if not os.path.isfile(fk):
    print('Cannot find '+fk)
    sys.exit()
else:
    K = readMx(fk)
    N = len(K)
    S = np.diag(K2M(K)) # \sigma_{ii}
    ms= np.mean(S)

    if not N == len(iw):
        print("Check shapes of %s and %s"%(fz, fk))
        sys.exit()

    r2 = np.zeros(N)
    if ns > 0:
        pn = 1 if ns < 100 else int(ns/100.)
        for s in range(0, ns):
            R = xyzSample(K)
            r2 += np.sum(R**2, axis=1)
            if s%pn == 0:
                print("K2cfg_r2: %6d / %6d"%(s, ns))
        r2 = r2/ns
    else:
        r2+= np.nan

    fw = open(fk+'.r2', 'w')
    fw.write("#N_cfg: %d\n#c i i_withMissing r^2_i \sigma_{ii} log(\sigma_{ii}/<\sigma_{ii}>)\n"%(ns))
    for i in range(0, N):
        lt  = "%2d %6d %6d "%(ic[i], i, iw[i])
        lt += "%11s "%('NaN') if np.isnan(r2[i]) else "%11.5e "%(r2[i])
        lt += "%11.5e %+12.5e " % (S[i], np.log(S[i]/ms))
        fw.write(lt+'\n')
    fw.close()


