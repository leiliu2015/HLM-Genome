import os
import sys
import time
import warnings
import cupy as cp
import numpy as np
from scipy.stats import pearsonr
cp.set_printoptions(precision=3, linewidth=200)
warnings.filterwarnings('ignore')

# An accelarated version of tino_nan.py with fdx = 0, ldx = 1, mdx = 4, and r^2_c = 1
# refs: PHi-C2 & core_nupy_potModified.py
#
if not len(sys.argv) == 2:
    print("usage:: python core_phic2_cupy.py normalized-HiC-Contact-Matrix")
    sys.exit()
fhic= str(sys.argv[1]) # HiC observation

# Initialize k_ij
def Init_K(K, N, INIT_K0):
    for i in range(1, N):
        j = i-1
        K[i,j] = K[j,i] = INIT_K0
    return K

# Onami's method. Note that here GAM = 1/(2*\gamma) = <r^2_{ij}>,
# which is different from \gamma used in Hyeon's paper in 2019.
def K2P(K):
    d = cp.sum(K, axis=0)
    L = cp.diag(d) - K

    Q = cp.linalg.inv(L[1:N,1:N])
    M = 0.5*(Q + cp.transpose(Q))
    A = cp.diag(M)

    G = cp.zeros((N, N))
    G[1:N,1:N] = -2*M + A + cp.reshape(A,(-1,1))
    G[0,1:N] = A
    G[1:N,0] = A
    P = (1.+3.*G)**(-1.5)
    return P

def Pdif2cost(P_dif):
    cost = cp.sqrt(cp.sum(P_dif**2)) / N
    return cost

def phic2(K, ETA=1.0e-4, ALPHA=1.0e-4, ITERATION_MAX=10000):
    stop_delta = ETA * ALPHA
    paras_fit = "%e\t%e\t%d\t" % (ETA, ALPHA, ITERATION_MAX)

    P_dif = K2P(K) - P_obs
    cost = Pdif2cost(P_dif)
    c_traj = cp.zeros((ITERATION_MAX+1, 2)); c_traj[0,0] = cost; c_traj[0,1] = time.time(); # {cost, time}

    iteration = 1
    while True:
        cost_bk = cost

        K = K - ETA*P_dif
        P_dif = K2P(K) - P_obs
        cost = Pdif2cost(P_dif)
        c_traj[iteration,0] = cost; c_traj[iteration,1] = time.time()

        cost_dif = cost_bk - cost
        if spq:
            print("%d\t%f\t%+12.5e" % (iteration, cost, -cost_dif))
        if (0 < cost_dif < stop_delta) or (iteration == ITERATION_MAX) or (cp.isnan(cost)):
            break
        iteration += 1
    c_traj = c_traj[:iteration+1]
    return [K, c_traj, paras_fit]

# Save an array of size n*m
def saveLg(fn, xy, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    n = len(xy)
    m = np.shape(xy)[1] if xy.ndim==2 else 0
    for i in range(n):
        if m == 0:
            lt = "%11s "%('NaN') if np.isnan(xy[i]) else "%11.5e"%(xy[i])
        else:
            lt = ''
            for v in xy[i]:
                lt += "%11s "%('NaN') if np.isnan(v) else "%11.5e "%(v)
        fw.write(lt+'\n')
    fw.close()

# Save a matrix of size n*n
def saveMx(fn, xy, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    n = len(xy)
    for i in range(n):
        lt = ''
        for v in xy[i]:
            lt += "%11s "%('NaN') if np.isnan(v) else "%11.5e "%(v)
        fw.write(lt+'\n')
    fw.close()

if True:
    # Read Hi-C
    if not os.path.isfile(fhic):
        print('Cannot find '+fhic)
        sys.exit()
    else:
        P_obs = []
        with open(fhic) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    P_obs.append( list(map(float, lt)) )
        P_obs = cp.array(P_obs)
        N = len(P_obs)
        cp.nan_to_num(P_obs, copy=False) # Replace NaN with 0
        P_obs = P_obs + cp.eye(N) # Set p_ii = 1

    # Minimization
    K_fit = cp.zeros((N, N))
    K_fit = Init_K(K_fit, N, INIT_K0=0.5)

    spq = False
    phic2_alpha = 1.0e-10
    K_fit, c_traj, paras_fit = phic2(K_fit, ETA=1.0e-4, ALPHA=phic2_alpha, ITERATION_MAX=1000000)
    P_fit = K2P(K_fit)

    # Save results
    c_traj = cp.asnumpy(c_traj)
    K_fit = cp.asnumpy(K_fit)
    P_fit = cp.asnumpy(P_fit)
    P_obs = cp.asnumpy(P_obs)

    dataDir = fhic[:fhic.rfind('.')] + "_phic2_a%7.1e_cupy"%(phic2_alpha)
    os.makedirs(dataDir, exist_ok=True)
    fo = "%s/N%d"%(dataDir, N)

    c_traj[:,1] = c_traj[:,1] - c_traj[0,1]
    saveLg(fo+'.log', c_traj, "#%s\n#cost systemTime\n"%(paras_fit))

    saveMx(fo+'.K_fit', K_fit, "#K_fit N %d min: %11.5e max: %11.5e\n"%(N, np.min(K_fit), np.max(K_fit)))

    triMask = np.where( np.triu(np.ones((N,N)),1)>0 ) # Matrix indices with j>i
    pijMask = np.where( np.triu(P_obs,1)>0 ) # Matrix indices with j>i and p_{ij}>0
    p1 = pearsonr(P_fit[triMask], P_obs[triMask])[0]
    p2 = pearsonr(P_fit[pijMask], P_obs[pijMask])[0]
    ct = "#P_fit N %d min: %11.5e max: %11.5e pearson: %11.5e %11.5e\n"%\
         (N, np.nanmin(P_fit), np.nanmax(P_fit), p1, p2)
    saveMx(fo+'.P_fit', P_fit, ct)


