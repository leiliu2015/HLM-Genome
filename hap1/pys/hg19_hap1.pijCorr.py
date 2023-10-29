import os
import sys
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import gaussian_kde

def pij2ps(mx):
    N = len(mx)
    ps = np.zeros(N)
    for s in range(0, N):
        y = []
        for i in range(0, N-s):
            j = i+s
            if np.isfinite(mx[i,j]):
                y.append( mx[i,j] )
        if len(y) == 0:
            ps[s] = np.nan
        else:
            ps[s] = np.mean( np.array(y) )
    return ps

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

# Compare P^{fit}_{ij} with P^{obs}_{ij}
#
if not len(sys.argv) == 3:
    print('usage:: python hg19_hap1.pijCorr.py P_hic P_fit')
    sys.exit()
fhic = str(sys.argv[1])
ffit = str(sys.argv[2])

if True:
    # Read genome info
    fx = fhic[:-4] + '.gnu-bc.apx'
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        with open(fx) as f:
            for line in f:
                if 'chroN[' in line:
                    lt = line.strip().split()[4:-1]
                    Nc = []
                    for x in lt:
                        if x[-1] == ',':
                            Nc.append( int(x[:-1]) )
                        else:
                            Nc.append( int(x) )
                    break
        Nc = np.array(Nc, dtype=int)
        Nx = len(Nc)
        Zg = np.array([sum(Nc[:i]) for i in range(Nx)], dtype=int)

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
        P_obs = np.array(P_obs)
        N = len(P_obs)
        if not np.sum(Nc) == N:
            print('Check genome size')
            sys.exit()
        mask_f = (~np.isnan(P_obs)).astype(int) # 0/1: P^{obs}_{ij}=nan/finite
        idx_p = np.where( (np.triu(np.ones((N,N)),1) * mask_f)>0 ) # Matrix indices with j>i and P^{obs}_{ij} is finite

        for i in range(0, N):
            P_obs[i,i] = 1
        mask_f = (~np.isnan(P_obs)).astype(int) # 0/1: P^{obs}_{ij}=nan/finite
        idx_q = np.where( mask_f>0 )

        mask_intra = np.zeros((N,N))
        for c in range(0, Nx):
            i = Zg[c]
            j = i + Nc[c]
            mask_intra[i:j, i:j] = 1
        idx_intra = np.where( (mask_intra*mask_f)>0 )
        idx_inter = np.where( ((1-mask_intra)*mask_f)>0 )


    # Read P_fit
    if not os.path.isfile(ffit):
        print('Cannot find '+ffit)
        sys.exit()
    else:
        P_fit = []
        with open(ffit) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    P_fit.append( list(map(float, lt)) )
        P_fit = np.array(P_fit)
        for i in range(0, N):
            P_fit[i,i] = 1

if True: # P(s) by chain
    ps_obs = []
    ps_fit = []
    for c in range(0, Nx):
        i = Zg[c]
        j = i + Nc[c]
        ps_obs.append( pij2ps(P_obs[i:j, i:j]) )
        ps_fit.append( pij2ps(P_fit[i:j, i:j]) )
    #
    fw = open(ffit+'.pijCorr1', 'w')
    for c in range(0, Nx):
        ct = "#chain: %d N:%d\n#s P^{hic}(s) P^{fit}(s)\n" % (c, Nc[c])
        fw.write(ct)
        for s in range(0, Nc[c]):
            lt = "%5d "%(s)
            lt += "%11s "%('NaN') if np.isnan(ps_obs[c][s]) else "%11.5e "%(ps_obs[c][s])
            lt += "%11s "%('NaN') if np.isnan(ps_fit[c][s]) else "%11.5e "%(ps_fit[c][s])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()

if True: # Histogram(log10(p))
    xa = -5
    xb = 0
    dx = 0.01
    xz = xa - dx/2.
    nx = int((xb-xa)/dx + 1)
    xm = np.zeros((nx,nx))
    ct = "#min(log10(p)): %f max(log10(p)):: %f delta(log10(p)): %f shape: %d"%(xa, xb, dx, nx)

    xu = np.arange(nx)*dx+xz
    xv = xu + dx
    xc = (xu+xv)/2.0

    us = np.log10(P_obs[idx_q])
    vs = np.log10(P_fit[idx_q])
    xs = ((us-xz)/dx).astype(int)
    ys = ((vs-xz)/dx).astype(int)
    for i,j in np.column_stack((xs,ys)):
        if (0<=i<nx) and (0<=j<nx):
            xm[i,j] += 1
    xm = xm/np.sum(xm)/dx**2

###    ker_x = gaussian_kde(us)
###    kde_x = ker_x(xc)
###    ker_y = gaussian_kde(vs)
###    kde_y = ker_y(xc)

    fw = open(ffit+'.pijCorr2', 'w') # hist(p^{hic}) vs. hist(p^{obs})
    fw.write('#p pdf[log10(p^{hic})] pdf[log10(p^{fit})]\n')
    # fw.write('#p pdf[log10(p^{hic})] pdf[log10(p^{fit})] kde[log10(p^{hic})] kde[log10(p^{fit})]\n')
    for i in range(0, nx):
        lt = "%11.5e %11.5e %11.5e " % (10**xc[i], np.sum(xm[i,:])*dx, np.sum(xm[:,i])*dx)
        # lt+= "%11.5e %11.5e" % (kde_x[i],kde_y[i])
        fw.write(lt+'\n')
    fw.close()
    saveMx(ffit+'.pijCorr3', xm, ct+'\n') # hist(p^{hic}_{ij}) vs. hist(p^{obs}_{ij})

###    ker_2 = gaussian_kde(np.vstack([us, vs])) # TOO SLOW!
###    x, y = np.mgrid[0:nx, 0:nx]/dx
###    grid = np.vstack([x.ravel(), y.ravel()])
###    kde_2 = np.reshape(ker_2(grid).T, x.shape)
###    saveMx(ffit+'.pijCorr4', kde_2, ct+'\n') # kde(p^{hic}_{ij}) vs.kde(p^{obs}_{ij})


if True: # Pearsosn
    p_odiag = pearsonr(P_fit[idx_p], P_obs[idx_p])[0] # Pearson correlation w. diagonal p_{ii}
    p_wdiag = pearsonr(P_fit[idx_q], P_obs[idx_q])[0] # Pearson correlation o. diagonal p_{ii}
    p_intra = pearsonr(P_fit[idx_intra], P_obs[idx_intra])[0]
    p_inter = pearsonr(P_fit[idx_inter], P_obs[idx_inter])[0]
    pcs_lin = np.array([p_odiag, p_wdiag, p_intra, p_inter])

    p_odiag = pearsonr(np.log10(P_fit[idx_p]), np.log10(P_obs[idx_p]))[0] # Pearson correlation w. diagonal p_{ii}
    p_wdiag = pearsonr(np.log10(P_fit[idx_q]), np.log10(P_obs[idx_q]))[0] # Pearson correlation o. diagonal p_{ii}
    p_intra = pearsonr(np.log10(P_fit[idx_intra]), np.log10(P_obs[idx_intra]))[0]
    p_inter = pearsonr(np.log10(P_fit[idx_inter]), np.log10(P_obs[idx_inter]))[0]
    pcs_log = np.array([p_odiag, p_wdiag, p_intra, p_inter])

    pcs_names = ['odiag', 'wdiag', 'intra', 'inter']
    ct = '#lin'
    for k in range(4):
        ct += "pc_%s: %11.5e " % (pcs_names[k], pcs_lin[k])
    ct+= '\n'
    ct+= '#log'
    for k in range(4):
        ct += "pc_%s: %11.5e " % (pcs_names[k], pcs_log[k])
    ct+= '\n'
    fw = open(ffit+'.pijCorr4', 'w')
    fw.write(ct)
    fw.close()

if True: # Histogram(p)
    xa = 0
    xb = 1
    dx = 0.001
    nx = int((xb-xa)/dx)
    xm = np.zeros((nx,nx))
    ct = "#min(p): %f max(p):: %f delta(p): %f shape: %d"%(xa, xb, dx, nx)

    xu = np.arange(nx)*dx
    xv = xu + dx
    xc = (xu+xv)/2.0

    us = P_obs[idx_q]
    vs = P_fit[idx_q]
    xs = ((us-xa)/dx).astype(int)
    ys = ((vs-xa)/dx).astype(int)
    for i,j in np.column_stack((xs,ys)):
        if (0<=i<nx) and (0<=j<nx):
            xm[i,j] += 1
    xm = xm/np.sum(xm)/dx**2

    saveMx(ffit+'.pijCorr5', xm, ct+'\n') # hist(p^{hic}_{ij}) vs. hist(p^{obs}_{ij})

