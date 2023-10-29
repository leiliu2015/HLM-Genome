import os
import sys
import h5py
from numpy import *

# distribution of asp
# ref: /home/kapok/Downloads/hap1/core/core_K2cfgs.py
#
if not len(sys.argv) == 3:
    print('usage:: python core_K2xyz-shape.py xxx.K_fit NSamples')
    sys.exit()
fk = str(sys.argv[1])
ncfgs = int(sys.argv[2])
rds= 1274; random.seed(rds)
fy = fk[:-6] + ".xyz-shape.asp"

# Read K-matrix
if True:
    if not os.path.isfile(fk):
        print('Cannot find '+fk)
        sys.exit()
    else:
        K_fit = []
        with open(fk) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    K_fit.append( list(map(float, lt)) )
        K_fit = array(K_fit)
        Ng = len(K_fit)

    # K to Laplacian matrix
    d = sum(K_fit, axis=0) + diag(K_fit)
    Lap = diag(d) - K_fit
    # Eigenvalues and eigenvectors
    lam, Qs = linalg.eigh(Lap)

    #fo = h5py.File(fk[:-6] + '.xyz-shape.h5', 'w')
    #R  = fo.create_dataset('xyz', (ncfgs, Ng, 3), dtype='f')
    sfs = zeros((ncfgs, 7)) # trQ, (trQ)^2, (trQ)^3, tr(\hat{Q}^2), det(\hat{Q}), asp, sf
    for c in range(0, ncfgs):
        # Collective coordinates X
        X = zeros((Ng, 3))
        for k in range(3):
            # X[1:,k] = sqrt(1.0/3/lam[1:])*random.randn(Ng-1)
            X[1:,k] = sqrt(1.0/lam[1:])*random.randn(Ng-1)
        # X -> 3D coordinates R
        xyz = zeros((Ng, 3))
        for k in range(3):
            xyz[:,k] = dot(Qs, X[:,k])
        xyz = xyz - mean(xyz, axis=0) # This step is necessary when diag(K) != 0
        #R[c,:,:] = xyz

        Q = zeros((3,3))
        for i in range(0, 3):
            for j in range(0, 3):
                Q[i,j] = mean(xyz[:,i]*xyz[:,j])
        trQ = trace(Q)
        QH = Q - trQ/3.0*eye(3)
        for k in range(0, 3):
            sfs[c,k] = trQ**(k+1)
        sfs[c,3] = trace(dot(QH,QH))
        sfs[c,4] = linalg.det(QH)
        sfs[c,5] = 1.5*sfs[c,3]/sfs[c,1]
        sfs[c,6] = 27.*sfs[c,4]/sfs[c,2]
    #fo.create_dataset('asp', data=sfs[:,5], dtype='f')
    #fo.create_dataset('spf', data=sfs[:,6], dtype='f')
    #fo.close()

    # Print
    asp = 1.5*mean(sfs[:,3])/mean(sfs[:,1])
    spf = 27.*mean(sfs[:,4])/mean(sfs[:,2])
    fw = open(fy, 'w')
    lt = "#Nsamples: %d " % (ncfgs)
    lt+= "asp: %11.5e " % (asp)
    lt+= "spf: %+12.5e " % (spf)
    lt+= "asp*: "
    for q in [5, 25, 50, 75, 95]:
        lt += "%11.5e " % (percentile(sfs[:,5], q))
    lt+= "spf*: "
    for q in [5, 25, 50, 75, 95]:
        lt += "%+12.5e " % (percentile(sfs[:,6], q))
    fw.write(lt+'\n')
    # P(\Delta)
    fw.write('#asp pdf\n')
    vs = sfs[:,5].copy()
    vd = 0.01
    v0 = 0
    v1 = 1
    vn = int((v1-v0)/vd)
    gx = array((vs - v0)/vd).astype(int)
    bx = bincount(gx)/len(gx)/vd
    pdf= zeros(vn)
    pdf[0:len(bx)] = bx
    for i in range(0, vn):
        lt = "%+12.5e %11.5e" % ((i+0.5)*vd+v0, pdf[i])
        fw.write(lt+'\n')
    fw.write('\n\n')
    # P(S)
    fw.write('#spf pdf\n')
    vs = sfs[:,6].copy()
    vd = 0.01
    v0 =-0.25
    v1 = 2.0
    vn = int((v1-v0)/vd)
    gx = array((vs - v0)/vd).astype(int)
    bx = bincount(gx)/len(gx)/vd
    pdf= zeros(vn)
    pdf[0:len(bx)] = bx
    for i in range(0, vn):
        lt = "%+12.5e %11.5e" % ((i+0.5)*vd+v0, pdf[i])
        fw.write(lt+'\n')
    fw.write('\n\n')
    fw.close()
    # P(\Delta) - P(S)
    vd = 0.01
    vnx = int((1.0 - 0.0)/vd)
    vny = int((2.0 +0.25)/vd)
    pdf = zeros((vnx, vny))
    for c in range(0, ncfgs):
        bx = int((sfs[c,5] - 0.0)/vd)
        by = int((sfs[c,6] +0.25)/vd)
        pdf[bx,by] += 1
    pdf = pdf/float(ncfgs*vd*vd)
    fw = open(fy+'2', 'w')
    fw.write("#N_Delta: %d N_S: %d\n"%(vnx, vny))
    for i in range(0, vnx):
        lt = ''
        for j in range(0, vny):
            lt += "%11.5e " % (pdf[i,j])
        fw.write(lt+'\n')
    fw.close()





