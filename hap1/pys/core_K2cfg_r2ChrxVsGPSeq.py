import os
import sys
from numpy import *

# concatenate xxx.r2Chrx, chromosome size, and GPSeq results
# ---
if not len(sys.argv) == 3:
    print('usage:: python core_K2cfg_r2ChrxVsGPSeq.py xxx.r2Chrx GPSeq.dat')
    sys.exit()
fx = str(sys.argv[1])
fq = str(sys.argv[2])

chrSizes = [249250621,243199373,198022430,191154276,180915260,\
            171115067,159138663,146364022,141213431,135534747,\
            135006516,133851895,115169878,107349540,102531392,\
            90354753,81195210,78077248,59128983,63025520,\
            48129895,51304566,155270560] # 1-22+X
L = len(chrSizes)

if not os.path.isfile(fx):
    print('Cannot find ' + fx)
    sys.exit()
else:
    r2s = []
    with open(fx) as f:
        for line in f:
            r2s.append( line[:-1] )
    if not L == len(r2s):
        print("Check the number of chromosomes in %s"%(fx))
        sys.exit()

#fq = 'GPSeq-Fig3d.txt'
if not os.path.isfile(fq):
    print('Cannot find ' + fq)
    sys.exit()
else:
    gps = []
    with open(fq) as f:
        for line in f:
            lt = line.strip().split()
            gps.append( array(list(map(float, lt[1:4]))) )
    gps = array(gps)
    if not L == len(gps):
        print("Check the number of chromosomes in %s"%(fq))
        sys.exit()

fw = open(fx+'VsGPSeq', 'w')
for c in range(0, L):
    lt = r2s[c]
    lt += " chrSize: %9d " % (chrSizes[c])
    lt += " GPSeq(25/50/75percentiles): "
    for p in range(0, 3):
        lt += "%11.5e " % (gps[c,p])
    fw.write(lt+'\n')
fw.close()


