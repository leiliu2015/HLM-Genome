#!/bin/bash

Q=erez2015_rs1000000_kq0_nq1_om
N=2869
X=${Q}_phic2_a1.0e-05_cupy_backup

#### This step requires GPU
###python pys/core_phic2_cupy.py ${Q}.txt

#-- Fig. 1B --
gnuplot -persist -c fig/hg19_hap1.pij.gw.gnu ${Q}.gnu-gw.apx ${X}/N${N}.P_fit 1
gnuplot -persist -c fig/hg19_hap1.pij.bc.gnu ${Q}.gnu-bc.apx ${X}/N${N}.P_fit 5 5 1
gnuplot -persist -c fig/hg19_hap1.pij.bc.gnu ${Q}.gnu-bc.apx ${X}/N${N}.P_fit 19 5 1

#-- Fig. 1C --
python pys/hg19_hap1.pij.interChain.py ${Q}.gnu-bc.apx ${Q}.txt
python pys/hg19_hap1.pij.interChain.py ${Q}.gnu-bc.apx ${X}/N${N}.P_fit
gnuplot -persist -c fig/hg19_hap1.pij.interChain.gnu ${Q}.txt.interChain ${X}/N${N}.P_fit.interChain

#-- Fig. 1D/E --
python pys/hg19_hap1.pijCorr.py ${Q}.txt ${X}/N${N}.P_fit
gnuplot -persist -c fig/hg19_hap1.pijCorr1.gnu ${X}/N${N}.P_fit.pijCorr1
gnuplot -persist -c fig/hg19_hap1.pijCorr3.gnu ${X}/N${N}.P_fit.pijCorr3

#-- Fig. 2B --
python pys/core_K2xyz-shape.py ${X}/N${N}.K_fit 10000
gnuplot -persist -c fig/shapeDistr.gnu ${X}/N${N}.xyz-shape.asp2

#-- Fig. 2E/F --
python pys/core_K2cfg_r2.py ${X}/N${N}.K_fit ${Q}.zi 0
python pys/core_K2cfg_r2Chrx.py ${X}/N${N}.K_fit.r2
python pys/core_K2cfg_r2ChrxVsGPSeq.py ${X}/N${N}.K_fit.r2Chrx exp/GPSeq.txt
gnuplot -persist -c fig/r2VsChroSize.gnu ${X}/N${N}.K_fit.r2ChrxVsGPSeq
gnuplot -persist -c fig/r2VsGPSeq.gnu ${X}/N${N}.K_fit.r2ChrxVsGPSeq

#-- Fig. 2G --
python pys/core_K2cfg_r2VsLamin.py ${X}/N${N}.K_fit.r2 exp/LB1-DamID.txt
gnuplot -persist -c fig/r2VsLamin.gnu ${X}/N${N}.K_fit.r2VsLamin

#-- Figs. 1A/2C/2D can be found in ./fig (see README.txt).

