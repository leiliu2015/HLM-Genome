## HLM-Genome

HLM-***Genome*** is a *many*-polymer model to reconstruct full 3D genome structures using Hi-C data as the *sole* input. We represent a genome by using *M* polymer chains consisting of total *N* monomers, each representing a chromatin segment at a prescribed resolution. The algorithm, while similar to its previous versions ([HLM](https://github.com/leiliu2015/HLM) and [HLM-Nbody](https://github.com/leiliu2015/HLM-Nbody)), is updated by considering [PHi-C2](https://github.com/soyashinkai/PHi-C2). What's more, we modified the optimization code with [CuPy](https://github.com/cupy/cupy) to use ***GPU***, in order to construct genome models composed of thousands of monomers. 

### System Requirements

Python, [Gnuplot](gnuplot.sourceforge.net), and [PyMOL](https://pymol.org/2/) scripts to reproduce the results about HAP1 genome in our recent [work]() are deposited here, which we tested on ubuntu 18.04 LTS. [Anaconda](https://www.anaconda.com/distribution/) is recommended to manage the Python environment (Python *3.8*) and required packages, such as CuPy, H5py, NumPy, and Scipy. Last but not least, all genome models were optimized by using a NVIDIA **A100** card. 

### File Description
- hap1/
  - [hap1.sh](hap1/hap1.sh) (A BASH script to perform all the modeling and analyses in this directory)
  - [erez2015_rs1000000_kq0_nq1_om.txt](hap1/erez2015_rs1000000_kq0_nq1_om.txt) (A genome-wide Hi-C matrix of HAP1 cells measured by [Sanborn *et al.*](https://www.pnas.org/doi/full/10.1073/pnas.1518552112))
  - [erez2015_rs1000000_kq0_nq1_om.zi](hap1/erez2015_rs1000000_kq0_nq1_om.zi) (A TXT file including genomic position of each monomer)
  - [erez2015_rs1000000_kq0_nq1_om.gnu-gw.apx](hap1/erez2015_rs1000000_kq0_nq1_om.gnu-gw.apx) (A Gnuplot header file)
  - [erez2015_rs1000000_kq0_nq1_om.gnu-bc.apx](hap1/erez2015_rs1000000_kq0_nq1_om.gnu-bc.apx) (A Gnuplot header file)
  - erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup/ (A directory including optimization results)
  - exp/ (A directory including GPSeq and lamin B1 DamID data)
  - fig/ (A directory including Gnuplot and PyMOL scripts)
  - pys/ (A directory including Python scripts for optimization and analyses)
- pals/ (A directory including color palettes)
- [clearAll.sh](clearAll.sh) (A BASH script to delete all outputs)


### User Guide

The only required input of HLM is a Hi-C contact probability matrix, in this case the file named [erez2015_rs1000000_kq0_nq1_om.txt](hap1/erez2015_rs1000000_kq0_nq1_om.txt). If a A100 GPU card is available and CuPy has been installed in your Python environment, the optimal parameters of HLM for this genome can be determined as follows:
```
$ cd ./hap1
$ python pys/core_phic2_cupy.py erez2015_rs1000000_kq0_nq1_om.txt
```
As shown in the following table, this optimization step takes us 1 to 37 hours for different genomes. 

| Species | Cells | *N* | Optimization Time /hour |
| ------- | ----- | --- | ----------------------- |
| Human   | HAP1  | 2869| 37.3                    |
| Mouse   | Thymocyte (WT) | 2490 | 7.7          |
| Mouse   | Thymocyte (Lbr-/-) | 2476 | 7.5      |
| Fruit Fly | Kc167 | 1205 | 3.1 |
| Fruit Fly | Early embryos | 2293 | 1.7 |
| Mosquito | *A. Aaegyti* | 2387 | 14.3 |
| Wheat | *T. aestivum* (shoots) | 3499 | 24.9 |
| Wheat | *T. aestivum* (leaves) | 2776 | 1.1 |
| Yeast | *S. cerevisiae* | 2309 | 13.9 |
| Peanut | *A. hypogaea* (leaves) | 2502 | 7.7 |

Once it is finished, the optimized parameters ([N2869.K_fit](hap1/erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup/N2869.K_fit)), the model contact probability matrix ([N2869.P_fit](hap1/erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup/N2869.P_fit)), and a log file ([N2869.log](hap1/erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup/N2869.log)) will be stored in a directory as being listed in [erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup](hap1/erez2015_rs1000000_kq0_nq1_om_phic2_a1.0e-05_cupy_backup). 

Next, to analyze the model, all you need to do is to run the BASH script in your terminal which takes less than 10 minutes on our desktop with a Intel® Core™ i5-9500 CPU. 
```
$ bash ./hap1.sh
```

Last, to plot representative model structures, launch PyMOL from the [fig](hap1/fig) subdirectory, and type the following line in PyMOL command input area. Other two PyMOL scripts can be used in the same way.
```
PyMOL>@hg19_hap1.cfgs.Fig1A.pml
```

All the output files can be deleted by typing `$ bash ./clearAll.sh` at the repository root. For further questions and possible applications about HLM-Genome, please contact Lei Liu (leiliu2015@163.com)

