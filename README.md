# DDmap
A Matlab Package for Solving the Double Digest Problems

It is the basic supplementary data of our paper entiled "DDmap: A Matlab package for double digest problem using multiple genetic operators" (preprint). With the purpose to keep the consistency of the results with that of in this paper, any modification from others should be at first send to wanglc2012@126.com, instead of editing this public project directly.

Note: Our all simulations in this project are conducted on a X1 Carbon laptop with Windows(TM) 8, Intel(R) Core(TM) i5-4300U CPU@1.90GHz/2.49GHz and 8GB RAM. Softwares: MATLAB R2014a, MAPLE 2018

The package DDmap consists of
13 Matlab algorithms:
–permGA.m, the MATLAB genetic algorithm for solving the DDP problem. This is the main algorithm, and its flowchart is depicted in Fig. 6. Note that this file also contains the definitions of five genetic operators — RWS, PCC, P4X, FLP, CSH and related MATLAB functions for calculating the fitness values.

–referIndexSort.m, the MATLAB algorithm for implementing the so-called referenced sorting (based on index).

–opPermCross.m, the MATLAB algorithm for implementing the RSC genetic operator.

–getInstance.m, the auxiliary MATLAB algorithm for outputting test DDP instances in [Sur-Kolay S. et al., 2005].

–randDDPinstance.m, the auxiliary MATLAB algorithm for producing a valid DDP instance according the given parameters.

–strABC.m, the auxiliary MATLAB algorithm for producing Maple commands for reading data before calling the Maple algorithm DDdraw.mws.

–simu1004.m, simu1004plots.m, simu1007.m, and simu1008.m, the auxiliary MATLAB algorithms for organizing simulations and producing the related figures.

–Trans.m, the Matlab algorithm for implementing Scaling-rounding-adjusting approach for Cases of .

–Plot1.m, the auxiliary Matlab algorithms for comparison of DDmap and algorithm in [3] and [13].

–Plot2.m, the auxiliary Matlab algorithms for comparison of DDmap and other algorithms under the condition of .

1 Maple algorithm, DDdraw.mws, is used for drawing the DDP solution in nested pie charts, with inputs A, B and C that are assigned by using Maple commands produced by strABC.m.

43 Data files: 42 of them are named as , where and , and the last is named as . These data files are in fact the running records of our simulations towards the 7 valid DDP test instances given in [13]. In these running records, many exact DDP solutions are provided.

References

[3].Ganjtabesh M., Ahrabian H., Nowari-Dalini A. & Moghadam Z.R.K. Genetic algorithm solution for double digest problem. Bioinformation, 2012, 8(10):453-456,.

[13].Sur-Kolay S , Banerjee S , Mukhopadhyaya S , et al. Genetic Algorithm for Double Digest Problem. International Conference on Pattern Recognition & Machine Intelligence. Springer Berlin Heidelberg, 2005.
