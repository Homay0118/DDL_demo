Matlab demo code accompanying the PR paper
[J. Yin, H. Zhu, D. Yuan, T. Xue, Sparse Representation over Discriminative Dictionary for Stereo Matching, PR]

Contributed by Hongmei Zhu (z_hongmei@buaa.edu.cn)


Important notes:

(1) The code is provided for academic use only. Use of the code in any commercial or industrial related activities is prohibited.

(2) If you use our codes we request that you cite the paper [J. Yin, H. Zhu, D. Yuan, T. Xue, Sparse Representation over Discriminative Dictionary for Stereo Matching, PR]

(3) Platform: Matlab 2010+ with C Compiler

Usage:

* Download stereo pairs from http://vision.middlebury.edu/stereo/data/ for training and test. Remember to change the data paths if necessary.

* Dictionary Learning is in the folder "training"
 - folders "all15","first8","last7" are used to save the corresponding training samples and learned dictionaries.
 - run CollectPatches.m to generate the training samples
 - run DictTraining.m to learn discriminative dictionary, pretrained dictionaries are already saved in folders "all15","first8","last7", respectively
 - mexLassoWeighted.m/mexLassoWeighted.mexw64, mexTrainDL.m/mexTrainDL.mexw64 are from SPAMS-Julien Toolbox (http://spams-devel.gforge.inria.fr/)

* Evaluation in the folder "Test"
 - run Middleburytest_2014.m to evaluate on Middlebury benchmark version 3
 - run Middleburytest_2006.m to evaluate on 30 stereo pairs from Middlebury benchmark version 2
 - CostComputingmex.c corresponds to Eq.(12).
 - SGMmex.c is the mex version to re-implement four directions cost aggregation like reference "H. Hirschmuller, Stereo processing by semiglobal matching and mutual information, PAMI".
 - PostProcmex.c is the mex version to re-implement post processing in reference "C. Rhemann, and A. Hosni, and M. Bleyer, and C. Rother, and M. Gelautz, Fast cost-volume filtering for visual correspondence and beyond, CVPR2011".
 - mexLasso.m/mexLasso.mexw64 is from SPAMS-Julien Toolbox (http://spams-devel.gforge.inria.fr/)
 - If you trained a new dictionary, you should change the dictionary name into "*_new.mat" in Middleburytest_2014.m and Middleburytest_2006.m

Enjoy it!

2017.05.18