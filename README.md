# Multiband-RF-pulse-Design
This is a MATLAB toolbox for RF pulse design used in Magnetic Resonance Imaging (MRI). It features design of 1D selective RF pulse with multiband magnitude profile, arbitrary phase profile and generalized flip angle. The implementation is built on Shinnar-Le Roux (SLR) algorithm and convex optimization. Spectral sparsity was exploited to further improve pulse performance, such as minimizing pulse duration, transition width or total energy with flexible trade-off. 

Example designs of RF pulses for bSSFP C-13 MRI at 14T and a dualband saturation pulse for H-1 MRS at 3T were developed.

It is hosted on this open-source, collaborative platform in order to encourage everyone in the MRI research community to contribute tools that will help our field rapidly progress.

Some public MRI toolboxes are used, including 
    A RF pulse design toolbox developed by John Pauly using the Shinnar-Le Roux algorithm. Reference: Pauly J, Le Roux P, Nishimra D, Macovski A. Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm. IEEE Tr Medical Imaging 1991; 10(1):53-65. http://rsl.stanford.edu/research/software.html
    A multiband spectral-spatial RF pulse design package developed by Adam B. Kerr and Peder E.Z. Larson. Reference: Kerr, A.B., Larson, P.E., Lustig, M., Cunningham, C.H., Chen, A.P., Vigneron, D.B., and Pauly, J.M. Multiband spectral-spatial design for high-field and hyperpolarized C-13 applications. In Proceedings of the 16th Annual Meeting of ISMRM (Toronto, 2008), p.226. http://rsl.stanford.edu/research/software.html
    A Bloch simulator written as a MEX function developed by Brian Hargreaves. http://mrsrl.stanford.edu/~brian/blochsim/
For readers' convenience, these public packages are downloaded and included here. The author acknowledges their work and making it public.

CVX is used to solve convex optimization problem in this work, which is a Matlab-based modeling system for specifying and solving convex programs. Reference: Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013. To use this toolbox, please install cvx first. More details can be found at http://cvxr.com/cvx/.
