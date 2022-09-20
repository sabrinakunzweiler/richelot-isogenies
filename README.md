# Efficient Computation of (2^n,2^n)-isogenies

This is the Sage/Magma code accompanying the paper

> S. Kunzweiler, Efficient Computation of (2^n,2^n)-Isogenies, Preprint, 2022

We provide an implementation of the algorithm described in that article and various examples. The implementation is available in Magma and Sagemath (see the respective folders). 
Furthermore, the repository contains code that can be used to verify the formulas developed in the manuscript.

Our algorithm can also be applied to speed-up the recent attack on SIDH by Castryck and Decru. A variant of the original implementation of this attack using the new algorithm for isogeny chains is contained in the folder SIDH-Attack-magma. For a (faster) variant in Sagemath go to https://github.com/sabrinakunzweiler/Castryck-Decru-SageMath/tree/new-algorithm. This builds on the improvements by various authors https://github.com/jack4818/Castryck-Decru-SageMath, but we replaced the isogeny chain computations with our algorithm. 
