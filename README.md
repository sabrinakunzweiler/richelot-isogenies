# Efficient Computation of (2^n,2^n)-isogenies

This is the Sage/Magma code accompanying the paper

> S. Kunzweiler, Efficient Computation of (2^n,2^n)-Isogenies, Preprint, 2022

It consists of three parts currently. 
- The folder "paper-2022" contains the algorithms and examples from the mentioned paper. The implementation is available in Magma and Sagemath (see the respective folders). 
Furthermore, we provide code that can be used to verify the formulas developed in the manuscript.

- The content of the folder "product-isogenies" was developed after our eprint appeared, and it complements that work by also providing an algorithm for computing product isogenies. Here, we provide two options: 
For the first method, we use the formulae from our paper (+ additional gluing and splitting maps)
For the second method, we use the set-up of the paper, but compute on the associated Kummer model of the Jacobian variety. This has the advantage of obtaining cleaner (and more efficient) formulae for isogeny evaluations. On the other hand, arithmetic on the Kummer surface is less efficient.

- The algorithm can also be applied to speed-up the attack on SIDH by Castryck and Decru. A variant of the original implementation of this attack using the new algorithm for isogeny chains is contained in the folder SIDH-Attack-magma. (Note that there are newer versions of the attack which are faster.)
