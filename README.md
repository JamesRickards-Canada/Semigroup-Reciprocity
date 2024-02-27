# Semigroup-Reciprocity

This repository is supplementary to the [paper](https://arxiv.org/abs/2401.01860)  "Reciprocity obstructions in semigroup orbits in $\text{SL}(2, \mathbb{Z})$", by James Rickards and Katherine E. Stange.

The goal of the package is to efficiently compute data about certain orbits of semigroups. Some methods have been optimized to be able to compute extremely large orbits. The files "paper.gp" and "cfracsearch.c" are aimed at verifying and recreating various claims in the paper.

See the final sections for information on how to install the package. You need to be running PARI/GP on a Linux based system. If you are running Windows, then you must use Windows Subsystem for Linux (WSL). See this [tutorial](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) for getting started with PARI/GP on WSL.

## Semigroups
For the purposes of this package, a semigroup is finitely generated, specified by a vector of $n\times n$ matrices of infinite order and non-zero determinant with nonnegative integer entries. This will typically be a subsemigroup of $\text{SL}(n, \mathbb{Z})^{\geq 0}$.

Let $\Gamma$ be such a semigroup, and $v\in\mathbb{Z}^n$ be a primitive vector with nonnegative entries. This gives an orbit $\Gamma\cdot v$. We have methods to compute the elements of $\Gamma\cdot v$ whose entries are bounded by $B$, the elements of $\Gamma$ that are bounded by $B$, as well as estimate the growth rate of $\Gamma\cdot v$. There is also a method to reduce a semigroup to a minimal set of generators.

Furthermore, there are methods for exploring the local-global conjecture for large values of $B$ ($10^{10}$ and beyond!). In particular, if we want to study the $k^{\text{th}}$ entry of the orbit $\Gamma\cdot v$, then we can compute the "missing curvatures" up to a bound $B$.

## Main methods
General semigroups:
* ```semigroup_mats```: computes all elements of the semigroup with bounded entries.
* ```semigroup_mgens```: computes a minimal generating set for the semigroup. Only works if the matrices all have determinant $\pm 1$.
* ```semigroup_growth```: estimates the growth rate of the semigroup orbit.

Semigroup orbits:
* ```semigroup_orbit```: computes a semigroup orbit up to a bound. This is a good method for initial testing and playing around (or if you want to do other operations with the orbit besides taking a single entry), but is not optimized for large-scale computations.
* ```semigroup_missing```: computes the missing integers in a semigroup orbit. The user must supply the relevant congruence obstructions, so this should not be used if they are not known (or, worse case scenario, no congruence obstructions are passed).
* ```semigroup_missing_parabolic```: same as ```semigroup_missing```, except optimized for memory usage, which is useful when one of the generating matrices being parabolic. This is roughly 80% slower, but should be used if you are computing beyond $10^9$, as ```semigroup_missing``` will take far too much memory.
* ```semigroup_missinglist```: similar to ```semigroup_missing```, except we supply a vector of integers that we want to check if it is completely missing. This is not optimized for the semigroup having parabolic elements.

Testing:
* ```runalltests```: run all the tests in paper.gp, which aims to computationally test many of the claims made in the paper, as well as provide supplemental computations.
* ```contfrac_tail_missing```: searches for positive integers that cannot be denominators of rational numbers of the form $[0;4a_1,4a_2,...,4a_n,a_{n+1},1,2]$.

Parallel:
* ```cfracsearch```: a separate C program used to compute positive integers that cannot be denominators of rational numbers of the form $[0;4a_1,4a_2,...,4a_n,a_{n+1},1,2]$. This is essentially ```contfrac_tail_missing```, except it is parallelized, and therefore is much more efficient.

## Installation
System requirements:
* Linux/Mac: none that I am aware of;
* Windows: you need to be running PARI/GP through WSL. See this [tutorial](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) for how to set this up.

This package has been tested on PARI/GP version 2.15.2 and 2.16.2, and should work on 2.15.0 and later.

### Main PARI/GP package:
You need to know where the version of PARI/GP you want to use is installed, in particular, the file ```pari.cfg```. The default location is ```/usr/local/lib/pari/pari.cfg```.
* If this is the correct location, call ```make``` to build the project.
* Otherwise, call ```make setup``` to search for the location of the file. By default the program searches in ```/usr```, but there is a chance it is not installed there (this sometimes happens on a server). If this is the case, you can supply an alternate location.
* If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to ```paricfg_loc.txt```. Once this step completes, a call to ```make``` will compile the project! Modifying the program (e.g. via ```git pull```) won't require redoing this setup, unless the version of PARI/GP you use changes.

Call ```gp semigroup``` to start gp and load the methods. ```?semigroup``` accesses the help.

Call ```make clean``` to clean up the object files (.o) created.

### cfracsearch

1. After cloning the repository, call ```make cfracsearch``` to compile it.
   
2. Call ```./cfracsearch 10 5000 8``` to compute the missing denominators between 10 and 5000 (inclusive), using 8 threads (with reference to Corollary 2.21 and Conjecture 2.22 of the [paper](https://arxiv.org/abs/2401.01860)).

## Troubleshooting installation

* **No library found**: the files were not found in the search location. Try asking to search in "/", which searches everywhere. This will be slow, but is guaranteed to find the correct files, if they exist.
* **Wrong version**: Maybe you found the libraries, but there were warnings with ```make```. The likely cause is the version of PARI/GP you found was too old. If there are multiple copies of PARI/GP on your computer, then perhaps you chose the wrong one! Check the shared object file created: it will be called "libsemigroup-X-Y-Z.so", where "X.Y.Z" is the version of PARI/GP in the libraries. If this does not match the version you are using, then you found the wrong one!
* **Miscellaneous**: when you compile with ```make```, object files (.o) are created. However, if the underlying code did not change, then nothing will happen. If you change the version of PARI/GP you are working with (or perhaps, trying different installations), then you should call ```make clean``` in between to clear out these object files. Otherwise, recompiling will do nothing!
* If you are still having issues with installation, please get in touch, and I will try to help sort it out!
 
