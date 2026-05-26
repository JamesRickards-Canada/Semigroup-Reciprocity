# Semigroups
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

Parallel:
* ```cfracsearch```: a separate C program used to compute positive integers that cannot be denominators of rational numbers of the form $[0;4a_1,4a_2,...,4a_n,a_{n+1},1,2]$. This is essentially ```contfrac_tail_missing```, except it is parallelized, and therefore is much more efficient. For example, call ```./cfracsearch 10 5000 8``` to compute the missing denominators between 10 and 5000 (inclusive), using 8 threads (with reference to Corollary 2.21 and Conjecture 2.22 of the [paper](https://doi.org/10.1215/00127094-2025-0017)/[Arxiv](https://arxiv.org/abs/2401.01860)).

## Testing
Load the file ```paper.gp```:
* ```runalltests```: run all the tests, which aims to computationally test many of the claims made in the paper, as well as provide supplemental computations.
* ```contfrac_tail_missing```: searches for positive integers that cannot be denominators of rational numbers of the form $[0;4a_1,4a_2,...,4a_n,a_{n+1},1,2]$.
