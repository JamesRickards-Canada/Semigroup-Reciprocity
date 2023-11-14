/*thin.c*/

/*SECTION 1: BASIC (SEMI)GROUP METHODS*/
GEN LRword(GEN M);
GEN semigroupgrowth(GEN mats, long binsize, long Nbins, GEN start, long prec);
GEN semigroupmats(GEN mats, long N);

/*SECTION 2: MISSING NUMBERS IN ORBITS*/
GEN semigroup_missing(GEN mats, GEN B, GEN start, GEN congs, long entry);
GEN semigroup_missing_parabolic(GEN mats, GEN B, GEN start, GEN congs, long entry);

/*SECTION 3: LINEAR REGRESSION*/
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN x, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);
