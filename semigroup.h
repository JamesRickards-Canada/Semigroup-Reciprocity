/*semigroup.c*/

/*SECTION 1: BASIC (SEMI)GROUP METHODS*/
GEN LRword(GEN M);
GEN semigroup_growth(GEN mats, long binsize, long Nbins, GEN start, long prec);
GEN semigroup_mats(GEN mats, long N);
GEN semigroup_mgens(GEN mats);

/*SECTION 2: SEMIGROUP ORBITS*/
GEN semigroup_missing(GEN mats, GEN B, GEN start, GEN congs, long entry, long load);
GEN semigroup_missing_parabolic(GEN mats, GEN B, GEN start, GEN congs, long entry, long load);
int semigroup_missinglist(GEN mats, GEN miss, GEN start, long entry);
GEN semigroup_orbit(GEN mats, long B, GEN start);

/*SECTION 3: PSI METHODS*/
GEN psi_mats(long N);
int psi_rep(long x, long y, long n, long entry);

/*SECTION 4: LINEAR REGRESSION*/
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN x, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);
