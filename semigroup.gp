print("\n\nType '?semigroup' for help.\n\n");
parigp_version=version();
semigroup_library=strprintf("./libsemigroup-%d-%d.so", parigp_version[1], parigp_version[2]);

/*semigroup.c*/
  addhelp(semigroup,"Basic (semi)group methods:\n\tLRword, semigroup_growth, semigroup_mats, semigroup_mgens.\n\nSemigroup orbits:\n\tsemigroup_missing, semigroup_missing_parabolic, semigroup_missinglist, semigroup_orbit.\n\nPsi methods:\n\tpsi_mats, psi_rep.\n\nLinear regression:\n\tOLS, OLS_nointercept, OLS_single, rsquared.\n\nPaper methods:\n\ttestkronaction, kronactioncorrect, psi_missingsquares, psi1_missingsquares, table1_line, table1_prediction, table1_iscorrect_psi, table1_iscorrect_psi1, table1_bigtest.");
  
  /*SECTION 1: BASIC (SEMI)GROUP METHODS*/
  install(LRword,"G",,semigroup_library);
  addhelp(LRword,"LRword(M): for M a SL(2, Z)^{>=0} hyperbolic matrix with positive entries, returns M as a word in L and R, where L=[1,1;0,1] and R=[1,0;1,1]. This is a Vecsmall with 0 representing L and 1 representing R. M must itself have entries that fit into the sice of a C long.");
  install(semigroup_growth,"GD100,L,D1000,L,DGp");
  addhelp(semigroup_growth,"semigroup_growth(mats, {binsize=100}, {Nbins=1000}, {start=[1, 1]}): estimates the growth rate of the orbits of the semigroup generated by mats, whose elements all have infinite order and nonnegative entries. The growth rate is of the form c*N^nu, and we estimate this by computing the orbit of the column vector in [binsize*Nbins, 2*binsize*Nbins], putting the sizes of the vectors into bins, and then running a linear regression. The return value is [c, nu, R^2], where R^2 is the R^2 value of the regression. The larger this value, the more confidence we have in (c, nu).");
  install(semigroup_mats,"GL");
  addhelp(semigroup_mats,"semigroup_mats(mats, N): assume the matrices in mats all have positive entries and infinite order. This returns the matrices in the semigroup they generate that have all entries at most N. If there are relations, the corresponding matrices will get counted multiple times.");
  install(semigroup_mgens,"G");
  addhelp(semigroup_mgens,"semigroup_mgens(mats): returns a minimal set of generators for the semigroup generated by mats. We assume that these matrices have all nonnegative entries, determinant +/-1, and remove the identity element if it occurs.");

  /*SECTION 2: SEMIGROUP ORBITS*/
  install(semigroup_missing,"GGGGLD1,L,");
  addhelp(semigroup_missing,"semigroup_missing(mats, B, start, congs, entry, {load=1}): prints to a file missing entries in a semigroup orbit, and al loads them if load=1. mats is the vector of matrices that generate the semigroup, necessarily all containing nonnegative entries. B is the bounds we wish to search, either a positive integer (for 1 to B), or a vector of two integers [B1, B2] for between B1 and B2. Note that we stop once any entry in the orbit is > B, which means that for some choices of matrices, some entries near B are declared missing but are not actually, as the first time they appear is with a larger entry. start is the starting vector for the orbit, which should be primitive and nonnegative. congs=[[r1, r2, ..., rk], modulus] is supplied as the congruence restictions: the entryth entry of a vector is the element we track, and the user must supply the possible congruence classes it can fall into, i.e. x mod(modulus) is one of r1, r2, ..., rk. If there are no congruence obstructions, you should supply [0, [1]].");
  install(semigroup_missing_parabolic,"GGGGLD1,L,");
  addhelp(semigroup_missing_parabolic,"semigroup_missing_parabolic(mats, B, start, congs, entry, {load=1}): same as semigroup_missing, except up to 80% slower but uses less memory. If one of the matrices is parabolic, the memory savings is extreme, and it is suggested to use this method when B gets large (10^9 or so). This method will work when none of the matrices are parabolic, though it is unlikely that the difference in memory will matter, and using semigroup_missing is suggested.");
  install(semigroup_missinglist,"iGGGL");
  addhelp(semigroup_missinglist,"semigroup_missinglist(mats, miss, start, entry): this returns 1 if the sorted vector of positive integers miss is entirely missed by the semigroup orbit of <mats>*start, where we look at the entryth entry of the orbit vectors. Calling this on a vector of square numbers gives evidence towards this being missing from an orbit. This is not optimized for parabolic generators, so should not be used for vectors miss which go past 10^9.");
  install(semigroup_orbit,"GLG");
  addhelp(semigroup_orbit,"semigroup_orbit(mats, B, start): for the semigroup given by mats (with all nonnegative entries and infinite order), this returns the orbit of the semigroup on the (column) vector start, where the entries are bounded by B. The output can have repeats, especially if there are relations among the matrices in mats. This is mainly useful for some initial testing, as the other semigroup_missing methods are significantly more efficient when searching for missing entries.");

  /*SECTION 3: PSI METHODS*/
  install(psi_mats,"L");
  addhelp(psi_mats,"psi_mats(N): Returns all the matrices in Psi with maximum entry N (Psi as defined in the paper: elements of Gamma_1(4) with nonnegative entries and Kronecker symbol of a row/column 1.");
  install(psi_rep,"iLLLD1,L,");
  addhelp(psi_rep,"psi_rep(x, y, n, {entry=1}): returns 1 if n is a numerator (entry=1) or denominator (entry=2) in the orbit Psi*[x,y]~. This is based on Lemma 9.2. All inputs must fit as C long's on your operating system, and (x, y) need to be nonnegative and coprime.");

  /*SECTION 4: LINEAR REGRESSION*/
  install(OLS,"GGD1,L,");
  addhelp(OLS,"OLS(X, y, {retrsqr=1}): performs ordinary least squares regression on the data, with the inputs being the n columns of the matrix X, and the outputs being the entries of the column vector y. We include a constant term, so the first row of X must be all 1's. If retrsqr=1, returns [params, R^2], and otherwise returns params, where params is the column vector of best fit parameters.");
  install(OLS_nointercept,"GGD1,L,");
  addhelp(OLS_nointercept,"OLS_nointercept(x, y, {retrsqr=1}): performs ordinary least squares regression assuming that y[i]=c*X[i] (where X and y are both (column) vectors), i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
  install(OLS_single,"GGD1,L,");
  addhelp(OLS_single,"OLS_single(x, y, {retrsqr=1}): performs linear regression for a single variable (essentially a macro for OLS with y=mx+b). Requires y to be a column vector and x is either a vector or column vector.");
  install(rsquared,"GGG");
  addhelp(rsquared,"rsquared(X, y, fit): returns the R^2 value for the proposed linear regression, where the input is X, output is y, and fit is the proposed parameters.");

/*paper.gp*/
  addhelp(testkronaction,"testkronaction(B, n, xymin, xymax): tests Proposition 3.2 by calling kronactioncorrect on n random matrices in SL(2, Z)^{>=0} with entries bounded by B. The values of x, y tried are all valid pairs with xymin<=x, y<=xymax. If the formula fails, we print the failing inputs, and raise an error. If no error occurs, all tests passed successfully.");
  addhelp(kronactioncorrect,"kronactioncorrect(M, xy): tests Proposition 3.2: M=[a, b;c, d] in SL(2, Z)^{>=0} a, b, c, d >= 0, xy=[x, y] with x, y>=0 coprime, gcd(x, d)=1, this proposition gives a formula for kron(ax+by/cx+dy). This function returns 1 if and only if the formula is correct for the given inputs. If the method returns anything other than 1, then the formula has failed.");
  addhelp(psi_missingsquares,"psi_missingsquares(xy, B, entry): tests if all the squares in the orbit of Psi*xy are missing up to B^2. Returns 1 if they are missing, 0 else.");
  addhelp(psi1_missingsquares,"psi1_missingsquares(xy, B, entry): tests if all the squares in the orbit of Psi_1*xy are missing up to B^2. Returns 1 if they are missing, 0 else.");
  addhelp(table1_line,"table1_line(xy): returns the line of Table 1 that our input corresponds to.");
  addhelp(table1_prediction,"table1_prediction(xy, entry): returns 1 if we know there are no squares in the orbit of Psi1*[x, y]~ due to a reciprocity obstruction, 0 if we expect there to be squares (based on Table 1), and -1 if there are no squares due to a congruence obstruction. entry=1 or 2 corresponding to numerator/denominator.");
  addhelp(table1_iscorrect_psi,"table1_iscorrect_psi(xy, B, entry): returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi orbit for [x, y]~. If squares are ruled out by congruence, we simply return 1.");
  addhelp(table1_iscorrect_psi1,"table1_iscorrect_psi1(xy, B, entry): returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi_1 orbit for [x, y]~. If squares are ruled out by congruence, we simply return 1.");
  addhelp(table1_bigtest,"table1_bigtest(xymax, B, {whichpsi = 1}): tests Table 1 for all 1<=x, y<=xymax and squares up to B, using the Psi orbit if whichpsi=0, and Psi_1 orbit else. Any failures of the table are output. If squares were predicted to not exist then this is a big issue! If squares were predicted to exist but did not, try increasing B for this example to see if you did not compute far enough. Returns the vector of counts of how many pairs of each line (1 to 22) were tested.");

read("paper.gp");
default(parisize, "4096M");/*Must come at the end*/