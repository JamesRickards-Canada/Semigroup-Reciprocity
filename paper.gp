/*Methods related to the paper, divided into three sections:

1. TESTING: tests the correctness of certain claims by computationally checking them in a number of cases.

2. SUPPLEMENTARY COMPUTATIONS: Recreate computations claimed in the paper.

3. SUPPORTING METHODS: self explanatory
*/


/*SECTION 1: TESTING*/

/*Run tests for all of the testing methods we wrote.*/
runalltests() = {
  printf("Let's computationally test various claims made in the paper.\n\n\n");
  printf("To begin, let's test that even continued fractions correspond to LR words.\n");
  test_evencontfrac(20000, 3000);
  printf("\n\nNext, let's check Lemma 3.1 that the Kronecker symbol of any row or column of an element of Gamma_1(4)^{>=0} is constant.\n");
  test_gamma14geq0_kronequal(20000, 3000);
  printf("\n\nNext, let's check Proposition 3.2 about the formula for kronecker(ax+by, cx+dy) in terms of kronecker(x, y) when [a,b;c,d] is a matrix in SL(2, Z)^{>=0}.\n");
  test_kronaction_many(2000, 3000, 1, 70);
  printf("\n\nNext, let's check that Psi is a semigroup.\n");
  test_psisemigroup(20000, 3000);
  printf("\n\nNext, let's check that Table 1 gives a correct and complete list of reciprocity obstructions for Psi, i.e. Theorem 2.6.\n");
  test_table1_psi_many();
  printf("\n\nNext, let's check that Table 1 gives a correct and complete list of reciprocity obstructions for Psi_1, i.e. Conjecture 2.13.\n");
  test_table1_psi1_many();
  printf("\n\nOne orbit did not work! Let's increase B to 400 and see if it works now.\n");
  if (test_table1_psi1([143, 190], 400, 1), printf("Square found!\n"), printf("No square found. Is B too small or is the conjecture false?"));
  printf("\n\nNext, let's check that the matrices M_k (from Lemma 7.7) must appear in any set of generators for Psi.\n");
  test_psioogens(60);
  printf("\n\nNext, let's check that the numerators in the orbit of Psi*[2, 3]~ are as claimed in Theorem 2.7.\n");
  test_psi_23orbit();
  printf("\n\nNext, let's estimate the Hausdorff dimension of Psi_1.\n");
  test_psi1_hdim();
  printf("\n\nFinally, let's check the numerators in the orbit of Psi_1*[2, 3]~ up to 10^7 to get Conjecture 2.14. This may take half a minute.\n");
  test_psi1_23orbit();
  printf("\n\nAll standard tests complete!! No critical errors found.\n");
}

/*Tests that the even continued fraction does correspond to orbits of [1,0]~ as described at the start of section 2. We pick random positive rational numbers x/y with 1<=x,y<=B, and do this test n times.*/
test_evencontfrac(n, B) = {
  my(x, c, M);
  printf("Testing %d random positive rational numbers with num/denom bounded by %d\n", n ,B);
  for (i = 1, n,
    x = (random(B) + 1)/(random(B) + 1);
    c = contfraceven(x);
    M = contfractoword(c);
    if (M[1, 1] != numerator(x) || M[2,1] != denominator(x), printf("Test failed with x=%Ps\n", x);error("Failed test."));
  );
  printf("All tests passed.\n");
}

/*Generates n elements of Gamma_1(4)^{>=0} whose coefficients are at most B, and checks Lemma 3.1: the Kronecker symbols of their rows and columns are all equal.*/
test_gamma14geq0_kronequal(n, B) = {
  my(M, k);
  printf("Testing %d random elements of Gamma_1(4)^{>=0} whose coefficients are at most %d:\n", n, B);
  for (i = 1, n,
    M = gamma14geq0_random(B);
    k = kronecker(M[1, 1], M[1, 2]);
    if (k != kronecker(M[1, 1], M[2, 1]), printf("Test failed with M = %Ps\n", M);error("Test failed."));
    if (k != kronecker(M[1, 2], M[2, 2]), printf("Test failed with M = %Ps\n", M);error("Test failed."));
    if (k != kronecker(M[2, 1], M[2, 2]), printf("Test failed with M = %Ps\n", M);error("Test failed."));
  );
  printf("All tests passed.\n");
}

/*Tests Proposition 3.2: M=[a, b;c, d] in SL(2, Z)^{>=0} a, b, c, d >= 0, x, y>=0 coprime, gcd(x, d)=1, this proposition gives a formula for kron(ax+by/cx+dy). This function returns 1 if and only if the formula is correct for the given inputs.*/
test_kronaction(M, xy) = {
  my(x, y, a, b, c, d, A, B, cxpdy, C, D, alpha, mu);
  x = xy[1]; y = xy[2];
  a = M[1, 1]; b = M[1, 2];
  c = M[2, 1]; d = M[2, 2];
  if (gcd(x, d) > 1 || gcd(x, y) > 1 || a < 0 || b < 0 || c < 0 || d < 0, return(1));/*Ignore this case*/
  A = ((oddpart(x) - 1) >> 1) % 2;
  B = ((oddpart(d) - 1) >> 1) % 2;
  cxpdy = c * x + d * y;
  C = ((oddpart(cxpdy) - 1) >> 1) % 2;
  D = ((oddpart(y) - 1) >> 1) % 2;
  alpha = (A * B + A * C + B * C + A * D) % 2;/*Formula for alpha*/
  if (valuation(x, 2) == 1 || valuation(d, 2) == 1,/*mu=mu_1 so far*/
    mu = kronecker(c*x*d*y + 1, 2);
  ,
    mu = 1;
  );
  if (valuation(cxpdy, 2) == 1, mu = mu * kronecker(b * x * cxpdy + 1, 2));/*Now mu=mu_1*mu_2*/
  return(kronecker(a * x + b * y, cxpdy) * (-1)^(alpha) * mu * kronecker(c, d) * kronecker(x, y));
}

/*Tests Proposition 3.2 by calling test_kronaction on n random matrices in SL(2, Z)^{>=0} with entries bounded by B. The values of x, y tried are all valid pairs with xymin<=x, y<=xymax*/
test_kronaction_many(n, B, xymin, xymax) = {
  my(M);
  printf("Testing %d random matrices with coefficients at most %d and all valid pairs [x, y]~ where %d<=x, y<=%d:\n", n, B, xymin, xymax);
  for (i = 1, n,
    M = sl2zgeq0_random(B);
    for (x = xymin, xymax,
      if (gcd(x, M[2, 2]) == 1,
        for (y = xymin, xymax,
          if (gcd(x, y) == 1,
            if (test_kronaction(M, [x, y]) != 1, printf("WRONG FORMULA: %Ps %d %d\n", M, x, y);error("FORMULA FAILED"));
          );
        );
      );
    );
    if (i % 500 == 0, printf("%d matrices tried\n", i));
  );
  printf("All tests passed.\n");
}

/*Tests that the matrices M_k=[12k+1,4k;12k+4,4k+1] are all in a minimal generating set up to k=n.*/
test_psioogens(n) = {
  my(v, mv, found, M);
  printf("Testing that the matrices M_k=[12k+1,4k;12k+4,4k+1] for 1<=k<=%d cannot be expressed as non-trivial words in Psi.\n", n);
  v = psi_mats(12 * n + 4);
  printf("%d elements of Psi found with coefficients bounded by %d\n", #v, 12 * n + 4);
  mv = semigroup_mgens(v);
  printf("Minimal generating set has %d elements\n", #mv);
  for (i = 1, n,
    found = 0;
    M = [12 * i + 1, 4 * i; 12 * i + 4, 4 * i + 1];
    for (j = 1, #mv, if (mv[j] == M, found = 1; break; ));
    if (!found, printf("MATRIX NOT FOUND: %Ps\n", M);error("Not a minimal generator."));
  );
  printf("All tests passed.\n");
}

/*We generate pairs of elements of Psi with coefficients of size at most B, and check that their product still lies in Psi.*/
test_psisemigroup(n, B) = {
  my(M, N, P);
  printf("Testing %d random pairs of elements of Psi:\n", n);
  for (i = 1, n,
    M = psi_random(B);
    N = psi_random(B);
    P = M * N;
    if (kronecker(P[1, 1], P[1, 2]) != 1, printf("Test failed with %Ps and %Ps\n", M, N); error("Test failed."));
  );
  printf("All tests passed.\n");
}

/*Returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi orbit for xy=[x, y]~.*/
test_table1_psi(xy, B, entry) = {
  my(pred);
  pred = psi_isreciprocity(xy, entry);
  if (pred == -1, return(1));/*Congruence obstruction, nothing to do.*/
  actual = psi_missingsquares(xy, B, entry);
  if (pred == 1 && actual == 0, printf("ERROR: squares found when there were not supposed to be with inputs xy=%Ps, B=%d, entry=%d\n", xy, B, entry);error("Squares found when not allowed."));
  if (pred == 0 && actual == 1, return(0));/*No obstruction but we didn't find any squares. We don't raise an error since perhaps B is not large enough.*/
  return(1);
}

/*Tests Table 1 for all 1<=x, y<=xymax and squares up to B^2. If B is too small this will of course fail. Returns the vector of counts of how many pairs of each line (1 to 9) were tested.*/
test_table1_psi_many(xymax = 200, B = 300) = {
  my(x, y, v);
  printf("Testing Table 1 on the orbits of Psi for 1<=x, y<=%d and squares up to %d^2:\n", xymax, B);
  v = vector(9);
  for (x = 1, xymax,
    for (y = 1, xymax,
      if (gcd(x, y) == 1,
        v[table1_line([x, y])]++;
        if (test_table1_psi([x, y], B, 1) != 1, printf("Squares expected but not found: (x, y, B, entry) = %d %d %d 1\n", x, y, B));
        if (test_table1_psi([x, y], B, 2) != 1, printf("Squares expected but not found: (x, y, B, entry) = %d %d %d 2\n", x, y, B));
      );
    );
  );
  printf("Done tests!\n");
  return(v);
}


/*SECTION 2: SUPPLEMENTARY COMPUTATIONS*/

/*Prints the missing non-square numerators in the orbit, defaulting to the value we need to check it up to.*/
test_psi_23orbit(n = 109120000) = {
  printf("Finding the non-square numerators of Psi[2, 3]~ between 1 and %d, printing those that are missing:\n", n);
  for (i = 1, n,
    if (!issquare(i) && !psi_rep(2, 3, i, 1), printf("%d\n", i));/*Missing non-square*/
  );
  printf("All numbers searched!\n");
}

/*Prints the missing non-square numerators in the orbit, defaulting to 10^7.*/
test_psi1_23orbit(n = 10000000) = {
  my(M, v);
  printf("Finding the non-square numerators of Psi_1[2, 3]~ between 1 and %d, printing those that are missing:\n", n);
  M = [[1, 1;0, 1], [1, 0;4, 1]];
  v = semigroup_missing(M, n, [2, 3]~, [[0], 1], 1, 1);
  for (i = 1, #v, if (!issquare(v[i]), printf("%d ", v[i])));
  printf("\nAll numbers searched!\n");
}

/*Estimates the Hausdorff dimension of Psi_1.*/
test_psi1_hdim() = {
  my(M, binsize, Nbins, v);
  binsize = 100;
  Nbins = 500;
  printf("Estimating the Hausdorff dimension of Psi_1 by running a least squares regression on the orbit of [1, 1]~ up to %d to estimate the critical exponent:\n", binsize * Nbins);
  M = [[1, 1;0, 1], [1, 0;4, 1]];
  v = semigroup_growth(M, binsize, Nbins, [1, 1]);
  printf("Hausdorff dimension: %.5f, R^2: %.5f\n", v[2]/2, v[3]);
}

/*Returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi_1 orbit for xy=[x, y]~.*/
test_table1_psi1(xy, B, entry) = {
  my(pred);
  pred = psi_isreciprocity(xy, entry);
  if (pred == -1, return(1));/*Congruence obstruction, nothing to do.*/
  actual = psi1_missingsquares(xy, B, entry);
  if (pred == 1 && actual == 0, printf("ERROR: squares found when there were not supposed to be with inputs xy=%Ps, B=%d, entry=%d\n", xy, B, entry);error("Squares found when not allowed."));
  if (pred == 0 && actual == 1, return(0));/*No obstruction but we didn't find any squares. We don't raise an error since perhaps B is not large enough (or, maybe Table 1 is wrong for Psi_1!).*/
  return(1);
}

/*Tests Table 1 for all 1<=x, y<=xymax and squares up to B^2. If B is too small this will of course fail. Returns the vector of counts of how many pairs of each line (1 to 9) were tested.*/
test_table1_psi1_many(xymax = 200, B = 300) = {
  my(x, y, v);
  printf("Testing Table 1 on the orbits of Psi_1 for 1<=x, y<=%d and squares up to %d^2:\n", xymax, B);
  v = vector(9);
  for (x = 1, xymax,
    for (y = 1, xymax,
      if (gcd(x, y) == 1,
        v[table1_line([x, y])]++;
        if (test_table1_psi1([x, y], B, 1) != 1, printf("Squares expected but not found: (x, y, B, entry) = %d %d %d 1\n", x, y, B));
        if (test_table1_psi1([x, y], B, 2) != 1, printf("Squares expected but not found: (x, y, B, entry) = %d %d %d 2\n", x, y, B));
      );
    );
  );
  printf("Done tests!\n");
  return(v);
}


/*SECTION 3: SUPPORTING METHODS*/

/*Returns a random element of Gamma_1(4)^{>=0}, where the size of the coefficients are bounded by n (not necessarily uniformly, but should be reasonably close).*/
gamma14geq0_random(n) = {
  my(bd1, bd2, a, b, c, d, bd, s);
  bd1 = floor((n - 1) >> 2);
  bd2 = floor(n >> 2);
  while (1 == 1,
    a = 4 * random(bd1) + 1;
    c = 4 * random(bd2);
    while (gcd(a, c) != 1,/*Pick random (a, c)'s until the gcd is 1*/
      a = 4 * random(bd1) + 1;
      c = 4 * random(bd2);
    );
    if (c == 0, return([1, random(n + 1);0, 1]));
    d = lift(1 / Mod(a, c));
    b = (a * d - 1) / c;/*Solve for minimal nonnegative b, d*/
    bd = floor((n - b) / a);
    bd = min(bd, floor((n - d) / c));/*See if we can shift them higher; pick a random shift. If either b,d>n, then we just redo it all.*/
    if (bd >= 0,/*All coeffs can be made to be at most n*/
      s = random(bd + 1);
      b = b + a * s;
      d = d + c * s;
      return([a, b;c, d]);
    );
  );
}

/*Returns a random element of Psi where the size of the coefficients are bounded by n (not necessarily uniformly, but should be reasonably close). Note that we use the convention of kronecker(a, c)=1, but by Lemma 3.1 (which we also test), this is equivalent.*/
psi_random(n) = {
  my(bd1, bd2, a, b, c, d, bd, s);
  bd1 = floor((n - 1) >> 2);
  bd2 = floor(n >> 2);
  while (1 == 1,
    a = 4 * random(bd1) + 1;
    c = 4 * random(bd2);
    while (kronecker(a, c) != 1,/*Pick random (a, c)'s until the Kronecker symbol is 1*/
      a = 4 * random(bd1) + 1;
      c = 4 * random(bd2);
    );
    if (c == 0, return([1, random(n + 1);0, 1]));
    d = lift(1 / Mod(a, c));
    b = (a * d - 1) / c;/*Solve for minimal nonnegative b, d*/
    bd = floor((n - b) / a);
    bd = min(bd, floor((n - d) / c));/*See if we can shift them higher; pick a random shift. If either b,d>n, then we just redo it all.*/
    if (bd >= 0,/*All coeffs can be made to be at most n*/
      s = random(bd + 1);
      b = b + a * s;
      d = d + c * s;
      return([a, b;c, d]);
    );
  );
}

/*Returns a random element of SL(2, Z)^{>=0}, where the size of the coefficients are bounded by n (not necessarily uniformly, but should be reasonably close).*/
sl2zgeq0_random(n) = {
  my(a, b, c, d, bd, s);
  while (1 == 1,
    a = random(n) + 1;
    c = random(n + 1);
    while (gcd(a, c) != 1,/*Pick random (a, c)'s until the gcd is 1*/
      a = random(n) + 1;
      c = random(n + 1);
    );
    if (c == 0, return([1, random(n + 1);0, 1]));
    d = lift(1 / Mod(a, c));
    b = (a * d - 1) / c;/*Solve for minimal nonnegative b, d*/
    bd = floor((n - b) / a);
    bd = min(bd, floor((n - d) / c));/*See if we can shift them higher; pick a random shift. If either b,d>n initially, then we just redo it all.*/
    if (bd >= 0,/*All coeffs can be made to be at most n*/
      s = random(bd + 1);
      b = b + a * s;
      d = d + c * s;
      return([a, b;c, d]);
    );
  );
}

/*Returns the odd part of an integer.*/
oddpart(w) = {
  my(v);
  v = valuation(w, 2);
  return(w >> v);
}

/*Returns 1 if we know there are no squares in the orbit of Psi1*xy due to a reciprocity obstruction, 0 if we expect there to be squares (based on Table 1), and -1 if there are no squares due to a congruence obstruction. */
psi_isreciprocity(xy, entry) = {
  my(line);
  line = table1_line(xy);
  if (line == 1, return(1));
  if (line == 2, return(0));
  if (line == 3, if (entry == 1, return(-1), return(1)));
  if (line == 4, if (entry == 1, return(-1), return(0)));
  if (line == 5, return(1));
  if (line == 6, return(0));
  if (line == 7, if (entry == 1, return(0), return(-1)));
  if (line == 8, if (entry == 1, return(1), return(-1)));
  if (line == 9, if (entry == 1, return(0), return(-1)));
  error("Invalid inputs");
}

/*Tests if all the squares in the orbit of Psi*[x, y]~ are missing up to B^2. Returns 1 if they are missing, 0 else.*/
psi_missingsquares(xy, B, entry) = {
  for (i = 1, B, if (psi_rep(xy[1], xy[2], i^2, entry), return(0)));
  return(1);
}

/*Tests if all the squares in the orbit of Psi_1*[x, y]~ are missing up to B^2. Returns 1 if they are missing, 0 else.*/
psi1_missingsquares(xy, B, entry) = {
  my(M, v);
  M = [[1, 1;0, 1], [1, 0;4, 1]];
  v = vector(B, i, i^2);
  return(semigroup_missinglist(M, v, xy, entry));
}

/*Given the pair (x, y), returns the line in Table 1 that they correspond to.*/
table1_line(xy) = {
  my(x, y);
  x = xy[1]; y = xy[2];
  if (gcd(x, y) > 1, return(0));
  if (y % 4 == 0,
    if (x % 4 == 1,/*(1, 0) mod 4: lines 1, 2*/
      if (kronecker(x, y) == -1, return(1), return(2));
    );
    /*(3, 0) mod 4: lines 3, 4*/
    if (kronecker(x, y) == -kronecker(-1, y), return(3), return (4)); 
  );
  if (y % 4 == 1,/*(*, 1) mod 4: lines 5, 6*/
    if (kronecker(x, y) == -1, return(5), return(6));
  );
  if (y % 4 == 2, return(7));/*(*, 2) mod 4: line 7*/
  if (kronecker(x, y) == -1, return(8), return(9));/*(*, 3) mod 4: lines 8, 9*/
}

