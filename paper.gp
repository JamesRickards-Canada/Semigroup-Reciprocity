/*Various methods to test results in the paper.*/

/*SECTION 1: TESTING*/

/*Tests that the even continued fraction does correspond to orbits of [1,0]~ as described at the start of section 2. We pick random positive rational numbers x/y with 1<=x,y<=B, and do this test n times.*/
test_evencontfrac(n, B) = {
  my (x, c, M);
  printf("Testing %d random positive rational numbers with num/denom bounded by %d\n", n ,B);
  for (i = 1, n,
    x = (random(B) + 1)/(random(B) + 1);
    c = contfraceven(x);
    M = contfractoword(c);
    if (M[1, 1] != numerator(x) || M[2,1] != denominator(x), printf("Test failed with x=%Ps\n", x);break);
  );
  printf("All tests passed.\n");
}


/*SECTION 2: SUPPORTING METHODS*/

/*Returns a random element of Gamma_1(4)^{>=0}, where the size of the coefficients are bounded by n (not necessarily uniformly, but should be reasonably close).*/
gamma14geq0_random(n) = {
  my (bd1, bd2, a, b, c, d, bd, s);
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
  my (bd1, bd2, a, b, c, d, bd, s);
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
  my (a, b, c, d, bd, s);
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





/*Tests Proposition 3.2 by calling kronactioncorrect on n random matrices in SL(2, Z)^{>=0} with entries bounded by B. The values of x, y tried are all valid pairs with xymin<=x, y<=xymax*/
testkronaction(B, n, xymin, xymax) = {
  my(pair, v, M);
  pair = [[1, 1;0, 1], [1, 0;1, 1]];
  v = semigroupmats(pair, B);
  for (i = 1, n,
    M = v[random(#v) + 1];
    for (x = xymin, xymax,
      if (gcd(x, M[2, 2]) == 1,
        for (y = xymin, xymax,
          if (gcd(x, y) == 1,
            if (kronactioncorrect(M, [x, y]) != 1, printf("WRONG FORMULA: %Ps %d %d\n", M, x, y);error("FORMULA FAILED"));
          );
        );
      );
    );
    if (i % 50 == 0, printf("%d matrices tried\n", i));
  );
}

/*Tests Proposition 3.2: M=[a, b;c, d] in SL(2, Z)^{>=0} a, b, c, d >= 0, x, y>=0 coprime, gcd(x, d)=1, this proposition gives a formula for kron(ax+by/cx+dy). This function returns 1 if and only if the formula is correct for the given inputs.*/
kronactioncorrect(M, xy) = {
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
  if (x % 4 == 0,
    if (y % 4 == 1,/*(0, 1) mod 4*/
      if (kronecker(x, y) == -1, return(1), return(2));
    );
    if (kronecker(x, y) == -1, return(3), return(4));/*(0, 3) mod 4*/
  );
  if (x % 4 == 1,
    if (y % 4 == 0,/*(1, 0) mod 4*/
      if (kronecker(x, y) == -1, return(5), return(6));
    );
    if (y % 4 == 1,/*(1, 1) mod 4*/
      if (kronecker(x, y) == -1, return(7), return(8));
    );
    if (y % 4 == 2, return(9));/*(1, 2) mod 4*/
    if (kronecker(x, y) == -1, return(10), return(11));/*(1, 3) mod 4*/
  );
  if (x % 4 == 2,
    if (y % 4 == 1,/*(2, 1) mod 4*/
      if (kronecker(x, y) == -1, return(12), return(13));
    );
    if (kronecker(x, y) == -1, return(14), return(15));/*(2, 3) mod 4*/
  );
  if (y % 4 == 0,
    if (kronecker(x, y) == -kronecker(-1, oddpart(y)), return(16), return(17));/*(3, 0) mod 4*/
  );
  if (y % 4 == 1,
    if (kronecker(x, y) == -1, return(18), return(19));/*(3, 1) mod 4*/
  );
  if (y % 4 == 2, return(20));/*(3, 2) mod 4*/
  if (kronecker(x, y) == -1, return(21), return(22));/*(3, 3) mod 4*/
}

/*Returns 1 if we know there are no squares in the orbit of Psi1*xy due to a reciprocity obstruction, 0 if we expect there to be squares (based on Table 1), and -1 if there are no squares due to a congruence obstruction. */
table1_prediction(xy, entry) = {
  my(line);
  line = table1_line(xy);
  if (line == 1,  return(1));
  if (line == 2,  return(0));
  if (line == 3,  if (entry == 1, return(1), return(-1)));
  if (line == 4,  if (entry == 1, return(0), return(-1)));
  if (line == 5,  return(1));
  if (line == 6,  return(0));
  if (line == 7,  return(1));
  if (line == 8,  return(0));
  if (line == 9,  if (entry == 1, return(0), return(-1)));
  if (line == 10, if (entry == 1, return(1), return(-1)));
  if (line == 11, if (entry == 1, return(0), return(-1)));
  if (line == 12, return(1));
  if (line == 13, return(0));
  if (line == 14, if (entry == 1, return(1), return(-1)));
  if (line == 15, if (entry == 1, return(0), return(-1)));
  if (line == 16, if (entry == 1, return(-1), return(1)));
  if (line == 17, if (entry == 1, return(-1), return(0)));
  if (line == 18, return(1));
  if (line == 19, return(0));
  if (line == 20, if (entry == 1, return(0), return(-1)));
  if (line == 21, if (entry == 1, return(1), return(-1)));
  if (line == 22, if (entry == 1, return(0), return(-1)));
  error("Invalid inputs");
}

/*Returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi orbit for xy=[x, y]~.*/
table1_iscorrect_psi(xy, B, entry) = {
  my(pred);
  pred = table1_prediction(xy, entry);
  if (pred == -1, return(1));
  actual = psi_missingsquares(xy, B, entry);
  if (pred != actual, return(0));
  return(1);
}

/*Returns 1 if Table 1 gives the correct prediction for the missing squares up to B^2 in the Psi_1 orbit for xy=[x, y]~.*/
table1_iscorrect_psi1(xy, B, entry) = {
  my(pred);
  pred = table1_prediction(xy, entry);
  if (pred == -1, return(1));
  actual = psi1_missingsquares(xy, B, entry);
  if (pred != actual, return(0));
  return(1);
}

/*Tests Table 1 for all 1<=x, y<=xymax and squares up to B. If B is too small this will of course fail. Returns the vector of counts of how many pairs of each line (1 to 22) were tested.*/
table1_bigtest(xymax, B, whichpsi = 1) = {
  my(x, y, v);
  v = vector(22);
  if (whichpsi,
    for (x = 1, xymax,
      for (y = 1, xymax,
        if (gcd(x, y) == 1,
          v[table1_line([x, y])]++;
          if (table1_iscorrect_psi1([x, y], B, 1) != 1,printf("WRONG: (x, y, B, entry) = %d %d %d 1\n", x, y, B));
          if (table1_iscorrect_psi1([x, y], B, 2) != 1,printf("WRONG: (x, y, B, entry) = %d %d %d 2\n", x, y, B));
        );
      );
    );
    return(v);
  );
  for (x = 1, xymax,
    for (y = 1, xymax,
      if (gcd(x, y) == 1,
        v[table1_line([x, y])]++;
        if (table1_iscorrect_psi([x, y], B, 1) != 1,printf("WRONG: (x, y, B, entry) = %d %d %d 1\n", x, y, B));
        if (table1_iscorrect_psi([x, y], B, 2) != 1,printf("WRONG: (x, y, B, entry) = %d %d %d 2\n", x, y, B));
      );
    );
  );
  return(v);
}


/*Supporting methods*/

/*Returns the odd part of an integer.*/
oddpart(w)={
  my(v);
  v = valuation(w, 2);
  return(w >> v);
}