/*Various methods to test results in the paper.*/

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
            if (kronactioncorrect(M, x, y) != 1, printf("WRONG FORMULA: %Ps %d %d\n", M, x, y);error("FORMULA FAILED"));
          );
        );
      );
    );
    if (i % 50 == 0, printf("%d matrices tried\n", i));
  );
}

/*Tests Proposition 3.2: M=[a, b;c, d] in SL(2, Z)^{>=0} a, b, c, d >= 0, x, y>=0 coprime, gcd(x, d)=1, this proposition gives a formula for kron(ax+by/cx+dy). This function returns 1 if and only if the formula is correct for the given inputs.*/
kronactioncorrect(M, x, y) = {
  my(a, b, c, d, A, B, cxpdy, C, D, alpha, mu);
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







/*Supporting methods*/

/*Returns the odd part of an integer.*/
oddpart(w)={
  my(v);
  v = valuation(w, 2);
  return(w >> v);
}