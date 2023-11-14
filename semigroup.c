/*Methods for thin (semi)-groups.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "semigroup.h"

/*STATIC DECLARATIONS*/

/*SECTION 1: BASIC (SEMI)GROUP METHODS*/
/*SECTION 2: SEMIGROUP ORBITS*/
static long** ZM_to_Cmatrix(GEN M);
static void findmissing(long Bmin, long Bmax, long ***mats, long nmats, long *start, long n, long *res, long nres, long modulus, long entry);
static void missing_tofile(long blocks, unsigned long **rclass, long Bmin, long Bmax, long *start, long n, long *res, long nres, long ubits, long modulus);
static void findmissing_parabolic(long Bmin, long Bmax, long ***mats, long ***matsinv, long nmats, long *start, long n, long *res, long nres, long modulus, long entry);
static int findmissinglist(hashtable *hmiss, long Bmin, long Bmax, long ***mats, long ***matsinv, long nmats, long *start, long n, long entry);
/*SECTION 3: LINEAR REGRESSION*/

/*MAIN BODY*/


/*SECTION 1: BASIC (SEMI)GROUP METHODS*/

/*For M a hyperbolic matrix with positive entries, returns M as a word in L and R, where L=[1,1;0,1] and R=[1,0;1,1]. This is a Vecsmall with 0 representing L and 1 representing R.*/
GEN
LRword(GEN M)
{
  pari_sp av = avma;
  long maxdepth = 100, i = 1;
  long a = itos(gcoeff(M, 1, 1)), b = itos(gcoeff(M, 1, 2));
  long c = itos(gcoeff(M, 2, 1)), d = itos(gcoeff(M, 2, 2));
  GEN v = cgetg(maxdepth + 1, t_VECSMALL);
  while(c || b) {
    if (a > c) {
	  v[i] = 0;/*L*/
	  a -= c;
	  b -= d;
	}
	else {
	  v[i] = 1;/*R*/
	  c -= a;
	  d -= b;
	}
	i++;
	if (i > maxdepth) {
	  maxdepth <<= 1;
	  v = vecsmall_lengthen(v, maxdepth);
	}
  }
  return gerepileupto(av, vecsmall_shorten(v, i - 1));
}

/*Returns the largest entry of the ZM M as a long.*/
INLINE long
ZM_largest_s(GEN M)
{
  long x = itos(gcoeff(M, 1, 1)), lM = lg(M), i, j;
  for (j = 2; j < lM; j++) x = maxss(x, itos(gcoeff(M, 1, j)));
  for (i = 2; i < lM; i++) {
    for (j = 1; j < lM; j++) x = maxss(x, itos(gcoeff(M, i, j)));
  }
  return x;
}

/*We estimate the growth, by computing all elements of the orbit (given by the column vector start) that are at most binsize * Nbins. We then group their sizes into bins, and runs a linear regression on the data to determine the growth rate. Returns [c, nu, R^2], where up to N we have c*N^nu, and the R^2 value of the regression is given.*/
GEN
semigroup_growth(GEN mats, long binsize, long Nbins, GEN start, long prec)
{
  pari_sp av = avma;
  if (!start) start = mkcol2s(1, 1);
  if (typ(start) == t_VEC) start = gtocol(start);
  GEN v = semigroup_mats(mats, binsize * Nbins);
  long lv, i;
  GEN images = cgetg_copy(v, &lv);
  for (i = 1; i < lv; i++) gel(images, i) = ZM_ZC_mul(gel(v, i), start);
  GEN images_uniq = vecsort0(images, NULL, 8);/*Sort and remove double counting.*/
  long lu = lg(images_uniq);
  GEN counts = const_vecsmall(Nbins, 0);
  for (i = 1; i < lu; i++) {
	GEN r;
	GEN SOS = addii(sqri(gmael(images_uniq, i, 1)), sqri(gmael(images_uniq, i, 2)));
	long rt = itos(sqrtremi(SOS, &r));
	if (!isintzero(r)) rt++;/*rt=ceil(size of the vector)*/
	long rem, quo;
	quo = sdivss_rem(rt, binsize, &rem);
	if (rem) quo++;/*quo=ceil(rt/binsize)*/
	if (quo <= Nbins) counts[quo]++;
  }
  counts = gerepilecopy(av, counts);/*Clear garbage, clone counts to the heap as this is all we need.*/
  GEN cumu = cgetg(Nbins + 1, t_COL);/*Cumulative counts*/
  gel(cumu, 1) = stoi(counts[1]);
  for (i = 2; i <= Nbins; i++) gel(cumu, i) = addis(gel(cumu, i - 1), counts[i]);
  for (i = 1; i <= Nbins; i++) gel(cumu, i) = glog(gel(cumu, i), prec);/*Take the log of all entries.*/
  GEN X = cgetg(Nbins + 1, t_MAT);
  for (i = 1; i <= Nbins; i++) {/*X is matrix with X[,i]=[1,log(binsize*i)].*/
	gel(X, i) = mkcol2(gen_1, glog(stoi(binsize * i), prec));
  }
  GEN reg = OLS(X, cumu, 1);/*cumu[i]=c*(binsize*i)^nu -> take logarithms.*/
  GEN c = gexp(gmael(reg, 1, 1), prec);
  return gerepilecopy(av, mkvec3(c, gmael(reg, 1, 2), gel(reg, 2)));
}

/*Assume the matrices in mats all have positive entries and infinite order. This returns the matrices in the semigroup they generate that have all entries at most N. If there are relations, the corresponding matrices will get counted multiple times.*/
GEN
semigroup_mats(GEN mats, long N)
{
  pari_sp av = avma;
  long lmats = lg(mats), nM = lg(gel(mats, 1)) - 1;
  long maxdepth = 100, maxfound = 10000, vind = 0;
  GEN v = cgetg(maxfound + 1, t_VEC);
  GEN depthseq = cgetg(maxdepth + 1, t_VEC);
  GEN swaps = const_vecsmall(maxdepth, 0);
  gel(depthseq, 1) = matid(nM);
  long ind = 2;
  while (ind > 1) {
    long cind = ++swaps[ind];
    if (cind == lmats) {/*Overflowed, go back.*/
      swaps[ind] = 0;
      ind--;
      continue;
    }
    GEN M = ZM_mul(gel(depthseq, ind - 1), gel(mats, cind));
    long maxM = ZM_largest_s(M);
    if (maxM > N) continue;/*Too big! Go back.*/
    vind++;
    if (vind > maxfound) {/*Double size.*/
      maxfound <<= 1;
      v = vec_lengthen(v, maxfound);
    }
    gel(v, vind) = M;
    gel(depthseq, ind) = M;
    ind++;
    if (ind > maxdepth) {/*Double depth sequence / swaps length.*/
      maxdepth <<= 1;
      depthseq = vec_lengthen(depthseq, maxdepth);
      swaps = vecsmall_lengthen(swaps, maxdepth);
      long i;
      for (i = ind; i <= maxdepth; i++) swaps[i] = 0;
    }
  }
  return gerepilecopy(av, vec_shorten(v, vind));
}


/*SECTION 2: SEMIGROUP ORBITS*/

/*Takes in a ZM M and returns a malloc'ed array representing M as a C array.*/
static long**
ZM_to_Cmatrix(GEN M)
{
  long ncols = lg(M) - 1, nrows = nbrows(M), i, j;
  long **A = (long **)pari_malloc(nrows * sizeof(long *));
  for (i = 0; i < nrows; i++) {
	A[i] = (long *)pari_malloc(ncols * sizeof(long));
	for (j = 0; j < ncols; j++) A[i][j] = itos(gcoeff(M, i + 1, j + 1));
  }
  return A;
}

/*Runs C code to find the missing positive entries up to the given bound, then saves them to a file. Must also supply the congruence conditions.*/
GEN
semigroup_missing(GEN mats, GEN B, GEN start, GEN congs, long entry)
{
  pari_sp av = avma;
  long Bmin = 0, Bmax = 0, t = typ(B);
  if (t == t_INT) { Bmin = 1; Bmax = itos(B); }
  else if (t == t_VEC || t == t_COL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 1)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 2)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = itou(gel(B, 1));
    Bmax = itou(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = B[1];
    Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  if (Bmin <= 0) Bmin = 1;
  long n = lg(start) - 1;
  if (entry < 0 || entry > n) pari_err_TYPE("entry must be an index from 1 to n, the dimension of the vectors computed", stoi(entry));
  long nmats = lg(mats) - 1, i;
  long ***Cmats = (long ***)pari_malloc(nmats * sizeof(long **));/*Stores the matrices as C arrays*/
  for (i = 0; i < nmats; i++) Cmats[i] = ZM_to_Cmatrix(gel(mats, i + 1));
  long *Cstart = (long *)pari_malloc(n * sizeof(long));/*Stores the starting tuple.*/
  for (i = 0; i < n; i++) Cstart[i] = itos(gel(start, i + 1));
  if (typ(congs) != t_VEC || lg(congs) < 3 || typ(gel(congs, 1)) != t_VEC) pari_err_TYPE("congs must be the relevant congruence conditions, formatted as [[r1, r2, ..., rk], n] for ri modulo n, which should be sorted.", congs);
  long modulus = itos(gel(congs, 2)), nres = lg(gel(congs, 1)) - 1;
  long *residues = (long *)pari_malloc(nres * sizeof(long));
  for (i = 0; i < nres; i++) residues[i] = itos(gmael(congs, 1, i + 1));
  findmissing(Bmin, Bmax, Cmats, nmats, Cstart, n, residues, nres, modulus, entry - 1);
  pari_free(residues);
  pari_free(Cstart);
  long j;
  for (i = 0; i < nmats; i++) {
	for (j = 0; j < n; j++) pari_free(Cmats[i][j]);
	pari_free(Cmats[i]);
  }
  pari_free(Cmats);
  set_avma(av);
  /*
  Now load them maybe!!
  */
  return gen_0;
}

/*Updates rclass with the new value. Assumes we are at most Bmax*/
INLINE void
missing_update(unsigned long **rclass, unsigned long *bitswap, long *maxmiss_block, int *maxmiss_bit, long Base, long *Bmax, long value, long *res, long nres, long Bmin, long ubits, long modulus)
{
  long shifted = value - Base;
  if (shifted <= 0) return;
  long b = shifted % modulus;
  long a = shifted / modulus;/*shifted=modulus*a+b. b gives the block to insert into, and we need to save "a"*/
  long v = a % ubits;
  long u = a / ubits;/*a=ubits * u+v. u gives the entry of the array, v gives the bit to swap.*/
  rclass[b][u] |= bitswap[v];
  if (u != maxmiss_block[b]) return;/*Not swapping in the last block.*/
  if (v != maxmiss_bit[b]) return;/*Not swapping the last bit.*/
  long i;
  for (i = v - 1; i >= 0; i--) {
    if (!(rclass[b][u] & bitswap[i])) {/*Zero, so we stop here.*/
      maxmiss_bit[b] = i;
      goto BMAXUPDATE;
    }
  }
  long j;/*We made it out of the block.*/
  for (j = u - 1; j >= 0; j--) {
    for (i = ubits - 1; i >= 0; i--) {
      if (!(rclass[b][j] & bitswap[i])) {/*Zero, so we stop here.*/
        maxmiss_block[b] = j;
        maxmiss_bit[b] = i;
        goto BMAXUPDATE;
      }
    }
  }
  maxmiss_block[b] = -1;/*Everything is gone!*/
  maxmiss_bit[b] = -1;
  BMAXUPDATE:;/*See if we update Bmax*/
  if (value != *Bmax) return;/*Did not remove the largest exception.*/
  long worst = res[0];
  for (i = 1; i < nres; i++) {
    if (maxmiss_block[res[i]] > maxmiss_block[worst]) {
      worst = res[i];
      continue;
    }
    if (maxmiss_block[res[i]] < maxmiss_block[worst]) continue;
    if (maxmiss_bit[res[i]] >= maxmiss_bit[worst]) worst = res[i];
  }
  *Bmax = Base + ((maxmiss_block[worst] * ubits + maxmiss_bit[worst]) * modulus ) + worst;/*Update Bmax*/
  if (Bmin > *Bmax) *Bmax = 0;/*All eliminated, let's quit early!*/
}

/*executes semigroup_missing*/
static void
findmissing(long Bmin, long Bmax, long ***mats, long nmats, long *start, long n, long *res, long nres, long modulus, long entry)
{
  long ubits = sizeof(unsigned long) << 3;/*How many bits in an unsigned long*/
  unsigned long *bitswap = (unsigned long*)pari_malloc(ubits * sizeof(unsigned long)), i;/*Used for swapping bits of unsigned longs.*/
  bitswap[0] = 1;
  for (i = 1; i < ubits; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long Base = Bmin - 1;
  Base -= (Base % modulus);/*We want to start at a multiple of modulus to not ruin the mod stuff.*/
  long classmax = (Bmax - Base)/ modulus + 1;/*Maximal number of entries found in each residue class.*/
  long blocks = ((classmax - 1) / ubits) + 1;/*This is the number of unsigned longs we need to store in each class.*/
  unsigned long **rclass = (unsigned long **)pari_malloc(modulus * sizeof(unsigned long *));/*Stores pointers to the individual classes.*/
  for (i = 0; i < nres; i++) {
    rclass[res[i]] = (unsigned long *)pari_calloc(blocks * sizeof(unsigned long));/*pari_calloc the classes we want, since we want them as 0 to start.*/
    if (!rclass[res[i]]) {
      printf("Insufficient memory to allocate to store the curvatures.\n");
      exit(1);
    }
  }
  long Bmaxoriginal = Bmax;/*Save for later in case we change it.*/
  long *maxmiss_block = (long *)pari_malloc(modulus * sizeof(long));/*Tracks the largest block in the residue class that still contains 0's*/
  int *maxmiss_bit = (int *)pari_malloc(modulus * sizeof(int));/*Tracks the largest bit of said class that is non-zero.*/
  int Bmaxmod = Bmax % modulus;/*We now will initialize it.*/
  long Bmaxbase = Bmax - Bmaxmod, j;
  int foundlargest = 0;
  for (i = 0; i < nres; i++) {
    long mc = Bmaxbase + res[i];
    if (res[i] > Bmaxmod) {
      if (!foundlargest) {
        foundlargest = 1;
        if (i) Bmax = Bmaxbase + res[i - 1];/*Update to the largest actually possible value.*/
        else Bmax = Bmaxbase + res[nres - 1] - modulus;
      }
      mc -= modulus;/*The last curvature of this type.*/
    }
    if (mc < Bmin) {
      maxmiss_bit[res[i]] = -1;
      maxmiss_block[res[i]] = -1;
      continue;
    }
    mc -= Base;/*Shift it back.*/
    long a = mc / modulus;/*Save a in block res[i]*/
    maxmiss_bit[res[i]] = a % ubits;
    maxmiss_block[res[i]] = a / ubits;
  }
  if (!foundlargest) Bmax = Bmaxbase + res[nres - 1];/*They all fit in with the +.*/
  for (i = 0; i <= nres; i++) {
    if (i == nres) {/*All were -1, so our range actually contains no residues.*/
      Bmax = 0;
    }
    else if (maxmiss_bit[res[i]] != -1) break;/*We have at least one valid residue.*/
  }
  long maxdepth = 100;/*Maximal depth, to start.*/
  long **depthseq = (long **)pari_malloc(maxdepth * sizeof(long *));/*Tracks the sequence of swaps.*/
  if (!depthseq) {
    printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps.*/
  if (!swaps) {
    printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) {
    depthseq[i] = (long *)pari_malloc(n * sizeof(long));
    swaps[i] = -1;/*Initialize to all -1's*/
  }
  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, start[entry], res, nres, Bmin, ubits, modulus);
  for (i = 0; i < n; i++) depthseq[0][i] = start[i];
  long ind = 1;/*Which depth we are working at.*/
  while (ind > 0) {/*We are coming in trying to swap this circle out.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == nmats) {/*Overflowed, go back.*/
      swaps[ind] = -1;
      ind--;
      continue;
    }
    long lastind = ind - 1;
	int toofar = 0;
	for (i = 0; i < n; i++) {/*Apply the matrix action.*/
	  long en = 0;
	  for (j = 0; j < n; j++) en += mats[cind][i][j] * depthseq[lastind][j];
	  if (en > Bmax) toofar = 1;
	  depthseq[ind][i] = en;
	}
	if (toofar) {/*We still want to update the entry on the off chance that we overflowed but are OK in our desired place*/
	  if (depthseq[ind][entry] <= Bmax) missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, depthseq[ind][entry], res, nres, Bmin, ubits, modulus);
	  continue;/*A single entry was too large.*/
	}
	missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, depthseq[ind][entry], res, nres, Bmin, ubits, modulus);
	/*Do the bitswap to update */
    ind++;
    if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
      long newdepth = maxdepth << 1;/*Double it.*/
      depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
      if (!depthseq) {
        printf("Insufficient memory to pari_reallocate the depth sequence.\n");
        exit(1);
      }
      swaps = pari_realloc(swaps, newdepth * sizeof(int));
      if (!swaps) {
        printf("Insufficient memory to pari_reallocate the swaps.\n");
        exit(1);
      }
      for (i = maxdepth; i < newdepth; i++) {
        depthseq[i] = (long *)pari_malloc(n * sizeof(long));
        swaps[i] = -1;
      }
      maxdepth = newdepth;
    }
  }
  missing_tofile(blocks, rclass, Bmin, Bmaxoriginal, start, n, res, nres, ubits, modulus);/*Print to file.*/
  /*Time to free all of the allocated memory.*/
  pari_free(swaps);
  for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
  pari_free(depthseq);
  pari_free(maxmiss_block);
  pari_free(maxmiss_bit);
  for (i = 0; i < nres; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);
  pari_free(bitswap);
  return;
}

/*Prints the found data to a file.*/
static void
missing_tofile(long blocks, unsigned long **rclass, long Bmin, long Bmax, long *start, long n, long *res, long nres, long ubits, long modulus)
{
  long Base = Bmin - 1, i;
  Base = Base - (Base % modulus);/*We started at a multiple of modulus to not ruin the mod stuff.*/
  if (!pari_is_dir("missing")) {
    int s = system("mkdir -p missing");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY missing");
  }
  
  char *fname = stack_sprintf("missing/%ld", start[0]);
  for (i = 1; i < n; i++) fname = stack_sprintf("%s_%ld", fname, start[i]);
  fname = stack_sprintf("%s_%ld-to-%ld.dat", fname, Bmin, Bmax);
  FILE *F = fopen(fname, "w");
  for (i = 0; i < nres; i++) {/*The ith residue*/
    long b = res[i];
    fprintf(F, "[");
    int found = 0;
    long u;
    for (u = 0; u < blocks; u++) {/*Figure out each block.*/
      unsigned long val = rclass[b][u];
      long v;
      for (v = 0; v < ubits; v++) {
        if (!(val & 1)) {/*A missing value*/
          long a = ubits * u + v;
          long num = modulus * a + b + Base;/*The correct one!*/
          if (num <= Bmax && num >= Bmin) {/*Correct range! We only check >=Bmin due to the initial shift by up to modulus.*/
            if (found) fprintf(F, ", %ld", num);/*Print it to the file.*/
            else { found = 1; fprintf(F, "%ld", num); }
          }
        }
        val >>= 1;/*Shift right one.*/
      }
    }
    fprintf(F, "]\n");
  }
  fclose(F);
}

/*semigroup_missing, except we store the sequence of swaps differently. This is slightly slower but suggested if a matrix is parabolic, as storing the whole depth sequence and swaps is extremely costly in terms of memory when B is large.*/
GEN
semigroup_missing_parabolic(GEN mats, GEN B, GEN start, GEN congs, long entry)
{
  pari_sp av = avma;
  long Bmin = 0, Bmax = 0, t = typ(B);
  if (t == t_INT) { Bmin = 1; Bmax = itos(B); }
  else if (t == t_VEC || t == t_COL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 1)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    if (typ(gel(B, 2)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = itou(gel(B, 1));
    Bmax = itou(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
    if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
    Bmin = B[1];
    Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  if (Bmin <= 0) Bmin = 1;
  long n = lg(start) - 1;
  if (entry < 0 || entry > n) pari_err_TYPE("entry must be an index from 1 to n, the dimension of the vectors computed", stoi(entry));
  long nmats = lg(mats) - 1, i;
  long ***Cmats = (long ***)pari_malloc(nmats * sizeof(long **));/*Stores the matrices as C arrays*/
  long ***Cmatsinv = (long ***)pari_malloc(nmats * sizeof(long **));/*Inverses*/
  for (i = 0; i < nmats; i++) {
	Cmats[i] = ZM_to_Cmatrix(gel(mats, i + 1));
	Cmatsinv[i] = ZM_to_Cmatrix(QM_inv(gel(mats, i + 1)));
  }
  long *Cstart = (long *)pari_malloc(n * sizeof(long));/*Stores the starting tuple.*/
  for (i = 0; i < n; i++) Cstart[i] = itos(gel(start, i + 1));
  if (typ(congs) != t_VEC || lg(congs) < 3 || typ(gel(congs, 1)) != t_VEC) pari_err_TYPE("congs must be the relevant congruence conditions, formatted as [[r1, r2, ..., rk], n] for ri modulo n, which should be sorted.", congs);
  long modulus = itos(gel(congs, 2)), nres = lg(gel(congs, 1)) - 1;
  long *residues = (long *)pari_malloc(nres * sizeof(long));
  for (i = 0; i < nres; i++) residues[i] = itos(gmael(congs, 1, i + 1));
  findmissing_parabolic(Bmin, Bmax, Cmats, Cmatsinv, nmats, Cstart, n, residues, nres, modulus, entry - 1);
  pari_free(residues);
  pari_free(Cstart);
  long j;
  for (i = 0; i < nmats; i++) {
	for (j = 0; j < n; j++) {
	  pari_free(Cmats[i][j]);
	  pari_free(Cmatsinv[i][j]);
	}
	pari_free(Cmats[i]);
	pari_free(Cmatsinv[i]);
  }
  pari_free(Cmats);
  pari_free(Cmatsinv);
  set_avma(av);
  /*
  Now load them maybe!!
  */
  return gen_0;
}

/*executes semigroup_missing_par*/
static void
findmissing_parabolic(long Bmin, long Bmax, long ***mats, long ***matsinv, long nmats, long *start, long n, long *res, long nres, long modulus, long entry)
{
  long ubits = sizeof(unsigned long) << 3;/*How many bits in an unsigned long*/
  unsigned long *bitswap = (unsigned long*)pari_malloc(ubits * sizeof(unsigned long)), i;/*Used for swapping bits of unsigned longs.*/
  bitswap[0] = 1;
  for (i = 1; i < ubits; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long Base = Bmin - 1;
  Base -= (Base % modulus);/*We want to start at a multiple of modulus to not ruin the mod stuff.*/
  long classmax = (Bmax - Base)/ modulus + 1;/*Maximal number of entries found in each residue class.*/
  long blocks = ((classmax - 1) / ubits) + 1;/*This is the number of unsigned longs we need to store in each class.*/
  unsigned long **rclass = (unsigned long **)pari_malloc(modulus * sizeof(unsigned long *));/*Stores pointers to the individual classes.*/
  for (i = 0; i < nres; i++) {
    rclass[res[i]] = (unsigned long *)pari_calloc(blocks * sizeof(unsigned long));/*pari_calloc the classes we want, since we want them as 0 to start.*/
    if (!rclass[res[i]]) {
      printf("Insufficient memory to allocate to store the curvatures.\n");
      exit(1);
    }
  }
  long Bmaxoriginal = Bmax;/*Save for later in case we change it.*/
  long *maxmiss_block = (long *)pari_malloc(modulus * sizeof(long));/*Tracks the largest block in the residue class that still contains 0's*/
  int *maxmiss_bit = (int *)pari_malloc(modulus * sizeof(int));/*Tracks the largest bit of said class that is non-zero.*/
  int Bmaxmod = Bmax % modulus;/*We now will initialize it.*/
  long Bmaxbase = Bmax - Bmaxmod, j;
  int foundlargest = 0;
  for (i = 0; i < nres; i++) {
    long mc = Bmaxbase + res[i];
    if (res[i] > Bmaxmod) {
      if (!foundlargest) {
        foundlargest = 1;
        if (i) Bmax = Bmaxbase + res[i - 1];/*Update to the largest actually possible value.*/
        else Bmax = Bmaxbase + res[nres - 1] - modulus;
      }
      mc -= modulus;/*The last curvature of this type.*/
    }
    if (mc < Bmin) {
      maxmiss_bit[res[i]] = -1;
      maxmiss_block[res[i]] = -1;
      continue;
    }
    mc -= Base;/*Shift it back.*/
    long a = mc / modulus;/*Save a in block res[i]*/
    maxmiss_bit[res[i]] = a % ubits;
    maxmiss_block[res[i]] = a / ubits;
  }
  if (!foundlargest) Bmax = Bmaxbase + res[nres - 1];/*They all fit in with the +.*/
  for (i = 0; i <= nres; i++) {
    if (i == nres) {/*All were -1, so our range actually contains no residues.*/
      Bmax = 0;
    }
    else if (maxmiss_bit[res[i]] != -1) break;/*We have at least one valid residue.*/
  }
  long maxdepth = 100;/*Maximal depth, to start.*/
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps.*/
  if (!swaps) {
    printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  unsigned long *swapfreq = (unsigned long*)pari_malloc(maxdepth * sizeof(unsigned long));/*Tracks the number of times each swap is executed in a row.*/
  for (i = 0; i < maxdepth; i++) swaps[i] = -1;/*Initialize to all -1's*/
  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, start[entry], res, nres, Bmin, ubits, modulus);/*Insert first entry.*/
  long v[n];/*Stores the current orbit entry.*/
  for (i = 0; i < n; i++) v[i] = start[i];
  long swapind = 1;/*Index of the swapping depth.*/
  while (swapind) {/*We are coming in trying to swap at this depth. We always pretend to insert the swaps at a new swapdepth, and double check if it actually is a new entry or we are just increasing the swapfreq.*/
    int cind = ++swaps[swapind];/*Increment the swapping index.*/
    if (cind == nmats) {/*Overflowed, go back.*/
	  swapind--;
	  if (!swapind) continue;/*Done!*/
	  cind = swaps[swapind];/*How to move backwards.*/
	  long vnew[n];
	  for (i = 0; i < n; i++) {/*Move v backwards.*/
		long en = 0;
		for (j = 0; j < n; j++) en += matsinv[cind][i][j] * v[j];
		vnew[i] = en;
	  }
	  for (i = 0; i < n; i++) v[i] = vnew[i];/*Update v to the previous.*/
	  if (--swapfreq[swapind]) {/*Decrease and check if non-zero*/
	    swapind++;
		swaps[swapind] = swaps[swapind - 1];/*Move on.*/
	  }
	  else swaps[swapind + 1] = -1;/*Reset the subsequent one.*/
      continue;
    }
	int toofar = 0;
	long vnew[n];
	for (i = 0; i < n; i++) {/*Apply the matrix action.*/
	  long en = 0;
	  for (j = 0; j < n; j++) en += mats[cind][i][j] * v[j];
	  if (en > Bmax) toofar = 1;
	  vnew[i] = en;
	}
	if (toofar) {/*We still want to update the entry on the off chance that we overflowed but are OK in our desired place*/
	  if (vnew[entry] <= Bmax) missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, vnew[entry], res, nres, Bmin, ubits, modulus);
	  continue;/*A single entry was too large.*/
	}
	for (i = 0; i < n; i++) v[i] = vnew[i];/*Update v,*/
	missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, v[entry], res, nres, Bmin, ubits, modulus);/*Do the bitswap to update */
	if (swapind > 1 && swaps[swapind - 1] == swaps[swapind]) {/*It gets added to the previous sequence.*/
	  swapfreq[swapind - 1]++;
	  swaps[swapind] = -1;
	}
	else {
	  swapfreq[swapind] = 1;
      swapind++;
      if (swapind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        if (!swaps) {
          printf("Insufficient memory to pari_reallocate the swaps.\n");
          exit(1);
        }
        for (i = maxdepth; i < newdepth; i++) swaps[i] = -1;
		swapfreq = pari_realloc(swapfreq, newdepth * sizeof(unsigned long));
        maxdepth = newdepth;
      }
	}
  }
  missing_tofile(blocks, rclass, Bmin, Bmaxoriginal, start, n, res, nres, ubits, modulus);/*Print to file.*/
  /*Time to free all of the allocated memory.*/
  pari_free(swaps);
  pari_free(swapfreq);
  pari_free(maxmiss_block);
  pari_free(maxmiss_bit);
  for (i = 0; i < nres; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);
  pari_free(bitswap);
  return;
}

/*Given a sorted list of positive integers miss, this returns 1 if none of them are in the given orbit, and 0 if they are.*/
int
semigroup_missinglist(GEN mats, GEN miss, GEN start, long entry)
{
  pari_sp av = avma;
  if (!RgV_is_ZVpos(miss)) pari_err_TYPE("miss must be a sorted vector of positive integers.", miss);
  long lmiss = lg(miss);
  if (lmiss == 1) return 1;/*Nothing to miss!*/
  long Bmin = itos(gel(miss, 1)), Bmax = itos(gel(miss, lmiss - 1));
  long n = lg(start) - 1;
  if (entry < 0 || entry > n) pari_err_TYPE("entry must be an index from 1 to n, the dimension of the vectors computed", stoi(entry));
  long nmats = lg(mats) - 1, i;
  long ***Cmats = (long ***)pari_malloc(nmats * sizeof(long **));/*Stores the matrices as C arrays*/
  long ***Cmatsinv = (long ***)pari_malloc(nmats * sizeof(long **));/*Inverses*/
  for (i = 0; i < nmats; i++) {
    Cmats[i] = ZM_to_Cmatrix(gel(mats, i + 1));
    Cmatsinv[i] = ZM_to_Cmatrix(QM_inv(gel(mats, i + 1)));
  }
  long *Cstart = (long *)pari_malloc(n * sizeof(long));/*Stores the starting tuple.*/
  for (i = 0; i < n; i++) Cstart[i] = itos(gel(start, i + 1));
  hashtable *hmiss = hash_create_ulong(lmiss, 1);/*Hashtable of the missing entries.*/
  for (i = 1; i < lmiss; i++) hash_insert(hmiss, (void *)itou(gel(miss, i)), NULL);/*Insert them.*/
  int result = findmissinglist(hmiss, Bmin, Bmax, Cmats, Cmatsinv, nmats, Cstart, n, entry - 1);
  pari_free(Cstart);
  long j;
  for (i = 0; i < nmats; i++) {
    for (j = 0; j < n; j++) {
      pari_free(Cmats[i][j]);
      pari_free(Cmatsinv[i][j]);
    }
    pari_free(Cmats[i]);
    pari_free(Cmatsinv[i]);
  }
  pari_free(Cmats);
  pari_free(Cmatsinv);
  hash_destroy(hmiss);
  return gc_int(av, result);
}

/*executes semigroup_missinglist*/
static int
findmissinglist(hashtable *hmiss, long Bmin, long Bmax, long ***mats, long ***matsinv, long nmats, long *start, long n, long entry)
{
  long maxdepth = 100, i;/*Maximal depth, to start.*/
  long **depthseq = (long **)pari_malloc(maxdepth * sizeof(long *));/*Tracks the sequence of moves.*/
  if (!depthseq) {
    printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps.*/
  if (!swaps) {
    printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) {
    depthseq[i] = (long *)pari_malloc(n * sizeof(long));
    swaps[i] = -1;/*Initialize to all -1's*/
  }
  if (hash_search(hmiss, (void *)start[entry])) {/*Initial entry was in our missing list.*/
    pari_free(swaps);
    for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
    pari_free(depthseq);
    return 0;
  }
  for (i = 0; i < n; i++) depthseq[0][i] = start[i];
  long ind = 1;/*Which depth we are working at.*/
  while (ind > 0) {/*We are coming in trying to swap this circle out.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == nmats) {/*Overflowed, go back.*/
      swaps[ind] = -1;
      ind--;
      continue;
    }
    long lastind = ind - 1, j;
    int toofar = 0;
    for (i = 0; i < n; i++) {/*Apply the matrix action.*/
      long en = 0;
      for (j = 0; j < n; j++) en += mats[cind][i][j] * depthseq[lastind][j];
      if (en > Bmax) toofar = 1;
      depthseq[ind][i] = en;
    }
    if (hash_search(hmiss, (void *)depthseq[ind][entry])) {/*Entry is in our missing list.*/
      pari_free(swaps);
      for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
      pari_free(depthseq);
      return 0;
    }
    if (toofar) continue;/*A single entry was too large.*/
    ind++;
    if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
      long newdepth = maxdepth << 1;/*Double it.*/
      depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
      if (!depthseq) {
        printf("Insufficient memory to pari_reallocate the depth sequence.\n");
        exit(1);
      }
      swaps = pari_realloc(swaps, newdepth * sizeof(int));
      if (!swaps) {
        printf("Insufficient memory to pari_reallocate the swaps.\n");
        exit(1);
      }
      for (i = maxdepth; i < newdepth; i++) {
        depthseq[i] = (long *)pari_malloc(n * sizeof(long));
        swaps[i] = -1;
      }
      maxdepth = newdepth;
    }
  }
  /*Time to free all of the allocated memory.*/
  pari_free(swaps);
  for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
  pari_free(depthseq);
  return 1;/*We made it!*/
}

/*Returns the largest entry of the ZC C as a long.*/
INLINE long
ZC_largest_s(GEN C)
{
  long x = itos(gel(C, 1)), lC = lg(C), i;
  for (i = 2; i < lC; i++) x = maxss(x, itos(gel(C, i)));
  return x;
}

/*Returns the orbit of mats*start up to max size of entries B. This is not as efficient as the other methods, so use this for intial observation finding, and use the faster methods for larger searches.*/
GEN
semigroup_orbit(GEN mats, long B, GEN start)
{
  pari_sp av = avma;
  start = ZC_copy(start);/*Ensure it is a column vector.*/
  long lmats = lg(mats);
  long maxdepth = 100, maxfound = 10000, vind = 0;
  GEN v = cgetg(maxfound + 1, t_VEC);
  GEN depthseq = cgetg(maxdepth + 1, t_VEC);
  GEN swaps = const_vecsmall(maxdepth, 0);
  gel(depthseq, 1) = start;
  long ind = 2;
  while (ind > 1) {
    long cind = ++swaps[ind];
    if (cind == lmats) {/*Overflowed, go back.*/
      swaps[ind] = 0;
      ind--;
      continue;
    }
    GEN M = ZM_ZC_mul(gel(mats, cind), gel(depthseq, ind - 1));
    long maxM = ZC_largest_s(M);
    if (maxM > B) continue;/*Too big! Go back.*/
    vind++;
    if (vind > maxfound) {/*Double size.*/
      maxfound <<= 1;
      v = vec_lengthen(v, maxfound);
    }
    gel(v, vind) = M;
    gel(depthseq, ind) = M;
    ind++;
    if (ind > maxdepth) {/*Double depth sequence / swaps length.*/
      maxdepth <<= 1;
      depthseq = vec_lengthen(depthseq, maxdepth);
      swaps = vecsmall_lengthen(swaps, maxdepth);
      long i;
      for (i = ind; i <= maxdepth; i++) swaps[i] = 0;
    }
  }
  return gerepilecopy(av, vec_shorten(v, vind));
}


/*SECTION 3: LINEAR REGRESSION*/

/*Perform ordinary least squares regression. X is a matrix whose columns are the parameters, and y is a column vector of results. Must include linear term as first variable of X. The formula is B=Bhat=(X*X^T)^(-1)Xy, for the ordinary least squares regression for y=X^T*B+error (formula differs to Wikipedia due to X here being the transpose of what they define there. Returns either best fit or [best fit, R^2]*/
GEN
OLS(GEN X, GEN y, int retrsqr)
{
  pari_sp av = avma;
  if (typ(y) != t_COL) y = gtocol(y);
  if (lg(y) != lg(X)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  GEN Xy = RgM_RgC_mul(X, y);
  GEN XTX = RgM_multosym(X, shallowtrans(X));/*X*X^T, which is symmetric*/
  GEN fit = RgM_solve(XTX, Xy);/*Best fit.*/
  if (!fit) pari_err(e_MISC, "Could not compute matrix inverse. Error with given data or with precision?");
  if (!retrsqr) return gerepileupto(av, fit);
  GEN rsqr = rsquared(X, y, fit);
  return gerepilecopy(av, mkvec2(fit, rsqr));
}

/*Performs OLS where we have one independant variable and assume the intercept is 0 (so y=ax). The formua is now sum(x*y)/sum(x^2).*/
GEN
OLS_nointercept(GEN x, GEN y, int retrsqr)
{
  pari_sp av = avma;
  GEN xysum = gen_0, xsqrsum = gen_0;
  long lx = lg(x), i;
  if (lx != lg(y)) pari_err_TYPE("The inputs must have the same length.", mkvec2(x, y));
  for (i = 1; i < lx; i++) {
    xysum = gadd(xysum, gmul(gel(x, i), gel(y, i)));
    xsqrsum = gadd(xsqrsum, gsqr(gel(x, i)));
  }
  GEN fit = gdiv(xysum, xsqrsum);
  if (!retrsqr) return gerepileuptoleaf(av, fit);
  GEN M = cgetg(lx, t_MAT);
  for (i = 1; i < lx; i++) gel(M, i) = mkcol2(gen_1, gel(x, i));
  GEN rsqr = rsquared(M, y, mkcol2(gen_0, fit));
  return gerepilecopy(av, mkvec2(fit, rsqr));
}

/*OLS, where there is only one input variable. This just puts it into a matrix form and calls OLS, and is included for convenience.*/
GEN
OLS_single(GEN x, GEN y, int retrsqr)
{
  pari_sp av = avma;
  long lgx = lg(x), i;
  GEN xmat = cgetg(lgx, t_MAT);
  for (i = 1; i < lgx; i++) gel(xmat, i) = mkcol2(gen_1, gel(x, i));
  return gerepileupto(av, OLS(xmat, y, retrsqr));
}

/*Given inputs for OLS and the proposed linear fit, this returns the R^2 value of the regression.*/
GEN
rsquared(GEN X, GEN y, GEN fit)
{
  pari_sp av = avma;
  long n = lg(y) - 1, i;/*Number of observations*/
  GEN predicted = RgV_RgM_mul(shallowtrans(fit), X);/*1xn matrix of the fitted values.*/
  GEN yavg = gen_0;
  for (i = 1; i <= n; i++) yavg = gadd(yavg, gel(y, i));
  yavg = gdivgs(yavg, n);/*Average value of y*/
  GEN sstot = gen_0;
  GEN ssres = gen_0;
  for (i = 1; i <= n; i++) {
    sstot = gadd(sstot, gsqr(gsub(gel(y, i), yavg)));
    ssres = gadd(ssres, gsqr(gsub(gel(y, i), gel(predicted, i))));
  }
  return gerepileupto(av, gsubsg(1, gdiv(ssres, sstot)));
}

