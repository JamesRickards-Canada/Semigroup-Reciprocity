/*
  Parallel implementation of contfrac_tail_missing in C.
  
  Compile with: cc cfracsearch.c -o cfracsearch -O3 -pthread
      Run with: ./cfracsearch Bmin Bmax Nthreads
*/

#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <pthread.h>
#include <time.h>

typedef struct _cfrac_t {/*[0;4a_n,...,4a_1,k,1,2]. The threads start with either [4a_1, k, 1, 2] or [k, 1, 2]*/
  int myid;
  long Bmin;
  long Bmax;
  int Nthreads;
  unsigned long *found;/*The found denominators. We actually do this locally to avoid mutexing a ton.*/
  /*We want to break up the first few values of k since they really dominate (we don't want 1 thread still running while the others are all done). Similarly, after a certain point, we should just assign each thread to look at that class modulo Nthreads, since the continual mutexing is a lot of wastage. */
  long Ninitial;/*The number of precomputed numerator and denominator pairs to run through.*/
  unsigned long *den0;/*The vector of initial den[0]'s*/
  unsigned long *den1;/*The vector of initial den[1]'s*/
  long *initialdentodo;/*The smallest den[i] task not yet done, -1 once done.*/
  long *kdentodo;/*The smallest value of (regular) k to do.*/
  long stoptodo;/*After this, we just run though k in an arithmetic progression modulo Nthreads.*/
} cfrac_t;

pthread_mutex_t mutex_cfrac;

static void missing_tofile(unsigned long *found, long nblocks, long Bmin, long Bmax);
static void *cfrac_par(void *args);


/*Updates found with N*/
inline void
found_update(unsigned long *found, unsigned long *bitswap, long Bmin, long N)
{
  long shifted = N - Bmin;
  if (shifted < 0) return;
  long v = shifted % 64;
  long u = shifted >> 6;/*shifted = 64*u+v. u gives the entry of the array, v gives the bit to swap.*/
  found[u] |= bitswap[v];
}

/*Computes all missing denominators between Bmin and Bmax (inclusive), storing them to a file.*/
int
main(int argc, char *argv[])
{
  struct timespec start, finish;/*For measuring the real time.*/
  clock_gettime(CLOCK_MONOTONIC, &start);
  if (argc < 4) return 1;/*Not enough inputs.*/
  long Bmin = atol(argv[1]);
  long Bmax = atol(argv[2]);
  int Nthreads = atoi(argv[3]);
  long Bdiff = Bmax - Bmin + 1, i;
  long nblocks = ((Bdiff - 1) >> 6) + 1;
  unsigned long *found = (unsigned long *)calloc(nblocks, sizeof(unsigned long));/*Stores the found denominators*/
  if (!found) {printf("Insufficient memory to store the denominators found.\n"); exit(1); }
  unsigned long *bitswap = (unsigned long *)malloc(sizeof(unsigned long) << 6);/*Used for swapping bits of unsigned longs. ASSUME 64 BIT SYSTEM.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long Ninitial = (Nthreads + 5);
  Ninitial = Ninitial * Ninitial ;/*Go Nthreads + 5 deep on the first 2.*/
  unsigned long *den0 = (unsigned long *)malloc(Ninitial * sizeof(unsigned long));
  unsigned long *den1 = (unsigned long *)malloc(Ninitial * sizeof(unsigned long));
  long k, a1, a2, ind = 0;
  if (3 <= Bmax) found_update(found, bitswap, Bmin, 3);/*Denominator of 3*/
  for (k = 1; k <= Nthreads + 5; k++) {/*Let's go a two levels deep for the first few entries.*/
    long d = 3 * k + 2;/*The denominator*/
    if (d <= Bmax) found_update(found, bitswap, Bmin, d);/*We add the first 1 only, no need for the second since that is done in the thread.*/
    for (a1 = 1; a1 <= Nthreads + 5; a1++) {
      long e = (a1 << 2) * d + 3;
      den0[ind] = d;/*Add 'em in!*/
      den1[ind] = e;
      ind++;
    }
  }
  free(bitswap);
  pthread_t thread_id[Nthreads];/*The thread ids*/
  cfrac_t data[Nthreads];
  long initialdentodo = 0, kdentodo = Nthreads + 6;/*We did up to Nthreads + 5 in the first three levels.*/
  long stoptodo = (Bmax >> 20) + Nthreads + 100;/*Don't think this matters too much as long as it's moderately small.*/
  for (i = 0; i < Nthreads; i++) {/*Initialize the structures holding our data.*/
    data[i].myid = i;
    data[i].Bmin = Bmin;
    data[i].Bmax = Bmax;
    data[i].Nthreads = Nthreads;
    data[i].found = found;
    data[i].Ninitial = Ninitial;
    data[i].den0 = den0;
    data[i].den1 = den1;
    data[i].initialdentodo = &initialdentodo;
    data[i].kdentodo = &kdentodo;
    data[i].stoptodo = stoptodo;
  }
  pthread_mutex_init(&mutex_cfrac, NULL);/*Make the mutex*/
  for (i = 0; i < Nthreads; i++) pthread_create(&thread_id[i], NULL, cfrac_par, (void *)&data[i]);/*Make the threads.*/
  for (i = 0; i < Nthreads; i++) pthread_join(thread_id[i], NULL);/*Wait for them to all finish.*/
  pthread_mutex_destroy(&mutex_cfrac);/*Eliminate the mutex*/
  free(den1);
  free(den0);
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Total computation time: %fs\n", elapsed);
  missing_tofile(found, nblocks, Bmin, Bmax);/*Save the data to file.*/
  free(found);/*Free the memory.*/
  return 0;
}

/*Prints the found data to a file.*/
static void
missing_tofile(unsigned long *found, long nblocks, long Bmin, long Bmax)
{
  char fname[200];
  int pos = 0;
  DIR* dir = opendir("missing");
  if (dir) {
    closedir(dir);/* Directory exists. */
    pos += sprintf(&fname[pos], "missing/");
  }
  else if (ENOENT == errno) {/* Directory does not exist. */
    mkdir("missing", 0777);
    pos += sprintf(&fname[pos], "missing/");
  }
  else {/* opendir() failed for some other reason. */
    printf("Could not create directory, saving to current folder");
  }
  pos += sprintf(&fname[pos], "contfrac_tail_par_%ld-to-%ld.dat", Bmin, Bmax);
  FILE *F;
  F = fopen(fname, "w");
  long u;
  for (u = 0; u < nblocks - 1; u++) {/*Figure out each block.*/
    unsigned long val = found[u];
    long v;
    for (v = 0; v < 64; v++) {
      if (!(val & 1)) {/*A missing value*/
        long n = (u << 6) + v + Bmin;/*The value.*/
        fprintf(F, "%ld\n", n);/*Print it to the file.*/
      }
      val >>= 1;/*Shift right one.*/
    }
  }
  unsigned long val = found[nblocks - 1];/*Final block: we need to make sure we don't go back Bmax here.*/
  long v;
  for (v = 0; v < 64; v++) {
    if (!(val & 1)) {/*A missing value*/
      long n = (u << 6) + v + Bmin;/*The value.*/
      if (n <= Bmax) fprintf(F, "%ld\n", n);/*Print it to the file.*/
      else break;/*All done!*/
    }
    val >>= 1;/*Shift right one.*/
  }
  fclose(F);
}

/*Runs the thread, where we search for the missing denominators of the form [0;4a_n,...,4a_1,k,1,2]. Each thread tackles a fixed k.*/
static void *
cfrac_par(void *args)
{
  cfrac_t *data = (cfrac_t *)args;/*Import the given data.*/
  long Bmin = data -> Bmin, Bmax = data -> Bmax;
  int Nthreads = data -> Nthreads, Ninitial = data -> Ninitial;
  unsigned long *bitswap = (unsigned long *)malloc(sizeof(unsigned long) << 6), i;/*Used for swapping bits of unsigned longs. ASSUME 64 BIT SYSTEM.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long nblocks = ((Bmax - Bmin) >> 6) + 1;
  unsigned long *localfound = (unsigned long *)calloc(nblocks, sizeof(unsigned long));/*Local version of found to avoid mutexing a lot.*/
  if (!localfound) {printf("Insufficient memory to store the denominators (in the thread) found.\n"); exit(1); }
  long maxdepth = 100;/*Maximal depth, to start.*/
  unsigned long *dens = (unsigned long *)malloc(maxdepth * sizeof(unsigned long));/*Tracks the sequence of denominators.*/
  if (!dens) { printf("Insufficient memory to store the denominator sequence.\n"); exit(1); }
  for (;;) {/*Where we execute going through den[0] and den[1]*/
    pthread_mutex_lock(&mutex_cfrac);/*Retrieve the next denominators*/
    int cont = 0;/*1 if we are at the last value through the initial k's, and should continue on with the rest of the a_1's.*/
    long *next = data -> initialdentodo;
    if (*next == -1) {/*Done the initial part.*/
      next = data -> kdentodo;
      if (*next > (data -> stoptodo)) break;/*Moving on to the arithmetic sequences of k's.*/
      dens[0] = 3;
      dens[1] = (*next) * 3 + 2;/*Set the initial two denominators.*/
      (*next)++;/*Increment.*/
    }
    else {/*In the initial part.*/
      dens[0] = (data -> den0)[*next];
      dens[1] = (data -> den1)[*next];
      (*next)++;/*Increment.*/
      if (*next == Ninitial) { *next = -1; cont = 1; }/*Done! Also, continue on the next time.*/
      else if ((data -> den0)[*next] > (data -> den0)[*next - 1]) cont = 1;/*Last one for this value of k, need to continue on.*/
    }
    pthread_mutex_unlock(&mutex_cfrac);/*Unlock it.*/
    if (dens[1] > Bmax) continue;/*Nothing to do here.*/
    found_update(localfound, bitswap, Bmin, dens[1]);/*Add the first denominator.*/
    long ind = 2;
    dens[2] = dens[0];/*Initially set to this.*/ 
    long baseind = cont? 0 : 1;/*If cont = 1, continue to change index 1. Else, only do index 2 and onward.*/
    while (ind > baseind) {
      unsigned long newden = dens[ind] + (dens[ind - 1] << 2);/*d_n=4a_n*d_{n-1}+d_{n-2}. We just need to add 4*d_{n-1} every time.*/
      if (newden > Bmax) {/*Too big! Go back.*/
        ind--;/*Back up.*/
        continue;
      }
      found_update(localfound, bitswap, Bmin, newden);/*Update the saved missing denominators.*/
      dens[ind] = newden;/*Update dens.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must reallocate the storage location.*/
        maxdepth <<= 1;/*Double it.*/
        dens = realloc(dens, maxdepth * sizeof(unsigned long));
        if (!dens) {
          printf("Insufficient memory to reallocate the denominator sequence.\n");
          exit(1);
        }
      }
      dens[ind] = dens[ind - 2];/*Initialize for the next loop.*/
    }
  }
  long *next = data -> kdentodo;/*Mutex still locked.*/
  dens[1] = (*next) * 3 + 2;/*Second denominator*/
  (*next)++;/*Increment*/
  pthread_mutex_unlock(&mutex_cfrac);/*Unlock it.*/
  dens[0] = 3;
  unsigned long threeNthreads = 3 * Nthreads;
  while (dens[1] <= Bmax) {/*Where we execute the search in arithmetic progressions of k modulo Nthreads.*/
    found_update(localfound, bitswap, Bmin, dens[1]);/*Add the first denominator.*/
    long ind = 2;
    dens[2] = dens[0];/*Initially set to this.*/ 
    while (ind > 1) {
      unsigned long newden = dens[ind] + (dens[ind - 1] << 2);/*d_n=4a_n*d_{n-1}+d_{n-2}. We just need to add 4*d_{n-1} every time.*/
      if (newden > Bmax) {/*Too big! Go back.*/
        ind--;/*Back up.*/
        continue;
      }
      found_update(localfound, bitswap, Bmin, newden);/*Update the saved missing denominators.*/
      dens[ind] = newden;/*Update dens.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must reallocate the storage location.*/
        maxdepth <<= 1;/*Double it.*/
        dens = realloc(dens, maxdepth * sizeof(unsigned long));
        if (!dens) {
          printf("Insufficient memory to reallocate the denominator sequence.\n");
          exit(1);
        }
      }
      dens[ind] = dens[ind - 2];/*Initialize for the next loop.*/
    }
    dens[1] += threeNthreads;/*Increment along the sequence.*/
  }
  pthread_mutex_lock(&mutex_cfrac);/*Update the counts*/
  for (i = 0; i < nblocks; i++) (data -> found)[i] |= localfound[i];
  pthread_mutex_unlock(&mutex_cfrac);
  free(dens);
  free(localfound);
  free(bitswap);
  return NULL;
}

