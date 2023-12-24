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

typedef struct _cfrac_t {
  int myid;
  long Bmin;
  long Bmax;
  unsigned long *found;/*The found denominators. We actually do this locally to avoid mutexing a ton.*/
  long stop;/*Largest value of a to try.*/
  long *todo;/*The smallest task not yet done, -1 if done.*/
} cfrac_t;

pthread_mutex_t mutex_cfrac;

static void missing_tofile(unsigned long *found, long nblocks, long Bmin, long Bmax);
static void *cfrac_par(void *args);


/*Computes all missing denominators between Bmin and Bmax (inclusive), storing them to a file.*/
int
main(int argc, char *argv[])
{
  if (argc < 4) return 1;/*Not enough inputs.*/
  long Bmin = atol(argv[1]);
  long Bmax = atol(argv[2]);
  int Nthreads = atoi(argv[3]);
  long Bdiff = Bmax - Bmin + 1, i;
  long nblocks = ((Bdiff - 1) >> 6) + 1;
  unsigned long *found = (unsigned long *)calloc(nblocks, sizeof(unsigned long));/*Stores the found denominators*/
  long stop = (Bmax + 1) / 3;/*One bigger than the largest value of a to consider in [a, 1, 2].*/
  pthread_t thread_id[Nthreads];/*The thread ids*/
  long todo = 1;
  cfrac_t data[Nthreads];
  for (i = 0; i < Nthreads; i++) {/*Initialize the structures holding our data.*/
    data[i].myid = i;
    data[i].Bmin = Bmin;
    data[i].Bmax = Bmax;
    data[i].found = found;
    data[i].stop = stop;
    data[i].todo = &todo;
  }
  pthread_mutex_init(&mutex_cfrac, NULL);/*Make the mutex*/
  for (i = 0; i < Nthreads; i++) pthread_create(&thread_id[i], NULL, cfrac_par, (void *)&data[i]);/*Make the threads.*/
  for (i = 0; i < Nthreads; i++) pthread_join(thread_id[i], NULL);/*Wait for them to all finish.*/
  pthread_mutex_destroy(&mutex_cfrac);/*Eliminate the mutex*/
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

/*Runs the thread, where we search for the missing denominators of the form [0;4a_n,...,4a_1,k,1,2]. Each thread tackles a fixed k.*/
static void *
cfrac_par(void *args)
{
  cfrac_t *data = (cfrac_t *)args;/*Import the given data.*/
  long Bmin = data -> Bmin, Bmax = data -> Bmax;
  unsigned long *bitswap = (unsigned long *)malloc(sizeof(unsigned long) << 6), i;/*Used for swapping bits of unsigned longs. ASSUME 64 BIT SYSTEM.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long nblocks = ((Bmax - Bmin) >> 6) + 1;
  unsigned long *localfound = (unsigned long *)calloc(nblocks, sizeof(unsigned long));/*Local version of found to avoid mutexing a lot.*/
  if (data -> myid == 0) {/*First thread will add in the first denominator of 3*/
    found_update(localfound, bitswap, Bmin, 3);
  }
  long maxdepth = 100;/*Maximal depth, to start.*/
  unsigned long *dens = (unsigned long *)malloc(maxdepth * sizeof(unsigned long));/*Tracks the sequence of denominators.*/
  dens[0] = 3;
  if (!dens) { printf("Insufficient memory to store the denominator sequence.\n"); exit(1); }
  for (;;) {/*Where we execute the threads.*/
    pthread_mutex_lock(&mutex_cfrac);/*Retrieve the next value of k.*/
    long *todoptr = data -> todo;
    if (*todoptr == -1) { pthread_mutex_unlock(&mutex_cfrac); break; }/*All tasks done or running. Compile data and exit. Don't forget to unlock!!*/
    dens[1] = (*todoptr) * 3 + 2;/*Second denominator.*/
    (*todoptr)++;/*Increment*/
    if (*todoptr == data -> stop) *todoptr = -1;/*Last one.*/
    pthread_mutex_unlock(&mutex_cfrac);/*Unlock it.*/
    /*Time to execute the depth first search.*/
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
  }
  pthread_mutex_lock(&mutex_cfrac);/*Update the counts*/
  for (i = 0; i < nblocks; i++) (data -> found)[i] |= localfound[i];
  pthread_mutex_unlock(&mutex_cfrac);
  free(dens);
  free(localfound);
  free(bitswap);
  return NULL;
}

