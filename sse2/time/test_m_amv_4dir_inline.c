/* time various approaches to adjoint su3_matrix * su3_vector */
#include "complex.h"
#include "su3.h"
#include "inline_sse.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <error.h>
#include <asm/msr.h>
#include <fcntl.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  su3_matrix a[4];
  su3_vector b, c[4], d[4];
  int i,j,k,iter=1;
  struct sched_param param={sched_priority:20};
  volatile unsigned long long timeA, timeMILC, timeSSE;
  unsigned int seed=1;
  int randomfd;
  
  if ((randomfd=open("/dev/urandom", O_RDONLY)) < 0)
    perror("Attempt to open /dev/urandom");
  if (read(randomfd, &seed, sizeof(seed)) < sizeof(seed))
    perror("Attempt to read /dev/urandom");
  close(randomfd);
  srand(seed);

  if (sched_setscheduler(0, SCHED_FIFO, &param) < 0)
    perror("Attempt to put in real time queue");

  if (argc > 1) sscanf(argv[1],"%d",&iter);

  for (k=0; k<4; k++) {
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	a[k].e[i][j].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
	a[k].e[i][j].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      }
    }
  }
  for (i=0; i<3; i++) {
    b.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
  }

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    mult_adj_su3_mat_vec_4dir(a, &b, c);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    _inline_sse_mult_adj_su3_mat_vec_4dir(a, &b, d);
  }
  rdtscll(timeSSE);
  timeSSE -= timeA;

  for (i=0; i<3; i++) {
      printf("%4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1f\n",
	   c[0].c[i].real, c[0].c[i].imag, d[0].c[i].real, d[0].c[i].imag,
	   c[1].c[i].real, c[1].c[i].imag, d[1].c[i].real, d[1].c[i].imag,
	   c[2].c[i].real, c[2].c[i].imag, d[2].c[i].real, d[2].c[i].imag,
	   c[3].c[i].real, c[3].c[i].imag, d[3].c[i].real, d[3].c[i].imag);
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}


    

