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
  su3_vector b, c0, c1, c2, c3, d0, d1, d2, d3;
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
    mult_adj_su3_mat_4vec(a, &b, &c0, &c1, &c2, &c3);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    _inline_sse_mult_adj_su3_mat_4vec(a, &b, &d0, &d1, &d2, &d3);
  }
  rdtscll(timeSSE);
  timeSSE -= timeA;

  for (i=0; i<3; i++) {
      printf("%4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1fi  %4.1f%+4.1f\n",
	   c0.c[i].real, c0.c[i].imag, d0.c[i].real, d0.c[i].imag,
	   c1.c[i].real, c1.c[i].imag, d1.c[i].real, d1.c[i].imag,
	   c2.c[i].real, c2.c[i].imag, d2.c[i].real, d2.c[i].imag,
	   c3.c[i].real, c3.c[i].imag, d3.c[i].real, d3.c[i].imag);
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}


    
