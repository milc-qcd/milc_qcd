/* time various approaches to su3_matrix * su3_vector */
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

void sse_mult_su3_mat_vec(su3_matrix *, su3_vector *, su3_vector *);

int main(int argc, char **argv)
{
  su3_matrix a[4];
  su3_vector b[4], c, d = {{{0.0},{0.0}}};
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
    for (i=0; i<3; i++) {
      b[k].c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      b[k].c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    }
  }
  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    mult_su3_mat_vec_sum_4dir(a, &b[0], &b[1], &b[2], &b[3], &c);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    _inline_sse_mult_su3_mat_vec_sum_4dir(a, &b[0], &b[1], &b[2], &b[3], &d);
  }
  rdtscll(timeSSE);
  timeSSE -= timeA;

  for (i=0; i<3; i++) {
    printf("%4.1f%+4.1fi    %4.1f%+4.1fi\n",
	   c.c[i].real, c.c[i].imag, d.c[i].real, d.c[i].imag);
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}


    
