/* time various approaches to su3_matrix + scalar * su3_matrix */
#include "complex.h"
#include "su3.h"
#include "../include/inline_sse.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <error.h>
#include <asm/msr.h>
#include <fcntl.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  su3_matrix a, b, c, d;
  Real s;
  int i,j,iter=1;
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
  
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      a.e[i][j].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      a.e[i][j].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      b.e[i][j].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      b.e[i][j].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    }
  }
  s = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    scalar_mult_add_su3_matrix(&a, &b, s, &c);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    _inline_sse_scalar_mult_add_su3_matrix(&a, &b, s, &d);
  }
  rdtscll(timeSSE);
  timeSSE -= timeA;

  printf("%4.1f\n", s);
  for (i=0; i<3; i++) {
    printf("%4.1f%+4.1fi %4.1f%+4.1fi %4.1f%+4.1fi | %4.1f%+4.1fi %4.1f%+4.1fi %4.1f%+4.1fi | %4.1f%+4.1fi %4.1f%+4.1fi %4.1f%+4.1fi | %4.1f%+4.1fi %4.1f%+4.1fi %4.1f%+4.1fi\n",
           a.e[i][0].real, a.e[i][0].imag, a.e[i][1].real, a.e[i][1].imag, a.e[i][2].real, a.e[i][2].imag,
           b.e[i][0].real, b.e[i][0].imag, b.e[i][1].real, b.e[i][1].imag, b.e[i][2].real, b.e[i][2].imag,
           c.e[i][0].real, c.e[i][0].imag, c.e[i][1].real, c.e[i][1].imag, c.e[i][2].real, c.e[i][2].imag,
           d.e[i][0].real, d.e[i][0].imag, d.e[i][1].real, d.e[i][1].imag, d.e[i][2].real, d.e[i][2].imag);
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}


    
