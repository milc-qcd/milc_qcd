/* time various approaches to adj su3_matrix * half-wilson vector */
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
  su3_matrix a;
  half_wilson_vector b, c, d;
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

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      a.e[i][j].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      a.e[i][j].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    }
    for (k=0; k<2; k++) {
      b.h[k].c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
      b.h[k].c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    }
  }

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    mult_adj_su3_mat_hwvec(&a, &b, &c);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    _inline_sse_mult_adj_su3_mat_hwvec(&a, &b, &d);
  }
  rdtscll(timeSSE);
  timeSSE -= timeA;

  for (i=0; i<3; i++) {
    printf("%4.1f%+4.1fi %4.1f%+4.1fi %4.1f%+4.1fi | %4.1f%+4.1fi | %4.1f%+4.1fi | %4.1f%+4.1fi | %4.1f%+4.1fi | %4.1f%+4.1fi | %4.1f%+4.1fi\n",
           a.e[i][0].real, a.e[i][0].imag, a.e[i][1].real, a.e[i][1].imag, a.e[i][2].real, a.e[i][2].imag,
	   b.h[0].c[i].real, b.h[0].c[i].imag, b.h[1].c[i].real, b.h[1].c[i].imag,
	   c.h[0].c[i].real, c.h[0].c[i].imag, c.h[1].c[i].real, c.h[1].c[i].imag,
	   d.h[0].c[i].real, d.h[0].c[i].imag, d.h[1].c[i].real, d.h[1].c[i].imag
	   );
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}
