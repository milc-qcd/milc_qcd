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

int main(int argc, char **argv)
{
  su3_vector a, b1, b2, b3, b4, c;
  int i,j,iter=1;
  struct sched_param param={sched_priority:20};
  volatile unsigned long long timeA, timeB, timeMILC, timeSSE;
  unsigned int seed=1;
  int randomfd;
  
  if ((randomfd=open("/dev/urandom", O_RDONLY)) < 0)
    perror("Attempt to open /dev/urandom");
  if (read(randomfd, &seed, sizeof(seed)) < sizeof(seed))
    perror("Attempt to read /dev/urandom");
  close(randomfd);
  //  srand(seed);

  if (sched_setscheduler(0, SCHED_FIFO, &param) < 0)
    perror("Attempt to put in real time queue");

  if (argc > 1) sscanf(argv[1],"%d",&iter);

  for (i=0; i<3; i++) {
    a.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    a.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b1.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b1.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b2.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b2.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b3.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b3.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b4.c[i].real = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
    b4.c[i].imag = (Real)(rand() - RAND_MAX/2)/(Real)(RAND_MAX/2)*2.0;
  }
  c = a;

  rdtscll(timeA);
  for (i=0; i<iter; i++) {
    sub_four_su3_vecs(&a, &b1, &b2, &b3, &b4);
  }
  rdtscll(timeMILC);
  timeMILC -= timeA;

  rdtscll(timeB);
  for (i=0; i<iter; i++) {
    _inline_sse_sub_four_su3_vecs(&c, &b1, &b2, &b3, &b4);
  }
  rdtscll(timeSSE);
  timeSSE -= timeB;

  for (i=0; i<3; i++) {
    printf("%4.1f%+4.1fi    %4.1f%+4.1fi\n",
           a.c[i].real, a.c[i].imag, c.c[i].real, c.c[i].imag);
  }

  printf("Time per iteration:\n  MILC: %Lu\n  SSE:  %Lu\n", timeMILC/(unsigned long long)iter,
	 timeSSE/(unsigned long long)iter);

  exit(0);
}


    
