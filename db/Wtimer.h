#ifndef _Wtimer_h_
#define _Wtimer_h_

#include <time.h>
#include <sys/time.h>

typedef struct {
  struct timeval start;
  struct timeval stop;
} Wtimer_t;

void   Wtimer_start ( Wtimer_t* );
double Wtimer_stop  ( Wtimer_t* );
double Wtimer_read  ( Wtimer_t* );

#endif
