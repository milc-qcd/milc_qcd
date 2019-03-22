#include <unistd.h>
#include "Wtimer.h"

void Wtimer_start ( Wtimer_t *tmr )
{
  gettimeofday ( &(tmr -> start), NULL );
}

double Wtimer_read ( Wtimer_t *tmr )
{
  long secs  = (tmr -> stop.tv_sec)  - (tmr -> start.tv_sec);
  long usecs = (tmr -> stop.tv_usec) - (tmr -> start.tv_usec);

  return ( (double) secs + 1.0e-6 * (double) usecs );
}

double Wtimer_stop ( Wtimer_t *tmr )
{
  gettimeofday ( &(tmr -> stop), NULL );
  return ( Wtimer_read ( tmr ) );
}

