/******* d_congrad5_fn_1sum_qdp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */


/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"


#if 0

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

void
print_mem(void)
{
  if(QDP_this_node==0) {
    int pid;
    char s[256];

    pid = getpid();
    //printf("vsize = ");
    fflush(stdout);
    snprintf(s, 256, "/usr/bin/ps -o vsz -p %i", pid);
    system(s);
  }
}
#endif

/* Standard MILC interface for the Asqtad inverter */

int
ks_congrad(field_offset f_src, field_offset f_dest, Real mass,
	   int niter, int nrestart, Real rsqmin, int prec, int parity, 
           Real *final_rsq_ptr, ferm_links_t *fn)
{

  int iteration;

  if(prec == 1)
    iteration = 
      ks_congrad_milc2qdp_F(f_src, f_dest, mass, niter, nrestart, 
			    rsqmin, parity, final_rsq_ptr, fn);
  else
    iteration = 
      ks_congrad_milc2qdp_D(f_src, f_dest, mass, niter, nrestart, 
			    rsqmin, parity, final_rsq_ptr, fn);
  
  //print_mem();
  return iteration;
}
