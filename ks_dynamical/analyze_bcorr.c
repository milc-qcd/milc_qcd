/* Project-specific measurements */
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */

/* All routines are called from control.c(main)
   
   setup_analyze: called immediately after calling setup()
   init_analyze:  called after each call to readin()
   analyze:       called at each measurement period
   end_analyze:   called at the end before returning to readin()

   */
#include "ks_dyn_includes.h"

int spect_iters,avspect_iters, avbcorr_iters;
complex plp_fuzzy;
int bar_corr();

void setup_analyze()
{
  if(this_node==0)printf("Fat Polyakov loop parameter %f\n",ALPHA_FUZZ);
  if(this_node==0)printf("bar_corr NBPRAND = %d\n",NBPRAND);
}

void init_analyze()
{
  avspect_iters = avbcorr_iters = 0;
}

void analyze(int meascount)
{
  plp_fuzzy = ploop_staple((Real)ALPHA_FUZZ);
  plp_fuzzy.real *= -1.0; plp_fuzzy.imag *= -1.0; /* KS phases! */
  
  /* call for various baryon density correlations */
  /* WARNING: This call must follow ploop and ploop_staple, 
     since it is assumed that the Polyakov loop link matrix 
     products are in the ploop and ploop_fuzz elements 
     of the site structure on even sites in the first two 
     time slices */
  avbcorr_iters += bar_corr();
  
  if(this_node==0)printf("GFUZZ %e %e\n",
			 (double)plp_fuzzy.real,(double)plp_fuzzy.imag);
}

void end_analyze(int meascount)
{
	    if(this_node==0)printf("average cg iters for spectrum = %e\n",
		(double)avbcorr_iters/meascount);
}


