/* Project-specific measurements */
/* MIMD version 6 */

/* All routines are called from control.c(main)
   
   setup_analyze: called immediately after calling setup()
   init_analyze:  called after each call to readin()
   analyze:       called at each measurement period
   end_analyze:   called at the end before returning to readin()

   */

#include "ks_dyn_includes.h"
int spectrum();

int spect_iters,avspect_iters, avbcorr_iters;

void setup_analyze()
{
  printf("With spectrum measurements\n");
}

void init_analyze()
{
  avspect_iters = 0;
}

void analyze(int meascount)
{
  spect_iters = spectrum();
  avspect_iters += spect_iters;
}

void end_analyze(int meascount)
{
  if(this_node==0)printf("average cg iters for spectrum = %e\n",
			 (double)avspect_iters/meascount);
}
