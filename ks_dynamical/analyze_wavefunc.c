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
int spectrum();
void setup_wavefunc_t();
void wavefunc_t();
void wf_pt_t();

int spect_iters,avspect_iters, avbcorr_iters;
complex plp_fuzzy;

void setup_analyze()
{
  node0_printf("With spectrum measurements\n");
  node0_printf("With wavefunction measurements\n");

  /* Set up gathers for wave function computations */
  setup_wavefunc_t();
}

void init_analyze()
{
  avspect_iters = 0;
}

void analyze(int meascount)
{
  spect_iters = spectrum();
  avspect_iters += spect_iters;
  wf_pt_t();
  wavefunc_t();
}

void end_analyze(int meascount)
{
	    if(this_node==0)printf("average cg iters for spectrum = %e\n",
		(double)avspect_iters/meascount);
}
