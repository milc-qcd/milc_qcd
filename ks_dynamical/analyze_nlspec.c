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
int mom_spec(); /* return the C.G. iteration number */

int spect_iters,avspect_iters, avbcorr_iters;

void setup_analyze()
{
  register int i;
  node0_printf("With spectrum measurements\n");
  node0_printf("With nonlocal spectrum measurements\n");
}

void init_analyze()
{
  avspect_iters = 0;
}

void analyze(int meascount)
{
  spect_iters = spectrum();
  spect_iters += mom_spec();
  source_inc = nt/2;
  n_sources = 2;
  source_start = 0;
  /* Fix TUP Coulomb gauge - gauge links only*/
  rephase( OFF );
  gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  rephase( ON );
  spect_iters += nl_spectrum( mass, F_OFFSET(phi), F_OFFSET(xxx), 
			      F_OFFSET(tempmat1), F_OFFSET(tempmat2));
  avspect_iters += spect_iters;
}

void end_analyze(int meascount)
{
	    if(this_node==0)printf("average cg iters for spectrum = %e\n",
		(double)avspect_iters/meascount);
}
