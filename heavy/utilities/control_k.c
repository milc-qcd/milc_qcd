/*****************control_k.c *****************************************/
/* set tabstop=2   for easy reading of this file */

/*
 * Main procedure for quenched heavy-light SU3 Wilson fermions --- summing up
 * hopping param expansion for heavies 	 
 */
/* MIMD version 6 */

/* THIS IS SCALAR CODE.  RUN ON A SINGLE NODE ONLY! */

#define CONTROL

#include "w_sum_includes.h"	/* global variables for lattice fields */



int main(int argc, char **argv)
{
  int prompt ; 
  double starttime, endtime;
  FILE *fp_m_in, *fp_k_out;
  int fb_m_in;
  int max_prop, file_nhop, do_hop;

  /***----------------------------------------****/

  initialize_machine(argc, argv);

  if (mynode() == 0)
  {				/* only node 0 does any work */
    /* set up */
    prompt = setup_k();

    /* loop over input sets */
    while (readin(prompt) == 0)
    {

      starttime = dclock();

      max_prop = 12;

      /* check if we should do hopping exp. or just use light propagator */
      if (fabs(kappa_h - kappa) < EPS)
      {
	do_hop = 0;
      } 
      else
      {
	do_hop = 1;
	if (kappa_h > 1.0)
	  kappa_h = kappa;	/* kappa_h > 1.0 is override to do degenerate
				 * meson but with hopping expansion for one
				 * light quark */
      }

      /* open files */



      /* open file for meson hopping input */
      if (startflag_m == RELOAD_ASCII)
      {
	fp_m_in = r_ascii_m_i(startfile_m, max_prop, &file_nhop);
	fb_m_in = -1;		/* i.e. file is NOT binary */
      }
      else if (startflag_m == RELOAD_BINARY)
      {
	fb_m_in = r_binary_m_i(startfile_m, max_prop, &file_nhop);
	fp_m_in = NULL;		/* i.e. file is NOT ascii */
      }
      if (file_nhop < nhop)
      {
	printf("asking for more hopping iters than exist in file\n");
	printf("file_nhop= %d,  program nhop= %d\n", file_nhop, nhop);
	terminate(1);
      }
      /* open file for output of current kappa_h sum */

      if (saveflag_k == SAVE_ASCII)
	fp_k_out = w_ascii_k_i(savefile_k, max_prop);

      if (do_hop)
      {
	sum(kappa_h, kappa_c, nhop, file_nhop, fp_m_in, fb_m_in, fp_k_out, writeflag);
      } else
      {
	sum_light(file_nhop, fp_m_in, fb_m_in, fp_k_out);
      }




      /* close files */

      if (startflag_m == RELOAD_ASCII)
	r_ascii_m_f(fp_m_in, startfile_m);
      else if (startflag_m == RELOAD_BINARY)
	r_binary_m_f(fb_m_in, startfile_m);

      if (saveflag_k == SAVE_ASCII)
	w_ascii_k_f(fp_k_out, savefile_k);




      printf("RUNNING COMPLETED\n");

      endtime = dclock();
      printf("Time = %e seconds\n", (double) (endtime - starttime));
      fflush(stdout);

    }				/* while(prompt) */
  }				/* node 0 */
  return 0;
}				/* main() */
