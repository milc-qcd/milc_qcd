/******* light_meson.c - find the contribution to the degenerate meson 
             propagators from the light quark with a given spin-color source */


/* MIMD version 7 */

/*
 * spin, color are just passed through to w_meson_hop to write in output file  
 */





#include "w_heavy_includes.h"




/* light_quark is an offset of a field of type wilson_vector */


void light_meson(field_offset light_quark, int color, int spin,
		   int wallflag, FILE * fp_m_out, int fb_m_out)
{
  double dtime;
  int channel, N_iter;
  double **meson_prop;
  wilson_vector *light_wall, *heavy_wall;


  /* Start Hopping */

  dtime = -dclock();

  N_iter = -1;			/* we treat this like the -1 iteration of the
				 * hopping expansion --- it goes before all
				 * the hopping results in the meson output
				 * file  */



  /* allocate memory for meson props */
  meson_prop = (double **) malloc(nchannels * sizeof(double *));
  for (channel = 0; channel < nchannels; channel++)
  {
    meson_prop[channel] = (double *) malloc(nt * sizeof(double));
  }

  /*
   * if wall sources, allocate memory for walls "heavy_wall" is needed,
   * despite the fact that the quarks are degenerate here, because
   * w_meson_hop is set up for the heavy-light case, which takes most of the
   * effort 
   */
  if (wallflag == CUTOFF_GAUSSIAN || wallflag == CUTOFF_GAUSSIAN_WEYL)
  {
    light_wall = (wilson_vector *) malloc(nt * sizeof(wilson_vector));
    heavy_wall = (wilson_vector *) malloc(nt * sizeof(wilson_vector));
  }
  /* now calculate and write out meson propagators */

  w_meson_hop(meson_prop, light_quark, light_quark, heavy_wall, light_wall,
	      EVENANDODD, spin, wallflag);

  if (saveflag_m == SAVE_MESON_ASCII)
  {
    w_ascii_m(fp_m_out, spin, color, N_iter, meson_prop);
  }
  if (saveflag_m == SAVE_MESON_BINARY)
  {
    w_binary_m(fb_m_out, spin, color, N_iter, meson_prop);
  }
  for (channel = 0; channel < nchannels; channel++)
  {
    free(meson_prop[channel]);
  }
  free(meson_prop);
  if (wallflag == CUTOFF_GAUSSIAN || wallflag == CUTOFF_GAUSSIAN_WEYL)
  {
    free(light_wall);
    free(heavy_wall);
  }
  dtime += dclock();
  if (this_node == 0)
    printf("LIGHT_MESON: time = %e = %e/site\n",
	   dtime, dtime / (volume));


  return;
}				/* light_meson  */
