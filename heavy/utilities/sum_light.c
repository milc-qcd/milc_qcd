/******* sum.c - summing up the degenerate light mesons */

/* MIMD version 6 */




#include "w_sum_includes.h"


void sum_light(int file_nhop, FILE * fp_m_in, int fb_m_in,
	         FILE * fp_k_out)
{
  register int t, channel;
  int N_iter, spin, color;
  double **meson_prop;
  double **meson_buf;
  double dtime ;



  if (mynode() == 0)
  {				/* only node 0 does any work */
    dtime = -dclock();



    /* allocate memory for meson props & buffer */
    meson_prop = (double **) malloc(nchannels * sizeof(double *));
    meson_buf = (double **) malloc(nchannels * sizeof(double *));

    for (channel = 0; channel < nchannels; channel++)
    {
      meson_prop[channel] = (double *) malloc(nt * sizeof(double));
      meson_buf[channel] = (double *) malloc(nt * sizeof(double));
    }

    /* clear meson_prop */
    for (channel = 0; channel < nchannels; channel++)
      for (t = 0; t < nt; t++)
      {
	meson_prop[channel][t] = 0.;
      }

    /* read in light-light results and sum over spin/color */
    for (spin = 0; spin < 4; spin++)
      for (color = 0; color < 3; color++)
      {
	for (N_iter = -1; N_iter < file_nhop; N_iter++)
	{
	  /*
	   * first read (N_iter == -1) just reads degenerate light meson;
	   * the rest is just ignored here 
	   */
	  if (startflag_m == RELOAD_ASCII)
	  {
	    r_ascii_m(fp_m_in, spin, color, N_iter, meson_buf);
	  }
	  if (startflag_m == RELOAD_BINARY)
	  {
	    r_binary_m(fb_m_in, spin, color, N_iter, meson_buf);
	  }
	  if (N_iter == -1)
	  {
	    for (channel = 0; channel < nchannels; channel++)
	      for (t = 0; t < nt; t++)
	      {
		meson_prop[channel][t] += meson_buf[channel][t];
	      }

	  }			/* N_iter == -1 */
	}			/* N_iter */
      }				/* spin color    --- now hopping results are
				 * all read in  */







    /* write out sum ; ascii is only option at this point */
    if (saveflag_k == SAVE_ASCII)
    {
      w_ascii_k(fp_k_out, -1, meson_prop);
    } else
    {
      printf(" Do not know how to write out results\n");
      terminate(1);
    }





    for (channel = 0; channel < nchannels; channel++)
    {
      free(meson_prop[channel]);
      free(meson_buf[channel]);
    }
    free(meson_prop);
    free(meson_buf);

    dtime += dclock();
    printf("SUM: time = %e\n", dtime);


  }				/* node 0 */
  return;
}				/* sum_light  */
