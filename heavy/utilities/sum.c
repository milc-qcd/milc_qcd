/******* sum.c - summing the hopping parameter expansion for heavy quarks */

/* MIMD version 6 */




#include "w_sum_includes.h"



void sum(Real kappa_h, Real kappa_c, int nhop, int file_nhop,
	   FILE * fp_m_in, int fb_m_in, FILE * fp_k_out, int writeflag)
{
  register int t, channel;
  int N_iter, spin, color;
  double **meson_hop;
  double **meson_prop;
  double **meson_buf;
  double kappa_power, kappa_ratio;
  int kappa_h_int;
  Real weight;
  double dtime ;



  if (mynode() == 0)
  {				/* only node 0 does any work */
    dtime = -dclock();


    /* criterion for whether we're to do static is kappa_h_int = 0 */
    kappa_h_int = (int) (1000 * kappa_h);

    if (kappa_h_int != 0)
    {				/* conventional */
      kappa_ratio = (double) (kappa_h / kappa_c);
    } else
    {				/* static */
      kappa_ratio = (double) (1.0 / (2.* kappa_c));
    }



    /* allocate memory for meson props & buffer */
    meson_hop = (double **) malloc(nchannels * sizeof(double *));
    meson_prop = (double **) malloc(nchannels * sizeof(double *));
    meson_buf = (double **) malloc(nchannels * sizeof(double *));
    for (channel = 0; channel < nchannels; channel++)
    {
      meson_hop[channel] = (double *) malloc(nhop * nt * sizeof(double));
      meson_prop[channel] = (double *) malloc(nt * sizeof(double));
      meson_buf[channel] = (double *) malloc(nt * sizeof(double));
    }

    /* clear meson_hop storage */
    for (channel = 0; channel < nchannels; channel++)
      for (t = 0; t < nt; t++)
      {
	for (N_iter = 0; N_iter < nhop; N_iter++)
	{
	  meson_hop[channel][t * nhop + N_iter] = 0.;
	}
      }

    /* read in hopping results and sum over spin/color only */
    /*
     * need to store things temporarily in meson_hop, since hopping results *
     * are stored with spin-color the slowest changes indices and N_iter the
     * fastest, * but we want to sum over spin-color before we sum over
     * iteration 
     */
    for (spin = 0; spin < 4; spin++)
      for (color = 0; color < 3; color++)
      {
	/******printf("%d %d\n",spin,color);******/
	for (N_iter = -1; N_iter < file_nhop; N_iter++)
	{
	  /*
	   * first read (N_iter=-1) just reads degenerate light meson; 
	   * result is ignored here 
	   */
	  if (startflag_m == RELOAD_ASCII)
	  {
	    r_ascii_m(fp_m_in, spin, color, N_iter, meson_buf);
	  }
	  if (startflag_m == RELOAD_BINARY)
	  {
	    r_binary_m(fb_m_in, spin, color, N_iter, meson_buf);
	  }
	  /* write out coef for this spin color */
	  /* write out coeffs to output */
	  /******t=40;
	  channel=0;
	  printf("%d   %.7E\n",N_iter,meson_buf[channel][t] );  ******/


	  if (N_iter > -1 && N_iter < nhop)
	  {			/* only include wanted iterations */
	    if (kappa_h_int != 0)
	    {			/* conventional case */
	      for (channel = 0; channel < nchannels; channel++)
		for (t = 0; t < nt; t++)
		{
		  meson_hop[channel][t * nhop + N_iter] += meson_buf[channel][t];
		}
	    } else
	    {			/* static case */
	      if (N_iter <= nt / 2)
	      {			/* halfway around the lattice for forward or
				 * backward movers */
		weight = 1.0;
		if (N_iter == 0 || (N_iter == nt / 2 && nt % 2 == 0))
		  weight = 0.5;
		/*
		 * average forward & backward movers at source slice and
		 * halfway across lattice (if nt even)  
		 */
		if (spin >= 2)
		{		/* these are the forward(!) moving spins */
		  t = source_t + N_iter;
		  if (t >= nt)
		    t -= nt;
		  for (channel = 0; channel < nchannels; channel++)
		  {
		    meson_hop[channel][t * nhop + N_iter] += weight * meson_buf[channel][t];
		  }
		}
		 /* forward spins */ 
		else
		{
		  /* these are the backward moving spins  */
		  t = source_t - N_iter;
		  if (t < 0)
		    t += nt;
		  for (channel = 0; channel < nchannels; channel++)
		  {
		    meson_hop[channel][t * nhop + N_iter] += weight * meson_buf[channel][t];
		  }
		}		/* backward spins */
	      }			/* N_iter < nt */
	    }			/* static case */

	  }			/* N_iter > -1 && N_iter < nhop */
	}			/* N_iter */
      }				/* spin color    --- now hopping results are
				 * all read in  */



    /* write out coeffs to output */
    /******    t=40;
        channel=0;
        for(N_iter=0;N_iter<nhop;N_iter++){
           printf("%d   %.7E\n",N_iter,meson_hop[channel][t*nhop + N_iter] );
        } ******/
    /* clear meson_prop */
    for (channel = 0; channel < nchannels; channel++)
      for (t = 0; t < nt; t++)
      {
	meson_prop[channel][t] = 0.;
      }

    kappa_power = 1.0;

    for (N_iter = 0; N_iter < nhop; N_iter++)
    {

      /* multiply by kappa_ratio and add on to meson_prop */
      for (channel = 0; channel < nchannels; channel++)
	for (t = 0; t < nt; t++)
	{
	  meson_prop[channel][t] += kappa_power * meson_hop[channel][t * nhop + N_iter];
	}

      /* write out current sum if desired; ascii is only option at this point */
      if (writeflag == WRITEALL || N_iter == nhop - 1)
      {
	if (saveflag_k == SAVE_ASCII)
	{
	  w_ascii_k(fp_k_out, N_iter, meson_prop);
	} else
	{
	  printf(" Don't know how to write out results\n");
	  terminate(1);
	}
      }
      kappa_power *= kappa_ratio;

    }				/* end of loop over N_iter  */


    for (channel = 0; channel < nchannels; channel++)
    {
      free(meson_hop[channel]);
      free(meson_prop[channel]);
      free(meson_buf[channel]);
    }
    free(meson_hop);
    free(meson_prop);
    free(meson_buf);


    dtime += dclock();
    printf("SUM: time = %e\n", dtime);


  }				/* node 0 */
  return;
}				/* sum  */
