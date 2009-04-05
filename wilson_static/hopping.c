/******* hopping.c - hopping parameter expansion for heavy quarks */
/* MIMD version 7 */
/* (Identical to heavy/hopping.c) */
/*
 * Assumes EVENFIRST and same number of even and odd sites on each node. 
   3/29/00 EVENFIRST is the rule now. CD.
 */

/* MIMD version 7 */

/*
 * The source vector is in the wilson vector field with offset "src", which
 * is overwritten at each iteration by kappa*dslash times itself, "temp" is
 * the offset of a temp wilson vector needed by kappa_dslash; "light_quark"
 * is the offset of the wilson vector field of the light quarks (assumed
 * already computed) --- it is needed by w_meson_hop to get the meson props.
 * src is also the residual vector (of the previous iteration to be precise)
 * because of the simplicity of jacobi.  nhops is the number of iterations to
 * do; kappa_c is an approximate value for kappa_critical (put in to make the
 * hopping parameter expansion numerically stable).  parity_of_source is
 * needed because the routine exploits the even-odd decomposition and has to
 * know which set of sites to use on a given iteration. spin, color are just
 * passed through to w_meson_hop to write in output file  
 */





#include "w_static_includes.h"



/*
 * src, temp and light_quark are offsets of fields of type wilson_vector 
 */


void hopping(field_offset src, field_offset temp,
	       field_offset light_quark,
	       int nhop, Real kappa_c, int parity_of_source,
	   int color, int spin, int wallflag, FILE * fp_m_out, int fb_m_out)
{
  double dtime ;
/**  double dtime1;  ****/
  int N_iter;
  register int i;
  register site *s;
  Real size_src, size_r;
  int old_parity, new_parity = 0x00, channel;
  double **meson_prop;
  wilson_vector *light_wall = NULL, *heavy_wall = NULL;



  /* Start Hopping */

  dtime = -dclock();


  /* Normalisation  */
  size_src = 0.0;
  FORSOMEPARITY(i, s, parity_of_source)
  {
    size_src += magsq_wvec(((wilson_vector *) F_PT(s, src)));
  }
  g_floatsum(&size_src);
  size_src = (Real) sqrt((double) size_src);

  if (this_node == 0)
    printf("beginning hopping--size_src=%e\n",
	   (double) size_src);


  /* allocate memory for meson props */

  meson_prop = (double **) malloc(nchannels * sizeof(double *));
  for (channel = 0; channel < nchannels; channel++)
  {
    meson_prop[channel] = (double *) malloc(nt * sizeof(double));
    if (meson_prop[channel] == NULL)
    {
      printf("NODE %d: no room for meson_prop, channel= %d\n", this_node, channel);
      terminate(1);
    }
  }


  /* if wall sources, allocate memory for walls */
  if (wallflag == CUTOFF_GAUSSIAN || wallflag == CUTOFF_GAUSSIAN_WEYL)
  {
    light_wall = (wilson_vector *) malloc(nt * sizeof(wilson_vector));
    heavy_wall = (wilson_vector *) malloc(nt * sizeof(wilson_vector));
  }
  old_parity = parity_of_source;
  for (N_iter = 0; N_iter < nhop; N_iter++)
  {

    /* dtime1 = -dclock();  */
    if (N_iter == 0)
    {
      new_parity = old_parity;	/* do nothing on zeroth pass */
    } else
    {
      /* find what new parity will be after dslash is applied */
      switch (old_parity)
      {
      case EVEN:
	new_parity = ODD;
	break;
      case ODD:
	new_parity = EVEN;
	break;
      case EVENANDODD:
	new_parity = EVENANDODD;
	break;
      }
      /* replace src by kappa*dslash(src)  */
      kappa_dslash(src, temp, kappa_c, new_parity);
    }


    /*
     * find residual, magnitude of current src --- [really this is the
     * residual of the previous iteration]  
     */
    size_r = 0.0;

    FORSOMEPARITY(i, s, new_parity)
    {
      size_r += magsq_wvec((wilson_vector *) F_PT(s, src));
    }

    g_floatsum(&size_r);
    size_r = (Real) sqrt((double) (size_r)) / size_src;

    if (this_node == 0 && ((N_iter % 100) == 0 || N_iter == nhop - 1))
      printf("iteration= %d, residue= %e\n", N_iter,
	     (double) (size_r));
    fflush(stdout);

    /* now calculate and write out meson propagators */

    w_meson_hop(meson_prop, src, light_quark, heavy_wall, light_wall,
		new_parity, spin, wallflag);

    if (saveflag_m == SAVE_MESON_ASCII)
    {

      w_ascii_m(fp_m_out, spin, color, N_iter, meson_prop);
    }
    if (saveflag_m == SAVE_MESON_BINARY)
    {
      w_binary_m(fb_m_out, spin, color, N_iter, meson_prop);
    }
    old_parity = new_parity;

    /*
     * dtime1 += dclock(); if(this_node==0) printf(" hopping:  time for this
     * iter= %e\n",dtime1); 
     */
  }				/* end of loop over N_iter  */

  for (channel = 0; channel < nchannels; channel++)
  {
    free(meson_prop[channel]);
  }
  free(meson_prop);
  if (wallflag == CUTOFF_GAUSSIAN || wallflag == CUTOFF_GAUSSIAN)
  {
    free(light_wall);
    free(heavy_wall);
  }
  dtime += dclock();
  if (this_node == 0)
    printf("HOPPING: time = %e = %e/site-iter\n",
	   dtime, dtime / (N_iter * volume));


  return;
}				/* hopping  */


void kappa_dslash(field_offset phi, field_offset temp, Real kappa, int parity)
{
  /*
   * phi  <-  kappa * dslash(phi)  ;  temp is offset of a temporary wilson
   * vector to hold dslash(phi)  
   */
  /* parity is parity of resultant phi */

  register int i;
  register site *s;

  dslash_w_site(phi, temp, PLUS, parity);
  FORSOMEPARITY(i, s, parity)
  {
    scalar_mult_wvec((wilson_vector *) F_PT(s, temp),
		     kappa, (wilson_vector *) F_PT(s, phi));
  }
  cleanup_dslash_wtemps();

}				/* kappa_dslash  */
