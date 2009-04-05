/******** w_meson_hop.c *************/
/* (Identical to heavy/w_meson_hop.c) */
/*** WARNING:  
  as of 9/22/92, this routine exposes certain optimizing compiler
  bugs on the ipsc with the portland group compiler:
  compile only with CFLAGS -> CFLAGS5 =  
  -O4 -Mquad -W0,-dalign -DPROTO -Di860 -I$( LIBDIR) $(<lattice.h>) 
  (note that -Mx,2,15 or -Mx,2,13 are removed) 
  *****************************************/

/* MIMD version 7 */

/*
 * computes & writes out meson propagator channels at each order in the
 * hopping param expansion 
 */

#include "w_static_includes.h"



/* light_quark and heavy_quark are offsets of fields of  type wilson_vector */

void w_meson_hop(double *prop[], field_offset heavy_quark, field_offset light_quark,
		   wilson_vector * heavy_wall, wilson_vector * light_wall,
		   int parity, int spin, int wallflag)
{
  /* pion and rho */

  register int i;
  register site *s;

  int channel, sign53 = 0;
  int t, tp, tb, te, cf, sf;
  complex g1, g2;
  int dx, dy, dz, x0, y0, z0, center_parity;
  static int calls = 0;


  wilson_vector localvec;	/* temporary storage for quark */
  wilson_vector localvec2;	/* temporary storage for quark */

  double dtime = 0. ;


  calls++;

  /* clear propagators */
  for (channel = 0; channel < nchannels; channel++)
    for (t = 0; t < nt; t++)
    {
      prop[channel][t] = 0.0;
    }


  /*
   * get signs corresponding to diagonal components in B&D conventions of
   * (-gamma5)*gamma3 (for the rho);  minus gamma5 because gamma5 transforms
   * with opposite sign from gamma_mu when switching conventions ---i.e. it
   * is not the identical product of gamma matrices in the 2 conventions  
   */
  /* (for (-5)3, component = sign53*I) */
  /*
   * quarks must be computed with a source in B&D gamma matrix conventions;
   * this is hard-wired in.   see bj_to_weyl.c for conventions   
   */
  switch (spin)
  {
  case 0:
    sign53 = -1;
    break;
  case 1:
    sign53 = 1;
    break;
  case 2:
    sign53 = 1;
    break;
  case 3:
    sign53 = -1;
    break;
  }




  /* now start various channels:    */

  /* first do channels that are the same for walls and points */

  /* gamma5-gamma5gamma0 channel  --- always local sink even if wall sources */
  channel = 1;

  FORSOMEPARITY(i, s, parity)
  {
    t = s->t;
    /* store even slices first, then odds; nt must be even */
    tp = (t % 2) * (nt / 2) + t / 2;
    mult_by_gamma_right_vec((wilson_vector *) F_PT(s, light_quark), &localvec, TUP);

    /* trace over propagators */
    for (sf = 0; sf < 4; sf++)
      for (cf = 0; cf < 3; cf++)
      {
	g2 = ((wilson_vector *) F_PT(s, heavy_quark))->d[sf].c[cf];
	g1 = localvec.d[sf].c[cf];
	prop[channel][tp] += g1.real * g2.real + g1.imag * g2.imag;
      }
  }
  /* sum propagator over nodes */
  g_vecdoublesum(prop[channel], nt);


  /* gamma3-gamma3 channel  --- point sink by definition */
  channel = 3;

  FORSOMEPARITY(i, s, parity)
  {
    t = s->t;
    /* store even slices first, then odds; nt must be even */
    tp = (t % 2) * (nt / 2) + t / 2;
    mult_by_gamma_right_vec((wilson_vector *) F_PT(s, light_quark), &localvec2, ZUP);
    mult_by_gamma_right_vec(&localvec2, &localvec, GAMMAFIVE);


    /* trace over propagators */
    for (sf = 0; sf < 4; sf++)
      for (cf = 0; cf < 3; cf++)
      {
	g2 = ((wilson_vector *) F_PT(s, heavy_quark))->d[sf].c[cf];
	g1 = localvec.d[sf].c[cf];

	/* this is (Real(sign53 * I *conj(g2) * g1))  */
	prop[channel][tp] += sign53 * (g2.imag * g1.real - g2.real * g1.imag);
      }
  }
  /* sum propagator over nodes */
  g_vecdoublesum(prop[channel], nt);



  /* now for channels that are different if wall or points */
  if (wallflag == POINT)
  {

    /* gamma5-gamma5 channel  --- different if wall or points */
    channel = 0;

    FORSOMEPARITY(i, s, parity)
    {
      t = s->t;
      /* store even slices first, then odds; nt must be even */
      tp = (t % 2) * (nt / 2) + t / 2;

      /* trace over propagators */
      for (sf = 0; sf < 4; sf++)
	for (cf = 0; cf < 3; cf++)
	{
	  g2 = ((wilson_vector *) F_PT(s, heavy_quark))->d[sf].c[cf];
	  g1 = ((wilson_vector *) F_PT(s, light_quark))->d[sf].c[cf];
	  prop[channel][tp] += g1.real * g2.real + g1.imag * g2.imag;
	}
    }
    /* sum propagator over nodes */
    g_vecdoublesum(prop[channel], nt);


    /* gamma3-gamma3 channel  --- different if wall or points */
    /* this is the same as channel 3 if POINT sources are run)  */
    channel = 2;

    prop[channel] = prop[3];

  }
   /* end of POINT sink channels  */ 
  else
  if (wallflag == CUTOFF_GAUSSIAN || wallflag == CUTOFF_GAUSSIAN_WEYL)
  {

    /* time this every 100th call */
    if (this_node == 0 && calls % 100 == 0)
      dtime -= dclock();


    /* loop over centers of wall sinks */


    for (dz = 0; dz < nz; dz += wall_separation)
      for (dy = 0; dy < ny; dy += wall_separation)
	for (dx = 0; dx < nx; dx += wall_separation)
	{
	  z0 = (wqs.z0 + dz) % nz;
	  y0 = (wqs.y0 + dy) % ny;
	  x0 = (wqs.x0 + dx) % nx;
	  center_parity = (x0 + y0 + z0) % 2;

	  /*
	   * find beginning and ending times for walls; remember even slices
	   * stored before odd slices 
	   */
	  if (parity == EVENANDODD)
	  {
	    tb = 0;
	    te = nt;
	  } else
	  {
	    tb = ((parity + center_parity) % 2) * (nt / 2);
	    te = tb + nt / 2;
	  }


	  /* find the walls */
	  wall_sink_h(heavy_quark, heavy_wall, parity, parity, x0, y0, z0);
	  /* heavy quark has only one non-zero parity at each order */
	  /* first "parity" argument of wall_sink_h tells which sites are */
	  /* non-zero; second "parity" argument is used to tell which time */
	  /* slices to skip */
	  wall_sink_h(light_quark, light_wall, EVENANDODD, parity, x0, y0, z0);
	  /*
	   * light quark has both parities, but skip slices where heavy wall
	   * vanishes 
	   */

	  /* gamma5-gamma5 channel  --- different if wall or points */
	  channel = 0;

	  for (tp = tb; tp < te; tp++)
	  {			/* trace over propagators */
	    /* don't bother if slice has wrong parity */

	    for (sf = 0; sf < 4; sf++)
	      for (cf = 0; cf < 3; cf++)
	      {
		g2 = heavy_wall[tp].d[sf].c[cf];
		g1 = light_wall[tp].d[sf].c[cf];
		prop[channel][tp] += g1.real * g2.real + g1.imag * g2.imag;
	      }
	  }
	  /*
	   * summing meson propagators over nodes is not done for wall sinks;
	   * quark wall propagators are already summed over nodes 
	   */


	  /* gamma3-gamma3 channel  --- different if wall or points */
	  channel = 2;

	  for (tp = tb; tp < te; tp++)
	  {			/* trace over propagators */
	    /* don't bother if slice has wrong parity */
	    mult_by_gamma_right_vec(&light_wall[tp], &localvec2, ZUP);
	    mult_by_gamma_right_vec(&localvec2, &localvec, GAMMAFIVE);

	    /* trace over propagators */
	    for (sf = 0; sf < 4; sf++)
	      for (cf = 0; cf < 3; cf++)
	      {
		g2 = heavy_wall[tp].d[sf].c[cf];
		g1 = localvec.d[sf].c[cf];

		/* this is (Real(sign53 * I *conj(g2) * g1))  */
		prop[channel][tp] += sign53 * (g2.imag * g1.real - g2.real * g1.imag);
	      }
	  }
	  /*
	   * summing meson propagators over nodes is not done for wall sinks;
	   * quark wall propagators are already summed over nodes 
	   */



	}			/* end of sum over wall sink centers */
    if (this_node == 0 && calls % 100 == 0)
    {
      dtime += dclock();
      printf("call %d to w_meson_hop: time for wall channels= %e\n", calls, dtime);
      fflush(stdout);
    }
    /*
     * if extra_sink chosen, loop over centers of wall sinks translated by
     * nx/4, ny/4, nz/4 
     */

    if (nchannels > NCHANNELS)
    {
      /* time this every 100th call */
      if (this_node == 0 && calls % 100 == 0)
	dtime -= dclock();

      for (dz = 0; dz < nz; dz += wall_separation)
	for (dy = 0; dy < ny; dy += wall_separation)
	  for (dx = 0; dx < nx; dx += wall_separation)
	  {
	    z0 = ((nz / 4) + wqs.z0 + dz) % nz;
	    y0 = ((ny / 4) + wqs.y0 + dy) % ny;
	    x0 = ((nx / 4) + wqs.x0 + dx) % nx;
	    center_parity = (x0 + y0 + z0) % 2;

	    /*
	     * find beginning and ending times for walls; remember even
	     * slices stored before odd slices 
	     */
	    if (parity == EVENANDODD)
	    {
	      tb = 0;
	      te = nt;
	    } else
	    {
	      tb = ((parity + center_parity) % 2) * (nt / 2);
	      te = tb + nt / 2;
	    }


	    /* find the walls */
	    wall_sink_h(heavy_quark, heavy_wall, parity, parity, x0, y0, z0);
	    /* heavy quark has only one non-zero parity at each order */
	    /* first "parity" argument of wall_sink_h tells which sites are */
	    /* non-zero; second "parity" argument is used to tell which time */
	    /* slices to skip */
	    wall_sink_h(light_quark, light_wall, EVENANDODD, parity, x0, y0, z0);
	    /*
	     * light quark has both parities, but skip slices where heavy
	     * wall vanishes 
	     */

	    /* gamma5-gamma5 channel  --- different if wall or points */
	    channel = NCHANNELS;

	    for (tp = tb; tp < te; tp++)
	    {			/* trace over propagators */
	      /* don't bother if slice has wrong parity */

	      for (sf = 0; sf < 4; sf++)
		for (cf = 0; cf < 3; cf++)
		{
		  g2 = heavy_wall[tp].d[sf].c[cf];
		  g1 = light_wall[tp].d[sf].c[cf];
		  prop[channel][tp] += g1.real * g2.real + g1.imag * g2.imag;
		}
	    }
	    /*
	     * summing meson propagators over nodes is not done for wall
	     * sinks; quark wall propagators are already summed over nodes 
	     */


	    /* gamma3-gamma3 channel  --- different if wall or points */
	    channel = NCHANNELS + 1;

	    for (tp = tb; tp < te; tp++)
	    {			/* trace over propagators */
	      /* don't bother if slice has wrong parity */
	      mult_by_gamma_right_vec(&light_wall[tp], &localvec2, ZUP);
	      mult_by_gamma_right_vec(&localvec2, &localvec, GAMMAFIVE);

	      /* trace over propagators */
	      for (sf = 0; sf < 4; sf++)
		for (cf = 0; cf < 3; cf++)
		{
		  g2 = heavy_wall[tp].d[sf].c[cf];
		  g1 = localvec.d[sf].c[cf];

		  /* this is (Real(sign53 * I *conj(g2) * g1))  */
		  prop[channel][tp] += sign53 * (g2.imag * g1.real - g2.real * g1.imag);
		}
	    }
	    /*
	     * summing meson propagators over nodes is not done for wall
	     * sinks; quark wall propagators are already summed over nodes 
	     */



	  }			/* end of sum over wall sink centers */
      if (this_node == 0 && (calls % 100) == 0)
      {
	dtime += dclock();
	printf("call %d to w_meson_hop: time for wall channels= %e\n", calls, dtime);
	fflush(stdout);
      }
    }				/* end of extra_sink stuff */
  }				/* end of WALL sink channels */
}				/* w_meson_hop */
