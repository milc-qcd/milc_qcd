/************************** wall_sink_h.c *****************************/
/* (Identical to heavy/wall_sink_h.c) */
/* MIMD version 7 */

/* CB wrote this */
/* 11/25/97 modifications for version 5 CD */
/*  2/14/98 corrected memory leak CD */

/* Wall sink for quark propagator in hopping expansion context */
/* sink is only non-zero on a fixed parity set of sites */

#include "w_static_includes.h"

static Real *sink_wall_template = NULL;

/* src is offset of field of type wilson_vector */

void wall_sink_h(field_offset src, wilson_vector * chi_out,
		   int parity, int keep_parity, int x0, int y0, int z0)
{
  register int i;
  register site *s;

  int t, tp, tb, te;

  int rx, ry, rz;
  Real scale;
  int center_parity;

/**  double dtime ;  **/
  
  /* Make template if it is not defined */
  if(sink_wall_template == NULL)
    sink_wall_template = make_template(1./(wqs.r0*wqs.r0),wqs.wall_cutoff);

  /* Gaussian sink centered on  spatial position x0,y0,z0;   */
  /* cutoff a rectangular distance wqs.wall_cutoff from center */

  /* clear quark propagators */
  for (t = 0; t < nt; t++)
  {
    clear_wvec(&(chi_out[t]));
  }
  center_parity = (x0 + y0 + z0) % 2;

  /*
   * find beginning and ending times for walls; remember even slices stored
   * before odd slices 
   */
  if (keep_parity == EVENANDODD)
  {
    tb = 0;
    te = nt;
  } else
  {
    /* on every time slice, heavy wall has spatial parity of source */
    /*
     * keep_parity may chosen different from parity, because you may not want
     * to bother doing a light wall on a time slice if it's just going to
     * multiplied by a heavy wall which is 0  
     */
    tb = ((keep_parity + center_parity) % 2) * (nt / 2);
    te = tb + nt / 2;
  }




  /* dtime = -dclock();  */
  FORSOMEPARITY(i, s, parity)
  {				/* propagator (i.e. src) is only non-zero for
				 * this parity; for light props, want
				 * parity=EVENANDODD */
    t = s->t;
    tp = (t % 2) * (nt / 2) + t / 2;	/* even time slices stored first in
					 * wall propagator */
    if (tp < tb || tp >= te)
      continue;

    rz = abs(s->z - z0);
    if (rz > (nz / 2))
      rz = nz - rz;
    if (rz > wqs.wall_cutoff)
      continue;			/* don't do anything outside wall_cutoff */
    ry = abs(s->y - y0);
    if (ry > (ny / 2))
      ry = ny - ry;
    if (ry > wqs.wall_cutoff)
      continue;			/* don't do anything outside wall_cutoff */
    rx = abs(s->x - x0);
    if (rx > (nx / 2))
      rx = nx - rx;
    if (rx > wqs.wall_cutoff)
      continue;			/* don't do anything outside wall_cutoff */
    scale = sink_wall_template[rx + (wqs.wall_cutoff + 1) * (ry + (wqs.wall_cutoff + 1) * rz)];
    scalar_mult_add_wvec(&(chi_out[tp]), (wilson_vector *) F_PT(s, src),
			 scale, &(chi_out[tp]));

  }				/* sites */
  /*
   * dtime += dclock(); if(this_node==0) printf(" wall_sink_h:  calc wall on
   * node; time= %e\n",dtime); dtime = -dclock();   
   */


  /*
   * do global sum of all relevant elements of chi_out at once; it takes much
   * too long to sum individual complex numbers 
   */
  g_veccomplexsum((complex *) (chi_out + tb), 12 * (te - tb));

  /*
   * dtime += dclock(); if(this_node==0) printf(" wall_sink_h: 
   * g_veccomplexsums; time= %e\n",dtime); 
   */


}				/* wall_sink_h */

/* Free wall template malloc */
void free_sink_template()
{
  if(sink_wall_template != NULL)free(sink_wall_template);
  sink_wall_template = NULL;
}
