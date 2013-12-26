/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
#include "wilson_flow_includes.h"
#include "../include/field_strength.h"

/* Computes the real trace of the su3 matrix product: ReTr(A*B) */
Real
real_trace_nn( su3_matrix *a, su3_matrix *b )
{
  register int i,j;
  register complex x;
  register Real sum;

  sum = 0.0;
  for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) {
    CMUL( a->e[i][j], b->e[j][i], x );
    sum += x.real;
  }

  return sum;
}

/* Computes the field strength components and topological charge */
void
fmunu_fmunu( double *time, double *space, double *charge )
{
  /* Site variables */
  register int i;
  register site *s;

  /* Temporary component storage */
  su3_matrix *ft, *fs;

  /* Initialize sums */
  *time = *space = *charge = 0;

  /* Compute 8*F_mu,nu at each site */
  make_field_strength( F_OFFSET(link), F_OFFSET(fieldstrength) );

  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s) {
    
    fs = &(s->fieldstrength[FS_XY]);
    ft = &(s->fieldstrength[FS_ZT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= real_trace_nn(fs, ft);

    fs = &(s->fieldstrength[FS_XZ]);
    ft = &(s->fieldstrength[FS_YT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= realtrace_su3(fs, ft); /* ReTr{ fs.dag * ft } */

    fs = &(s->fieldstrength[FS_YZ]);
    ft = &(s->fieldstrength[FS_XT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= real_trace_nn(fs, ft);
  }

  /* Sum over all nodes */
  g_doublesum(time);
  g_doublesum(space);
  g_doublesum(charge);

  /* Normalizations */
  *time /= (volume*64.0);
  *space /= (volume*64.0);
  *charge *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
}
