/************* make_gammas.c (in su3.a) **************************/

#include "gammatypes.h"

/* Starting from gamma_x, gamma_y, gamma_z, gamma_t, form the rest of the
   16 gamma matrices 
   Gamma matrices are represented via a "gamma_matrix" structure */

/* Does g3 = phase*g1 * g2  - used only to initialize the gamma matrices */
void mult_gamma(int phase, gamma_matrix *g1, gamma_matrix *g2, gamma_matrix *g3)
{
  int r,s;
  for(r=0;r<4;r++)
    {
      s = g1->row[r].column;
      g3->row[r].column = g2->row[s].column;
      g3->row[r].phase  = (g1->row[r].phase + g2->row[s].phase + phase) % 4;
    }
}

void make_gammas(gamma_matrix *gamma)
{
  int r;

  /* gamma_yz = i * gamma_y * gamma_z = sigma_{yz}, etc */
  mult_gamma(1,&gamma[GY ],&gamma[GZ ],&gamma[GYZ]);
  mult_gamma(1,&gamma[GZ ],&gamma[GX ],&gamma[GZX]);
  mult_gamma(1,&gamma[GX ],&gamma[GY ],&gamma[GXY]);

  mult_gamma(1,&gamma[GX ],&gamma[GT ],&gamma[GXT]);
  mult_gamma(1,&gamma[GY ],&gamma[GT ],&gamma[GYT]);
  mult_gamma(1,&gamma[GZ ],&gamma[GT ],&gamma[GZT]);

  /* gamma_5t = gamma_5 * gamma_t = gamma_x * gamma_y * gamma_z */
  /* phase 3 -> -i compensates for the i in gamma_xy */
  mult_gamma(3,&gamma[GX ],&gamma[GYZ],&gamma[G5T]);
  /* gamma_5 = gamma_x * gamma_y * gamma_z * gamma_t */
  mult_gamma(0,&gamma[G5T],&gamma[GT ],&gamma[G5 ]);
  mult_gamma(0,&gamma[G5 ],&gamma[GX ],&gamma[G5X]);
  mult_gamma(0,&gamma[G5 ],&gamma[GY ],&gamma[G5Y]);
  mult_gamma(0,&gamma[G5 ],&gamma[GZ ],&gamma[G5Z]);
 
  /* G1 is the unit matrix */
  for(r=0;r<4;r++)
    {
      gamma[G1].row[r].column = r;
      gamma[G1].row[r].phase = 0;
    }
  gamma_initialized = 1;
}
