/***************** u1plaq.c *************************************/

/* MIMD version 7 */
/* Calculate average of exp(iA) for spatial and temporal links  */
/* ************************************************************	*/

#include "generic_u1_includes.h"

void u1avlink(double *sLink, double *tLink, double charge)
{
  int i;

  *sLink = *tLink = 0.;

  FORALLFIELDSITES(i){
    int dir;
    FORALLUPDIRBUT(TUP,dir){
      *sLink += cos(charge*u1_A[4*i+dir]);
    }
    *tLink += cos(charge*u1_A[4*i+TUP]);
  }
  *sLink /= 3.;
  g_doublesum(sLink);
  g_doublesum(tLink);
  *sLink /= volume;
  *tLink /= volume;

} /* end of u1avlink */

/* ************************************************************	*/

