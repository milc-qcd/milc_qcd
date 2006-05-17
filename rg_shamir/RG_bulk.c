/************************* RG_bulk.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

/* Map hypercube vector v to index i in the rvs table */
int find_corner(int v[RG_Nd])
{
  int i,vs[RG_Nd];
  
  for (i=0; i<RG_Ncn; i++)
    {
      if((v[0]==rvs[i][0]) && (v[1]==rvs[i][1]) &&
	 (v[2]==rvs[i][2]) && (v[3]==rvs[i][3]))
	{
	  return i;
	  break;
	}
    }
  
  return ;
}


/* For hypercube displacement cmp = S, T, etap, parallel transport
   src[eta(cmp)] to the origin s of the hypercube, using the gauge
   connection wlink */
   
void RG_decompose(QDP_ColorVector *dest, 
		  QDP_ColorVector *src[RG_Ncn],
		  int cmp[4], 
		  QDP_ColorMatrix *wlink[RG_Ncn],
		  QDP_Sub_Block s)
{
  int etap,eta,sign,r[RG_Nd];
  QDP_Shift offset;
  QDP_ColorVector *phi;
  QLA_Real trace;
  

   phi = QDP_create_V();
   QDP_V_eq_zero(phi,QDP_all);

   /* Gamma trace is a sign and a hypercube displacement vector r */
   RG_trace(&sign,r,cmp);
   trace = (QLA_Real) sign;
   eta = find_corner(r);

   cmp[3] = eta;
   etap = cmp[2];

   /* Parallel transport src at hypercube site etap to the
      hypercube origin giving phi */
   offset = QDP_create_shift(rvs[etap]);
#ifndef NOTRANS
   /* Shift with parallel transport */
   SQDP_V_eq_M_times_sV(phi,wlink[etap],src[eta],offset,QDP_forward,s);
#else
   /* If no parallel transport, just shift */
   SQDP_V_eq_sV(phi,src[eta],offset,QDP_forward,s);
#endif
   /* dest = trace * phi */
   SQDP_V_eq_r_times_V(dest,&trace,phi,s);

   QDP_destroy_shift(offset);
   QDP_destroy_V(phi);

return;

}


