/************************* RG_fields.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

/*****************************************************************************/
/*    RG_initialize_field : set N_f dirac spinors to be equal to 1           */
/*    RG_create_field : builds RG_Ncn phi_c(2x) = W(2x,2x+r)*phi(2x),where  */
/*                      r runs over the RG_Ncn corners of the hypercube      */
/*                      ( Q^{\dagger)phi  Shamir notation )                  */
/*                                                                           */
/*****************************************************************************/


#include <stdio.h>
#include <qdp.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"


void point_V(QLA_ColorVector *s, int coords[])
{
  int i;

  for(i=0; i<QLA_Nc; i++) 
  {
    if((coords[0]==0)&&(coords[1]==0)&&(coords[2]==0)&&(coords[3]==0)) 
    {
      if (i == 0){ 
       QLA_c_eq_r(QLA_elem_V(*s,i), 1);
       }
      else
      { 
       QLA_c_eq_r(QLA_elem_V(*s,i), 0);
       }
    }
    else 
    {
    QLA_c_eq_r(QLA_elem_V(*s,i), 0);
    }
  }
}

void point_2V(QLA_ColorVector *s, int coords[])
{
  int i;
  for(i=0; i<QLA_Nc; i++) 
  {
    if((coords[0]==4)&&(coords[1]==4)&&(coords[2]==4)&&(coords[3]==4)) 
    {
      if (i == 0){ 
       QLA_c_eq_r(QLA_elem_V(*s,i), -1);
       }
      else
      { 
       QLA_c_eq_r(QLA_elem_V(*s,i), 0);
       }
    }
/* 
    else 
    {
    QLA_c_eq_r(QLA_elem_V(*s,i), 0);
    }
*/
  }

}

/* Input: 

   wlink[i] -- the gauge connection from the hypercube site i to the
   hypercube origin.

   phi -- the field at the hypercube origin.

   s -- the subset of hypercube origins

   Output:

   phi_c[i] -- wlink[i] * phi.  Lives at the hypercube origin. Still
   needs to be shifted to the hypercube site.

*/

/* First step in parallel transport phi on hypercube origin s to
   phi_c[i] for each hypercube displacement i.  This step just
   multiplies by the gauge connection.  Must follow up with
   RG_transf_field to shift the result to the hypercube location. */

void RG_create_field(QDP_ColorVector *phi_c[RG_Ncn], 
		     QDP_ColorVector *phi, 
		     QDP_ColorMatrix *wlink[RG_Ncn],QDP_Sub_Block s)
{
  int i,j;
  QLA_Real norm = 1.0/16.0;
  //QLA_Real norm = 1.0;
  QDP_ColorVector *phi_1[RG_Ncn]; 
  
  
  for(j=0; j<RG_Ncn; ++j)
    phi_1[j] = QDP_create_V();
  
  for (i=0;i<RG_Ncn;i++)
    {
#ifdef NOTRANS
      /* Just copy */
      SQDP_V_eq_V(phi_1[i],phi,s);
#else
      /* Parallel transport */
      SQDP_V_eq_Ma_times_V(phi_1[i],wlink[i],phi,s);
#endif
      SQDP_V_eq_r_times_V(phi_c[i],&norm,phi_1[i],s);
    }
  
  for(j=0; j<RG_Ncn; ++j)
    QDP_destroy_V(phi_1[j]);
  
  return;
  
}

