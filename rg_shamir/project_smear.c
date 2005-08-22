/************************* project_smear.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

#define TOL 1e-5
#define MAXCOUNT 100

void project_qdp(QDP_ColorMatrix *link_src[],QDP_ColorMatrix *link_dest[],int *space_only)
{
  int i,j,n;
  site *s;
  QLA_ColorMatrix *temp_dest, *temp_src;
  su3_matrix  c,h;

/* Expose QDP field */
   if ( *space_only == 1 ) n = 3;
   if ( *space_only == 0 ) n = 4;
   if ( *space_only == RG_Ncn) n = RG_Ncn;

  printf("node %d: projecting %d......... \n",this_node,n); fflush(stdout);
   for (i=0; i< n; i++)
   {
   temp_src = QDP_expose_M(link_src[i]);
   temp_dest = QDP_expose_M(link_dest[i]);
   FORALLSITES(j,s)
    {
     /* Copy temp[j] to c */
     memcpy(&c,(void *)&temp_src[j],sizeof(QLA_ColorMatrix));
     h = c ;
     reunit_su3(&h);
     project_su3(&h,&c,3*MAXCOUNT,TOL);
     memcpy((void *)&temp_dest[j],&h,sizeof(QLA_ColorMatrix));
    }
     QDP_reset_M(link_src[i]);
     QDP_reset_M(link_dest[i]);
    
   }
//   printf(".........done \n"); fflush(stdout);


  return ;
}
