/************************* RG_subsests.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
/**********************************************************************/
/*
   Create subsets s[i] = { x: x % pow(2,NRG-i) }
   with lattice spacing a = pow(2,NRG-i) x a'
   where a' is the lattice spacing of the finer lattice
*/
/***********************************************************************/

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

void RG_check_subset(QDP_Sub_Block QDP_block[NRG+1])
{
int i,j,len;
QDP_ColorMatrix *link_qdp[RG_Nd],*prova[RG_Nd];
QLA_Complex unit;

  for(i=0; i< RG_Nd; ++i)
   {
    link_qdp[i] = QDP_create_M();
    prova[i] = QDP_create_M();
   }

  QLA_c_eq_r(unit,1.0)


  for(i=0; i<RG_Nd; ++i)
  SQDP_M_eq_c(link_qdp[i],&unit,QDP_block[nrg]);
  
//  printf("Created!!! this node %d\n",this_node); fflush(stdout);

  for(i=0; i<RG_Nd; ++i)
   SQDP_M_eq_sM(prova[i], link_qdp[i], QDP_neighbor[i], QDP_forward, QDP_block[nrg-1]);

//  printf("I am out!!! this node %d\n",this_node); fflush(stdout);
//  SQDP_M_eq_func(prova[0],print_gl,QDP_block[nrg-1]);
//  printf("I have printed!!! this node %d\n",this_node); fflush(stdout);
 
 
  for(i=0; i< RG_Nd; ++i)
   {
    QDP_destroy_M(link_qdp[i]);
    QDP_destroy_M(prova[i]);
   }

  printf("I have destroyed every thing!!! this node %d\n",this_node); fflush(stdout);

return;
}
