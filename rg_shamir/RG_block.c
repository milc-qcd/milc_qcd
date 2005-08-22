/************************* RG_block.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )

int RG_check_for_setup()
{
int min,size,status=0;

 if(QDP_ndim() != RG_Nd )
  {
   fprintf(stderr,"Wrong dimensions\n");
   fprintf(stderr,"wanted %d found %d\n",QDP_ndim(),RG_Nd);
   status = 1;
  }

/* Find shortest lattice extent */
  min = nx;
  min = MIN(min,ny);
  min = MIN(min,nz);
  min = MIN(min,nt);

  size = min;
  node0_printf("min size %d\n",size);

/* The shortest extent must have NRG powers of two */
  if (( size % intpow(2,NRG)) != 0 )
  {
   node0_printf("Wrong number of RG blockings\n");
   node0_printf("%d is not divisible by %d\n",size,intpow(2,NRG));
   fflush(stdout);
   status = 1;
  }

return status;
}


void RG_setup(QDP_Sub_Block QDP_block[NRG+1],QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
int i,j,len;
QDP_ColorMatrix *link_qdp[RG_Nd],*rg_link[NRG][RG_Nd];

  node0_printf("Setup................\n"); fflush(stdout);
  
  rephase(OFF);
#ifdef FIX_GAUGE
   node0_printf("Fix the gauge\n"); fflush(stdout);
/* Lorentz gauge */
  gaugefix(8,(float)1.8,600,(float)GAUGE_FIX_TOL,F_OFFSET(stapleg),F_OFFSET(tempvecg),0,NULL,NULL,0,NULL,NULL);
#endif


  for(i=0; i< RG_Nd; ++i)
   {
    link_qdp[i] = QDP_create_M();
    set_M_from_field(link_qdp[i],F_OFFSET(link[i]));
   }
  rephase(ON);

  for(j=0; j<NRG; ++j)
  for(i=0; i<RG_Nd; ++i)
   {
    rg_link[j][i] = QDP_create_M();
    if (rg_link[j][i] == NULL) terminate(-1);
    SQDP_M_eq_zero(rg_link[j][i],QDP_block[j+1]);
   }


  RG_gauge(rg_link,link_qdp,QDP_block);
  
  for(i=0; i<RG_Nd; ++i)
    QDP_destroy_M(link_qdp[i]);

/* Create paths W(2x,2x+r) that connect the 16 corners of the hypercube */
  for (i=0; i < NRG; i++)
  {
  len = intpow(2,NRG-i-1);

  printf("node %d:paths from lattice %d to lattice %d\n",this_node,2*len,len);
  fflush(stdout);

   for (j = 0; j < RG_Ncn; j++)
    wlink[i][j] = QDP_create_M();
   
  RG_create_path(wlink[i],rg_link[i],QDP_block[i],len);
  //RG_check_hermicity(wlink[i],QDP_block[i]);
  }


  for(j=0; j<NRG; ++j)
  for(i=0; i<RG_Nd; ++i)
    QDP_destroy_M(rg_link[j][i]);

  
return;

}

void RG_coarse_to_fine(QDP_ColorVector *dest,QDP_Sub_Block QDP_block[NRG+1],QDP_ColorVector *src,QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
int i,len;
QDP_ColorVector *phi_c[RG_Ncn],*phi_f;

  phi_f = QDP_create_V();
  for(i=0; i<RG_Ncn; ++i)
   phi_c[i] = QDP_create_V();

  SQDP_V_eq_V(dest,src,QDP_block[0]);


  for (i=0; i < NRG; i++)
  {
   len = intpow(2,NRG-i-1);

   //printf("RG %d from lattice %d to lattice %d:this node %d\n",i,2*len,len,this_node);
   //fflush(stdout);
   
   RG_create_field(phi_c,dest,wlink[i],QDP_block[i]);
   RG_transf_field(phi_f,phi_c,len);

   QDP_V_eq_zero(dest,QDP_all);
   SQDP_V_eq_V(dest,phi_f,QDP_block[i+1]);

  }

    for(i=0; i<RG_Ncn; ++i)
     QDP_destroy_V(phi_c[i]);

     QDP_destroy_V(phi_f);

return;
}

void RG_fine_to_coarse(QDP_ColorVector *dest,QDP_Sub_Block QDP_block[NRG+1],QDP_ColorVector *src,QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
int i,j,len;
QDP_ColorVector *phi_f;

  phi_f = QDP_create_V();
  SQDP_V_eq_V(dest,src,QDP_block[NRG-1]);


  for (i=0; i < NRG-1; i++)
  {

   len = intpow(2,i+1);

   //printf("RG %d from lattice %d to lattice %d:this node %d\n",NRG-i-2,len,2*len,this_node);
   //fflush(stdout);
   
   RG_inv_transf_field(phi_f,dest,wlink[NRG-i-2],QDP_block[NRG-i-2],len);

   QDP_V_eq_zero(dest,QDP_all);
   SQDP_V_eq_V(dest,phi_f,QDP_block[NRG-i-2]);

  }


  QDP_destroy_V(phi_f);


return;
}
