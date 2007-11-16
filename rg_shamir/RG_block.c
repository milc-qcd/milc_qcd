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
  if (( size % intpow(2,nrg)) != 0 )
  {
   node0_printf("Wrong number of RG blockings\n");
   node0_printf("%d is not divisible by %d\n",size,intpow(2,nrg));
   fflush(stdout);
   status = 1;
  }

return status;
}

/* Fix Lorentz gauge and then construct the parallel transports wlink
   needed for all the RG blocking steps. */
void RG_setup(QDP_Sub_Block QDP_block[NRG+1],
	      QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
  int i,j,len;
  QDP_ColorMatrix *link_qdp[RG_Nd],*rg_link[NRG][RG_Nd];
  
  node0_printf("Setup................\n"); fflush(stdout);

  /* Fix gauge and convert from MILC to QDP format */
  rephase(OFF);
#ifdef FIX_GAUGE
  node0_printf("Fix the gauge\n"); fflush(stdout);
  /* Lorentz gauge */
  gaugefix(8,(float)1.8,600,(float)GAUGE_FIX_TOL);
#endif
  

#ifndef NOTRANS  
  /* The field link_qdp is set to a copy of the gauge-fixed gauge field
     without the KS phases */

  for(i=0; i< RG_Nd; ++i)
    {
      link_qdp[i] = QDP_create_M();
      set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
    }
#endif

  rephase(ON);

  /* Skip the rest if we aren't doing parallel transport */  
#ifndef NOTRANS
  /* Initialize intermediate rg_link matrices */
  for(j=0; j<nrg; ++j)
    for(i=0; i<RG_Nd; ++i)
      {
	rg_link[j][i] = QDP_create_M();
	if (rg_link[j][i] == NULL) terminate(-1);
	SQDP_M_eq_zero(rg_link[j][i],QDP_block[j+1]);
      }
  
  /* Construct the on-axis links rg_link for each level of coarseness */
  RG_gauge(rg_link,link_qdp,QDP_block);

  for(i=0; i<RG_Nd; ++i)
    QDP_destroy_M(link_qdp[i]);

  /* Create paths W(2x,2x+r) that connect the 16 corners of the
     hypercube to the hypercube origin for each level of coarseness */
  for (i=0; i < nrg; i++)
    {
      len = intpow(2,nrg-i-1);
      
      printf("node %d:paths from lattice %d to lattice %d\n",this_node,2*len,len);
      fflush(stdout);
      
      for (j = 0; j < RG_Ncn; j++)
	wlink[i][j] = QDP_create_M();
      
      RG_create_path(wlink[i],rg_link[i],QDP_block[i],len);
      //RG_check_hermicity(wlink[i],QDP_block[i]);
    }
  
  /* Destroy the on-axis links */
  for(j=0; j<nrg; ++j)
    for(i=0; i<RG_Nd; ++i)
      QDP_destroy_M(rg_link[j][i]);

#endif /* NOTRANS */
  
  return;
}

/* Starting from a color field defined on the coarsest lattice do
   a sequence of finer and finer inverse RG transformations to 
   distribute the coarse field over the fine lattice */
void RG_coarse_to_fine(QDP_ColorVector *dest,
		       QDP_Sub_Block QDP_block[NRG+1],
		       QDP_ColorVector *src,
		       QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
  int i,len;
  QDP_ColorVector *phi_c[RG_Ncn],*phi_f;
  
  phi_f = QDP_create_V();
  for(i=0; i<RG_Ncn; ++i)
    phi_c[i] = QDP_create_V();
  
  /* First "dest" is the src field */
  SQDP_V_eq_V(dest,src,QDP_block[0]);
  
  
  /* Iterate over levels of RG transformations */
  for (i=0; i < nrg; i++)
    {
      /* The lattice constant at the next fine level */
      len = intpow(2,nrg-i-1);
      /* The current coarse lattice has lattice constant 2*len */
      
      //printf("RG %d from lattice %d to lattice %d:this node %d\n",i,2*len,len,this_node);
      //fflush(stdout);
      
      /* Distribute the coarse field over the hypercube of size len^4
	 whose origin is the coarse lattice site */
      RG_create_field(phi_c,dest,wlink[i],QDP_block[i]);
      RG_transf_field(phi_f,phi_c,len);
      
      QDP_V_eq_zero(dest,QDP_all);
      /* Finally dest is the result of inverse blocking */
      SQDP_V_eq_V(dest,phi_f,QDP_block[i+1]);
    }
  
  for(i=0; i<RG_Ncn; ++i)
    QDP_destroy_V(phi_c[i]);
  
  QDP_destroy_V(phi_f);
  
  return;
}

/* RG block the field src from fine to coarse lattice.  Result in dest */
void RG_fine_to_coarse(QDP_ColorVector *dest,
		       QDP_Sub_Block QDP_block[NRG+1],
		       QDP_ColorVector *src,
		       QDP_ColorMatrix *wlink[NRG][RG_Ncn])
{
  int i,j,len;
  QDP_ColorVector *phi_f;
  
  phi_f = QDP_create_V();

  /* First dest = src on the fine hypercube origins */
  SQDP_V_eq_V(dest,src,QDP_block[nrg-1]);
  
  /* Run through all RG blockings from fine to coarse */
  for (i=0; i < nrg-1; i++)
    {
      
      len = intpow(2,i+1);
      
      //printf("RG %d from lattice %d to lattice %d:this node %d\n",nrg-i-2,len,2*len,this_node);
      //fflush(stdout);
      
      /* Do RG blocking on the hypercubes with side len */
      RG_inv_transf_field(phi_f,dest,wlink[nrg-i-2],QDP_block[nrg-i-2],len);
      
      QDP_V_eq_zero(dest,QDP_all);
      SQDP_V_eq_V(dest,phi_f,QDP_block[nrg-i-2]);
      
    }
  
  
  QDP_destroy_V(phi_f);


return;
}
