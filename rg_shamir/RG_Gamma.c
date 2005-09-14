/************************* RG_Gamma.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include "RG_include.h"
#include <math.h>

int len2;
int r[RG_Nd];


 
void RG_Gamma_Id()
{
 int len,i,j,k,isf,iss,is;
 QLA_Real fact2 = 1.0/2.0;
 QDP_DiracFermion *phi_check[RG_Nf],*phi[RG_Nf],*phi_old[RG_Nf],*phi_new[RG_Nf];
 QDP_DiracFermion *res_old,*res_new,*tempD[RG_Nf],*phi_fin[RG_Nf];
 QDP_ColorVector *res,*chi,*tempV;
 QDP_Subset QDP_block[NRG];



  fprintf(stderr,"-------START CHECK: Gamma_r* Gamma_r = 1 --------------\n");

  len = (int) size/pow(2,1);
  fprintf(stderr,"Check RG trans from lattice %d to lattice %d\n",2*len,len);

  chi = QDP_create_V();
  res = QDP_create_V();
  tempV = QDP_create_V();
  res_old = QDP_create_D();
  res_new = QDP_create_D();
  for(i=0; i<RG_Nf; ++i)
  {
    phi[i] = QDP_create_D();
    phi_old[i] = QDP_create_D();
    phi_new[i] = QDP_create_D();
    tempD[i] = QDP_create_D();
    phi_fin[i] = QDP_create_D();
    phi_check[i] = QDP_create_D();
  }

  RG_create_block(&QDP_block[0],nrg);
  RG_create_block(&QDP_block[1],nrg-1);
  

  RG_initialize_field(phi,0,0,QDP_block[0]);

  for(iss=0;iss<RG_Nf;iss++)
   QDP_D_eq_zero(phi_fin[iss],QDP_block[0]);
    

  for (j=0;j<RG_Ncn;j++)
  {
   QDP_V_eq_zero(res,QDP_block[0]);

   for(isf=0; isf<RG_Nf; ++isf)
   {

   QDP_D_eq_conj_D(phi_old[isf],phi[isf],QDP_block[0]);
  
   for(i=0;i<RG_Nd;i++)
    if (rvs[j][RG_Nd-i-1] != 0 )
     {
     QDP_D_eq_gamma_times_D(phi_new[isf],phi_old[isf],RG_Nd-i-1,QDP_block[0]);
     QDP_D_eq_D(phi_old[isf],phi_new[isf],QDP_block[0]);
     }
   
    QDP_D_eq_conj_D(phi_new[isf],phi_old[isf],QDP_block[0]);
    QDP_D_eq_r_times_D(phi_old[isf],&fact2,phi_new[isf],QDP_block[0]);
    
    QDP_V_eq_colorvec_D(chi,phi_old[isf],isf,QDP_block[0]); 
    QDP_V_peq_V(res,chi,QDP_block[0]); 

   }
    
    for(iss=0;iss<RG_Nf;iss++)
    QDP_D_eq_zero(tempD[iss],QDP_block[0]);

    for(is=0; is<RG_Ns; ++is)
    {
    QDP_D_eq_zero(res_old,QDP_block[0]);
    QDP_D_eq_colorvec_V(res_old,res,is,QDP_block[0]);

    for(i=0;i<RG_Nd;i++)
    if (rvs[j][RG_Nd-i-1] != 0 )
     {
     QDP_D_eq_gamma_times_D(res_new,res_old,RG_Nd-i-1,QDP_block[0]);
     QDP_D_eq_D(res_old,res_new,QDP_block[0]);
     }

     for(iss=0;iss<RG_Nf;iss++)
     {
     QDP_V_eq_colorvec_D(tempV,res_old,iss,QDP_block[0]);
     QDP_D_eq_colorvec_V(tempD[iss],tempV,is,QDP_block[0]);
     }

    
    }
    for(iss=0;iss<RG_Nf;iss++)
     QDP_D_peq_r_times_D(phi_fin[iss],&fact2,tempD[iss],QDP_block[0]);

    }

    for(i=0;i<RG_Nf;i++)
     {
     QDP_D_eq_D_minus_D(phi_check[i],phi_fin[i],phi[i],QDP_block[0]);
     QDP_D_eq_func(phi_check[i], check_df, QDP_block[0]);
     }

fprintf(stderr,"-------END CHECK: Gamma_r* Gamma_r = 1 --------------\n");

  QDP_destroy_V(chi);
  QDP_destroy_V(res);
  QDP_destroy_D(res_old);
  QDP_destroy_D(res_new);
  QDP_destroy_V(tempV);
  for(i=0; i<RG_Nf; ++i)
  {
    QDP_destroy_D(phi[i]);
    QDP_destroy_D(phi_new[i]);
    QDP_destroy_D(phi_old[i]);
    QDP_destroy_D(tempD[i]);
    QDP_destroy_D(phi_fin[i]);
    QDP_destroy_D(phi_check[i]);
  }

return;
}
