/************************* RG_smear.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
//#include <stdlib.h>
#include <stdio.h>
#include <qdp.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

void RG_smear_dir (QDP_ColorMatrix *sm_link, 
		   QDP_ColorMatrix *link[], 
		   QLA_Real w_l, QLA_Real w_s, 
		   QLA_Int dir, QDP_Sub_Block s, int len)
{
  int i,v[RG_Nd],n;
  QLA_Int nu,mu=dir;
  QDP_Subset sub;
  QDP_Shift offset[RG_Nd];
  QLA_Complex unit;
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4, *temp5, *temp6;
  
  
  temp1 = QDP_create_M();
  temp2 = QDP_create_M();
  temp3 = QDP_create_M();
  temp4 = QDP_create_M();
  temp5 = QDP_create_M();
  temp6 = QDP_create_M();
  
  for(nu=0; nu < RG_Nd ; nu++)
    {
      for(i=0; i<RG_Nd;i++) v[i] = 0;
      v[nu] = len;
      offset[nu] = QDP_create_shift(v);
    }
  
  SQDP_M_eq_r_times_M(temp6,&w_l,link[mu],s);
  
  /* Set temp4 to zero */
  SQDP_M_eq_zero(temp4,s);
  
  n = RG_Nd;
#ifdef CHECK_SMEAR_QDP_MILC
  n = 3;
#endif
  
  /* Sum on staples */
  for(nu=0; nu < n ; nu++)if(nu != mu)
    {
      
      /* For forward staples */
      SQDP_M_eq_sM(temp1, link[mu], offset[nu], QDP_forward, s);
      SQDP_M_eq_sM(temp2, link[nu], offset[mu], QDP_forward, s);
      SQDP_M_eq_M_times_Ma(temp3, temp1, temp2, s);
      SQDP_M_peq_M_times_M(temp4, link[nu], temp3, s);
      
      /* For backward staples */
      SQDP_M_eq_M_times_M(temp3, link[mu], temp2, s);
      SQDP_M_eq_Ma_times_M(temp1, link[nu], temp3, s);
      SQDP_M_eq_sM(temp5, temp1, offset[nu], QDP_backward, s);
      SQDP_M_peq_M(temp4, temp5, s);
      
    }
  
  /* U_smeared = w_l * U + w_s * U_staple */
  SQDP_M_eq_r_times_M_plus_M(sm_link,&w_s,temp4,temp6,s);
  
  
  QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);
  QDP_destroy_M(temp5);
  QDP_destroy_M(temp6);
  
  for(nu=0; nu < RG_Nd ; nu++)
    QDP_destroy_shift(offset[nu]);
  
  return ;
  
}

void RG_smearing_qdp (QDP_ColorMatrix *sm_link[], 
		      QDP_ColorMatrix *link[], 
		      QLA_Real *s_w, QLA_Real *l0, 
		      QDP_Sub_Block s, int len)
{
  QLA_Int nstaples,dir,n;
  QLA_Real w_l,w_s,norm_fact;
  
  nstaples = 6;
  n=RG_Nd;
  
#ifdef CHECK_SMEAR_QDP_MILC
  nstaples = 4;
  n=3;
#endif
  //node0_printf("Smearing......... \n"); fflush(stdout);
  norm_fact = 1.0/(1.0 + nstaples*(*l0)*(*l0)*(*s_w));
  w_l = norm_fact;
  w_s = norm_fact*(*s_w);
  
  for ( dir = 0; dir < n; dir++)
    RG_smear_dir(sm_link[dir],link,w_l,w_s,dir,s,len);
  
  //node0_printf(".........done \n"); fflush(stdout);
  
  return;
}


