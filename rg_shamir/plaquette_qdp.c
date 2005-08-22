/************************* plaquette_qdp.c *******************************/
/* MIMD version 7 */
/* J. Osborn F. Maresca Jul 2005 */
//#include <stdlib.h>
#include <stdio.h>
#include <qdp.h>
#include "RG_include.h"
#include "RG_Shamir_includes.h"
//#define LOCAL_SUM

QLA_Real plaquette_qdp(QDP_ColorMatrix *link_qdp[],QDP_Sub_Block s,int n)
{
  int vol,len,v[RG_Nd];
  int mu, nu;
  QLA_Real plaq;
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4;
  QDP_Shift offset[RG_Nd];
#ifdef LOCAL_SUM
  QDP_Real *treal1, *treal2;
#else
  QLA_Real tplaq;
#endif

  len = intpow(2,n);

  for(nu=0; nu<RG_Nd; ++nu) 
  {
   for(mu=0; mu<RG_Nd; ++mu)  v[mu] = 0;

   v[nu] = len;
   offset[nu] = QDP_create_shift(v);
  }

  
  
#ifdef LOCAL_SUM
  treal1 = QDP_create_R();
  treal2 = QDP_create_R();
  QDP_R_eq_zero(treal2, QDP_all);
#else
  plaq = 0;
#endif

  vol= nx*ny*nz*nt;
  temp1 = QDP_create_M();
  temp2 = QDP_create_M();
  temp3 = QDP_create_M();
  temp4 = QDP_create_M();

  for(mu=0; mu<3; ++mu) {
    for(nu=mu+1; nu<4; ++nu) {

      SQDP_M_eq_sM(temp1, link_qdp[nu], offset[mu], QDP_forward, s);
      SQDP_M_eq_sM(temp2, link_qdp[mu], offset[nu], QDP_forward, s);

      SQDP_M_eq_Ma_times_M(temp3, link_qdp[nu], link_qdp[mu], s);

      SQDP_M_eq_M_times_M(temp4, temp3, temp1, s);
      QDP_discard_M(temp1);

#ifdef LOCAL_SUM
      SQDP_R_eq_re_M_dot_M(treal1, temp2, temp4, s);
      QDP_discard_M(temp2);
      QDP_R_peq_R(treal2, treal1, QDP_all);
#else
      SQDP_r_eq_re_M_dot_M(&tplaq, temp2, temp4, s);
      QDP_discard_M(temp2);
      plaq += tplaq;
#endif

    }
  }

#ifdef LOCAL_SUM
  QDP_r_eq_sum_R(&plaq, treal2, s);
  QDP_destroy_R(treal1);
  QDP_destroy_R(treal2);
#endif

  QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);

  return plaq/(3*vol);
}

