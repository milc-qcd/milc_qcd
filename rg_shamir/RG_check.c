/************************* RG_check.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"
#define TOL 10e-6

void RG_check_hermicity(QDP_ColorMatrix *w[RG_Ncn],QDP_Sub_Block s)
{
int i,j;
QLA_Real c =1.0;
QLA_Complex one;
QDP_ColorMatrix *unit,*unit1,*unit2,*check1,*check2;

  QLA_C_eq_R(&one,&c);

   unit = QDP_create_M();
   unit1 = QDP_create_M();
   unit2 = QDP_create_M();
   check1 = QDP_create_M();
   check2 = QDP_create_M();

  for (j = 0; j < RG_Ncn; j++)
  {

   printf("Check hermicity of link %d\n",j);
   SQDP_M_eq_zero(unit,s);
   SQDP_M_eq_zero(unit1,s);
   SQDP_M_eq_zero(unit2,s);
   SQDP_M_eq_zero(check1,s);
   SQDP_M_eq_zero(check2,s);

   SQDP_M_eq_c(unit, &one, s);

   SQDP_M_eq_M_times_Ma(unit1,w[j],w[j],s);
   SQDP_M_eq_M_minus_M(check1,unit1,unit,s);
   SQDP_M_eq_func(check1,check_gl,s);
   SQDP_M_eq_Ma_times_M(unit2,w[j],w[j],s);
   SQDP_M_eq_M_minus_M(check2,unit2,unit,s);
   SQDP_M_eq_func(check2,check_gl,s);

   }

  QDP_destroy_M(unit);
  QDP_destroy_M(unit1);
  QDP_destroy_M(unit2);
  QDP_destroy_M(check1);
  QDP_destroy_M(check2);

return;
}

void RG_transformation( QDP_Sub_Block QDP_block[NRG+1])
{
  int i,j,len;
  QLA_Real c, fact_n,norm=1.0/16.0;
  QDP_ColorVector *phi_c,*phi;
  QDP_ColorVector *phi_1,*phi_f,*res;
  QDP_ColorMatrix *wlink[NRG][RG_Ncn];



  for (i = 0; i < nrg; i++)
  for (j = 0; j < RG_Ncn; j++)
   wlink[i][j] = QDP_create_M();

  RG_setup(QDP_block,wlink);
  node0_printf("-------START CHECK: hermicity of the paths--------------\n");
  node0_printf("Expected a zero matrix\n"); fflush(stdout);
  for (i = 0; i < nrg; i++)
   RG_check_hermicity(wlink[i],QDP_block[i]);
  node0_printf("-------END CHECK: hermicity of the paths --------------\n");

  node0_printf("-------START CHECK: Wilson<->Wilson --------------------\n");
  fflush(stdout);

  res = QDP_create_V();
  phi = QDP_create_V();
  phi_c = QDP_create_V();
  phi_f = QDP_create_V();
  phi_1 = QDP_create_V();


  SQDP_V_eq_func(phi, point_V, QDP_block[0]);
  
  printf("Color source set to one at the origin\n"); fflush(stdout);
  SQDP_V_eq_func(phi, print_cv, QDP_block[0]);
 
  RG_coarse_to_fine(phi_c,QDP_block,phi,wlink);
  
  //printf("RG %d from lattice %d to lattice %d\n",nrg-1,1,2); fflush(stdout);
 
  RG_inv_transf_field(phi_1,phi_c,wlink[nrg-1],QDP_block[nrg-1],1); 

  RG_fine_to_coarse(phi_f, QDP_block, phi_1, wlink);

  fact_n = pow(norm,(double)nrg);

  SQDP_V_eq_r_times_V_minus_V(res,&fact_n,phi,phi_f,QDP_block[0]);
 
  node0_printf("Expected a zero vector\n"); fflush(stdout);
  SQDP_V_eq_func(res, check_cv, QDP_block[0]);

 QDP_destroy_V(res);
 QDP_destroy_V(phi);
 QDP_destroy_V(phi_c);
 QDP_destroy_V(phi_f);
 QDP_destroy_V(phi_1);
  node0_printf("---------------END CHECK: Wilson<->Wilson --------------\n");
  fflush(stdout);

 for (i = 0; i < nrg; i++)
  for(j=0; j<RG_Ncn; ++j)
   QDP_destroy_M(wlink[i][j]);

return;

}

int RG_check(QDP_Sub_Block s[NRG+1])
{
int status = 0;


/* Check DeGrand smearing */
#ifdef CHECK_LINK
  if(RG_check_smear(s) != 0 ) status = 1;
#endif

#ifdef CHECK_TRANS
   RG_transformation(s);
#endif

return status;
}

int RG_check_smear(QDP_Sub_Block s[NRG+1])
{
 int i,j,len,status = 0;
 QLA_Real dssplaq[NRG],dssplaq_g[NRG];
 QLA_Real dssplaq1,dssplaq_g1,dssplaqn;
 QDP_ColorMatrix *link_qdp[RG_Nd],*link_qdp_g[RG_Nd];
 QDP_ColorMatrix *rg_link[NRG][RG_Nd],*rg_link_g[NRG][RG_Nd];
 QDP_ColorMatrix *pr_link[RG_Nd],*pr_link_g[RG_Nd];


#ifdef CHECK_PLAQ 
 printf("Check plaquette eavluated with QDP\n");

 for(i=0; i< RG_Nd; ++i)
  {
   link_qdp[i] = QDP_create_M();
   set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
  }

  dssplaq1=plaquette_qdp(link_qdp,s[nrg],0);
  rephase(OFF);
  rand_gauge(F_OFFSET(rgt));
  rephase(ON);
  
  for(i=0; i< RG_Nd; ++i)
   {
    link_qdp_g[i] = QDP_create_M();
    pr_link_g[i] = QDP_create_M();
    set_M_from_site(link_qdp_g[i],F_OFFSET(link[i]),EVENANDODD);
   }
 
  dssplaq_g1=plaquette_qdp(link_qdp_g,s[nrg],0);


  fprintf(stderr,"CHECK PLAQ: gauge diff: %e\n",dssplaq_g1-dssplaq1);
  if (dssplaq_g1-dssplaq1 > TOL )   
    {
     fprintf(stderr,"Error: plaq not gauge invariant \n");
     fprintf(stderr,"diff: %e\n",dssplaq_g1-dssplaq1);
     status = 1;
    }

#endif

#ifdef CHECK_SMEAR_QDP_MILC
 printf("Check difference QDP and MILC in smearing \n");

 rephase(OFF);
 for(i=0; i< RG_Nd; ++i)
  {
   pr_link[i] = QDP_create_M();
   link_qdp[i] = QDP_create_M();
   set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
  }


 RG_smearing(pr_link,link_qdp,s[nrg],1);
 SQDP_M_eq_M(pr_link[3],link_qdp[3],s[nrg]);
 dssplaq1=plaquette_qdp(pr_link,s[nrg],0);
 smearing();
 for(i=0; i< RG_Nd; ++i)
  set_M_from_site(pr_link[i],F_OFFSET(link[i]),EVENANDODD);
 dssplaqn=plaquette_qdp(pr_link,s[nrg],0);
 
 printf("QDP/MILC %e/%e\n",dssplaq1,dssplaqn);
 printf("diff: %e\n",dssplaqn-dssplaq1);
  if (fabs(dssplaq1-dssplaqn) > TOL )   
    {
     printf("Error: QDP not equal to MILC \n");
     printf("diff: %e\n",dssplaqn-dssplaq1);
     status = 1;
    }

#endif
 
#ifdef CHECK_SMEAR_MILC_2
 printf("Check difference in MILC, s(l) = sg(l)\n");
 rephase(OFF);
 for(i=0; i< RG_Nd; ++i)
  {
   pr_link[i] = QDP_create_M();
   set_M_from_site(pr_link[i],F_OFFSET(link[i]),EVENANDODD);
  }

 smearing();
 for(i=0; i< 3; ++i)
  set_M_from_site(pr_link[i],F_OFFSET(sm_link[i]),EVENANDODD);
 dssplaq1=plaquette_qdp(pr_link,s[nrg],0);

 rand_gauge(F_OFFSET(rgt));

 for(i=0; i< RG_Nd; ++i)
  {
   pr_link_g[i] = QDP_create_M();
   set_M_from_site(pr_link_g[i],F_OFFSET(link[i]),EVENANDODD);
  }

 smearing();
 for(i=0; i< 3; ++i)
  set_M_from_site(pr_link_g[i],F_OFFSET(sm_link[i]),EVENANDODD);
 dssplaq_g1=plaquette_qdp(pr_link_g,s[nrg],0);
 
 printf("SMEAR-GAUGED TYPE 2 MILC %e/%e\n",dssplaq1,dssplaq_g1);
 
  if (fabs(dssplaq1-dssplaq_g1) > TOL )   
    {
     printf("Error: MILC SMEAR not gauge inv. \n");
     printf("diff: %e\n",dssplaq_g1-dssplaq1);
     status = 1;
    }

#endif


#ifdef CHECK_SMEAR_GAUGE
 printf("Check difference in QDP s(l) = gs(l)\n"); fflush(stdout);
 rephase(OFF);
 for(i=0; i< RG_Nd; ++i)
  {
   link_qdp[i] = QDP_create_M();
   pr_link[i] = QDP_create_M();
   pr_link_g[i] = QDP_create_M();
   set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
  }

 
 RG_smearing(pr_link,link_qdp,s[nrg],1);
 dssplaq1=plaquette_qdp(pr_link,s[nrg],0);
 for(i=0; i< RG_Nd; ++i)
   set_site_from_M(F_OFFSET(link[i]),pr_link[i],EVENANDODD);


 rand_gauge(F_OFFSET(rgt));
 for(i=0; i< RG_Nd; ++i)
   set_M_from_site(pr_link_g[i],F_OFFSET(link[i]),EVENANDODD);
 dssplaq_g1=plaquette_qdp(pr_link_g,s[nrg],0);
 
 printf("SMEAR-GAUGED %e/%e\n",dssplaq1,dssplaq_g1);
  if (fabs(dssplaq1-dssplaq_g1) > TOL )   
    {
     printf("Error: QDP-SMEAR not gauge invariant \n");
     printf("diff: %e\n",dssplaq1-dssplaq_g1);
     status = 1;
    }

#endif

#ifdef CHECK_SMEAR_GAUGE_2
printf("Check difference in QDP s(l) = sg(l)\n");
 rephase(OFF);
 for(i=0; i< RG_Nd; ++i)
  {
   link_qdp[i] = QDP_create_M();
   link_qdp_g[i] = QDP_create_M();
   pr_link[i] = QDP_create_M();
   pr_link_g[i] = QDP_create_M();
   set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
  }

 
 RG_smearing(pr_link,link_qdp,s[nrg],1);
 SQDP_M_eq_M(pr_link[3],link_qdp[3],s[nrg]);
 dssplaq1=plaquette_qdp(pr_link,s[nrg],0);

 rand_gauge(F_OFFSET(rgt));
 for(i=0; i< RG_Nd; ++i)
   set_M_from_site(link_qdp_g[i],F_OFFSET(link[i]),EVENANDODD);

 RG_smearing(pr_link_g,link_qdp_g,s[nrg],1);
 SQDP_M_eq_M(pr_link_g[3],link_qdp_g[3],s[nrg]);
 dssplaq_g1=plaquette_qdp(pr_link_g,s[nrg],0);
 
 printf("SMEAR-GAUGED TYPE 2 %e/%e\n",dssplaq1,dssplaq_g1);
  if (fabs(dssplaq1-dssplaq_g1) > TOL )   
    {
     printf("Error: QDP-SMEAR not gauge invariant \n");
     printf("diff: %e\n",dssplaq1-dssplaq_g1);
     status = 1;
    }

#endif


#ifdef CHECK_DEGRAND
printf("Check DeGrand trick\n");

 rephase(OFF);
 for(i=0; i< RG_Nd; ++i)
  {
   for(j=0; j<nrg; ++j)
   {
    rg_link_g[j][i] = QDP_create_M();
    rg_link[j][i] = QDP_create_M();
   }
   link_qdp[i] = QDP_create_M();
   set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
  }


  RG_gauge(rg_link,link_qdp,s);

  for(j=0; j<nrg; ++j)
   dssplaq[j]=plaquette_qdp(rg_link[j],s[j+1],nrg-j-1);

  rand_gauge(F_OFFSET(rgt));
  for(i=0; i< RG_Nd; ++i)
   {
    link_qdp_g[i] = QDP_create_M();
    set_M_from_site(link_qdp_g[i],F_OFFSET(link[i]),EVENANDODD);
   }

  RG_gauge(rg_link_g,link_qdp_g,s);
  for(j=0; j<nrg; ++j)
   dssplaq_g[j]=plaquette_qdp(rg_link_g[j],s[j+1],nrg-j-1);


  for(j=0; j<nrg; ++j)
  {
   printf("DEGRAND/GAUGED %e/%e\n",dssplaq[j],dssplaq_g[j]);
   if (fabs(dssplaq[j]-dssplaq_g[j]) > TOL )   
    {
    printf("Error: QDP-SMEAR not gauge invariant \n");
    printf("diff for len %d: %e\n",intpow(2,nrg-j-1),fabs(dssplaq[j]-dssplaq_g[j]));
    status = 1;
    }
  }

#endif

return status;

}


