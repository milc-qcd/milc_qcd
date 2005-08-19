#include <stdio.h>
#include <qdp.h>
#include "ks_imp_includes_qdp.h"
#include "RG_include.h"
#define tol  10e-5

int len2;
int r[4];

int hyp_func ( int x[], void *arg)
{
  int i,c[4];
  
  for(i=0;i<QDP_ndim();i++)
    c[i] = (x[i]%len2 - r[i]);

  if ((c[0]==0) && (c[1]==0) && (c[2]==0) && (c[3]==0) ) return 0;
   else return 1; 


}


void RG_transf_field(QDP_ColorVector *dest, QDP_ColorVector *src[RG_Ncn], int len)
{
  int i,j;
  QDP_Subset *QDP_hyp; 
  QDP_Shift offset;


  
  for(i=0;i<RG_Ncn;i++)
  {

   for(j=0;j<RG_Nd;j++)
    r[j] = len*rvs[i][j];

    offset = QDP_create_shift(r);
 
    len2 = 2*len;
    QDP_hyp = QDP_create_subset(hyp_func,NULL,1);
    QDP_V_eq_sV(dest,src[i],offset,QDP_backward,*QDP_hyp); 

    QDP_destroy_shift(offset);
    QDP_destroy_subset(QDP_hyp);

  }
 
return;
}
 
void RG_inv_transf_field(QDP_ColorVector *dest, QDP_ColorVector *src, QDP_ColorMatrix *wlink[RG_Ncn],QDP_Sub_Block s,int len)
{
  int i,j;
  QDP_Shift offset;
  QLA_Real norm = 1.0/16.0;
//  QLA_Real norm = 1.0;
  int rv[4]; 
  QDP_ColorVector *phi;

  
  phi = QDP_create_V();
  QDP_V_eq_zero(dest,QDP_all);

  for (i=0; i<RG_Ncn; i++)
  {

    for(j=0;j<RG_Nd;j++)
     rv[j] = len*rvs[i][j];

    offset = QDP_create_shift(rv);
 
    SQDP_V_eq_M_times_sV(phi,wlink[i],src,offset,QDP_forward,s);
    SQDP_V_peq_r_times_V(dest,&norm,phi,s);

    QDP_destroy_shift(offset);

  }
  
  QDP_destroy_V(phi);

return;
}

void RG_M_inv(QDP_ColorVector *dest,QDP_ColorVector *src)
{
int iter;
QLA_Real qmass,qmass2,qrsqmin,qfinal_rsq;
Real final_rsq; 
QDP_ColorVector *phi_s,*phi_d;
QDP_ColorVector *phi_check,*phi_check1;

  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  phi_s = QDP_create_V();
  phi_d = QDP_create_V();
  
  qmass = (QLA_Real) mass;
  qrsqmin = (QLA_Real) rsqmin;

  QDP_V_eq_zero(phi_s,QDP_all);
  QDP_V_eq_zero(phi_d,QDP_all);


#ifdef CHECK_TRACE
  set_field_from_V(F_OFFSET(ttt),src);
#else
  iter=ks_congrad_qdp(src, phi_s, qmass, niter, qrsqmin, QDP_all,&qfinal_rsq);
  set_field_from_V(F_OFFSET(ttt),phi_s);
#endif 

  dslash_fn_site( F_OFFSET(ttt), F_OFFSET(phi2), EVENANDODD);

#ifdef CHECK_TRACE
  set_V_from_field(dest,F_OFFSET(phi2));
#else
  set_V_from_field(phi_d,F_OFFSET(phi2));
  qmass2 = 2.0*qmass;
  QDP_V_eq_r_times_V_minus_V(dest, &qmass2, phi_s, phi_d, QDP_all);
#endif
  
 
#ifdef CHECK_INV

  //fprintf(stderr,"QDP\n");
  //QDP_V_eq_func(dest,print_cv,QDP_all);

  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  set_field_from_V(F_OFFSET(phi2),dest);
  set_field_from_V(F_OFFSET(ttt),src);
  fprintf(stderr,"Check inversion for QDP\n");
  check_invert( F_OFFSET(phi2), F_OFFSET(ttt), mass, tol);

  phi_check = QDP_create_V();
  phi_check1 = QDP_create_V();

  clear_latvec(F_OFFSET(cg_p) , EVENANDODD );
  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  set_field_from_V(F_OFFSET(ttt),src);
  ks_congrad(F_OFFSET(ttt),F_OFFSET(cg_p),mass,niter,rsqmin,EVENANDODD,&final_rsq );

  dslash_site( F_OFFSET(cg_p), F_OFFSET(phi2), EVENANDODD);
  scalar_mult_add_latvec( F_OFFSET(phi2),F_OFFSET(cg_p),-2.0*mass,F_OFFSET(phi2), EVENANDODD);

  scalar_mult_latvec( F_OFFSET(phi2), -1.0, F_OFFSET(phi2), EVENANDODD );
  set_V_from_field(phi_check,F_OFFSET(phi2));

  //fprintf(stderr,"MILC\n");
  //QDP_V_eq_func(phi_check,print_cv,QDP_all);

  fprintf(stderr,"Check inversion for MILC\n");
  check_invert( F_OFFSET(phi2), F_OFFSET(ttt), mass, tol);


  fprintf(stderr," Check diff: QDP-MILC\n");
  QDP_V_eq_V_minus_V(phi_check1,phi_check,dest,QDP_all);
  QDP_V_eq_func(phi_check1,check_cv,QDP_all);
  
  QDP_destroy_V(phi_check1);
  QDP_destroy_V(phi_check);
 
#endif
  QDP_destroy_V(phi_s);
  QDP_destroy_V(phi_d);


return;
}

void RG_bulk(QDP_ColorVector *dest[RG_Ncn],QDP_ColorVector *src)
{
int j,i,n ;
QDP_ColorVector *phi,*phi_dest;
QDP_Subset *QDP_hyp;


  phi = QDP_create_V();
  phi_dest = QDP_create_V();


 for(i=0;i<RG_Ncn;i++)
 {

   QDP_V_eq_zero(phi,QDP_all);
   QDP_V_eq_zero(phi_dest,QDP_all);

   for(j=0;j<RG_Nd;j++)
    r[j] = rvs[i][j];

   len2 = 2;
   QDP_hyp = QDP_create_subset(hyp_func,NULL,1);
#ifdef CHECK_TRACE_UNIT
   QDP_V_eq_V(dest[i],src,*QDP_hyp);
#else

   QDP_V_eq_V(phi,src,*QDP_hyp);
   QDP_destroy_subset(QDP_hyp);

   RG_M_inv(phi_dest,phi);
   
   QDP_V_eq_V(dest[i],phi_dest,QDP_all);
#endif
  }


  QDP_destroy_V(phi);

return;
}

void RG_check_inversion(QDP_ColorVector *src,QDP_ColorVector *chi[RG_Ncn])
{
int i,iter;
QLA_Real qmass,qmass2,qrsqmin,qfinal_rsq;
QDP_ColorVector *dest,*phi_s,*phi_check,*phi_check1,*phi_d;

  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  phi_s = QDP_create_V();
  phi_d = QDP_create_V();
  dest = QDP_create_V();
  phi_check = QDP_create_V();
  phi_check1 = QDP_create_V();

  QDP_V_eq_zero(phi_check,QDP_all);
  QDP_V_eq_zero(phi_check1,QDP_all);
  QDP_V_eq_zero(phi_s,QDP_all);
  QDP_V_eq_zero(phi_d,QDP_all);
  QDP_V_eq_zero(dest,QDP_all);

  qmass = (QLA_Real) mass;
  qrsqmin = (QLA_Real) rsqmin;

  iter=ks_congrad_qdp(src, phi_s, qmass, niter, qrsqmin, QDP_all,&qfinal_rsq);
  
  set_field_from_V(F_OFFSET(ttt),phi_s);
  dslash_site( F_OFFSET(ttt), F_OFFSET(phi2), EVENANDODD);
  set_V_from_field(phi_d,F_OFFSET(phi2));

  qmass2 = 2.0*qmass;
  QDP_V_eq_r_times_V_minus_V(dest, &qmass2, phi_s, phi_d,QDP_all);

  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  set_field_from_V(F_OFFSET(phi2),dest);
  set_field_from_V(F_OFFSET(ttt),src);
  fprintf(stderr,"Check inversion for QDP entire field\n");
  check_invert( F_OFFSET(phi2), F_OFFSET(ttt), mass, tol);

  for (i=0;i<RG_Ncn;i++)
   QDP_V_peq_V(phi_check,chi[i],QDP_all);

  QDP_V_eq_V_minus_V(phi_check1,dest,phi_check,QDP_all);

  node0_printf("Check split in 16 region: Expected a zero vector\n");
  QDP_V_eq_func(phi_check1,check_cv,QDP_all);
  node0_printf("Check completed\n");

  QDP_destroy_V(phi_s);
  QDP_destroy_V(phi_d);
  QDP_destroy_V(phi_check1);
  QDP_destroy_V(phi_check);

return;
}
