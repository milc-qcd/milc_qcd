/************************* RG_transf.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
#include <stdio.h>
#include <qdp.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"
#define tol  10e-5

int len2;
int r[4];

/* On a lattice with lattice constant len2/2, select the subset of
   sites in each hypercube with hypercube offset r, where the
   components of r are 0 or len2/2.
*/
int hyp_func ( int x[], void *arg)
{
  int i,c[4];
  
  for(i=0;i<QDP_ndim();i++)
    c[i] = (x[i]%len2 - r[i]);

  if ((c[0]==0) && (c[1]==0) && (c[2]==0) && (c[3]==0) ) return 0;
   else return 1; 


}

/* On a lattice with lattice constant len, completes the parallel
   transport of the field from the origin of the hypercube to
   the hypercube site with offset r */
void RG_transf_field(QDP_ColorVector *dest, QDP_ColorVector *src[RG_Ncn], int len)
{
  int i,j;
  QDP_Subset *QDP_hyp; 
  QDP_Shift offset;


  
  for(i=0;i<RG_Ncn;i++)
  {

    /* Shift displacement */
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
 
/* For a lattice with spacing len, parallel transport the hypercube
   values in src to the hypercube origin and average them */
void RG_inv_transf_field(QDP_ColorVector *dest, 
			 QDP_ColorVector *src, 
			 QDP_ColorMatrix *wlink[RG_Ncn],
			 QDP_Sub_Block s,
			 int len)
{
  int i,j;
  QDP_Shift offset;
  QLA_Real norm = 1.0/16.0;
  //  QLA_Real norm = 1.0;
  int rv[4]; 
  QDP_ColorVector *phi;
  
  
  phi = QDP_create_V();
  QDP_V_eq_zero(dest,QDP_all);
  
  /* For all sites in hypercube with side len */
  for (i=0; i<RG_Ncn; i++)
    {
      
      for(j=0;j<RG_Nd;j++)
	rv[j] = len*rvs[i][j];
      
      offset = QDP_create_shift(rv);

      /* Parallel transport the src on the hypercube to the 
	 hypercube origin */
#ifndef NOTRANS
      /* Do parallel transport */
      SQDP_V_eq_M_times_sV(phi,wlink[i],src,offset,QDP_forward,s);
#else
      /* If no parallel transport, just shift */
      SQDP_V_eq_sV(phi,src,offset,QDP_forward,s);
#endif
      /* Sum result in dest */
      SQDP_V_peq_r_times_V(dest,&norm,phi,s);
      
      QDP_destroy_shift(offset);
      
    }
  
  QDP_destroy_V(phi);
  
  return;
}

/* Compute dest = M^{-1} src */
void RG_M_inv(QDP_ColorVector *dest, QDP_ColorVector *src)
{
  int iter;
  QLA_Real qmass,qmass2,qrsqmin,qfinal_rsq;
  Real final_rsq; 
  QDP_ColorVector *phi_s,*phi_d;
  QDP_ColorVector *phi_check,*phi_check1;
  quark_invert_control qic;
  
  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );
  
  phi_s = QDP_create_V();
  phi_d = QDP_create_V();
  
  qmass = (QLA_Real) mass;
  qrsqmin = (QLA_Real) rsqprop;
  
  QDP_V_eq_zero(phi_s,QDP_all);
  QDP_V_eq_zero(phi_d,QDP_all);
  
  
#ifdef CHECK_TRACE
  set_site_from_V(F_OFFSET(ttt),src,EVENANDODD);
#else
  /* Do phi_s = (M^\dagger M)^{-1} src */
  /* Pack structure */
  qic.prec      = PRECISION;
  qic.min       = 0;
  qic.max       = niter;
  qic.nrestart  = 5;
  qic.parity    = EVENANDODD;
  qic.start_flag = 0;
  qic.nsrc = 1;
  qic.resid     = sqrt(qrsqmin);
  qic.relresid  = 0;     /* Suppresses this test */
  load_ferm_links(&fn_links);
#if ( QDP_Precision == 'F' )
  qic.prec      = 1;
  iter=ks_congrad_qdp_F(src, phi_s, &qic, qmass, &fn_links);
#else
  qic.prec      = 2;
  iter=ks_congrad_qdp_D(src, phi_s, &qic, qmass, &fn_links);
#endif
  /* Then do dest = M^\dagger phi_s */
  set_site_from_V(F_OFFSET(ttt),phi_s,EVENANDODD);
#endif 
  
  load_ferm_links(&fn_links);
  dslash_fn_site( F_OFFSET(ttt), F_OFFSET(phi2), EVENANDODD, &fn_links);
  
#ifdef CHECK_TRACE
  set_V_from_site(dest,F_OFFSET(phi2),EVENANDODD);
#else
  set_V_from_site(phi_d,F_OFFSET(phi2),EVENANDODD);
  qmass2 = 2.0*qmass;
  QDP_V_eq_r_times_V_minus_V(dest, &qmass2, phi_s, phi_d, QDP_all);
#endif
  
  
#ifdef CHECK_INV
  
  //fprintf(stderr,"QDP\n");
  //QDP_V_eq_func(dest,print_cv,QDP_all);
  
  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );
  
  set_site_from_V(F_OFFSET(phi2),dest,EVENANDODD);
  set_site_from_V(F_OFFSET(ttt),src,EVENANDODD);
  fprintf(stderr,"Check inversion for QDP\n");
  check_invert( F_OFFSET(phi2), F_OFFSET(ttt), mass, tol);
  
  phi_check = QDP_create_V();
  phi_check1 = QDP_create_V();
  
  clear_latvec(F_OFFSET(cg_p) , EVENANDODD );
  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );
  
  set_site_from_V(F_OFFSET(ttt),src,EVENANDODD);
  load_ferm_links(&fn_links);
  ks_congrad(F_OFFSET(ttt),F_OFFSET(cg_p),mass,niter,rsqprop,PRECISION,
	     EVENANDODD,&final_rsq, &fn_links );
  
  load_ferm_links(&fn_links);
  dslash_site( F_OFFSET(cg_p), F_OFFSET(phi2), EVENANDODD, &fn_links);
  scalar_mult_add_latvec( F_OFFSET(phi2),F_OFFSET(cg_p),-2.0*mass,F_OFFSET(phi2), EVENANDODD);
  
  scalar_mult_latvec( F_OFFSET(phi2), -1.0, F_OFFSET(phi2), EVENANDODD );
  set_V_from_site(phi_check,F_OFFSET(phi2),EVENANDODD);
  
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

/* Starting from the src on the fine lattice, restricted to each hypercube
   offset i, multiply by M^{-1} to get dest[i] */
void RG_bulk(QDP_ColorVector *dest[RG_Ncn], QDP_ColorVector *src)
{
  int j,i,n ;
  QDP_ColorVector *phi,*phi_dest;
  QDP_Subset *QDP_hyp;
  
  
  phi = QDP_create_V();
  phi_dest = QDP_create_V();
  
  
  /* Iterate over hypercube offsets */
  for(i=0;i<RG_Ncn;i++)
    {
      
      QDP_V_eq_zero(phi,QDP_all);
      QDP_V_eq_zero(phi_dest,QDP_all);
      
      /* Define the subset of sites with hypercube offset r[i] 
	 and let phi be the restriction of the source src on
	 those offsets */
      for(j=0;j<RG_Nd;j++)
	r[j] = rvs[i][j];
      
      len2 = 2;
      QDP_hyp = QDP_create_subset(hyp_func,NULL,1);
#ifdef CHECK_TRACE_UNIT
      QDP_V_eq_V(dest[i],src,*QDP_hyp);
#else
      
      QDP_V_eq_V(phi,src,*QDP_hyp);
      QDP_destroy_subset(QDP_hyp);
      
      /* Compute phi_dest = M^{-1} phi */
      RG_M_inv(phi_dest,phi);
      
      /* Copy the result to dest[i] */
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
 quark_invert_control qic;

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
  qrsqmin = (QLA_Real) rsqprop;

  /* Pack structure */
  qic.parity    = EVENANDODD;
  qic.max       = niter;
  qic.nrestart  = 5;
  qic.resid     = sqrt(qrsqmin);
  qic.relresid  = 0;     /* Suppresses this test */
  load_ferm_links(&fn_links);
#if ( QDP_Precision == 'F' )
  qic.prec      = 1;
  iter=ks_congrad_qdp_F(src, phi_s, &qic, qmass, &fn_links );
#else
  qic.prec      = 2;
  iter=ks_congrad_qdp_D(src, phi_s, &qic, qmass, &fn_links );
#endif
  
  set_site_from_V(F_OFFSET(ttt),phi_s,EVENANDODD);
  load_ferm_links(&fn_links);
  dslash_site( F_OFFSET(ttt), F_OFFSET(phi2), EVENANDODD, &fn_links);
  set_V_from_site(phi_d,F_OFFSET(phi2),EVENANDODD);

  qmass2 = 2.0*qmass;
  QDP_V_eq_r_times_V_minus_V(dest, &qmass2, phi_s, phi_d,QDP_all);

  clear_latvec(F_OFFSET(phi2) , EVENANDODD );
  clear_latvec(F_OFFSET(ttt) , EVENANDODD );

  set_site_from_V(F_OFFSET(phi2),dest,EVENANDODD);
  set_site_from_V(F_OFFSET(ttt),src,EVENANDODD);
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
