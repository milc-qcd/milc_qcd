/****** fermion_links_asqtad_qdp.c  ************************************/
/* MIMD version 7 */
/* Compute fat links and long links using QDP
 *
 * Entry points
 *
 * init_ferm_links
 * load_ferm_links
 * load_ferm_links_dmdu0
 * invalidate_all_ferm_links
 *
 * takes gauge field from site structure "links"
 * If DBLSTORE_FN is defined, also sets double-stored links in the fn_links
 * parameter
 *
 */
/* J. Osborn 2005 */
/* CD 10/06 Wrapped for QDP */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static void
compute_gen_staple(QDP_ColorMatrix *staple, int mu, int nu,
		   QDP_ColorMatrix *link, double dcoef,
		   QDP_ColorMatrix *gauge[], QDP_ColorMatrix *fl[]);

static void 
create_fn_links_qdp(QDP_ColorMatrix *fl[], QDP_ColorMatrix *ll[],
		    QDP_ColorMatrix *gf[], asqtad_path_coeff *coeffs)
{
  
  int i, dir;
  QDP_ColorMatrix *staple, *tempmat1;
  int  nu,rho,sig ;
  QLA_Real one_link;
#ifdef LLTIME
  double nflopfl = 61632;
  double nflopll = 1804;
#endif
  double dtimefl,dtimell;

  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    ll[i] = QDP_create_M();
  }
  staple = QDP_create_M();
  tempmat1 = QDP_create_M();

  dtimefl = -dclock();

  /* to fix up the Lepage term, included by a trick below */
  one_link = coeffs->one_link - 6.0*coeffs->lepage;

  for(dir=0; dir<4; dir++) {
    QDP_M_eq_r_times_M(fl[dir], &one_link, gf[dir], QDP_all);
    for(nu=0; nu<4; nu++) if(nu!=dir) {
      compute_gen_staple(staple, dir, nu, gf[dir],
			 (double)coeffs->three_staple, gf, fl);
      compute_gen_staple(NULL, dir, nu, staple, coeffs->lepage, gf, fl);
      for(rho=0; rho<4; rho++) if((rho!=dir)&&(rho!=nu)) {
	compute_gen_staple(tempmat1, dir, rho, staple,
			   (double)coeffs->five_staple, gf, fl);
	for(sig=0; sig<4; sig++) {
	  if((sig!=dir)&&(sig!=nu)&&(sig!=rho)) {
	    compute_gen_staple(NULL, dir, sig, tempmat1,
			       (double)coeffs->seven_staple, gf, fl);
	  }
	} /* sig */
      } /* rho */
    } /* nu */
  } /* dir */

  dtimell = -dclock();
  dtimefl -= dtimell;
#ifdef LLTIME
  node0_printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtimefl,
         (Real)nflopfl*volume/(1e6*dtimefl*numnodes()) );
#endif

  /* long links */
  for(dir=0; dir<4; dir++) {
    QLA_Real naik = coeffs->naik;
    QDP_M_eq_sM(staple, gf[dir], QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(tempmat1, gf[dir], staple, QDP_all);
    QDP_discard_M(staple);
    QDP_M_eq_sM(staple, tempmat1, QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(ll[dir], gf[dir], staple, QDP_all);
    QDP_M_eq_r_times_M(ll[dir], &naik, ll[dir], QDP_all);
  }
  
  dtimell += dclock();
#ifdef LLTIME
  node0_printf("LLTIME(long): time = %e (Asqtad opt) mflops = %e\n",dtimell,
         (Real)nflopll*volume/(1e6*dtimell*numnodes()) );
#endif

  QDP_destroy_M(staple);
  QDP_destroy_M(tempmat1);
}

/* Computes the staple :
                 mu
              +-------+
        nu    |       |
              |       |
              X       X
  Where the mu link can be any su3_matrix. The result is saved in staple.
  if staple==NULL then the result is not saved.
  It also adds the computed staple to the fatlink[mu] with weight coef.
*/
static void
compute_gen_staple(QDP_ColorMatrix *staple, int mu, int nu,
		   QDP_ColorMatrix *link, double dcoef,
		   QDP_ColorMatrix *gauge[], QDP_ColorMatrix *fl[])
{
  QLA_Real coef = dcoef;
  QDP_ColorMatrix *ts0, *ts1;
  QDP_ColorMatrix *tmat1, *tmat2;
  QDP_ColorMatrix *tempmat;

  ts0 = QDP_create_M();
  ts1 = QDP_create_M();
  tmat1 = QDP_create_M();
  tmat2 = QDP_create_M();
  tempmat = QDP_create_M();

  /* Upper staple */
  QDP_M_eq_sM(ts0, link, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);

  if(staple!=NULL) {  /* Save the staple */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(staple, gauge[nu], tmat1, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(tmat2, gauge[nu], tmat1, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, tmat2, QDP_all);
  }

  /* lower staple */
  QDP_M_eq_sM(ts0, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(tmat1, gauge[nu], link, QDP_all);
  QDP_M_eq_M_times_M(tempmat, tmat1, ts0, QDP_all);
  QDP_M_eq_sM(ts0, tempmat, QDP_neighbor[nu], QDP_backward, QDP_all);

  if(staple!=NULL) { /* Save the staple */
    QDP_M_peq_M(staple, ts0, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, staple, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_peq_r_times_M(fl[mu], &coef, ts0, QDP_all);
  }

  QDP_destroy_M(ts0);
  QDP_destroy_M(ts1);
  QDP_destroy_M(tmat1);
  QDP_destroy_M(tmat2);
  QDP_destroy_M(tempmat);
} /* compute_gen_staple */

void load_asqtad_links(int both, ferm_links_t *fn, ks_action_paths *ap) {

  su3_matrix **t_fl = &fn->fat;
  su3_matrix **t_ll = &fn->lng;
  Real *act_path_coeff = ap->act_path_coeff;

  QDP_ColorMatrix *fl[4];
  QDP_ColorMatrix *ll[4];
  QDP_ColorMatrix *gf[4];
  int dir;
  double remaptime = -dclock();
  char myname[] = "load_asqtad_links";
  
  asqtad_path_coeff c;

  if( phases_in != 1){
    node0_printf("%s: BOTCH: needs phases in\n",myname);
    terminate(1);
  }

  /* Create QDP fields for fat links, long links, and temp for gauge field */
  FORALLUPDIR(dir){
    fl[dir] = QDP_create_M();
    ll[dir] = QDP_create_M();
    gf[dir] = QDP_create_M();
  }

  /* Map gauge links to QDP */
  set4_M_from_site(gf, F_OFFSET(link), EVENANDODD);

  /* Load Asqtad path coefficients from table */
  c.one_link     = act_path_coeff[0]; 
  c.naik         = act_path_coeff[1];
  c.three_staple = act_path_coeff[2];
  c.five_staple  = act_path_coeff[3];
  c.seven_staple = act_path_coeff[4];
  c.lepage       = act_path_coeff[5];

  /* Compute fat and long links as QDP fields */
  remaptime += dclock();
  create_fn_links_qdp(fl, ll, gf, &c);
  remaptime -= dclock();

  /* Clean up */
  FORALLUPDIR(dir){
    QDP_destroy_M(gf[dir]);
  }

  /* Allocate space for t_fatlink if NULL */
  if(*t_fl == NULL){
    *t_fl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fl==NULL){
      printf("%s(%d): no room for t_fatlink\n",myname,this_node);
      terminate(1);
    }
  }
  
  /* Allocate space for t_longlink if NULL and we are doing both fat
     and long */
  if(*t_ll == NULL && both){
    *t_ll = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_ll==NULL){
      printf("%s(%d): no room for t_longlink\n",myname,this_node);
      terminate(1);
    }
  }

  /* Map QDP fields to MILC order */
  set4_field_from_M(*t_fl, fl, EVENANDODD);
  if(both)set4_field_from_M(*t_ll, ll, EVENANDODD);

  /* Clean up */
  FORALLUPDIR(dir){
    QDP_destroy_M(fl[dir]);
    QDP_destroy_M(ll[dir]);
  }
  
  fn->valid = 1;

  remaptime += dclock();
#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif
}


/* Wrappers for MILC call to QDP */
void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap){

  if(fn->valid == 1)return;
#ifdef FN
  load_asqtad_links(1, fn, ap);
#endif
  fn->ap = ap;

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif

  fn->valid = 1;
}

#ifdef DM_DU0
/* Wrappers for MILC call to QDP */
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap){
  su3_matrix *null = NULL;

  if(fn->valid == 1)return;
#ifdef FN
  load_asqtad_links(0, fn, ap);
#endif
  fn->ap = ap;
  fn->valid = 1;
}
#endif

void
invalidate_all_ferm_links(ferm_links_t *fn)
{
  fn->valid = 0;
}

/* For compatibility with fermion_links_hisq.c */
void
invalidate_fn_links(ferm_links_t *fn)
{
  invalidate_all_ferm_links(fn);
}

void 
init_ferm_links(ferm_links_t *fn){
  fn->valid = 0;
  fn->fat = NULL;
  fn->lng = NULL;
  fn->fatback = NULL;
  fn->lngback = NULL;
  fn->ap = NULL;
}
