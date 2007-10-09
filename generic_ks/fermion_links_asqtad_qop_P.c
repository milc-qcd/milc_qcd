/****************** fermion_links_asqtad_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */

/* Note: This is an include file for fermion_links_asqtad_qop_F.c and
   fermion_links_asqtad_qop_D.c.  Any edits must be consistent with
   this use. */

/* Entry points (must be defined for the correct precision) */

/* 

   LOAD_FN_LINKS
   LOAD_FN_LINKS_DMDU0
   CREATE_QOP_ASQTAD_FERMION_LINKS
   DESTROY_QOP_ASQTAD_FERMION_LINKS
   INVALIDATE_FN_LINKS

*/

/* Redefinitions according to selected precision */

#if ( QOP_Precision == 1 )

#define LOAD_QOP_ASQTAD_COEFFS load_qop_F_asqtad_coeffs
#define LOAD_FN_LINKS load_fn_links_F
#define LOAD_FN_LINKS_DMDU0 load_fn_links_dmdu0_F
#define CREATE_L_FROM_FIELDS create_F_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_F_L_from_site_gauge
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_F_asqtad_fermion_links
#define DESTROY_QOP_ASQTAD_FERMION_LINKS destroy_qop_F_asqtad_fermion_links
#define INVALIDATE_FN_LINKS invalidate_fn_links_F
#define UNLOAD_L_TO_FIELDS unload_F_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_F_G_from_site
#define DESTROY_RAW4_G destroy_raw4_F_G

#else

#define LOAD_QOP_ASQTAD_COEFFS load_qop_D_asqtad_coeffs
#define LOAD_FN_LINKS load_fn_links_D
#define LOAD_FN_LINKS_DMDU0 load_fn_links_dmdu0_D
#define CREATE_L_FROM_FIELDS create_D_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_D_L_from_site_gauge
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_D_asqtad_fermion_links
#define DESTROY_QOP_ASQTAD_FERMION_LINKS destroy_qop_D_asqtad_fermion_links
#define INVALIDATE_FN_LINKS invalidate_fn_links_D
#define UNLOAD_L_TO_FIELDS unload_D_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_D_G_from_site
#define DESTROY_RAW4_G destroy_raw4_D_G

#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static QOP_FermionLinksAsqtad *qop_links = NULL;

#ifdef HAVE_NO_CREATE_L_FROM_G

/***********************************************************************/
/* Create QOP and MILC fat links from MILC gauge field using MILC code */
/***********************************************************************/

/* The fat and long links are created in BOTH the structure fn and in
   the global object qop_links. In the structure fn the links are
   created in the prevailing precision.  In qop_links they are created
   in the precision defined by QOP_Precision.

   We have two ways to create the fat and long links, depending on the
   level of QOP support.  If QOP is able to create its links from our
   gauge field, we ask it to create them and then unload them into fn.
   If QOP can't, then we use our own fattening routines to create them
   in fn and then pack them from there into the QOP link object. */
   
static QOP_FermionLinksAsqtad *
create_asqtad_links(int both, fn_links_t *fn, ks_action_paths *ap) {
  Real *act_path_coeff = ap->act_path_coeff;

  double remaptime;
  QOP_FermionLinksAsqtad *qop_l;
  char myname[] = "create_asqtad_links";

  if( phases_in != 1){
    node0_printf("create_asqtad_links: BOTCH: needs phases in\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  /* Use MILC link fattening routines */
  load_fatlinks(fn, ap);
  load_longlinks(fn, ap);

  /* Map to MILC fat and long links to QOP including possible change
     of precision */
  remaptime = -dclock();
  qop_l = CREATE_L_FROM_FIELDS(fn->fat, fn->lng, EVENANDODD);
  remaptime += dclock();

#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif

  return qop_l;
}

#else

/**********************************************************************/
/* Create QOP and MILC fat links from MILC gauge field using QOP code */
/**********************************************************************/

static QOP_FermionLinksAsqtad *
CREATE_L_FROM_SITE_GAUGE( QOP_info_t *info,
    QOP_asqtad_coeffs_t *coeffs, field_offset src, int parity)
{
  su3_matrix **raw;
  QOP_FermionLinksAsqtad *qop;
  QOP_GaugeField *gauge;
  raw = CREATE_RAW4_G_FROM_SITE(src, parity);
  if(raw == NULL)terminate(1);
  gauge = QOP_create_G_from_raw((Real **)raw, milc2qop_parity(parity));
  DESTROY_RAW4_G(raw); raw = NULL;
  qop = QOP_asqtad_create_L_from_G(info, coeffs, gauge);
  QOP_destroy_G(gauge); gauge = NULL;
  return qop;
}

static QOP_FermionLinksAsqtad *
create_asqtad_links(int both, fn_links_t *fn, ks_action_paths *ap) {

  Real *act_path_coeff = ap->act_path_coeff;
  su3_matrix **t_fl = &fn->fat;
  su3_matrix **t_ll = &fn->lng;
  char myname[] = "create_asqtad_links";
  QOP_FermionLinksAsqtad *qop_l;
  QOP_info_t info;
  QOP_asqtad_coeffs_t coeffs;
#ifdef LLTIME
  double nflopfl = 61632;
  double nflopll = 1804;
  double nflop = nflopfl + nflopll;
#endif
  double dtime;
  double remaptime = -dclock();

  LOAD_QOP_ASQTAD_COEFFS(&coeffs, 0.5, act_path_coeff);

  if( phases_in != 1){
    node0_printf("load_fermion_links_fn: BOTCH: needs phases in\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  remaptime += dclock();
  dtime  = -dclock();
  qop_l = CREATE_L_FROM_SITE_GAUGE( &info, &coeffs, F_OFFSET(link), 
					EVENANDODD );
  dtime += dclock();
  remaptime -= dclock();

  /* Allocate space for t_fl if NULL */
  if(*t_fl == NULL){
    *t_fl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fl==NULL){
      printf("%s(%d): no room for t_fl\n",myname,this_node);
      terminate(1);
    }
  }

  /* Allocate space for t_ll if NULL and we are doing both fat and long */
  if(*t_ll == NULL && both){
    *t_ll = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_ll==NULL){
      printf("%s(%d): no room for t_ll\n",myname,this_node);
      terminate(1);
    }
  }

  UNLOAD_L_TO_FIELDS( *t_fl, *t_ll, qop_l, EVENANDODD );

  remaptime += dclock();

#ifdef LLTIME
  node0_printf("LLTIME(total): time = %e (Asqtad opt) mflops = %e\n",dtime,
         (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif

  return qop_l;
}

#endif

/*********************************************************************/
/* Create fat and long links and qop_links                           */
/*********************************************************************/
/* Wrappers for MILC call to QOP */
  // &t_fatlink, get_quark_path_coeff(), get_q_paths());
void 
LOAD_FN_LINKS( fn_links_t *fn, ks_action_paths *ap ){
  if(fn->valid == 1)return;
  if(qop_links != NULL)
    DESTROY_QOP_ASQTAD_FERMION_LINKS();
  qop_links = create_asqtad_links(1, fn, ap);

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif
  fn->valid = 1;
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void LOAD_FN_LINKS_DMDU0( fn_links_t *fn, ks_action_paths *ap ){
  QOP_FermionLinksAsqtad *qop_l;

  if(fn->valid == 1)return;
  qop_l = create_asqtad_links(0, fn, ap);

  /* We don't use the QOP version of links based on dMdu0*/
  QOP_asqtad_destroy_L(qop_l);
  fn->valid = 1;
}
#endif

/*********************************************************************/
/* Create QOP_FermionLinksAsqtad object from MILC fat and long links */
/*********************************************************************/
QOP_FermionLinksAsqtad *
CREATE_QOP_ASQTAD_FERMION_LINKS( fn_links_t *fn, ks_action_paths *ap )
{
  /* If qop_links are unavailable, we have to rebuild them */
  /* They are restored only when we rebuild fat and long links */
  if(qop_links == NULL) 
    INVALIDATE_FN_LINKS(fn);
  /* Create fat and long links if necessary */
  LOAD_FN_LINKS(fn, ap);
  return qop_links;
}

void
DESTROY_QOP_ASQTAD_FERMION_LINKS( void )
{
  QOP_asqtad_destroy_L(qop_links);
  qop_links = NULL;
}

void
INVALIDATE_FN_LINKS( fn_links_t *fn )
{
  fn->valid = 0;
  /* Necessary, since fn->valid does not distinguish QOP_Precision,
     but qop_links does. */
  if(qop_links != NULL)
    DESTROY_QOP_ASQTAD_FERMION_LINKS();
}
