/****************** fermion_links_asqtad_qop_P.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */

/* Note: This is an include file for fermion_links_asqtad_qop_F.c and
   fermion_links_asqtad_qop_D.c.  Any edits must be consistent with
   this use. */

/* Entry points (must be defined for the correct precision) */

/* 

   LOAD_FERM_LINKS
   LOAD_FERM_LINKS_DMDU0
   CREATE_QOP_ASQTAD_FERMION_LINKS
   INVALIDATE_ALL_FERM_LINKS

*/

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_qopmilc.h"
#include "../include/generic_ks_qop.h"

/* Redefinitions according to selected precision */

#if ( QOP_Precision == 1 )

#define MYSU3MATRIX fsu3_matrix
#define MYREAL float
#define LOAD_QOP_ASQTAD_COEFFS load_qop_F_asqtad_coeffs
#define LOAD_FERM_LINKS load_ferm_links_F
#define LOAD_FERM_LINKS_DMDU0 load_ferm_links_dmdu0_F
#define CREATE_L_FROM_FIELDS create_F_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_F_L_from_site_gauge
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_F_asqtad_fermion_links
#define DESTROY_QOP_ASQTAD_FERMION_LINKS destroy_qop_F_asqtad_fermion_links
#define INVALIDATE_ALL_FERM_LINKS invalidate_all_ferm_links_F
#define UNLOAD_L_TO_FIELDS unload_F_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_F_G_from_site
#define DESTROY_RAW4_G destroy_raw4_F_G
#define QOP_L qop_F_l
#define VALID_QOP valid_qop_F
static void destroy_qop_F_asqtad_fermion_links( ferm_links_t *fn );

#else

#define MYSU3MATRIX dsu3_matrix
#define MYREAL double
#define LOAD_QOP_ASQTAD_COEFFS load_qop_D_asqtad_coeffs
#define LOAD_FERM_LINKS load_ferm_links_D
#define LOAD_FERM_LINKS_DMDU0 load_ferm_links_dmdu0_D
#define CREATE_L_FROM_FIELDS create_D_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_D_L_from_site_gauge
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_D_asqtad_fermion_links
#define DESTROY_QOP_ASQTAD_FERMION_LINKS destroy_qop_D_asqtad_fermion_links
#define INVALIDATE_ALL_FERM_LINKS invalidate_all_ferm_links_D
#define UNLOAD_L_TO_FIELDS unload_D_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_D_G_from_site
#define DESTROY_RAW4_G destroy_raw4_D_G
#define QOP_L qop_D_l
#define VALID_QOP valid_qop_D
static void destroy_qop_D_asqtad_fermion_links( ferm_links_t *fn );

#endif

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/***********************************************************************/
/* Create QOP links from MILC fat and long links                       */
/***********************************************************************/
static void
create_qop_links_from_milc_fn(ferm_links_t *fn)
{
  double remaptime;
  char myname[] = "create_qop_links_from_milc";

  remaptime = -dclock();

  DESTROY_QOP_ASQTAD_FERMION_LINKS(fn);
  fn->QOP_L = CREATE_L_FROM_FIELDS(fn->fat, fn->lng, EVENANDODD);
  remaptime += dclock();

#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif
}

#ifdef HAVE_NO_CREATE_L_FROM_G

/***********************************************************************/
/* Create QOP and MILC fat links from MILC gauge field using MILC code */
/***********************************************************************/

/* The fat and long links are created in the prevailing MILC precision
   and their pointers are stored in the structure fn.  The qop links
   are created only in the precision defined by QOP_Precision and the
   pointer is also stored in the structure fn.

   We have two ways to create the fat and long links, depending on the
   level of QOP support.  If QOP is able to create its links from our
   gauge field, we ask it to create them and then unload them into fn.
   If QOP can't, then we use the MILC fattening routines to create
   them in fn and then map them from there into the QOP link
   object. */
   
static void
create_asqtad_links(int both, ferm_links_t *fn, ks_action_paths *ap) {
  Real *act_path_coeff = ap->act_path_coeff;

  double remaptime;
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
  create_qop_links_from_milc_fn(fn);
}

#else

/**********************************************************************/
/* Create QOP and MILC fat links from MILC gauge field using QOP code */
/**********************************************************************/

static QOP_FermionLinksAsqtad *
CREATE_L_FROM_SITE_GAUGE( QOP_info_t *info,
    QOP_asqtad_coeffs_t *coeffs, field_offset src, int parity)
{
  MYSU3MATRIX **raw;
  QOP_FermionLinksAsqtad *qop;
  QOP_GaugeField *gauge;
  raw = CREATE_RAW4_G_FROM_SITE(src, parity);
  if(raw == NULL)terminate(1);
  gauge = QOP_create_G_from_raw((MYREAL **)raw, milc2qop_parity(parity));
  DESTROY_RAW4_G(raw); raw = NULL;
  qop = QOP_asqtad_create_L_from_G(info, coeffs, gauge);
  QOP_destroy_G(gauge); gauge = NULL;
  return qop;
}

static void
create_asqtad_links(int both, ferm_links_t *fn, ks_action_paths *ap) {

  Real *act_path_coeff = ap->act_path_coeff;
  su3_matrix **t_fl = &fn->fat;
  su3_matrix **t_ll = &fn->lng;
  char myname[] = "create_asqtad_links";
  QOP_info_t info = {0., 0., 0, 0, 0};
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

  DESTROY_QOP_ASQTAD_FERMION_LINKS(fn);
  fn->QOP_L = CREATE_L_FROM_SITE_GAUGE( &info, &coeffs, F_OFFSET(link), 
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

  UNLOAD_L_TO_FIELDS( *t_fl, *t_ll, fn->QOP_L, EVENANDODD );

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
}

#endif

/*********************************************************************/
/* Create fat and long links and qop links                           */
/*********************************************************************/
/* Wrappers for MILC call to QOP */
void 
LOAD_FERM_LINKS( ferm_links_t *fn, ks_action_paths *ap ){
  if(fn->valid == 1 && fn->VALID_QOP == 1)return;
  create_asqtad_links(1, fn, ap);

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif
  fn->ap = ap;
  fn->valid = 1;
  fn->VALID_QOP = 1;
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void LOAD_FERM_LINKS_DMDU0( ferm_links_t *fn, ks_action_paths *ap ){

  if(fn->valid == 1)return;
  create_asqtad_links(0, fn, ap);

  /* Currently we don't use the QOP version of links based on dMdu0*/
  DESTROY_QOP_ASQTAD_FERMION_LINKS( fn );
  fn->ap = ap;
  fn->valid = 1;
}
#endif

/*********************************************************************/
/* Create QOP_FermionLinksAsqtad object from MILC fat and long links */
/*********************************************************************/
/* We assume that before calling, the fat and long links in fn have
   already been properly constructed */
void
CREATE_QOP_ASQTAD_FERMION_LINKS( ferm_links_t *fn )
{
  if(fn->valid != 1){
    printf("create_qop_asqtad_fermion_links: invalid fn links\n");
    terminate(1);
  }
  /* No need to create them if they are already valid */
  if(fn->VALID_QOP == 1 && fn->QOP_L != NULL){
    /* Debug */
    //node0_printf("create_qop_asqtad_fermion_links: Reusing qop links\n");
    return;
  }

  create_qop_links_from_milc_fn(fn);
}

static void
DESTROY_QOP_ASQTAD_FERMION_LINKS( ferm_links_t *fn )
{
  if(fn->QOP_L != NULL)
    QOP_asqtad_destroy_L(fn->QOP_L);
  fn->QOP_L = NULL;
  fn->VALID_QOP = 0;
}

void
INVALIDATE_ALL_FERM_LINKS( ferm_links_t *fn )
{
  fn->valid = 0;
  fn->valid_qop_F = 0;
  fn->valid_qop_D = 0;
  if(fn->QOP_L != NULL)
    DESTROY_QOP_ASQTAD_FERMION_LINKS(fn);
}
