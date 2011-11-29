OBSOLETE

/****************** fermion_links_hisq_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */

/* Note: This is an include file for fermion_links_hisq_qop_F.c and
   fermion_links_hisq_qop_D.c.  Any edits must be consistent with
   this use. */

/* Entry points (must be defined for the correct precision) */

/* 

   LOAD_FERM_LINKS
   LOAD_FERM_LINKS_DMDU0
   CREATE_QOP_HISQ_FERMION_LINKS
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
#define LOAD_QOP_HISQ_COEFFS load_qop_F_hisq_coeffs
#define LOAD_FERM_LINKS load_ferm_links_F
#define LOAD_FERM_LINKS_DMDU0 load_ferm_links_dmdu0_F
#define CREATE_L_FROM_FIELDS create_F_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_F_L_from_site_gauge
#define CREATE_QOP_HISQ_FERMION_LINKS create_qop_F_hisq_fermion_links
#define DESTROY_QOP_HISQ_FERMION_LINKS destroy_qop_F_hisq_fermion_links
#define INVALIDATE_ALL_FERM_LINKS invalidate_all_ferm_links_F
#define UNLOAD_HISQ_L_TO_FIELDS unload_F_hisq_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_F_G_from_site
#define DESTROY_RAW4_G destroy_raw4_F_G
#define QOP_L qop_F_l
#define VALID_QOP valid_qop_F
static void destroy_qop_F_hisq_fermion_links( ferm_links_t *fn );
#define LOOK_AT_LINK look_at_link_F

#else

#define MYSU3MATRIX dsu3_matrix
#define MYREAL double
#define LOAD_QOP_HISQ_COEFFS load_qop_D_hisq_coeffs
#define LOAD_FERM_LINKS load_ferm_links_D
#define LOAD_FERM_LINKS_DMDU0 load_ferm_links_dmdu0_D
#define CREATE_L_FROM_FIELDS create_D_L_from_fields
#define CREATE_L_FROM_SITE_GAUGE create_D_L_from_site_gauge
#define CREATE_QOP_HISQ_FERMION_LINKS create_qop_D_hisq_fermion_links
#define DESTROY_QOP_HISQ_FERMION_LINKS destroy_qop_D_hisq_fermion_links
#define INVALIDATE_ALL_FERM_LINKS invalidate_all_ferm_links_D
#define UNLOAD_HISQ_L_TO_FIELDS unload_D_hisq_L_to_fields
#define CREATE_RAW4_G_FROM_SITE create_raw4_D_G_from_site
#define DESTROY_RAW4_G destroy_raw4_D_G
#define QOP_L qop_D_l
#define VALID_QOP valid_qop_D
static void destroy_qop_D_hisq_fermion_links( ferm_links_t *fn );
#define LOOK_AT_LINK look_at_link_D

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
  char myname[] = "create_qop_links_from_milc_fn";

  //AB THIS ROUTINE IS NOT SUPPOSED TO BE CALLED SO FAR
  printf("%s is not supposed to be called\n",myname);
  terminate(1);

  remaptime = -dclock();

  //AB DESTROY_QOP_HISQ_FERMION_LINKS(fn);
  //AB fn->QOP_L = CREATE_L_FROM_FIELDS(fn->fat, fn->lng, EVENANDODD);
  remaptime += dclock();

#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif
}

// #ifdef HAVE_NO_CREATE_L_FROM_G
// 
// /***********************************************************************/
// /* Create QOP and MILC fat links from MILC gauge field using MILC code */
// /***********************************************************************/
// 
// /* The fat and long links are created in the prevailing MILC precision
//    and their pointers are stored in the structure fn.  The qop links
//    are created only in the precision defined by QOP_Precision and the
//    pointer is also stored in the structure fn.
// 
//    We have two ways to create the fat and long links, depending on the
//    level of QOP support.  If QOP is able to create its links from our
//    gauge field, we ask it to create them and then unload them into fn.
//    If QOP can't, then we use the MILC fattening routines to create
//    them in fn and then map them from there into the QOP link
//    object. */
//    
// static void
// create_hisq_links(int both, ferm_links_t *fn) {
//   ks_action_paths_hisq *ap = fn->ap;
//   ks_component_paths p1 = ap->p1;
// 
//   //AB NO SUPPORT FOR SUCH SITUATION YET
//   node0_printf("HAVE_NO_CREATE_L_FROM_G macro is not supported for HISQ\n");
//   terminate(1);
// 
//   double remaptime;
//   char myname[] = "create_hisq_links";
// 
//   if( phases_in != 1){
//     node0_printf("create_hisq_links: BOTCH: needs phases in\n");
//     terminate(1);
//   }
// 
//   /* Initialize QOP */
//   if(initialize_qop() != QOP_SUCCESS){
//     printf("%s(%d): Error initializing QOP\n",myname,this_node);
//     terminate(1);
//   }
// 
//   /* Use MILC link fattening routines */
//   load_fatlinks(fn, p1);
//   load_longlinks(fn, ap);
// 
//   /* Map to MILC fat and long links to QOP including possible change
//      of precision */
//   create_qop_links_from_milc_fn(fn);
// }
// 
// #else

/**********************************************************************/
/* Create QOP and MILC fat links from MILC gauge field using QOP code */
/**********************************************************************/

int
CREATE_L_FROM_SITE_GAUGE( QOP_info_t *info,
    QOP_hisq_coeffs_t *coeffs, field_offset src, int parity,
    QOP_FermionLinksHisq *qop, int new_naik_set)
{
  int hisq_flag;
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter CREATE_L_FROM_SITE_GAUGE in fermion_links_hisq_qop_P.c\n");
#endif
  //AB A PROBLEM HERE IS THAT IF inaik IS SUCH THAT WE DO NOT NEED
  //   TO CALCULATE ANYTHING BUT SIMPLY NEED TO SET POINTERS
  //   TO DIFFERENT ARRAYS IT WOULD BE A WASTE OF TIME TO DO ALL THIS
  //   CONVERSIONS FROM SITE TO RAW
  MYSU3MATRIX **raw;
  QOP_GaugeField *gauge;
  raw = CREATE_RAW4_G_FROM_SITE(src, parity);
  if(raw == NULL)terminate(1);
  gauge = QOP_create_G_from_raw((MYREAL **)raw, milc2qop_parity(parity));
  DESTROY_RAW4_G(raw); raw = NULL;

#ifdef HISQ_SVD_COUNTER
  QOP_info_hisq_svd_counter(info) = hisq_svd_counter;
#endif

  hisq_flag = QOP_hisq_create_L_from_G(info, coeffs, gauge, qop, new_naik_set);

#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = QOP_info_hisq_svd_counter(info);
#endif

  QOP_destroy_G(gauge); gauge = NULL;

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  CREATE_L_FROM_SITE_GAUGE in fermion_links_hisq_qop_P.c\n");
#endif
  return hisq_flag;
}

static void
create_hisq_links(int both, ferm_links_t *fn) {
  ks_action_paths *ap = fn->ap;
  su3_matrix **t_fl = &fn->fl.fat;
  su3_matrix **t_ll = &fn->fl.lng;
  char myname[] = "create_hisq_links";
  QOP_info_t info;
  QOP_hisq_coeffs_t coeffs;
  int hisq_flag;
#ifdef LLTIME
//  double nflopfl = 61632;
//  double nflopll = 1804;
  double nflopfl = 0; //AB NEED TO COUNT THE FLOPS EVENTUALLY
  double nflopll = 0;
  double nflop = nflopfl + nflopll;
#endif
  double dtime;
  double remaptime = -dclock();
  MYREAL epsilon_naik[MAX_NAIK];
  int i;

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter create_hisq_links in fermion_links_hisq_qop_P.c\n");fflush(stdout);
#endif

  LOAD_QOP_HISQ_COEFFS(&coeffs, 0.5, ap);

  if( phases_in != 1){
    node0_printf("%s: BOTCH: needs phases in\n",myname);
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  remaptime += dclock();
  dtime  = -dclock();

  //AB Allocate the FermionLinksHisq structure if it is NULL.
  //   Presently it is never destroyed since flags inside this
  //   structure are needed for fat/long links routines to work
  if( fn->ql.QOP_L==NULL ) {
    // remap array with epsilons to QOP precision
    for( i=0; i<(fn->hl).n_naiks; i++) epsilon_naik[i]=(fn->hl).eps_naik[i];
    // allocate FermionLinksHisq structure
    fn->ql.QOP_L=QOP_allocate_hisq_fermion_links( (fn->hl).n_naiks, epsilon_naik );
  }

  //AB SO FAR THIS FLAG IS SET TO VALID AND NEVER CHANGED
  fn->ql.VALID_QOP=1;

  //AB CURRENTLY WE DO NOT DESTROY ANYTHING
  //AB DESTROY_QOP_HISQ_FERMION_LINKS(fn);


  hisq_flag=CREATE_L_FROM_SITE_GAUGE( &info, &coeffs, F_OFFSET(link), 
    EVENANDODD, fn->ql.QOP_L, (fn->hl).current_X_set );


  //AB AT THIS STAGE QOP EITHER CREATED NEW LINKS OR RESET POINTERS
  //   TO THE NEEDED SET; AT THE MOMENT WE WANT MILC FAT/LONG LINKS
  //   TO BE SYNCHRONIZED WITH QOP LINKS, THUS, WE FOLLOW THE SAME
  //   PROCEDURE AS QOP: CHECK THE FLAGS AND IF NEEDED UNLOAD QOP
  //   LINKS INTO MILC, OR SIMPLY CHANGE POINTERS
  //   FOR THE MILC CODE REORDERING (FOR THE INVERTER) AND STORAGE
  //   OF BACK LINKS ARE STILL REQUIRED

#ifdef AB_NEEDED_FOR_DEBUG
  switch(hisq_flag) {
    case QOP_HISQ_NEW_LINKS: // new links created, need to unload
      printf("**** Need to unload into milc links\n");
      break;
    case QOP_HISQ_CHANGED_POINTER: // keep links but change pointers
      printf("**** Need to change pointers\n");
      break;
    case QOP_HISQ_DONE_NOTHING: // do nothing
    default:
      ;
      printf("**** Do nothing with fat/long links\n");
  }
#endif


  dtime += dclock();
  remaptime -= dclock();

#ifdef THIS_SHOULD_NOT_BE_NEEDED_SINCE_ALL_SPACE_IS_ALLOCATED_IN_load_ferm_links_WRAPPER
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
#endif

  //AB NEED TO SORT OUT HOW TO UNLOAD FROM QOP TO MILC
  UNLOAD_HISQ_L_TO_FIELDS( *t_fl, *t_ll, fn->ql.QOP_L, EVENANDODD );
//  QOP_hisq_extract_L_to_raw(fn->hl.U_link, fn->hl.V_link, fn->hl.Y_unitlink, 
//                            fn->hl.W_unitlink, fn->hl.XX_fatlink, fn->hl.XX_longlink,
//			    fn->ql.QOP_L, EVENANDODD);

  remaptime += dclock();

#ifdef LLTIME
  node0_printf("LLTIME(total): time = %e (Hisq opt) mflops = %e\n",dtime,
         (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

#ifdef LLTIME
#ifdef REMAP
  node0_printf("LLREMAP:  time = %e\n",remaptime);
#endif
#endif

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  create_hisq_links in fermion_links_hisq_qop_P.c\n");
#endif
}

//#endif

/*********************************************************************/
/* Create fat and long links and qop links                           */
/*********************************************************************/
/* Wrappers for MILC call to QOP */
void 
LOAD_FERM_LINKS( ferm_links_t *fn ){
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter LOAD_FERM_LINKS in fermion_links_hisq_qop_P.c\n");
#endif

  //AB NEED MORE SOPHISTICATED CHECKS THAN THESE
//  if(fn->fl.valid == 1 && fn->ql.VALID_QOP == 1)return;

  create_hisq_links(1, fn);

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif
  fn->fl.valid = 1;
  fn->ql.VALID_QOP = 1;

//AB FOR DEBUGGING OUTPUT GIVEN LINK
//{
//  int x[4]={0,0,0,0};
//  LOOK_AT_LINK(fn,x,0);
//}

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  LOAD_FERM_LINKS in fermion_links_hisq_qop_P.c\n");
#endif
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void LOAD_FERM_LINKS_DMDU0( ferm_links_t *fn ){

  if(fn->fl.valid == 1)return;
  create_hisq_links(0, fn);

  /* Currently we don't use the QOP version of links based on dMdu0*/
  DESTROY_QOP_HISQ_FERMION_LINKS( fn );
  fn->fl.valid = 0;
}
#endif

/*********************************************************************/
/* Create QOP_FermionLinksHisq object from MILC fat and long links */
/*********************************************************************/
/* We assume that before calling, the fat and long links in fn have
   already been properly constructed */
void
CREATE_QOP_HISQ_FERMION_LINKS( ferm_links_t *fn )
{
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter CREATE_QOP_HISQ_FERMION_LINKS in fermion_links_hisq_qop_P.c\n");
#endif
  if(fn->fl.valid != 1){
    printf("create_qop_hisq_fermion_links: invalid fn links\n");
    terminate(1);
  }
  /* No need to create them if they are already valid */
  if(fn->ql.VALID_QOP == 1 && fn->ql.QOP_L != NULL){
    /* Debug */
    //node0_printf("create_qop_hisq_fermion_links: Reusing qop links\n");
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  CREATE_QOP_HISQ_FERMION_LINKS in fermion_links_hisq_qop_P.c  -- NO CHANGE\n");
#endif
    return;
  }

  printf("*** ERROR: Program flow should not enter create_qop_links_from_milc_fn in CREATE_QOP_HISQ_FERMION_LINKS in fermion_links_hisq_qop_P.c\n");
  create_qop_links_from_milc_fn(fn);
}

static void
DESTROY_QOP_HISQ_FERMION_LINKS( ferm_links_t *fn )
{
  if(fn->ql.QOP_L != NULL)
    QOP_hisq_destroy_L(fn->ql.QOP_L);
  fn->ql.QOP_L = NULL;
  fn->ql.VALID_QOP = 0;
}

void
INVALIDATE_ALL_FERM_LINKS( ferm_links_t *fn )
{
  fn->fl.valid = 0;
  fn->ql.valid_qop_F = 0;
  fn->ql.valid_qop_D = 0;
  if(fn->ql.QOP_L != NULL)
    DESTROY_QOP_HISQ_FERMION_LINKS(fn);
}

void
LOOK_AT_LINK( ferm_links_t *fn, int *x, int dir )
{
  QOP_hisq_look_at_link( fn->ql.QOP_L, x, dir );
}
