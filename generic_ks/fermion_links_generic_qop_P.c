/****************** fermion_links_generic_qop_P.c *********************/
/* MIMD version 7 */

/* This file collects QOP fermion links utilities shared by asqtad and
   HISQ applications */

/* Note: This is an include file for fermion_links_generic_qop_F.c and
   fermion_links_generic_qop_D.c.  Any edits must be consistent with
   this use. */

/* Entry points (must be defined for the correct precision) */

/* 

   CREATE_QOP_LINKS_FROM_MILC_FN
   CREATE_QOP_FERMION_LINKS
   INVALIDATE_QOP_FERM_LINKS

*/

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_qopmilc.h"
#include "../include/generic_ks_qop.h"

/* Redefinitions according to selected precision */

#if ( QOP_Precision == 1 )

#define MYSU3MATRIX fsu3_matrix
#define MYREAL float
#define CREATE_L_FROM_FIELDS create_F_L_from_fields
#define CREATE_QOP_LINKS_FROM_MILC_FN create_qop_F_links_from_milc_fn
#define CREATE_QOP_FERMION_LINKS create_qop_F_fermion_links
#define DESTROY_QOP_FERMION_LINKS destroy_qop_F_fermion_links
#define INVALIDATE_QOP_FERM_LINKS invalidate_qop_ferm_links_F
#define QOP_L qop_F_l
#define VALID_QOP valid_qop_F
static void destroy_qop_F_fermion_links( qop_fn_links_t *ql );

#else

#define MYSU3MATRIX dsu3_matrix
#define MYREAL double
#define CREATE_L_FROM_FIELDS create_D_L_from_fields
#define CREATE_QOP_LINKS_FROM_MILC_FN create_qop_D_links_from_milc_fn
#define CREATE_QOP_FERMION_LINKS create_qop_D_fermion_links
#define DESTROY_QOP_FERMION_LINKS destroy_qop_D_fermion_links
#define INVALIDATE_QOP_FERM_LINKS invalidate_qop_ferm_links_D
#define QOP_L qop_D_l
#define VALID_QOP valid_qop_D
static void destroy_qop_D_fermion_links( qop_fn_links_t *ql );

#endif

#if FERM_ACTION == HISQ

static void
DESTROY_QOP_FERMION_LINKS( qop_fn_links_t *ql )
{
  if(ql->QOP_L != NULL)
    QOP_hisq_destroy_L(ql->QOP_L);
  ql->QOP_L = NULL;
  ql->VALID_QOP = 0;
}

#else

/***********************************************************************/
/* Create QOP links from MILC fat and long links                       */
/***********************************************************************/
void
CREATE_QOP_LINKS_FROM_MILC_FN(qop_fn_links_t *ql, fn_links_t *fl)
{
  double remaptime;
  //  char myname[] = "create_qop_links_from_milc_fn";

  remaptime = -dclock();

  DESTROY_QOP_FERMION_LINKS(ql);
  ql->QOP_L = CREATE_L_FROM_FIELDS(fl->fat, fl->lng, EVENANDODD);
  ql->VALID_QOP = 1;
  remaptime += dclock();

#ifdef FLTIME
#ifdef REMAP
  node0_printf("FLREMAP:  time = %e\n",remaptime);
#endif
#endif
}

/*********************************************************************/
/* Create QOP_FermionLinksAsqtad object from MILC fat and long links */
/*********************************************************************/
/* We assume that before calling, the fat and long links in fn have
   already been properly constructed */
void
CREATE_QOP_FERMION_LINKS( qop_fn_links_t *ql, fn_links_t *fl )
{
  if(fl->valid != 1){
    printf("create_qop__fermion_links: invalid fl links\n");
    terminate(1);
  }
  /* No need to create them if they are already valid */
  if(ql->VALID_QOP == 1 && ql->QOP_L != NULL){
    /* Debug */
    //node0_printf("create_qop_fermion_links: Reusing qop links\n");
    return;
  }

  CREATE_QOP_LINKS_FROM_MILC_FN(ql, fl);
}

static void
DESTROY_QOP_FERMION_LINKS( qop_fn_links_t *ql )
{
  if(ql->QOP_L != NULL)
    QOP_asqtad_destroy_L(ql->QOP_L);
  ql->QOP_L = NULL;
  ql->VALID_QOP = 0;
}

#endif

void
INVALIDATE_QOP_FERM_LINKS( qop_fn_links_t *ql )
{
  ql->VALID_QOP = 0;
}


