/******************** dslash_fn_qdp_P.c ****************************/
/* MILC Version 7 */

/* Note: This is an include file for dslash_fn_qdp_F.c and
   dslash_fn_qdp_D.c so any edits must be consistent with this purpose */

/* WARNING: The FATLINKS and LONGLINKS here appear to require
   a call to ks_congrad_qdp for initialization, which makes it
   misleading to think this dslash can be used independently. */

/* 12/07/06 C. DeTar */

#if ( QDP_Precision == 'F' )

#define SETUP_DSLASH setup_dslash_F
#define DSLASH_QDP_FN_SPECIAL dslash_qdp_F_fn_special
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_F_fn_special2
#define DSLASH_QDP_FN dslash_qdp_F_fn
#define BCKLINK bcklink_F
#define FATLINKS          fatlinks_F
#define LONGLINKS         longlinks_F
#define IMPLINKS          implinks_F

#else

#define SETUP_DSLASH setup_dslash_D
#define DSLASH_QDP_FN_SPECIAL dslash_qdp_D_fn_special
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_D_fn_special2
#define DSLASH_QDP_FN dslash_qdp_D_fn
#define BCKLINK bcklink_D
#define FATLINKS          fatlinks_D
#define LONGLINKS         longlinks_D
#define IMPLINKS          implinks_D

#endif

#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include <lattice_qdp.h>

static int dslash_setup=0;
//static QDP_ColorVector *qsrc, *qdest;
//static QDP_ColorVector *temp1[16], *temp2[16];
static QDP_ColorVector *tempvec[8];
static QDP_ColorMatrix *tlinks[16];
static QDP_Shift shift[16];
static QDP_ShiftDir sd[16];
QDP_ColorMatrix *BCKLINK[8];
QDP_ColorMatrix **FATLINKS, **LONGLINKS, *IMPLINKS[8];

void
SETUP_DSLASH(void)
{
  if(dslash_setup==0) {
    int i;

    for(i=0; i<8; i++) {
      BCKLINK[i] = QDP_create_M();
    }
    for(i=0; i<8; i++) {
      IMPLINKS[i] = QDP_create_M();
    }
    FATLINKS = IMPLINKS;
    LONGLINKS = IMPLINKS + 4;
    //qsrc = QDP_create_V();
    //qdest = QDP_create_V();
#if 0
    for(i=0; i<8; i++) {
      tempvec[i] = QDP_create_V();
    }
    for(i=0; i<16; i++) {
      temp1[i] = QDP_create_V();
      temp2[i] = QDP_create_V();
    }
#endif
    dslash_setup = 1;
    for(i=0; i<4; i++) {
      tlinks[4*i]   = IMPLINKS[i];
      tlinks[4*i+1] = IMPLINKS[i+4];
      tlinks[4*i+2] = BCKLINK[i];
      tlinks[4*i+3] = BCKLINK[i+4];
      shift[4*i]   = QDP_neighbor[i];
      shift[4*i+1] = neighbor3[i];
      shift[4*i+2] = QDP_neighbor[i];
      shift[4*i+3] = neighbor3[i];
      sd[4*i]   = QDP_forward;
      sd[4*i+1] = QDP_forward;
      sd[4*i+2] = QDP_backward;
      sd[4*i+3] = QDP_backward;
    }
  }
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void
DSLASH_QDP_FN_SPECIAL2(QDP_ColorVector *src, QDP_ColorVector *dest,
		       QDP_Subset parity, QDP_ColorVector *temp[])
{
  int i;
  QDP_ColorVector *srcvec[16], *destvec[16];

  for(i=0; i<16; i++) {
    srcvec[i] = src;
    destvec[i] = dest;
  }

  /* Start gathers all directions */
  //QDP_V_veq_sV(temp, srcvec, shiftdirs, shiftfwd, parity, 8);
  //QDP_V_veq_sV(temp, srcvec, shift, sd, parity, 16);
  //for(i=0; i<4; i++) {
  //  QDP_V_veq_sV(&temp[4*i], &srcvec[4*i], &shift[4*i], &sd[4*i], parity, 4);
  //}
  for(i=0; i<8; i++) {
    QDP_V_veq_sV(&temp[2*i], &srcvec[2*i], &shift[2*i], &sd[2*i], parity, 2);
  }
  //for(i=0; i<16; i++) {
  //  QDP_V_eq_sV(temp[i], srcvec[i], shift[i], sd[i], parity);
  //}

  /* multiply by matrix and accumulate */
  //QDP_V_eq_M_times_V(destvec[0], tlinks[0], temp[0], parity);
  //QDP_V_vpeq_M_times_V(destvec+1, tlinks+1, temp+1, parity, 15);

  //QDP_V_eq_zero(dest, parity);
  //QDP_V_vpeq_M_times_V(destvec, tlinks, temp, parity, 8);
  //QDP_V_vpeq_M_times_V(destvec+8, tlinks+8, temp+8, parity, 8);

  QDP_V_eq_zero(dest, parity);
  QDP_V_vpeq_M_times_V(destvec, tlinks, temp, parity, 16);

  for(i=0; i<16; ++i) {
    QDP_discard_V(temp[i]);
  }
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void
DSLASH_QDP_FN_SPECIAL(QDP_ColorVector *src, QDP_ColorVector *dest,
		      QDP_Subset parity, QDP_ColorVector *temp[])
{
  int i;
  QDP_Subset otherparity;
  QDP_ColorVector *srcvec[8], *destvec[8];

  for(i=0; i<8; i++) {
    srcvec[i] = src;
    destvec[i] = dest;
  }

  if(parity==QDP_even) {
    otherparity = QDP_odd;
  } else if(parity==QDP_odd) {
    otherparity = QDP_even;
  } else {
    otherparity = QDP_all;
  }

  /* Start gathers from positive directions */
  QDP_V_veq_sV(temp, srcvec, shiftdirs, shiftfwd, parity, 8);
  //for(i=0; i<4; i++) {
  //QDP_V_eq_sV(temp[i], src, QDP_neighbor[i], QDP_forward, parity);
  //QDP_V_eq_sV(temp[i+4], src, neighbor3[i], QDP_forward, parity);
  //}

  /* Multiply by adjoint matrix at other sites */
  QDP_V_veq_Ma_times_V(tempvec, IMPLINKS, srcvec, otherparity, 8);
  //for(i=0; i<4; i++) {
  //QDP_V_eq_Ma_times_V(tempvec[i], FATLINKS[i], src, otherparity);
  //QDP_V_eq_Ma_times_V(tempvec3[i], LONGLINKS[i], src, otherparity);
  //}

  /* Start gathers from negative directions */
  QDP_V_veq_sV(temp+8, tempvec, shiftdirs, shiftbck, parity, 8);
  //for(i=0; i<4; i++) {
  //QDP_V_eq_sV(temp[i+8], tempvec[i], QDP_neighbor[i], QDP_backward, parity);
  //QDP_V_eq_sV(temp[i+12], tempvec[i+4], neighbor3[i], QDP_backward, parity);
  //}

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  QDP_V_eq_zero(dest, parity);
  QDP_V_vpeq_M_times_V(destvec, IMPLINKS, temp, parity, 8);
  for(i=0; i<8; ++i) {
    QDP_discard_V(temp[i]);
  }

  /*------------------------------------------------------------*/
  /* Wait gathers from negative directions, accumulate (negative) */

  QDP_V_vmeq_V(destvec, temp+8, parity, 8);
  for(i=0; i<8; ++i) {
    QDP_discard_V(temp[i+8]);
  }
}

void
DSLASH_QDP_FN(QDP_ColorVector *src, QDP_ColorVector *dest, QDP_Subset parity)
{
  QDP_ColorVector *vtemp[16];
  int i;
  for(i=0; i<8; i++) {
    tempvec[i] = QDP_create_V();
  }
  for(i=0; i<16; i++) {
    vtemp[i] = QDP_create_V();
  }
  DSLASH_QDP_FN_SPECIAL(src, dest, parity, vtemp);
#if 0
  if(parity==QDP_odd)
    DSLASH_QDP_FN_SPECIAL(src, dest, parity, temp2);
  else
    DSLASH_QDP_FN_SPECIAL(src, dest, parity, temp1);
#endif
  for(i=0; i<16; i++) {
    QDP_destroy_V(vtemp[i]);  vtemp[i] = NULL;
  }
  for(i=0; i<8; i++) {
    QDP_destroy_V(tempvec[i]); tempvec[i] = NULL;
  }
}

/*********************/
/* old MILC routines */
/*********************/

/* commented out for now since they require too much overhead to be useful */
/* DON'T REVIVE WITHOUT MAKING ALL ENTRY POINTS PRECISION-SPECIFIC! */
#if 0

void
cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
}

void
cleanup_dslash_temps(void)
{
}

void
dslash_fn_site(field_offset src, field_offset dest, int parity, ferm_links_t *fn)
{
  QDP_Subset subset;
  //printf("dslash_fn %i\n", parity);
  if(!dslash_setup) SETUP_DSLASH();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  set_V_from_site(qsrc, src,EVENANDODD;
  set_V_from_site(qdest, dest,EVENANDODD;
  if(!fn->fl.valid){
    printf("dslash_fn_site: invalid fn links!\n");
    terminate(1);
  }
  set4_M_from_field(FATLINKS, fn->fl.fat,EVENANDODD, EVENANDODD);
  set4_M_from_field(LONGLINKS, fn->fl.long,EVENANDODD, EVENANDODD);
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
}

/* Special dslash for use by congrad.  Uses restart_gather_site() when
   possible. Last argument is an array of message tags, to be set
   if this is the first use, otherwise reused. If start=1, use
   start_gather_site, otherwise use restart_gather_site.
   The calling program must clean up the gathers! */
void
dslash_fn_site_special(field_offset src, field_offset dest,
		       int parity, msg_tag **tag, int start,
		       ferm_links_t *fn)
{
  QDP_Subset subset;
  //printf("dslash_fn_sp %i\n", parity);
  if(!dslash_setup) SETUP_DSLASH();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  set_V_from_site(qsrc, src,EVENANDODD);
  set_V_from_site(qdest, dest,EVENANDODD);
  if(start) {
    if(!fn->fl.valid){
      printf("dslash_fn_site_special: invalid fn links!\n");
      terminate(1);
    }
    set4_M_from_field(FATLINKS, fn->fl.fat,EVENANDODD);
    set4_M_from_field(LONGLINKS, fn->fl.long,EVENANDODD);
  }
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
}

void
dslash_fn_field(su3_vector *src, su3_vector *dest, int parity,
		ferm_links_t *fn)
{
  QDP_Subset subset;
  //printf("dslash_fn_t %i\n", parity);
  if(!dslash_setup) SETUP_DSLASH();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  set_V_from_field(qsrc, src,EVENANDODD);
  set_V_from_field(qdest, dest,EVENANDODD);
  if(!fn->fl.valid){
    printf("dslash_fn_field: invalid fn links!\n");
    terminate(1);
  }
  set4_M_from_field(FATLINKS, fn->fl.fat,EVENANDODD);
  set4_M_from_field(LONGLINKS, fn->fl.lng,EVENANDODD);
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
   possible. Next to last argument is an array of message tags, to be set
   if this is the first use, otherwise reused. If start=1,use
   start_gather_field, otherwise use restart_gather_field.
   The calling program must clean up the gathers and temps! */
void
dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			int parity, msg_tag **tag, int start,
			ferm_links_t *fn)
{
  QDP_Subset qparity;
  //printf("dslash_fn_t_sp %i\n", parity);
  if(!dslash_setup) SETUP_DSLASH();
  if(parity==EVEN) qparity = QDP_even;
  else if(parity==ODD) qparity = QDP_odd;
  else qparity = QDP_all;
  set_V_from_field(qsrc, src,EVENANDODD);
  set_V_from_field(qdest, dest,EVENANDODD);
  if(start) {
    if(!fn->fl.valid){
      printf("dslash_fn_field_special: invalid fn links!\n");
      terminate(1);
    }
    set4_M_from_field(FATLINKS, fn->fl.fat,EVENANDODD);
    set4_M_from_field(LONGLINKS, fn->fl.long,EVENANDODD);
  }
  dslash_qdp_fn(qsrc, qdest, qparity);
  set_site_from_V(dest, qdest,EVENANDODD);
}

#endif
