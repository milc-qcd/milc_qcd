/* NOT MAINTAINED! */
#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include <lattice_qdp.h>
#include "mvmult.c"

static int dslash_setup=0;
static QDP_ColorVector *qsrc, *qdest;
static QDP_ColorVector *temp1[16];//, *temp2[16];
static QDP_ColorVector *tempvec[8];
static QDP_ColorMatrix *tlinks[16];
static QDP_Shift shift[16];
static QDP_ShiftDir sd[16];
QDP_ColorMatrix *bcklink[8];

void
setup_dslash(void)
{
  int i;
  //printf("setup dslash = %i\n", dslash_setup);
  if(dslash_setup==0) {
    dslash_setup = 1;
    for(i=0; i<8; i++) {
      bcklink[i] = QDP_create_M();
    }
#if 0
    qsrc = QDP_create_V();
    qdest = QDP_create_V();
#endif
#if 0
    for(i=0; i<8; i++) {
      tempvec[i] = QDP_create_V();
    }
#endif
#if 0
    for(i=0; i<16; i++) {
      temp1[i] = QDP_create_V();
      temp2[i] = QDP_create_V();
    }
#endif
  }
  for(i=0; i<4; i++) {
    tlinks[4*i]   = implinks[i];
    tlinks[4*i+1] = implinks[i+4];
    tlinks[4*i+2] = bcklink[i];
    tlinks[4*i+3] = bcklink[i+4];
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

void
unset_dslash(void)
{
  //printf("setup dslash = %i\n", dslash_setup);
  if(dslash_setup) {
    int i;
    dslash_setup = 0;
    for(i=0; i<8; i++) {
      QDP_destroy_M(bcklink[i]); bcklink[i] = NULL;
    }
  }
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void
dslash_qdp_fn_special2(QDP_ColorVector *src, QDP_ColorVector *dest,
		       QDP_Subset parity, QDP_ColorVector *temp[])
{
  int i;
  QDP_ColorVector *srcvec[16], *destvec[16];

  for(i=0; i<16; i++) {
    srcvec[i] = src;
    destvec[i] = dest;
  }

  /* Start gathers all directions */
  QDP_V_veq_sV(temp, srcvec, shift, sd, parity, 16);
  //QDP_V_veq_sV(temp, srcvec, shiftdirs, shiftfwd, parity, 8);
  //for(i=0; i<4; i++) {
    //QDP_V_veq_sV(&temp[4*i], &srcvec[4*i], &shift[4*i], &sd[4*i], parity, 4);
  //}
  //for(i=0; i<8; i++) {
  //QDP_V_veq_sV(&temp[2*i], &srcvec[2*i], &shift[2*i], &sd[2*i], parity, 2);
  //}
  //for(i=0; i<16; i++) {
  //QDP_V_eq_sV(temp[i], srcvec[i], shift[i], sd[i], parity);
  //}

  /* multiply by matrix and accumulate */
  QDP_V_eq_M_times_V(destvec[0], tlinks[0], temp[0], parity);
  QDP_V_vpeq_M_times_V(destvec+1, tlinks+1, temp+1, parity, 15);

  //QDP_V_eq_zero(dest, parity);
  //mvmult(destvec, tlinks, temp, parity, 16);

  //QDP_V_eq_zero(dest, parity);
  //QDP_V_vpeq_M_times_V(destvec, tlinks, temp, parity, 8);
  //QDP_V_vpeq_M_times_V(destvec+8, tlinks+8, temp+8, parity, 8);

  //QDP_V_eq_zero(dest, parity);
  //QDP_V_vpeq_M_times_V(destvec, tlinks, temp, parity, 16);

  for(i=0; i<16; ++i) {
    QDP_discard_V(temp[i]);
  }
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void
dslash_qdp_fn_special(QDP_ColorVector *src, QDP_ColorVector *dest,
		      QDP_Subset parity, QDP_ColorVector *temp[])
{
  int i;
  QDP_Subset otherparity;
  QDP_ColorVector *srcvec[8], *destvec[8];

  for(i=0; i<8; i++) {
    tempvec[i] = QDP_create_V();
  }

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
  //for(i=0; i<8; i++) {
  //QDP_V_eq_sV(temp[i], src, shiftdirs[i], QDP_forward, parity);
  //}
  //for(i=0; i<4; i++) {
  //QDP_V_eq_sV(temp[i], src, QDP_neighbor[i], QDP_forward, parity);
  //QDP_V_eq_sV(temp[i+4], src, neighbor3[i], QDP_forward, parity);
  //}

  /* Multiply by adjoint matrix at other sites */
  QDP_V_veq_Ma_times_V(tempvec, implinks, srcvec, otherparity, 8);
  //for(i=0; i<4; i++) {
  //QDP_V_eq_Ma_times_V(tempvec[i], fatlinks[i], src, otherparity);
  //QDP_V_eq_Ma_times_V(tempvec3[i], longlinks[i], src, otherparity);
  //}

  /* Start gathers from negative directions */
  QDP_V_veq_sV(temp+8, tempvec, shiftdirs, shiftbck, parity, 8);
  //for(i=0; i<8; i++) {
  //QDP_V_eq_sV(temp[i+8], tempvec[i], shiftdirs[i], QDP_backward, parity);
  //}
  //for(i=0; i<4; i++) {
  //QDP_V_eq_sV(temp[i+8], tempvec[i], QDP_neighbor[i], QDP_backward, parity);
  //QDP_V_eq_sV(temp[i+12], tempvec[i+4], neighbor3[i], QDP_backward, parity);
  //}

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  QDP_V_eq_zero(dest, parity);
  QDP_V_vpeq_M_times_V(destvec, implinks, temp, parity, 8);
  for(i=0; i<8; ++i) {
    QDP_discard_V(temp[i]);
  }

  /*------------------------------------------------------------*/
  /* Wait gathers from negative directions, accumulate (negative) */

  QDP_V_vmeq_V(destvec, temp+8, parity, 8);
  for(i=0; i<8; ++i) {
    QDP_discard_V(temp[i+8]);
  }
  for(i=0; i<8; i++) {
    QDP_destroy_V(tempvec[i]); tempvec[i] = NULL;
  }
}

void
dslash_qdp_fn(QDP_ColorVector *src, QDP_ColorVector *dest, QDP_Subset parity)
{
  int i;
  for(i=0; i<16; i++) {
    temp1[i] = QDP_create_V();
  }
  //if(parity==QDP_odd)
  //dslash_qdp_fn_special(src, dest, parity, temp2);
  //else
  dslash_qdp_fn_special(src, dest, parity, temp1);
  for(i=0; i<16; i++) {
    QDP_destroy_V(temp1[i]); temp1[i] = NULL;
  }
}

/*********************/
/* old MILC routines */
/*********************/

void
cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
}

void
cleanup_dslash_temps(void)
{
}

void
dslash_fn(field_offset src, field_offset dest, int parity)
{
  QDP_Subset subset;
  //printf("dslash_fn %i\n", parity);
  //setup_dslash();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  qsrc = QDP_create_V();
  qdest = QDP_create_V();
  set_V_from_site(qsrc, src,EVENANDODD);
  set_V_from_site(qdest, dest,EVENANDODD);
  set4_M_from_field(fatlinks, t_fatlink,EVENANDODD);
  set4_M_from_field(longlinks, t_longlink,EVENANDODD);
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
  QDP_destroy_V(qdest); qdest = NULL;
  QDP_destroy_V(qsrc);  qsrc = NULL;
}

/* Special dslash for use by congrad.  Uses restart_gather_site() when
   possible. Last argument is an array of message tags, to be set
   if this is the first use, otherwise reused. If start=1, use
   start_gather_site, otherwise use restart_gather_site.
   The calling program must clean up the gathers! */
void
dslash_fn_special(field_offset src, field_offset dest,
		  int parity, msg_tag **tag, int start)
{
  QDP_Subset subset;
  //printf("dslash_fn_sp %i\n", parity);
  //setup_dslash();
  load_ferm_links();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  qsrc = QDP_create_V();
  qdest = QDP_create_V();
  set_V_from_site(qsrc, src,EVENANDODD);
  set_V_from_site(qdest, dest,EVENANDODD);
  if(start) {
    set4_M_from_field(fatlinks, t_fatlink,EVENANDODD);
    set4_M_from_field(longlinks, t_longlink,EVENANDODD);
  }
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
  QDP_destroy_V(qdest); qdest = NULL;
  QDP_destroy_V(qsrc);  qsrc = NULL;
}

void
dslash_fn_on_temp(su3_vector *src, su3_vector *dest, int parity)
{
  QDP_Subset subset;
  //printf("dslash_fn_t %i\n", parity);
  //setup_dslash();
  load_ferm_links();
  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;
  qsrc = QDP_create_V();
  qdest = QDP_create_V();
  set_V_from_field(qsrc, src,EVENANDODD);
  set_V_from_field(qdest, dest,EVENANDODD);
  set4_M_from_field(fatlinks, t_fatlink,EVENANDODD);
  set4_M_from_field(longlinks, t_longlink,EVENANDODD);
  dslash_qdp_fn(qsrc, qdest, subset);
  set_site_from_V(dest, qdest,EVENANDODD);
  QDP_destroy_V(qdest); qdest = NULL;
  QDP_destroy_V(qsrc);  qsrc = NULL;
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
   possible. Next to last argument is an array of message tags, to be set
   if this is the first use, otherwise reused. If start=1,use
   start_gather_field, otherwise use restart_gather_field.
   The calling program must clean up the gathers and temps! */
void
dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			  int parity, msg_tag **tag, int start)
{
  QDP_Subset qparity;
  //printf("dslash_fn_t_sp %i\n", parity);
  //setup_dslash();
  load_ferm_links();
  if(parity==EVEN) qparity = QDP_even;
  else if(parity==ODD) qparity = QDP_odd;
  else qparity = QDP_all;
  qsrc = QDP_create_V();
  qdest = QDP_create_V();
  set_V_from_field(qsrc, src,EVENANDODD);
  set_V_from_field(qdest, dest,EVENANDODD);
  if(start) {
    set4_M_from_field(fatlinks, t_fatlink,EVENANDODD);
    set4_M_from_field(longlinks, t_longlink,EVENANDODD);
  }
  dslash_qdp_fn(qsrc, qdest, qparity);
  set_site_from_V(dest, qdest,EVENANDODD);
  QDP_destroy_V(qdest); qdest = NULL;
  QDP_destroy_V(qsrc);  qsrc = NULL;
}
