/**************** fn_links_milc.c *****************************/
/* MILC Version 7 */
/* Methods for the FN links "class"  */

#include "generic_ks_includes.h"
#include "../include/fn_links.h"
#include <string.h>

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/*-------------------------------------------------------------------*/
/* Special memory allocation for field with 4*su3_matrix per site */
/*-------------------------------------------------------------------*/

static su3_matrix *
create_G_special(void){
  char myname[] = "create_G_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  memset(m, '\0', sites_on_node*4*sizeof(su3_matrix));
  return m;
}

static void
destroy_G_special(su3_matrix *m){
  if(m == NULL)return;
  special_free(m);
}
/*-------------------------------------------------------------------*/
/* Special memory allocations for field with one su3_matrix per site */
/*-------------------------------------------------------------------*/

static su3_matrix *
create_m_special(void){
  char myname[] = "create_m_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));

  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  return m;
}

static void
destroy_m_special(su3_matrix *m){
  special_free(m);
}

/*-------------------------------------------------------------------*/
/* Create/destroy long links                                         */
/*-------------------------------------------------------------------*/
su3_matrix *
create_lnglinks(void) {

  su3_matrix *lng;
  lng = create_G_special();
  return lng;
}
  
/*-------------------------------------------------------------------*/

void 
destroy_lnglinks(su3_matrix *lng){
  if(lng == NULL)return;
  destroy_G_special(lng);
}

/*-------------------------------------------------------------------*/
/* Create/destroy fat links                                         */
/*-------------------------------------------------------------------*/
su3_matrix *
create_fatlinks(void) {

  su3_matrix *fat;
  fat = create_G_special();
  return fat;
}
  
/*-------------------------------------------------------------------*/

void 
destroy_fatlinks(su3_matrix *fat){
  if(fat == NULL)return;
  destroy_G_special(fat);
}

/*-------------------------------------------------------------------*/
/* Create/destroy lngback links                                      */
/*-------------------------------------------------------------------*/

/* Create backward long links from forward long links */

static su3_matrix *
create_lngbacklinks(su3_matrix *lng){
  su3_matrix *lngback;
  register int i;
  int dir;
  su3_matrix *tempmat1;
  msg_tag *tag[4];
  //  char myname[] = "create_lngbacklinks";

  lngback = create_G_special();
  tempmat1 = create_m_special();

  /* gather backwards longlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLFIELDSITES(i){
      tempmat1[i] = lng[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_3_DIR(DIR3(dir)), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLFIELDSITES(i){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i], lngback + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  destroy_m_special(tempmat1);

  return lngback;
}

/*-------------------------------------------------------------------*/

static void 
destroy_lngbacklinks(su3_matrix *lngback){
  if(lngback == NULL)return;
  destroy_G_special(lngback);
}

/*-------------------------------------------------------------------*/
/* Create/destroy fatback links                                      */
/*-------------------------------------------------------------------*/

/* Move up the backward fatlinks.  Result in t_fbl */
static su3_matrix *
create_fatbacklinks(su3_matrix *fat){
  su3_matrix *fatback;
  register int i;
  int dir;
  su3_matrix *tempmat1;
  msg_tag *tag[4];
  // char myname[] = "create_fatbacklinks";

  /* Allocate space for fatback if NULL */
  fatback = create_G_special();
  tempmat1 = create_m_special();

  /* gather backwards fatlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLFIELDSITES(i){
      tempmat1[i] = fat[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLFIELDSITES(i){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i],
      fatback + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  destroy_m_special(tempmat1);

  return fatback;
}

/*-------------------------------------------------------------------*/
static void 
destroy_fatbacklinks(su3_matrix *fatback){
  if(fatback == NULL)return;
  destroy_G_special(fatback);
}

/*-------------------------------------------------------------------*/
/* Load back links                                                   */
/*-------------------------------------------------------------------*/

void destroy_fn_backlinks(fn_links_t *fn){
  if(fn == NULL)return;
  destroy_lngbacklinks(fn->lngback);
  destroy_fatbacklinks(fn->fatback);
  fn->lngback = NULL;
  fn->fatback = NULL;
}


void
load_fn_backlinks(fn_links_t *fn){
  if(fn == NULL)return;
  
  destroy_fn_backlinks(fn);
  fn->lngback = create_lngbacklinks(fn->lng);
  fn->fatback = create_fatbacklinks(fn->fat);
}

/*-------------------------------------------------------------------*/
/* Create/destroy fn links                                           */
/*-------------------------------------------------------------------*/

/* The fat/long members are not created */

fn_links_t *
create_fn_links(void){

  fn_links_t *fn;
  char myname[] = "create_fn_links";
  
  /* Create the structure and allocate the fat and lng members, but
     not the back links.  To create the back links, use
     load_fn_backlinks */

  fn = (fn_links_t *)malloc(sizeof(fn_links_t));
  if(fn == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  fn->phase = create_link_phase_info();
  fn->fat = create_fatlinks();
  fn->lng = create_lnglinks();
  fn->fatback = NULL;
  fn->lngback = NULL;

  return fn;
}

/*-------------------------------------------------------------------*/
void 
destroy_fn_links(fn_links_t *fn){
  if(fn == NULL)return;

  destroy_link_phase_info(fn->phase);
  destroy_fatlinks(fn->fat);
  destroy_lnglinks(fn->lng);
  destroy_fn_backlinks(fn);
  free(fn);
}

/*-------------------------------------------------------------------*/
/* Accessors                                                         */
/*-------------------------------------------------------------------*/

su3_matrix *get_lnglinks(fn_links_t *fn){
  return fn->lng;
}

su3_matrix *get_fatlinks(fn_links_t *fn){
  return fn->fat;
}

su3_matrix *get_lngbacklinks(fn_links_t *fn){
  return fn->lngback;
}

su3_matrix *get_fatbacklinks(fn_links_t *fn){
  return fn->fatback;
}

/*-------------------------------------------------------------------*/
/* Some methods                                                      */
/*-------------------------------------------------------------------*/

/* Copy: fndst = fnsrc */

void 
copy_fn(fn_links_t *fn_src, fn_links_t *fn_dst){
  int i, dir;

  su3_matrix *fatsrc = get_fatlinks(fn_src);
  su3_matrix *lngsrc = get_lnglinks(fn_src);
  su3_matrix *fatbacksrc = get_fatbacklinks(fn_src);
  su3_matrix *lngbacksrc = get_lngbacklinks(fn_src);

  su3_matrix *fatdst = get_fatlinks(fn_dst);
  su3_matrix *lngdst = get_lnglinks(fn_dst);
  su3_matrix *fatbackdst = get_fatbacklinks(fn_dst);
  su3_matrix *lngbackdst = get_lngbacklinks(fn_dst);

  if(fn_src == fn_dst)return;

  FORALLFIELDSITES(i) {
    for(dir=XUP;dir<=TUP;dir++) {
      fatdst[4*i + dir] = fatsrc[4*i + dir];
      lngdst[4*i + dir] = lngsrc[4*i + dir];
      if(fatbackdst != NULL && fatbacksrc != NULL)
	fatbackdst[4*i + dir] = fatbacksrc[4*i + dir];
      if(lngbackdst != NULL && lngbacksrc != NULL)
	lngbackdst[4*i + dir] = lngbacksrc[4*i + dir];
    }
  }
}

/* Multipy by scalar: fndst = fnsrc * s.  OK to do this in place. */

void 
scalar_mult_fn(fn_links_t *fn_src, Real s, fn_links_t *fn_dst){
  int i, dir;

  su3_matrix *fatsrc = get_fatlinks(fn_src);
  su3_matrix *lngsrc = get_lnglinks(fn_src);
  su3_matrix *fatbacksrc = get_fatbacklinks(fn_src);
  su3_matrix *lngbacksrc = get_lngbacklinks(fn_src);

  su3_matrix *fatdst = get_fatlinks(fn_dst);
  su3_matrix *lngdst = get_lnglinks(fn_dst);
  su3_matrix *fatbackdst = get_fatbacklinks(fn_dst);
  su3_matrix *lngbackdst = get_lngbacklinks(fn_dst);

  FORALLFIELDSITES(i) {
    for(dir=XUP;dir<=TUP;dir++) {
      scalar_mult_su3_matrix( fatsrc + 4*i + dir, s, fatdst + 4*i + dir );
      scalar_mult_su3_matrix( lngsrc + 4*i + dir, s, lngdst + 4*i + dir );
      if(fatbacksrc != NULL && fatbackdst != NULL)
	scalar_mult_su3_matrix( fatbacksrc + 4*i + dir, s, fatbackdst + 4*i + dir );
      if(lngbacksrc != NULL && lngbackdst != NULL)
	scalar_mult_su3_matrix( lngbacksrc + 4*i + dir, s, lngbackdst + 4*i + dir );
    }
  }
}

/* Add fnC = fnA + fnB */

void 
add_fn(fn_links_t *fn_A, fn_links_t *fn_B, fn_links_t *fn_C){
  int i, dir;

  su3_matrix *fatA = get_fatlinks(fn_A);
  su3_matrix *lngA = get_lnglinks(fn_A);
  su3_matrix *fatbackA = get_fatbacklinks(fn_A);
  su3_matrix *lngbackA = get_lngbacklinks(fn_A);

  su3_matrix *fatB = get_fatlinks(fn_B);
  su3_matrix *lngB = get_lnglinks(fn_B);
  su3_matrix *fatbackB = get_fatbacklinks(fn_B);
  su3_matrix *lngbackB = get_lngbacklinks(fn_B);

  su3_matrix *fatC = get_fatlinks(fn_C);
  su3_matrix *lngC = get_lnglinks(fn_C);
  su3_matrix *fatbackC = get_fatbacklinks(fn_C);
  su3_matrix *lngbackC = get_lngbacklinks(fn_C);

  FORALLFIELDSITES(i) {
    for(dir=XUP;dir<=TUP;dir++) {
      add_su3_matrix( fatA + 4*i + dir, fatB + 4*i + dir, fatC + 4*i + dir );
      add_su3_matrix( lngA + 4*i + dir, lngB + 4*i + dir, lngC + 4*i + dir );
      if(fatbackA != NULL && fatbackC != NULL)
	add_su3_matrix( fatbackA + 4*i + dir, fatbackB + 4*i + dir, fatbackC + 4*i + dir );
      if(lngbackA != NULL && lngbackC != NULL)
	add_su3_matrix( lngbackA + 4*i + dir, lngbackB + 4*i + dir, lngbackC + 4*i + dir );
    }
  }
}

