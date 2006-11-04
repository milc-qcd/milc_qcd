/**************** fermion_links_helpers.c *****************************/
/* MILC Version 7 */
/* Routines used by all of the fermion_links*.c files */

#include "generic_ks_includes.h"	/* definitions files and prototypes */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/* Move up the backward longlinks.  Result in t_lbl */
void 
load_longbacklinks(su3_matrix **t_lbl, su3_matrix *t_ll){
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_longbacklinks";

  /* Allocate space for t_lbl if NULL */
  if(*t_lbl == NULL){
    *t_lbl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_lbl==NULL){
      printf("%s(%d): no room for t_lbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards longlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_ll[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_3_DIR(DIR3(dir)), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i],
      (*t_lbl) + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

/* Move up the backward fatlinks.  Result in t_fbl */
void 
load_fatbacklinks(su3_matrix **t_fbl, su3_matrix *t_fl){
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_fatbacklinks";

  /* Allocate space for t_fbl if NULL */
  if(*t_fbl == NULL){
    *t_fbl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fbl==NULL){
      printf("%s(%d): no room for t_fbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards fatlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_fl[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i],
      (*t_fbl) + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

static void 
free_t_links(su3_matrix **t_l){
  if(*t_l != NULL) special_free(*t_l);
  *t_l = NULL;
}

/* Wrappers for MILC call to QOP */
void 
free_fn_links(){
  free_t_links(&t_fatlink);
  free_t_links(&t_longlink);
#ifdef DBLSTORE_FN
  free_t_links(&t_fatbacklink);
  free_t_links(&t_longbacklink);
#endif

  valid_fn_links = 0;
}

#ifdef DM_DU0
/* Routines for dDslash/du0 */

void free_fn_links_dmdu0(){
  free_t_links(&t_dfatlink_du0);
  valid_fn_links_dmdu0 = 0;
}
#endif

