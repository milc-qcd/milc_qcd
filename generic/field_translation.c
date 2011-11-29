/************************** field_translation.c *****************************/
/* MIMD version 7 */

/* Utilities for translating fields by a fixed displacement */
/* No parallel transport here */

/* Translate the gauge field by a fixed coordinate displacement rshift */

#include "generic_includes.h"
#include "../include/generic_wilson.h"

#ifndef NO_GAUGE_FIELD

/*------------------------------------------------------------------*/
void shift_gauge(int rshift[]){

  msg_tag *tag;
  int dir, doshift;
  int ndim[4];
  int i;
  site *s;
  su3_matrix *old_link;

  ndim[XUP] = nx; ndim[YUP] = ny; ndim[ZUP] = nz; ndim[TUP] = nt;

  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = rshift[dir] % ndim[dir];
    if(rshift[dir] != 0)doshift = 1;
  }

  /* Don't do a null shift */
  if(!doshift)return;

  old_link = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(old_link == NULL){
    printf("shift_gauge(%d): No room for old_link\n",this_node);
  }

  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      su3mat_copy( &s->link[dir], old_link+i);
    }

    tag = start_general_gather_field( old_link, sizeof(su3_matrix), 
				      rshift, EVENANDODD, gen_pt[0]);
    wait_general_gather(tag);

    FORALLSITES(i,s){
      su3mat_copy((su3_matrix *)gen_pt[0][i], &s->link[dir]);
    }

    cleanup_general_gather(tag);
  }

  free(old_link);
}

#endif

/*------------------------------------------------------------------*/
/* Translate a complex field by a fixed coordinate displacement rshift */

void shift_complex(complex *src, int rshift[]){

  msg_tag *tag;
  int dir, doshift;
  int ndim[4];
  int i;
  site *s;
  complex *tmp;
  complex *c;

  ndim[XUP] = nx; ndim[YUP] = ny; ndim[ZUP] = nz; ndim[TUP] = nt;

  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = rshift[dir] % ndim[dir];
    if(rshift[dir] != 0)doshift = 1;
  }

  /* Don't do a null shift */
  if(!doshift)return;

  tmp = create_c_field();
  copy_c_field(tmp, src);

  tag = start_general_gather_field( tmp, sizeof(complex), 
				    rshift, EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);

  FORALLSITES(i,s){
    c = (complex *)(gen_pt[0][i]);
    src[i].real = c->real;
    src[i].imag = c->imag;
  }

  cleanup_general_gather(tag);
  free(tmp);

}

/*------------------------------------------------------------------*/
/* Translate an su3_vector field by a fixed coordinate displacement rshift */

void shift_su3_vector(su3_vector *src, int rshift[]){

  msg_tag *tag;
  int dir, doshift;
  int ndim[4];
  int i;
  site *s;
  su3_vector *tmp;

  ndim[XUP] = nx; ndim[YUP] = ny; ndim[ZUP] = nz; ndim[TUP] = nt;

  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = rshift[dir] % ndim[dir];
    if(rshift[dir] != 0)doshift = 1;
  }

  /* Don't do a null shift */
  if(!doshift)return;

  tmp = create_v_field();
  copy_v_field(tmp, src);

  tag = start_general_gather_field( tmp, sizeof(su3_vector), 
				    rshift, EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);

  FORALLSITES(i,s){
    su3vec_copy((su3_vector *)(gen_pt[0][i]), src+i);
  }
  cleanup_general_gather(tag);
  free(tmp);

}

/*------------------------------------------------------------------*/
/* Translate a wilson_vector field by a fixed coordinate displacement rshift */

void shift_wilson_vector(wilson_vector *src, int rshift[]){

  msg_tag *tag;
  int dir, doshift;
  int ndim[4];
  int i;
  site *s;
  wilson_vector *tmp;

  ndim[XUP] = nx; ndim[YUP] = ny; ndim[ZUP] = nz; ndim[TUP] = nt;

  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = rshift[dir] % ndim[dir];
    if(rshift[dir] != 0)doshift = 1;
  }

  /* Don't do a null shift */
  if(!doshift)return;

  tmp = create_wv_field();
  copy_wv_field(tmp, src);

  tag = start_general_gather_field( tmp, sizeof(wilson_vector), 
				    rshift, EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);

  FORALLSITES(i,s){
    copy_wvec((wilson_vector *)(gen_pt[0][i]), src+i);
  }
  cleanup_general_gather(tag);
  free(tmp);

}

