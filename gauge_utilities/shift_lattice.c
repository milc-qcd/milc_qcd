/************************  shift_lattice.c  -- ******************/
/* MIMD version 7 */

/* Translate the gauge field by a fixed coordinate displacement rshift */

#include "gauge_utilities_includes.h"

void shift_gauge(int rshift[]){

  msg_tag *tag;
  int dir, doshift;
  int ndim[4];
  int i;
  site *s;

  ndim[XUP] = nx; ndim[YUP] = ny; ndim[ZUP] = nz; ndim[TUP] = nt;

  /* The general gather pulls back by an amount rshift but we want to
     push foward by the specified amount, so we reverse the rshift
     here */
  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = (ndim[dir] - rshift[dir]) % ndim[dir];
    if(rshift[dir] != 0)doshift = 1;
  }

  /* Don't do a null shift */
  if(!doshift)return;

  gauge_field_copy(F_OFFSET(link), F_OFFSET(old_link));

  tag = start_general_gather_site( F_OFFSET(old_link), 4*sizeof(su3_matrix), 
				   rshift, EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);

  FORALLSITES(i,s){
    su3_matrix *shifted_link = (su3_matrix *)gen_pt[0][i];
    FORALLUPDIR(dir){
      su3mat_copy(shifted_link + dir, &s->link[dir]);
    }
  }
  cleanup_general_gather(tag);
}
