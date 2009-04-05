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

  doshift = 0;
  FORALLUPDIR(dir){
    rshift[dir] = rshift[dir] % ndim[dir];
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

/* copy a gauge field - an array of four su3_matrices */
void gauge_field_copy(field_offset src,field_offset dest){
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    su3mat_copy( (su3_matrix *)F_PT(s,src2),
		(su3_matrix *)F_PT(s,dest2) );
	    src2 += sizeof(su3_matrix);
	    dest2 += sizeof(su3_matrix);
	}
    }
#ifdef FN
  invalidate_all_ferm_links(&fn_links);
#endif
}
