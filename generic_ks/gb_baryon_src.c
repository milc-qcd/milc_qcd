#include "generic_ks_includes.h"
#include <math.h>


/*------------------------------------------------------------------*/
/**
    Almost identical to the definitions of sym_shift_field in :    
    flavor_ops2.c, spin_taste_ops.c                     
    Modify source by symmetric shifting in both directions.
    Rephases new source such that it has same phase on each color
    as initial source.
 */
static void
one_link_sym_shift_source(int dir, su3_vector *dest, su3_vector *src,
                          su3_matrix *links, int use_links){
  register int i;
  register site *s;
  msg_tag *tag[2];
  su3_vector *tvec0 = create_v_field();

  tag[0] = start_gather_field(src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0]);
  /* With ONE_SIDED_SHIFT_GB defined, the shift is asymmetric */
#ifndef ONE_SIDED_SHIFT_GB
  
  tag[1] = start_gather_field(src, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1]);
#endif
  wait_gather(tag[0]);

  FORALLSITES(i,s){
     /* gen_pt -> dest */
    su3vec_copy((su3_vector *)gen_pt[0][i], dest+i );
  }
  //}
  cleanup_gather(tag[0]);
#ifndef ONE_SIDED_SHIFT_GB
  wait_gather(tag[1]);
  FORALLSITES(i,s){
    add_su3_vector(dest+i, (su3_vector*)gen_pt[1][i], dest+i) ;
  }
  /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
  FORALLSITES(i,s){
    scalar_mult_su3_vector(dest+i, .5, dest+i);
  }
  cleanup_gather(tag[1]);
#endif
  destroy_v_field(tvec0);
}

/*------------------------------------------------------------------*/
void
apply_par_xport_src_v(su3_vector *dest, su3_vector *src,
      quark_source_sink_op *qss_op, su3_matrix *links){
  double dtime = start_timing();

  int i,c;
  int n = qss_op->disp;
  int use_links = 0; // hack, do not use links
  int dir[3] = {0};
  su3_vector *tvec0 = create_v_field();
  su3_vector *tvec1 = create_v_field();
  su3_vector *tsrc  = create_v_field();
  /* save to a temporary vector to prevent overwriting during computation */
  copy_v_field(tsrc,src); 

  if(n==0){
    copy_v_field(dest,tsrc); 
    

  }
  if(n==1){
    one_link_sym_shift_source(qss_op->dir1,dest,tsrc,links,use_links); 
  }
  else if(n==2){
    one_link_sym_shift_source(qss_op->dir1,tvec0,tsrc,links,use_links); 
    one_link_sym_shift_source(qss_op->dir2,tvec1,tvec0,links,use_links); 
    one_link_sym_shift_source(qss_op->dir2,tvec0,tsrc,links,use_links); 
    one_link_sym_shift_source(qss_op->dir1,dest,tvec0,links,use_links); 
    FORALLFIELDSITES(i){
      add_su3_vector( dest+i, tvec1+i, dest+i ); 
      scalar_mult_su3_vector(dest+i,.5,dest+i);
    }

  }
  else if(n==3){
    for(c=0;c<6;c++){
      switch(c){
        case 0: dir[0]=XUP; dir[1]=YUP; dir[2]=ZUP; break;
        case 1: dir[0]=YUP; dir[1]=ZUP; dir[2]=XUP; break;
        case 2: dir[0]=ZUP; dir[1]=XUP; dir[2]=YUP; break;
        case 3: dir[0]=XUP; dir[1]=ZUP; dir[2]=YUP; break;
        case 4: dir[0]=YUP; dir[1]=XUP; dir[2]=ZUP; break;
        case 5: dir[0]=ZUP; dir[1]=YUP; dir[2]=XUP; break;
      }
      one_link_sym_shift_source(dir[0],tvec0,tsrc,links,use_links); 
      one_link_sym_shift_source(dir[1],tvec1,tvec0,links,use_links); 
      one_link_sym_shift_source(dir[2],tvec0,tvec1,links,use_links); 
      if(c==0) copy_v_field(dest,tvec0);
      else FORALLFIELDSITES(i){ add_su3_vector(dest+i,tvec0+i,dest+i); }
    }
    FORALLFIELDSITES(i){ scalar_mult_su3_vector(dest+i,1./6.,dest+i); }
  }

  destroy_v_field(tvec0);
  destroy_v_field(tvec1);
  destroy_v_field(tsrc);
 
  print_timing(dtime, "parallel transporting source");
}

