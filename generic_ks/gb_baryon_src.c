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
  
  /* Never use links for source for staggered baryon */
  //if(use_links){
  //  FORALLSITES(i,s)
  //    {
  //      mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec0+i );
  //    }
  // } else {
    /* src -> tvec0 */
  //}
  tag[1] = start_gather_field(src, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1]);
#endif
  wait_gather(tag[0]);

  /* Never use links for source for staggered baryon */
  //if(use_links){
  //  FORALLSITES(i,s)
  //    {
  //      mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i );
  //    }
  //} else {
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
    
    /*
    // DEBUG
    site *s;
    g_sync();
    node0_printf("===========================================\n");
    node0_printf("symmetric shift\n");
    node0_printf("===========================================\n");
    g_sync();
    FORALLSITES(i,s){
      if(s->x==4 && s->y==20 && s->z==14 && s->t==0){
        printf("(4,20,14,0): %.5f\n", (dest+i)->c[0].real);
      }
      if(s->x==5 && s->y==20 && s->z==14 && s->t==0){
        printf("(5,20,14,0): %.5f\n", (dest+i)->c[0].real);
      }
      if(s->x==6 && s->y==20 && s->z==14 && s->t==0){
        printf("(6,20,14,0): %.5f\n", (dest+i)->c[0].real);
      }
      if(s->x==7 && s->y==20 && s->z==14 && s->t==0){
        printf("(7,20,14,0): %.5f\n", (dest+i)->c[0].real);
      }
      if(s->x==8 && s->y==20 && s->z==14 && s->t==0){
        printf("(8,20,14,0): %.5f\n", (dest+i)->c[0].real);
      }
    }

    g_sync();
    node0_printf("===========================================\n");
    g_sync();
    */

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
 
  /*
  double global_sum_src = 0.0;
  int icolor;
  FORALLFIELDSITES(i){ 
    for (icolor=0; icolor<3; icolor++){
      global_sum_src += ((dest+i)->c[icolor].real)*((dest+i)->c[icolor].real); 
    }
  }
  g_doublesum( &global_sum_src );
  node0_printf("apply_xport: source_norm = %e\n", (double)global_sum_src);
  */
  //FORALLFIELDSITES(i){ scalar_mult_su3_vector(dest+i, 1./sqrt(global_sum_src),dest+i); }
  
  /* not necessary if phases are not in links */
  //if(n == 1){
  //  /* one link */
  //  d[0][0] = qss_op->dir1;
  //  apply_sym_shift_src_v(n,d[0],r0,dest,src,links);
  //}
  //else if (n == 2){
  //  /* two link */
  //  d[0][0] = qss_op->dir1; d[0][1] = qss_op->dir2;
  //  d[1][1] = qss_op->dir1; d[1][0] = qss_op->dir2;
  //  apply_sym_shift_src_v(n,d[0],r0,tvec0,src,links);
  //  apply_sym_shift_src_v(n,d[1],r0,tvec1,src,links);
  //  FORALLSITES(i,s){
  //    add_su3_vector( tvec0+i, tvec1+i, dest+i );
  //    scalar_mult_su3_vector( dest+i, 0.5, dest+i );
  //  }
  //}
  //else if (n == 3){
  //  /* three link */
  //  /* use the d given */
  //  for(j=0;j<6;j++){
  //    apply_sym_shift_src_v(n,d[j],r0,tvec0,src,links);
  //    if(j==0) copy_v_field(tvec1,tvec0);
  //    else add_v_fields(tvec1,tvec1,tvec0);
  //  }
  //  FORALLSITES(i,s){
  //    scalar_mult_su3_vector( tvec1+i, 1./6., dest+i );
  //  }
  //}
}

