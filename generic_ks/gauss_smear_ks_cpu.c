/************************** gauss_smear_ks.c ********************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source, Chroma style
 * 
 * CD 4/07 Stolen from Chroma. MILC version.
 */

#include "generic_ks_includes.h"
#include <string.h>
#include "../include/openmp_defs.h"

static su3_vector *wtmp[8] ;
static const char *prec_label[2] = {"F", "D"};

/*--------------------------------------------------------------*/
/* Procedures for allowing reuse of the 2-link gauge connection */

su3_matrix *twolink = NULL;
static int recompute_2link = 1;

void
gauss_smear_reuse_2link_cpu( int flag )
{
  recompute_2link = ! flag;
}

/* Create saved two-link. */
static void
gauss_smear_create_2link_cpu(){
  twolink  = create_G();
}

/* Delete saved two-link. */
void
gauss_smear_delete_2link_cpu()
{
  free(twolink);
  twolink = NULL;
  gauss_smear_reuse_2link_cpu(0);
}


/* Compute two-link -----------------------*/
static void
gauss_smear_compute_twolink(su3_matrix *t_links){

  char myname[] = "gauss_smear_compute_twolink";
  msg_tag *tag;
  su3_matrix * link1 = create_m_field();

  int dir;

  if(twolink == NULL) gauss_smear_create_2link_cpu();

  FORALLUPDIRBUT(TUP,dir){

    tag = declare_strided_gather( t_links+dir, sizeof(su3_matrix)*4, sizeof(su3_matrix),
                                  dir, EVENANDODD, gen_pt[dir] );
    prepare_gather( tag );
    do_gather( tag );   
    wait_gather( tag );

    int i;
    FORALLFIELDSITES_OMP(i,){
      su3mat_copy( (su3_matrix *)gen_pt[dir][i], link1 + i );
    } END_LOOP_OMP;

    cleanup_gather( tag );
    
    FORALLFIELDSITES_OMP(i,){
      mult_su3_nn( t_links + 4*i + dir, link1 + i, twolink + 4*i + dir ); 
    } END_LOOP_OMP;
  }

  destroy_m_field(link1);
}

/*------------------------------------------------------------*/
static void 
malloc_kg_temps(){
  int dir;
  
  for(dir=0;dir<8;dir++)wtmp[dir] = NULL;

  FORALLUPDIRBUT(TUP,dir){
    wtmp[dir] =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    if(wtmp[dir] == NULL){
      printf("node %d can't malloc wtmp[%d]\n",this_node,dir);
      terminate(1);
    }
    memset(wtmp[dir],'\0',sites_on_node*sizeof(su3_vector));
    
    wtmp[OPP_DIR(dir)] =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    if(wtmp[OPP_DIR(dir)] == NULL){
      printf("node %d can't malloc wtmp[%d]\n",this_node,OPP_DIR(dir));
      terminate(1);
    }
    memset(wtmp[OPP_DIR(dir)],'\0',sites_on_node*sizeof(su3_vector));
  }
}

/*------------------------------------------------------------*/
static void 
cleanup_kg_temps(){
  int i ;
  for(i=0;i<8;i++){
    if(wtmp[i] != NULL){
      free(wtmp[i]); 
      wtmp[i] = NULL;
    }
  }
}

/*------------------------------------------------------------*/
/* Double forward parallel transport the quick and dirty way.
   Result in dest */
/* Original code */

static void 
forward2(int dir, su3_vector *dest, su3_vector *src, 
	 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  su3_vector *tmp = create_v_field();
  
  /* start parallel transport of src from up dir */
  tag = start_gather_field( src, sizeof(su3_vector),
			    dir, EVENANDODD, gen_pt[dir] );
  wait_gather(tag);
  
  /* tmp <- U(up,dir) shift(up,dir)(src) */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      mult_su3_mat_vec( t_links + 4*i + dir,  
			(su3_vector * )(gen_pt[dir][i]), 
			tmp + i ); 
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
	mult_su3_mat_vec( t_links + 4*i + dir,  
			  (su3_vector * )(gen_pt[dir][i]), 
			  tmp + i ); 
      }
    }
  }

  cleanup_gather(tag);

  /* start parallel transport of tmp from up dir */
  tag = start_gather_field( tmp, sizeof(su3_vector),
			    dir, EVENANDODD, gen_pt[dir] );
  wait_gather(tag);

  /* dest <- U(up,dir) shift(up,2dir)(src) */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      mult_su3_mat_vec( t_links + 4*i + dir,  
			(su3_vector * )(gen_pt[dir][i]), 
			dest + i ); 
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
	mult_su3_mat_vec( t_links + 4*i + dir,  
			  (su3_vector * )(gen_pt[dir][i]), 
			  dest + i ); 
      }
    }
  }

  cleanup_gather(tag);
  destroy_v_field(tmp);
}

/*------------------------------------------------------------*/
/* Double forward parallel transport using second-neighbor gathers
   Result in dest */
/* Hwancheol Jeong 4/2024 */

static void 
forward2_twolink(int dir, su3_vector *dest, su3_vector *src, 
                 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;

  tag = start_gather_field( src, sizeof(su3_vector),
                            DIR2(dir), EVENANDODD, gen_pt[dir] );
  wait_gather(tag);

  /* dest <- U(up,dir) shift(up,2dir)(src) */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      mult_su3_mat_vec( t_links + 4*i + dir, (su3_vector * )(gen_pt[dir][i]),  
                        dest + i ); 
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
        mult_su3_mat_vec( t_links + 4*i + dir, (su3_vector * )(gen_pt[dir][i]),  
                          dest + i ); 
      }
    }
  }

  cleanup_gather(tag);
}


/*------------------------------------------------------------*/
/* Double backward parallel transport the quick and dirty way.
   Result in dest */
/* Original code */

static void 
backward2(int dir, su3_vector *dest, su3_vector *src, 
	 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  su3_vector *tmp = create_v_field();

  /* prepare parallel transport of psi from down dir */

  /* dest <- U^dagger(down,dir) src */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      /* Work only on the specified time slice(s) */
      mult_adj_su3_mat_vec( t_links +  4*i + dir, src + i, 
			    dest + i );
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      /* Work only on the specified time slice(s) */
      if(s->t == t0){
	mult_adj_su3_mat_vec( t_links +  4*i + dir, src + i, 
			      dest + i );
      }
    }
  }
  
  /* gen_pt_array <- shift(down,dir)(dest) */
  tag = start_gather_field(dest, 
			   sizeof(su3_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );
  wait_gather(tag);
  
  /* tmp <- U^dagger(down,dir) gen_pt_array */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      mult_adj_su3_mat_vec( t_links + 4*i + dir,  
			    (su3_vector * )(gen_pt[OPP_DIR(dir)][i]), 
			    tmp + i ); 
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
	mult_adj_su3_mat_vec( t_links + 4*i + dir,  
			      (su3_vector * )(gen_pt[OPP_DIR(dir)][i]), 
			      tmp + i ); 
      }
    }
  }
  
  cleanup_gather(tag);

  /* gen_pt_array <- shift(down,dir)(tmp) */
  tag = start_gather_field(tmp, 
			   sizeof(su3_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );

  wait_gather(tag);

  /* dest <- gen_pt_array */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      su3vec_copy( (su3_vector *)gen_pt[OPP_DIR(dir)][i], dest + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
	su3vec_copy( (su3_vector *)gen_pt[OPP_DIR(dir)][i], dest + i);
      }
    }
  }

  cleanup_gather(tag);
  
  destroy_v_field(tmp);
}

/*------------------------------------------------------------*/
/* Double backward parallel transport using second-neighgor gathers
   Result in dest */
/* Hwancheol Jeong 4/2024 */

static void 
backward2_twolink(int dir, su3_vector *dest, su3_vector *src, 
                  su3_matrix *twolink, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  su3_vector *tmp = create_v_field();

  /* prepare parallel transport of psi from down dir */

  /* tmp <- U^dagger(down,dir) src */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      mult_adj_su3_mat_vec( twolink +  4*i + dir, src + i, 
                            tmp + i );
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
        mult_adj_su3_mat_vec( twolink +  4*i + dir, src + i, 
                              tmp + i );
      }
    }
  }
  
  /* gen_pt_array <- shift(down,dir)(dest) */
  tag = start_gather_field(tmp, 
                           sizeof(su3_vector), DIR2(OPP_DIR(dir)),
                           EVENANDODD, gen_pt[OPP_DIR(dir)] );
  wait_gather(tag);

  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      su3vec_copy( (su3_vector *)gen_pt[OPP_DIR(dir)][i], dest + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0){
        su3vec_copy( (su3_vector *)gen_pt[OPP_DIR(dir)][i], dest + i);
      }
    }
  } 
  
  cleanup_gather(tag);
  
  destroy_v_field(tmp);
}


/*------------------------------------------------------------*/
/* For staggered fermions we compute the Laplacian on sites displaced
   by 2 lattice units */
/* Compute chi <- msq * psi - Lapl_3d psi
   where Lapl_3d psi(r) = -6 psi + sum_{dir=1}^3 [psi(r+2*dir) + psi(r-2*dir)] */
/* Original code */

void 
klein_gord_field(su3_vector *psi, su3_vector *chi, 
		 su3_matrix *t_links, Real msq, int t0)
{
  Real ftmp = 6 + msq;  /* for 3D */
  int i, dir;
  site *s;

  malloc_kg_temps();

  /* chi = psi * ftmp; */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      scalar_mult_su3_vector(psi + i, ftmp, chi + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0)
	scalar_mult_su3_vector(psi + i, ftmp, chi + i);
    }
  }

  /* do 2-link parallel transport of psi in all dirs */
  FORALLUPDIRBUT(TUP,dir){
    forward2(dir, wtmp[dir], psi, t_links, t0);
    backward2(dir, wtmp[OPP_DIR(dir)], psi, t_links, t0);
  }
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) -
     sum_dir shift(down,dir) U^\dagger(down,dir)(psi) */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      sub_su3_vector( chi + i, wtmp[XUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[YUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[ZUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[XDOWN] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[YDOWN] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[ZDOWN] + i, chi + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(t0 == ALL_T_SLICES || s->t == t0){
	sub_su3_vector( chi + i, wtmp[XUP] + i, chi + i);
	sub_su3_vector( chi + i, wtmp[YUP] + i, chi + i);
	sub_su3_vector( chi + i, wtmp[ZUP] + i, chi + i);
	sub_su3_vector( chi + i, wtmp[XDOWN] + i, chi + i);
	sub_su3_vector( chi + i, wtmp[YDOWN] + i, chi + i);
	sub_su3_vector( chi + i, wtmp[ZDOWN] + i, chi + i);
      }
    }
  }

  cleanup_kg_temps();
}

/*------------------------------------------------------------*/
/* For staggered fermions we compute the Laplacian on sites displaced
   by 2 lattice units.  */
/* Compute chi <- msq * psi - Lapl_3d psi
   where Lapl_3d psi(r) = -6 psi + sum_{dir=1}^3 [psi(r+2*dir) + psi(r-2*dir)] */
/* Uses second-neighbor gathers  */
/* t_links is ignored */
/* Hwancheol Jeong 4/2024 */

void 
klein_gord_field_twolink(su3_vector *psi, su3_vector *chi, 
                         su3_matrix *t_links, Real msq, int t0)
{
  Real ftmp = 6 + msq;  /* for 3D */
  int i, dir;
  site *s;

  malloc_kg_temps();

  /* chi = psi * ftmp; */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      scalar_mult_su3_vector(psi + i, ftmp, chi + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(s->t == t0)
        scalar_mult_su3_vector(psi + i, ftmp, chi + i);
    }
  }

  /* do 2-link parallel transport of psi in all dirs */
  FORALLUPDIRBUT(TUP,dir){
    forward2_twolink(dir, wtmp[dir], psi, twolink, t0);
    backward2_twolink(dir, wtmp[OPP_DIR(dir)], psi, twolink, t0);
  }
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) -
     sum_dir shift(down,dir) U^\dagger(down,dir)(psi) */
  if(t0 == ALL_T_SLICES){
    FORALLFIELDSITES_OMP(i,){
      sub_su3_vector( chi + i, wtmp[XUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[YUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[ZUP] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[XDOWN] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[YDOWN] + i, chi + i);
      sub_su3_vector( chi + i, wtmp[ZDOWN] + i, chi + i);
    } END_LOOP_OMP;
  } else {
    FORALLSITES(i,s){
      if(t0 == ALL_T_SLICES || s->t == t0){
        sub_su3_vector( chi + i, wtmp[XUP] + i, chi + i);
        sub_su3_vector( chi + i, wtmp[YUP] + i, chi + i);
        sub_su3_vector( chi + i, wtmp[ZUP] + i, chi + i);
        sub_su3_vector( chi + i, wtmp[XDOWN] + i, chi + i);
        sub_su3_vector( chi + i, wtmp[YDOWN] + i, chi + i);
        sub_su3_vector( chi + i, wtmp[ZDOWN] + i, chi + i);
      }
    }
  }

  cleanup_kg_temps();
}

/*------------------------------------------------------------*/

/* Computes 
   src <- exp[-width^2/4 * Lapl_3d] src 
   by approximating exp(a) as (1 + a/iters)^iters 
   and Lap_3d is the discrete three dimensional Laplacian
*/
/* Original code */

void 
gauss_smear_v_field_cpu(su3_vector *src, su3_matrix *t_links,
			Real width, int iters, int t0)
{
  su3_vector *tmp;
  Real ftmp = -(width*width)/(4*iters*4);  /* Extra 4 to compensate for stride 2 */
  Real ftmpinv = 1. / ftmp;
  int i, j;
  site *s;

  if(t_links == NULL){
    printf("gauss_smear_v_field(%d): NULL t_links\n",this_node);
    terminate(1);
  }

  double dtime = -dclock();

  tmp = create_v_field();
  
  /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */

  for(j = 0; j < iters; j++)
    {
      if(t0 == ALL_T_SLICES){
	FORALLFIELDSITES_OMP(i,){
	  /* tmp = src * ftmp; */
	  scalar_mult_su3_vector(src+i, ftmp, tmp+i);
	} END_LOOP_OMP;
      } else {
	FORALLSITES(i,s){
	  /* tmp = src * ftmp; */
	  if(s->t == t0)
	    scalar_mult_su3_vector(src+i, ftmp, tmp+i);
	}
      }
      klein_gord_field(tmp, src, t_links, ftmpinv, t0);
    }

  dtime += dclock();

  if(this_node==0){
    printf("GSMEAR: time = %e (fn %s) iters = %d\n",
	   dtime, prec_label[MILC_PRECISION-1], iters);
    fflush(stdout);
  }
  
  destroy_v_field(tmp);
}

/*------------------------------------------------------------*/

/* Computes 
   src <- exp[-width^2/4 * Lapl_3d] src 
   by approximating exp(a) as (1 + a/iters)^iters 
   and Lap_3d is the discrete three dimensional Laplacian
*/
/* Hwancheol Jeong 4/2024 */

void 
gauss_smear_v_field_cpu_twolink(su3_vector *src, su3_matrix *t_links,
                                Real width, int iters, int t0)
{
  char myname[] = "gauss_smear_v_field_cpu_twolink";
  su3_vector *tmp;
  Real ftmp = -(width*width)/(4*iters*4);  /* Extra 4 to compensate for stride 2 */
  Real ftmpinv = 1. / ftmp;
  int i, j;
  site *s;
  static su3_matrix *t_links_last = NULL;

  if(t_links == NULL){
    printf("%s(%d): NULL t_links\n", __func__, this_node);
    terminate(1);
  }

  double dtime = -dclock();

  if ( t_links != t_links_last){
    /* If the t_link pointer changed, refresh the two-link */
    node0_printf("%s: [Warning] Input field for two-links changed.  Will recompute them\n", myname);
    gauss_smear_delete_2link_cpu();
    t_links_last = t_links;
  }

  if(recompute_2link || twolink == NULL){
    gauss_smear_compute_twolink(t_links);
    gauss_smear_reuse_2link_cpu(1);
  }    

  /*-----------------------------------------*/
  
  tmp = create_v_field();
  
  /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */

  for(j = 0; j < iters; j++)
  {
    if(t0 == ALL_T_SLICES){
      FORALLFIELDSITES_OMP(i,){
        /* tmp = src * ftmp; */
        scalar_mult_su3_vector(src+i, ftmp, tmp+i);
      } END_LOOP_OMP;
    } else {
      FORALLSITES(i,s){
        /* tmp = src * ftmp; */
        if(s->t == t0)
          scalar_mult_su3_vector(src+i, ftmp, tmp+i);
      }
    }
    klein_gord_field_twolink(tmp, src, t_links, ftmpinv, t0);
  }

  dtime += dclock();

  if(this_node==0){
    printf("GSMEAR: time = %e (fn %s) iters = %d\n",
           dtime, prec_label[MILC_PRECISION-1], iters);
    fflush(stdout);
  }
  
  destroy_v_field(tmp);
}

