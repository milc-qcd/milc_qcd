/************************** gauss_smear_ks.c ********************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source, Chroma style
 * 
 * CD 4/07 Stolen from Chroma. MILC version.
 */

#include "generic_ks_includes.h"
#include <string.h>

static su3_vector *wtmp[8] ;

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
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_su3_mat_vec( t_links + 4*i + dir,  
			(su3_vector * )(gen_pt[dir][i]), 
			tmp + i ); 
    }
  }

  cleanup_gather(tag);

  /* start parallel transport of tmp from up dir */
  tag = start_gather_field( tmp, sizeof(su3_vector),
			    dir, EVENANDODD, gen_pt[dir] );
  wait_gather(tag);

  /* dest <- U(up,dir) shift(up,2dir)(src) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_su3_mat_vec( t_links + 4*i + dir,  
			(su3_vector * )(gen_pt[dir][i]), 
			dest + i ); 
    }
  }

  cleanup_gather(tag);
  destroy_v_field(tmp);
}

/*------------------------------------------------------------*/
/* Double backward parallel transport the quick and dirty way.
   Result in dest */

static void 
backward2(int dir, su3_vector *dest, su3_vector *src, 
	 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  su3_vector *tmp = create_v_field();

  /* prepare parallel transport of psi from down dir */
  FORALLSITES(i,s){
    /* Work only on the specified time slice(s) */
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_adj_su3_mat_vec( t_links +  4*i + dir, src + i, 
			    dest + i );
    }
  }
  
  tag = start_gather_field(dest, 
			   sizeof(su3_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );
  wait_gather(tag);
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_adj_su3_mat_vec( t_links + 4*i + dir,  
			    (su3_vector * )(gen_pt[dir][i]), 
			    tmp + i ); 
    }
  }
  
  cleanup_gather(tag);

  tag = start_gather_field(tmp, 
			   sizeof(su3_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );

  wait_gather(tag);

  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      su3vec_copy( (su3_vector *)gen_pt[dir][i], dest + i);
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

static void 
klein_gord_field(su3_vector *psi, su3_vector *chi, 
		 su3_matrix *t_links, Real msq, int t0)
{
  Real ftmp = 6 + msq;  /* for 3D */
  int i, dir;
  site *s;

  malloc_kg_temps();

  /* chi = psi * ftmp; */
  FORALLSITES(i,s){
    scalar_mult_su3_vector(psi + i, ftmp, chi + i);
  }

  /* do 2-link parallel transport of psi in all dirs */
  FORALLUPDIRBUT(TUP,dir){
    forward2(dir, wtmp[dir], psi, t_links, t0);
    backward2(dir, wtmp[OPP_DIR(dir)], psi, t_links, t0);
  }
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) -
     sum_dir shift(down,dir) U^\dagger(down,dir)(psi) */
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

  cleanup_kg_temps();
}

/*------------------------------------------------------------*/

/* Computes 
   src <- exp[-width^2/4 * Lapl_3d] src 
   by approximating exp(a) as (1 + a/iters)^iters 
   and Lap_3d is the discrete three dimensional Laplacian

*/

void 
gauss_smear_v_field(su3_vector *src, su3_matrix *t_links,
		    Real width, int iters, int t0)
{
  su3_vector *tmp;
  Real ftmp = -(width*width)/(4*iters);
  Real ftmpinv = 1. / ftmp;
  int i, j;
  site *s;

  if(t_links == NULL){
    printf("gauss_smear_v_field(%d): NULL t_links\n",this_node);
    terminate(1);
  }

  tmp = create_v_field();
  
  /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */

  for(j = 0; j < iters; j++)
    FORALLSITES(i,s){
      /* tmp = src * ftmp; */
      scalar_mult_su3_vector(src+i, ftmp, tmp+i);
    }
  klein_gord_field(tmp, src, t_links, ftmpinv, t0);
  
  destroy_v_field(tmp);
}

void 
gauss_smear_ks_prop_field(ks_prop_field *src, su3_matrix *t_links,
			  Real width, int iters, int t0)
{
  int color;

  for(color = 0; color < src->nc; color++)
    gauss_smear_v_field(src->v[color], t_links, width, iters, t0);

}
