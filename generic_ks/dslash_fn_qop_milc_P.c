/**** dslash_fn_qop_milc_P.c - dslash for improved KS fermions ****/
/* MIMD version 7 */

/* This version supports the MILC QOP API d_congrad5_fn_qop_milc.c
 * It is intended for comparing results with other QOP routines.
 * It is based on dslash_fn2.c */

/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */

/* This version waits for gathers from both positive and negative
   directions before computing, thereby combining two lattice loops in
   an attempt to gain prefetching time for sub_four_su3_vecs */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced the versions */

#if ( QOP_Precision == 1 )

#define CLEANUP_GATHERS_QOP_MILC cleanup_gathers_qop_milc_F
#define CLEANUP_DSLASH_QOP_MILC_TEMPS cleanup_dslash_qop_milc_temps_F
#define DSLASH_FN_QOP_MILC dslash_fn_qop_milc_F
#define DSLASH_FN_QOP_MILC_FIELD_SPECIAL dslash_fn_qop_milc_field_special_F

#else

#define CLEANUP_GATHERS_QOP_MILC cleanup_gathers_qop_milc_D
#define CLEANUP_DSLASH_QOP_MILC_TEMPS cleanup_dslash_qop_milc_temps_D
#define DSLASH_FN_QOP_MILC dslash_fn_qop_milc_D
#define DSLASH_FN_QOP_MILC_FIELD_SPECIAL dslash_fn_qop_milc_field_special_D

#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */

/* Temporary work space for dslash_fn_field_special */ 
static su3_vector *temp[9] ;
/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

void CLEANUP_GATHERS_QOP_MILC(msg_tag *tags1[], msg_tag *tags2[])
{
  int i;

  for(i=XUP;i<=TUP;i++){
    cleanup_gather( tags1[i] );
    cleanup_gather( tags1[OPP_DIR(i)] );
    cleanup_gather( tags2[i] );
    cleanup_gather( tags2[OPP_DIR(i)] );
  }

  for(i=X3UP;i<=T3UP;i++){
    cleanup_gather( tags1[i] );
    cleanup_gather( tags1[OPP_3_DIR(i)] );
    cleanup_gather( tags2[i] );
    cleanup_gather( tags2[OPP_3_DIR(i)] );
  }
}

void CLEANUP_DSLASH_QOP_MILC_TEMPS(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<9;i++) {
      free(temp[i]) ; 
    }
  temp_not_allocated=1 ;
}


void DSLASH_FN_QOP_MILC( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity)
{
   register int dir;
   msg_tag *tag[16];

   DSLASH_FN_QOP_MILC_FIELD_SPECIAL(fatlinks, longlinks, 
				    src, dest, parity, tag, 1);
   
   /* free up the buffers */
   for(dir=XUP; dir<=TUP; dir++){
     cleanup_gather(tag[dir]);
     cleanup_gather(tag[OPP_DIR(dir)]);
   }
   
   for(dir=X3UP; dir<=T3UP; dir++){
     cleanup_gather(tag[dir]);
     cleanup_gather(tag[OPP_3_DIR(dir)]);
   }
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers and temps! */
void DSLASH_FN_QOP_MILC_FIELD_SPECIAL(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start)
{
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register su3_matrix *fat, *lng;
  
  /* allocate temporary work space only if not already allocated */
  if(temp_not_allocated)
    {
      for( dir=XUP; dir<=TUP; dir++ ){
	temp[dir]  =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
	temp[dir+4]=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      }
      temp[8]=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      temp_not_allocated = 0 ;
    }
  
  switch(parity)
    {
    case EVEN:	otherparity=ODD; break;
    case ODD:	otherparity=EVEN; break;
    case EVENANDODD:	otherparity=EVENANDODD; break;
    }
  
  /* Start gathers from positive directions */
  /* And start the 3-step gather too */
  for( dir=XUP; dir<=TUP; dir++ ){
    if(start==1)
      {
	tag[dir] = start_gather_field( src, sizeof(su3_vector), 
					   dir, parity,gen_pt[dir] );
	tag[DIR3(dir)] = start_gather_field(src, sizeof(su3_vector),
						DIR3(dir),parity, 
						gen_pt[DIR3(dir)] );
      }
    else
      {
	restart_gather_field( src, sizeof(su3_vector), 
				  dir, parity,gen_pt[dir], tag[dir]);
	restart_gather_field(src, sizeof(su3_vector), DIR3(dir), parity, 
				 gen_pt[DIR3(dir)], tag[DIR3(dir)]);
      }
  }
  
  /* Multiply by adjoint matrix at other sites */
  /* Use fat link for single link transport */
  FORSOMEPARITY( i, s, otherparity ){
    FORALLUPDIR(dir){
      fat  = fatlinks + dir*sites_on_node + i;
      lng = longlinks + dir*sites_on_node + i;
      mult_adj_su3_mat_vec( fat, &(src[i]), &(temp[dir][i]));
      /* multiply by 3-link matrices too */
      mult_adj_su3_mat_vec( lng, &(src[i]),&(temp[4+dir][i]));
    }
  } END_LOOP
      
  /* Start gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
      if (start==1) tag[OPP_DIR(dir)] = start_gather_field( temp[dir],
	   sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
      else restart_gather_field( temp[dir], sizeof(su3_vector), 
	   OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
   }

  /* Start 3-neighbour gathers from negative directions */
  for( dir=X3UP; dir <= T3UP; dir++){
      if (start==1) tag[OPP_3_DIR(dir)]=start_gather_field(
             temp[INDEX_3RD(dir)+4], sizeof(su3_vector), 
	     OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
      else restart_gather_field(temp[INDEX_3RD(dir)+4], 
	    sizeof(su3_vector), OPP_3_DIR( dir),parity, 
	    gen_pt[OPP_3_DIR(dir)], tag[OPP_3_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
	wait_gather(tag[DIR3(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	fat = fatlinks + i;
	lng = longlinks + i;

	mult_su3_mat_vec( fat, (su3_vector *)gen_pt[XUP][i], &(dest[i]) );
	mult_su3_mat_vec_sum( fat+sites_on_node, 
			      (su3_vector *)gen_pt[YUP][i], &(dest[i]) );
	mult_su3_mat_vec_sum( fat+sites_on_node*2, 
			      (su3_vector *)gen_pt[ZUP][i], &(dest[i]) );
	mult_su3_mat_vec_sum( fat+sites_on_node*3, 
			      (su3_vector *)gen_pt[TUP][i], &(dest[i]) );

	mult_su3_mat_vec( lng, (su3_vector *)gen_pt[X3UP][i], &(temp[8][i]) );
	mult_su3_mat_vec_sum( lng+sites_on_node, 
			      (su3_vector *)gen_pt[Y3UP][i], &(temp[8][i]) );
	mult_su3_mat_vec_sum( lng+sites_on_node*2, 
			      (su3_vector *)gen_pt[Z3UP][i], &(temp[8][i]) );
	mult_su3_mat_vec_sum( lng+sites_on_node*3, 
			      (su3_vector *)gen_pt[T3UP][i], &(temp[8][i]) );
    } END_LOOP
   
    /* Wait gathers from negative directions, accumulate (negative) */
    /* and the same for the negative 3-rd neighbours */
    for(dir=XUP; dir<=TUP; dir++){
      wait_gather(tag[OPP_DIR(dir)]);
    }
    for(dir=X3UP; dir<=T3UP; dir++){
      wait_gather(tag[OPP_3_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVVV( 
		       &(dest[i+FETCH_UP]),
		       (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
	prefetch_VVVVV( 
		       &(temp[8][i+FETCH_UP]), 
		       (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
		       (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
      }

      sub_four_su3_vecs( &(dest[i]),
			 (su3_vector *)(gen_pt[XDOWN][i]),
			 (su3_vector *)(gen_pt[YDOWN][i]),
			 (su3_vector *)(gen_pt[ZDOWN][i]),
			 (su3_vector *)(gen_pt[TDOWN][i]) );
      sub_four_su3_vecs( &(temp[8][i]), 
			 (su3_vector *)(gen_pt[X3DOWN][i]),
			 (su3_vector *)(gen_pt[Y3DOWN][i]),
			 (su3_vector *)(gen_pt[Z3DOWN][i]),
			 (su3_vector *)(gen_pt[T3DOWN][i]) );
      /* Now need to add these things together */
      add_su3_vector(&(dest[i]), &(temp[8][i]),&(dest[i]));
    } END_LOOP 

}
