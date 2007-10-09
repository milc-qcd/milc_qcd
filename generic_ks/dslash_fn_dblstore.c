// dslash_fn.dstore.c - double store long and fat links

/******* dslash_fn.c - dslash for improved KS fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */

/* This version overlaps computation and gathers from negative
   directions, and has an extra lattice loop devoted to exclusively to
   sub_four_vectors (traditional algorithm) */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced versions */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */

#ifndef DBLSTORE_FN
BOMB.  Requires compilation with -DDBLSTORE_FN
#endif

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
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

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_fn_site( field_offset src, field_offset dest, int parity,
		     fn_links_t *fn, ks_action_paths *ap) {
    register int dir;
    msg_tag *tag[16];

    dslash_fn_site_special(src, dest, parity, tag, 1, fn, ap );

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

/* Special dslash for use by congrad.  Uses restart_gather_site() when
  possible. Last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_site, otherwise use restart_gather_site. 
  The calling program must clean up the gathers! */
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     fn_links_t *fn, ks_action_paths *ap)
{
    register int i;
    register site *s;
    register int dir,otherparity=0;
    register su3_matrix *fat4, *long4;
    su3_matrix *t_fatlink;
    su3_matrix *t_longlink;

    load_fn_links(fn, ap);
    t_fatlink = fn->fat;
    t_longlink = fn->lng;
    
    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

    /* Start gathers from positive directions */
    for(dir=XUP; dir<=TUP; dir++){
/**printf("dslash_fn_site_special: up gathers, start=%d\n",start);**/
	if(start==1) tag[dir] = start_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] );
	else restart_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] , tag[dir] ); 
    }

    /* and start the 3rd neighbor gather */
    for(dir=X3UP; dir<=T3UP; dir++){
        if(start==1) tag[dir] = start_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] );
	else restart_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] , tag[dir] ); 
    }

    /* Multiply by adjoint matrix at other sites */
    FORSOMEPARITY(i,s,otherparity){

      fat4 = &(t_fatlink[4*i]);
      long4 = &(t_longlink[4*i]);
	mult_adj_su3_mat_vec_4dir( fat4,
	    (su3_vector *)F_PT(s,src), s->tempvec );
	/* multiply by 3-link matrices too */
	mult_adj_su3_mat_vec_4dir( long4,
	    (su3_vector *)F_PT(s,src), s->templongvec );
    } END_LOOP

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
/**printf("dslash_fn_site_special: down gathers, start=%d\n",start);**/
	if (start==1) tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather_site( F_OFFSET(tempvec[dir]), sizeof(su3_vector),
	    OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] , tag[OPP_DIR(dir)] );
    }

    /* and 3rd neighbours */
    for( dir=X3UP; dir <= T3UP; dir++){
/**printf("dslash_fn_site_special: down gathers, start=%d\n",start);**/
	if (start==1) tag[OPP_3_DIR(dir)] = 
	  start_gather_site( F_OFFSET(templongvec[INDEX_3RD(dir)]),
	  sizeof(su3_vector), OPP_3_DIR(dir), parity, gen_pt[OPP_3_DIR(dir)] );
	else restart_gather_site( F_OFFSET(templongvec[INDEX_3RD(dir)]),
	  sizeof(su3_vector), OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)],
	  tag[OPP_3_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
    }

    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=X3UP; dir<=T3UP; dir++){
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
      fat4 = &(t_fatlink[4*i]);
      long4 = &(t_longlink[4*i]);
      mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    (su3_vector *)F_PT(s,dest));
      mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    (su3_vector *) &(s->templongv1));
    } END_LOOP

    /* Wait gathers from negative directions, accumulate (negative) */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    } 

    /* and the same for the negative 3-rd neighbours */

    for(dir=X3UP; dir<=T3UP; dir++){
	wait_gather(tag[OPP_3_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
	sub_four_su3_vecs( & (s->templongv1), 
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
        /*** Now need to add these things together ***/
        add_su3_vector((su3_vector *)F_PT(s,dest), &(s->templongv1),
				(su3_vector *)F_PT(s,dest));
    } END_LOOP

}

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity,
		      fn_links_t *fn, ks_action_paths *ap) {
    register int dir;
    msg_tag *tag[16];

    dslash_fn_field_special(src, dest, parity, tag, 1, fn, ap );

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
void dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			     int parity, msg_tag **tag, int start,
			     fn_links_t *fn, ks_action_paths *ap){
  register int i;
  register site *s;
  register int dir,otherparity;
  register su3_matrix *fat4, *long4, *fatback4, *longback4;
  su3_vector tvec;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;
  su3_matrix *t_fatbacklink;
  su3_matrix *t_longbacklink;
#ifdef D_FN_GATHER13
  int coords[4]; /* for avoiding gathers */
  static int d_fn_g13_checked = 0;
#endif

  /* load fatlinks and longlinks */
  load_fn_links(fn, ap);
  t_fatlink = fn->fat;
  t_longlink = fn->lng;
  t_fatbacklink = fn->fatback;
  t_longbacklink = fn->lngback;

  switch(parity)
    {
    case EVEN:	otherparity=ODD; break;
    case ODD:	otherparity=EVEN; break;
    case EVENANDODD:	otherparity=EVENANDODD; break;
    }
  
#ifdef D_FN_GATHER13
  /* Start gathers from positive directions, first and third
     neighbors */
  /* Then change pointers from first neighbor to point into third
      neighbor results, so we don't have to restart those gathers,
      provided third neighbors are done */
  for( dir=XUP; dir<=TUP; dir++ ){
    if(start==1) 
      {
	tag[DIR3(dir)] = start_gather_field(src, sizeof(su3_vector),
					    DIR3(dir),parity, 
					    gen_pt[DIR3(dir)] );
	tag[dir] = start_gather_field( src, sizeof(su3_vector), 
				       dir, parity,gen_pt[dir] );
	/* if the first neighbor came from another node, we should be able
	   to find it in the third neighbor list, layout permitting --
	   it's the third neighbor of the site two back from us */
	FORSOMEPARITY(i,s,parity){
	  if( gen_pt[dir][i] < (char *)src 
	    || gen_pt[dir][i] >= (char *)(src+sites_on_node) ){
	      coords[XUP]=s->x; coords[YUP]=s->y; coords[ZUP]=s->z;
	      coords[TUP]=s->t; coords[dir]-=2;
	      gen_pt[dir][i]=gen_pt[DIR3(dir)][node_index(coords[XUP],
		coords[YUP],coords[ZUP], coords[TUP])];
	      if(d_fn_g13_checked == 0)
		if(node_number(coords[XUP],coords[YUP],
			       coords[ZUP],coords[TUP]) != this_node){
		  printf("node %d Can't use D_FN_GATHER13 with this layout!\n",
			 this_node);
		  terminate(1);
		}
	  }
	} END_LOOP
    }
    else {
      restart_gather_field(src, sizeof(su3_vector), DIR3(dir), parity, 
			   gen_pt[DIR3(dir)], tag[DIR3(dir)]);
      //First nearest neighbor gather doesn't need restarting - pointer are OK
    }
  }

  d_fn_g13_checked = 1;
  
  /* Start gathers from negative directions */
  /* Start 3-neighbour gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
      if (start==1){
	tag[OPP_3_DIR(DIR3(dir))]=start_gather_field( src,
	  sizeof(su3_vector), OPP_3_DIR(DIR3(dir)), parity, gen_pt[OPP_3_DIR(DIR3(dir))] );
        tag[OPP_DIR(dir)] = start_gather_field( src,
	   sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
	FORSOMEPARITY(i,s,parity){
	  if( gen_pt[OPP_DIR(dir)][i] < (char *)src 
	    || gen_pt[OPP_DIR(dir)][i] >= (char *)(src+sites_on_node) ){
	      coords[XUP]=s->x; coords[YUP]=s->y; coords[ZUP]=s->z;
	      coords[TUP]=s->t; coords[dir]+=2;
	      gen_pt[OPP_DIR(dir)][i]=gen_pt[OPP_3_DIR(DIR3(dir))][node_index(coords[XUP],
		coords[YUP],coords[ZUP], coords[TUP])];
	    /* don't need the check here - if it works forward it
	       should work backward */
	  }
	} END_LOOP
      }
      else{
	restart_gather_field( src, sizeof(su3_vector),
         OPP_3_DIR(DIR3(dir)),parity, gen_pt[OPP_3_DIR(DIR3(dir))], 
			      tag[OPP_3_DIR(DIR3(dir))] );
        //Don't restart first neighbor
      }
   }

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  /* wait for the 3-neighbours from positive directions, multiply */
  for(dir=XUP; dir<=TUP; dir++){
    if(start==1) wait_gather(tag[dir]);
    wait_gather(tag[DIR3(dir)]);
  }

#else
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
  
  /* Start gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
      if (start==1) tag[OPP_DIR(dir)] = start_gather_field( src,
	   sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
      else restart_gather_field( src, sizeof(su3_vector), 
	   OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
   }

  /* Start 3-neighbour gathers from negative directions */
  for( dir=X3UP; dir <= T3UP; dir++){
      if (start==1) tag[OPP_3_DIR(dir)]=start_gather_field(
        src, sizeof(su3_vector), OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
      else restart_gather_field( src, sizeof(su3_vector),
        OPP_3_DIR( dir),parity, gen_pt[OPP_3_DIR(dir)], tag[OPP_3_DIR(dir)] );
  }

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  /* wait for the 3-neighbours from positive directions, multiply */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[dir]);
    wait_gather(tag[DIR3(dir)]);
  }
#endif
  
  FORSOMEPARITY(i,s,parity){
      fat4 = &(t_fatlink[4*i]);
      long4 = &(t_longlink[4*i]);
      mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    &(dest[i]) );

      mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    &tvec );
      add_su3_vector(&(dest[i]), &tvec, &(dest[i]) );
  } END_LOOP
   
#ifdef D_FN_GATHER13
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  for(dir=XUP; dir<=TUP; dir++){
      if( start==1 )wait_gather(tag[OPP_DIR(dir)]);
      wait_gather(tag[OPP_3_DIR(DIR3(dir))]);
  }
#else
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  for(dir=XUP; dir<=TUP; dir++){
      wait_gather(tag[OPP_DIR(dir)]);
    }
    for(dir=X3UP; dir<=T3UP; dir++){
      wait_gather(tag[OPP_3_DIR(dir)]);
    }
#endif

  FORSOMEPARITY(i,s,parity){
      fatback4 = &(t_fatbacklink[4*i]);
      longback4 = &(t_longbacklink[4*i]);
      mult_su3_mat_vec_sum_4dir( fatback4,
	    (su3_vector *)gen_pt[XDOWN][i], (su3_vector *)gen_pt[YDOWN][i],
	    (su3_vector *)gen_pt[ZDOWN][i], (su3_vector *)gen_pt[TDOWN][i],
	    &tvec );
      sub_su3_vector(&(dest[i]), &tvec, &(dest[i]) );
      mult_su3_mat_vec_sum_4dir( longback4,
	    (su3_vector *)gen_pt[X3DOWN][i], (su3_vector *)gen_pt[Y3DOWN][i],
	    (su3_vector *)gen_pt[Z3DOWN][i], (su3_vector *)gen_pt[T3DOWN][i],
	    &tvec );
      sub_su3_vector(&(dest[i]), &tvec, &(dest[i]) );
  } END_LOOP 

}

/* We don't need temps, but d_congrad5_fn thinks we do.  Do nothing */
void cleanup_dslash_temps(){
}


#ifdef DM_DU0

/* d(D_slash)/d(u0) routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */
void ddslash_fn_du0_site( field_offset src, field_offset dest, int parity,
			  fn_links_t *fn, ks_action_paths *ap,
			  fn_links_t *fn_dmdu0, ks_action_paths *ap_dmdu0) {
   register int i;
   register site *s;
   register int dir,otherparity=0;
   register su3_matrix *fat4, *long4;
   msg_tag *tag[16];
   su3_matrix *t_dfatlink_du0;
   su3_matrix *t_longlink;

   load_fn_links(fn, ap);
   t_longlink = fn->lng;
   load_fn_links_dmdu0(fn_dmdu0, ap_dmdu0);
   t_dfatlink_du0 = fn_dmdu0->fat;
    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

    /* Start gathers from positive directions */
    /* And start the 3-step gather too */
    for( dir=XUP; dir<=TUP; dir++ ){
	tag[dir] = start_gather_site( src, sizeof(su3_vector), dir, parity,
	    gen_pt[dir] );
	tag[DIR3(dir)] = start_gather_site( src, sizeof(su3_vector), DIR3(dir),
	    parity, gen_pt[DIR3(dir)] );
    }

    /* Multiply by adjoint matrix at other sites */
    /* Use fat link for single link transport */
    FORSOMEPARITY( i, s, otherparity ){

      fat4 = &(t_dfatlink_du0[4*i]);
      long4 = &(t_longlink[4*i]);
	mult_adj_su3_mat_vec_4dir( fat4,
	    (su3_vector *)F_PT(s,src), s->tempvec );
	/* multiply by 3-link matrices too */
	mult_adj_su3_mat_vec_4dir( long4,
	    (su3_vector *)F_PT(s,src), s->templongvec );
	for( dir=XUP; dir<=TUP; dir++ )
	  scalar_mult_su3_vector( &(s->templongvec[dir]), -2.0/u0,
				  &(s->templongvec[dir]) );
    } END_LOOP

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
	tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity,
	    gen_pt[OPP_DIR(dir)] );
    }

    /* Start 3-neighbour gathers from negative directions */
    for( dir=X3UP; dir <= T3UP; dir++){
	tag[OPP_3_DIR(dir)] 
           = start_gather_site( F_OFFSET(templongvec[INDEX_3RD(dir)]),
			   sizeof(su3_vector), OPP_3_DIR( dir), parity,
			   gen_pt[OPP_3_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
	wait_gather(tag[DIR3(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
      fat4 = &(t_dfatlink_du0[4*i]);
      long4 = &(t_longlink[4*i]);

      mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    (su3_vector *)F_PT(s,dest));

      mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    (su3_vector *) &(s->templongv1));
      scalar_mult_su3_vector( (su3_vector *) &(s->templongv1), -2.0/u0,
			      (su3_vector *) &(s->templongv1) );

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

        sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
        sub_four_su3_vecs( &(s->templongv1), 
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
        /* Now need to add these things together */
        add_su3_vector((su3_vector *)F_PT(s,dest), &(s->templongv1),
		       (su3_vector *)F_PT(s,dest));
    } END_LOOP

    /* free up the buffers */
    for(dir=XUP; dir<=TUP; dir++){
	cleanup_gather(tag[dir]);
	cleanup_gather(tag[OPP_DIR(dir)]);
    }
    for(dir=X3UP; dir<=T3UP; dir++){
	cleanup_gather(tag[dir]);
	cleanup_gather(tag[OPP_3_DIR(dir)]);
    }
} /* end ddslash_fn_du0_site */


void ddslash_fn_du0_field( su3_vector *src, su3_vector *dest, int parity,
			  fn_links_t *fn, ks_action_paths *ap,
			  fn_links_t *fn_dmdu0, ks_action_paths *ap_dmdu0) 
{
   register int i;
   register site *s;
   register int dir,otherparity=0;
   msg_tag *tag[16];
   su3_vector *tempvec[4], *templongvec[4], *templongv1 ;
   register su3_matrix *fat4, *long4;
   su3_matrix *t_dfatlink_du0;
   su3_matrix *t_longlink;
    
   for( dir=XUP; dir<=TUP; dir++ )
     {
       tempvec[dir]    =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
       templongvec[dir]=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
     }
   templongv1=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

   load_fn_links(fn, ap);
   t_longlink = fn->lng;
   load_fn_links_dmdu0(fn_dmdu0, ap_dmdu0);
   t_dfatlink_du0 = fn_dmdu0->fat;

   switch(parity)
     {
     case EVEN:	otherparity=ODD; break;
     case ODD:	otherparity=EVEN; break;
     case EVENANDODD:	otherparity=EVENANDODD; break;
     }
   
   /* Start gathers from positive directions */
   /* And start the 3-step gather too */
   for( dir=XUP; dir<=TUP; dir++ ){
     tag[dir] = start_gather_field( src, sizeof(su3_vector), dir, parity,
					gen_pt[dir] );
     tag[DIR3(dir)] = start_gather_field( src, sizeof(su3_vector), 
					      DIR3(dir),parity, 
					      gen_pt[DIR3(dir)] );
   }

   /* Multiply by adjoint matrix at other sites */
   /* Use fat link for single link transport */
   FORSOMEPARITY( i, s, otherparity ){
     fat4 = &(t_dfatlink_du0[4*i]);
     long4 = &(t_longlink[4*i]);

     mult_adj_su3_mat_4vec( fat4, &(src[i]), &(tempvec[0][i]),
			    &(tempvec[1][i]), &(tempvec[2][i]), 
			    &(tempvec[3][i]) );
     /* multiply by 3-link matrices too */
     mult_adj_su3_mat_4vec( long4, &(src[i]),&(templongvec[0][i]),
			    &(templongvec[1][i]), &(templongvec[2][i]), 
			    &(templongvec[3][i]) );
     for( dir=XUP; dir<=TUP; dir++ )
       scalar_mult_su3_vector( &(templongvec[dir][i]), -2.0/u0,
			       &(templongvec[dir][i]) );
   } END_LOOP

   /* Start gathers from negative directions */
   for( dir=XUP; dir <= TUP; dir++){
     tag[OPP_DIR(dir)] = start_gather_field( tempvec[dir],
	   sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
   }

  /* Start 3-neighbour gathers from negative directions */
    for( dir=X3UP; dir <= T3UP; dir++){
      tag[OPP_3_DIR(dir)]=start_gather_field(templongvec[INDEX_3RD(dir)],
	sizeof(su3_vector), OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
	wait_gather(tag[DIR3(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
     fat4 = &(t_dfatlink_du0[4*i]);
     long4 = &(t_longlink[4*i]);

     mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    &(dest[i]) );

     mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    &(templongv1[i]));
     scalar_mult_su3_vector( &(templongv1[i]), -2.0/u0, &(templongv1[i]) );
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
      sub_four_su3_vecs( &(dest[i]),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
      sub_four_su3_vecs( &(templongv1[i]), 
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
      /* Now need to add these things together */
      add_su3_vector(&(dest[i]), &(templongv1[i]),&(dest[i]));
    } END_LOOP 

    /* free up the buffers */
    for(dir=XUP; dir<=TUP; dir++){
      cleanup_gather(tag[dir]);
      cleanup_gather(tag[OPP_DIR(dir)]);
    }
    
    for(dir=X3UP; dir<=T3UP; dir++){
      cleanup_gather(tag[dir]);
      cleanup_gather(tag[OPP_3_DIR(dir)]);
    }
    
    for( dir=XUP; dir<=TUP; dir++ ){
      free(tempvec[dir]);
      free(templongvec[dir]);
    }
    free(templongv1);
} /* end ddslash_fn_du0_field */

#endif /* DM_DU0 */
