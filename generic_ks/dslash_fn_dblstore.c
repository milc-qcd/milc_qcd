/******* dslash_fn_dblstore.c - dslash for improved KS fermions ****/
/* MIMD version 7 */

/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */

/* This version double stores long and fat links  */


/* This version overlaps computation and gathers from negative
   directions, and has an extra lattice loop devoted to exclusively to
   sub_four_vectors (traditional algorithm) */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced versions */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#include "../include/fn_links.h"

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */

#ifndef DBLSTORE_FN
BOMB.  Requires compilation with -DDBLSTORE_FN
#endif

static void 
cleanup_one_gather_set(msg_tag *tags[])
{
  int i;

  for(i=XUP;i<=TUP;i++){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_DIR(i)] );
  }

  for(i=X3UP;i<=T3UP;i++){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_3_DIR(i)] );
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_fn_site( field_offset src, field_offset dest, int parity,
		     fn_links_t *fn ) 
{
  msg_tag *tag[16];

  dslash_fn_site_special(src, dest, parity, tag, 1, fn);
  cleanup_one_gather_set(tag);
}

/* Special dslash_site for use by congrad.  Uses restart_gather_site() when
  possible. Third to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_site, otherwise use restart_gather_site. 
  The calling program must clean up the gathers! */
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     fn_links_t *fn)
{
    register int i;
    register site *s;
    register int dir,otherparity=0;
    register su3_matrix *fat4, *long4;
    su3_matrix *t_fatlink;
    su3_matrix *t_longlink;
    su3_vector *tempvec,*templongvec, *templongv1;
    char myname[] = "dslash_fn_site_special";

    //    if(!fn->fl.valid){
    if(fn==NULL){
      printf("dslash_fn_site_special: invalid fn links!\n");
      terminate(1);
    }
    t_fatlink = get_fatlinks(fn);
    t_longlink = get_lnglinks(fn);

    tempvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
    if(tempvec == NULL){
      printf("%s(%d)No room for temporary\n",myname, this_node);
      terminate(1);
    }

    templongvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
    if(templongvec == NULL){
      printf("%s(%d)No room for temporary\n",myname, this_node);
      terminate(1);
    }
    
    templongv1 = create_v_field();
    
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
	    (su3_vector *)F_PT(s,src), (tempvec+4*i) );
	/* multiply by 3-link matrices too */
	mult_adj_su3_mat_vec_4dir( long4,
	    (su3_vector *)F_PT(s,src), (templongvec+4*i) );
    } END_LOOP

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
/**printf("dslash_fn_site_special: down gathers, start=%d\n",start);**/
      if (start==1){
	/* We need the strided gather so we can pick off one of a
	   group of four vectors in tempvec */
	tag[OPP_DIR(dir)] = 
	  declare_strided_gather( (char *)(tempvec+dir), 4*sizeof(su3_vector), 
				  sizeof(su3_vector), OPP_DIR( dir), parity, 
				  gen_pt[OPP_DIR(dir)] );
	prepare_gather(tag[OPP_DIR(dir)]);
	do_gather(tag[OPP_DIR(dir)]);
      }	else {
	do_gather(tag[OPP_DIR(dir)]);
      }
    }
    
    /* and 3rd neighbours */
    for( dir=X3UP; dir <= T3UP; dir++){
/**printf("dslash_fn_site_special: down gathers, start=%d\n",start);**/
      if (start==1){
	tag[OPP_3_DIR(dir)] = 
	  declare_strided_gather( (char *)(templongvec+INDEX_3RD(dir)), 
				  4*sizeof(su3_vector), 
				  sizeof(su3_vector), OPP_3_DIR(dir), 
				  parity, gen_pt[OPP_3_DIR(dir)] );
	prepare_gather(tag[OPP_3_DIR(dir)]);
	do_gather(tag[OPP_3_DIR(dir)]);
      }	else {
	do_gather(tag[OPP_3_DIR(dir)]);
      }
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
	    templongv1+i);
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
	sub_four_su3_vecs( templongv1+i, 
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
        /*** Now need to add these things together ***/
        add_su3_vector((su3_vector *)F_PT(s,dest), templongv1+i,
				(su3_vector *)F_PT(s,dest));
    } END_LOOP

	destroy_v_field(templongv1);
    free(templongvec);
    free(tempvec);

}

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity,
		      fn_links_t *fn) {
    register int dir;
    msg_tag *tag[16];

    dslash_fn_field_special(src, dest, parity, tag, 1, fn );

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
			     fn_links_t *fn ){
  register int i;
  register site *s;
  register int dir;
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
  if(fn == NULL){
    printf("dslash_fn_field_special: invalid fn links!\n");
    terminate(1);
  }
  t_fatlink = get_fatlinks(fn);
  t_longlink = get_lnglinks(fn);
  t_fatbacklink = get_fatbacklinks(fn);
  t_longbacklink = get_lngbacklinks(fn);

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

/* Apply a Fat-Naik dslash-type shift in direction "dir", either
   forward or backward and with specified weights for the fat and Naik
   components.  Do this for dest sites of the specified
   parity. Accumulate result in dest */


void 
dslash_fn_dir(su3_vector *src, su3_vector *dest, int parity,
	      fn_links_t *fn, int dir, int fb, 
	      Real wtfat, Real wtlong)
{
  register int i ;
  msg_tag *tag[2];
  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
  su3_matrix *fatback = get_fatbacklinks(fn);
  su3_matrix *lngback = get_lngbacklinks(fn);
  su3_vector tmp;
  int do_long = (lng != NULL) && (wtlong != 0.);
  char myname[] = "fn_shift";
  
  if(fat == NULL) {
    printf("%s(%d): fat or lng member is null\n", myname, this_node);
    terminate(1);
  }
  
  if(fb > 0){

    /* Shift from forward direction */

    tag[0] = start_gather_field( src, sizeof(su3_vector), dir, 
				 parity, gen_pt[0] );
    if(do_long)
      tag[1] = start_gather_field( src, sizeof(su3_vector), DIR3(dir), 
				   parity, gen_pt[1] );
    wait_gather(tag[0]);
    if(do_long)
      wait_gather(tag[1]);
  
    FORSOMEFIELDPARITY(i,parity)
      {
	mult_su3_mat_vec( fat+4*i+dir, (su3_vector *)gen_pt[0][i], &tmp );
	scalar_mult_add_su3_vector( dest+i, &tmp, wtfat, dest+i ) ;    
	if(do_long){
	  mult_su3_mat_vec( lng+4*i+dir, (su3_vector *)gen_pt[1][i], &tmp );
	  scalar_mult_add_su3_vector( dest+i, &tmp, wtlong, dest+i ) ;    
	}
      } END_LOOP

  } else {

    /* Shift from backward direction */


    tag[0] = start_gather_field( src, sizeof(su3_vector), OPP_DIR(dir),
				 parity, gen_pt[0] );
    if(do_long)
      tag[1] = start_gather_field( src, sizeof(su3_vector), OPP_3_DIR(DIR3(dir)),
				   parity, gen_pt[1] );

    wait_gather(tag[0]);
    if(do_long)
      wait_gather(tag[1]);
  
    /* NOTE minus sign convention! */

    FORSOMEFIELDPARITY(i,parity)
      {
	mult_su3_mat_vec( fatback+4*i+dir, (su3_vector *)gen_pt[0][i], &tmp );
	scalar_mult_add_su3_vector( dest+i, &tmp, -wtfat, dest+i ) ;    
	if(do_long){
	  mult_su3_mat_vec( lngback+4*i+dir, (su3_vector *)gen_pt[1][i], &tmp );
	  scalar_mult_add_su3_vector( dest+i, &tmp, -wtlong, dest+i ) ;    
	}
      } END_LOOP
  }

  cleanup_gather(tag[0]);

  if(do_long)cleanup_gather(tag[1]);

}
	 
