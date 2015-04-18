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
#include "../include/openmp_defs.h"
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#include <quark_action.h>

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */

/* Temporary work space for dslash_fn_field_special */ 
static su3_vector *temp[9] ;
/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

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

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<9;i++) {
      free(temp[i]) ; 
    }
  temp_not_allocated=1 ;
}


/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_fn_site( field_offset src, field_offset dest, int parity,
		     fn_links_t *fn )
{
  msg_tag *tag[16];
  
  dslash_fn_site_special(src, dest, parity, tag, 1, fn );
  cleanup_one_gather_set(tag);
}

/* Special dslash_site for use by congrad.  Uses restart_gather_site() when
  possible. Third to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_site, otherwise use restart_gather_site. 
  The calling program must clean up the gathers! */
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     fn_links_t *fn){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register su3_matrix *fat4;
  su3_matrix *t_fatlink;
#ifndef NO_LONG_LINKS
  register su3_matrix *long4;
  su3_matrix *t_longlink;
  su3_vector *templongvec, *templongv1;
#endif
  su3_vector *tempvec;
  char myname[] = "dslash_fn_site_special";
  
  if(fn == NULL){
    printf("dslash_fn_site_special: invalid fn links!\n");
    terminate(1);
  }
#ifndef NO_LONG_LINKS
  t_longlink = get_lnglinks(fn);
#endif
  t_fatlink = get_fatlinks(fn);

  tempvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
  if(tempvec == NULL){
    printf("%s(%d)No room for temporary\n",myname, this_node);
    terminate(1);
  }
  
#ifndef NO_LONG_LINKS
  templongvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
  if(templongvec == NULL){
    printf("%s(%d)No room for temporary\n",myname, this_node);
    terminate(1);
  }
  
  templongv1 = create_v_field();
#endif
  
  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Start gathers from positive directions */
  for(dir=XUP; dir<=TUP; dir++){
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
  FORSOMEPARITYDOMAIN_OMP(i,s,otherparity,private(fat4,long4)){
    if( i < loopend-FETCH_UP ){
      fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      prefetch_4MV4V( 
		     fat4,
		     (su3_vector *)F_PT(s+FETCH_UP,src),
		     tempvec+4*i+FETCH_UP );
#ifndef NO_LONG_LINKS
      long4 = &(t_longlink[4*(i+FETCH_UP)]);
      prefetch_4MV4V(
		     long4,
		     (su3_vector *)F_PT(s+FETCH_UP,src),
		     templongvec+4*i+FETCH_UP );
#endif
    }
    
    fat4 = &(t_fatlink[4*i]);
#ifndef NO_LONG_LINKS
    long4 = &(t_longlink[4*i]);
#endif
    mult_adj_su3_mat_vec_4dir( fat4,
			       (su3_vector *)F_PT(s,src), (tempvec+4*i) );
    /* multiply by 3-link matrices too */
#ifndef NO_LONG_LINKS
    mult_adj_su3_mat_vec_4dir( long4,
			       (su3_vector *)F_PT(s,src), (templongvec+4*i) );
#endif
  } END_LOOP_OMP
      
  /* Start gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
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
  
#ifndef NO_LONG_LINKS
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
#endif
  
  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[dir]);
  }
  
  /* wait for the 3-neighbours from positive directions, multiply */
  for(dir=X3UP; dir<=T3UP; dir++){
    wait_gather(tag[dir]);
  }
  FORSOMEPARITYDOMAIN_OMP(i,s,parity, private(fat4,long4) ){
    if( i < loopend-FETCH_UP ){
      fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      prefetch_4MVVVV( 
		      fat4,
		      (su3_vector *)gen_pt[XUP][i+FETCH_UP],
		      (su3_vector *)gen_pt[YUP][i+FETCH_UP],
		      (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
		      (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
#ifndef NO_LONG_LINKS
      long4 = &(t_longlink[4*(i+FETCH_UP)]);
      prefetch_VV(
		  (su3_vector *)F_PT(s+FETCH_UP,dest),
		  templongv1+i+FETCH_UP);
      prefetch_4MVVVV( 
		      long4,
		      (su3_vector *)gen_pt[X3UP][i+FETCH_UP],
		      (su3_vector *)gen_pt[Y3UP][i+FETCH_UP],
		      (su3_vector *)gen_pt[Z3UP][i+FETCH_UP],
		      (su3_vector *)gen_pt[T3UP][i+FETCH_UP] );
#endif
    }
    fat4 = &(t_fatlink[4*i]);
    mult_su3_mat_vec_sum_4dir( fat4,
	       (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	       (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	       (su3_vector *)F_PT(s,dest));
#ifndef NO_LONG_LINKS
    long4 = &(t_longlink[4*i]);
    mult_su3_mat_vec_sum_4dir( long4,
	       (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	       (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	       templongv1+i);
#endif
  } END_LOOP_OMP
      
  /* Wait gathers from negative directions, accumulate (negative) */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[OPP_DIR(dir)]);
  } 
  
  /* and the same for the negative 3-rd neighbours */
  
  for(dir=X3UP; dir<=T3UP; dir++){
    wait_gather(tag[OPP_3_DIR(dir)]);
  }
  
  FORSOMEPARITYDOMAIN_OMP(i,s,parity, ){
    if( i < loopend-FETCH_UP ){
#ifndef NO_LONG_LINKS
      prefetch_VV(
		  (su3_vector *)F_PT(s+FETCH_UP,dest),
		  templongv1+i+FETCH_UP);
#endif
      prefetch_VVVV( 
		    (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
      prefetch_VVVV( 
		    (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
    }
    sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
		       (su3_vector *)(gen_pt[XDOWN][i]),
		       (su3_vector *)(gen_pt[YDOWN][i]),
		       (su3_vector *)(gen_pt[ZDOWN][i]),
		       (su3_vector *)(gen_pt[TDOWN][i]) );
#ifndef NO_LONG_LINKS
    sub_four_su3_vecs( templongv1+i,
		       (su3_vector *)(gen_pt[X3DOWN][i]),
		       (su3_vector *)(gen_pt[Y3DOWN][i]),
		       (su3_vector *)(gen_pt[Z3DOWN][i]),
		       (su3_vector *)(gen_pt[T3DOWN][i]) );
    /*** Now need to add these things together ***/
    add_su3_vector((su3_vector *)F_PT(s,dest), templongv1+i,
		   (su3_vector *)F_PT(s,dest));
#endif
  } END_LOOP_OMP
      
}

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity,
		      fn_links_t *fn) {
  msg_tag *tag[16];
    
   dslash_fn_field_special(src, dest, parity, tag, 1, fn);
   cleanup_one_gather_set(tag);
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers and temps! */
void dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			     int parity, msg_tag **tag, int start,
			     fn_links_t *fn){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register su3_matrix *fat4;
  su3_matrix *t_fatlink;
#ifndef NO_LONG_LINKS
  register su3_matrix *long4;
  su3_matrix *t_longlink;
#endif
  
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
  
  /* load fatlinks and longlinks */
  if(fn == NULL){
    printf("dslash_fn_field_special: invalid fn links!\n");
    terminate(1);
  }
#ifndef NO_LONG_LINKS
  t_longlink = get_lnglinks(fn);
#endif
  t_fatlink = get_fatlinks(fn);

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
#ifndef NO_LONG_LINKS
	tag[DIR3(dir)] = start_gather_field(src, sizeof(su3_vector),
						DIR3(dir),parity, 
						gen_pt[DIR3(dir)] );
#endif
      }
    else
      {
	restart_gather_field( src, sizeof(su3_vector), 
				  dir, parity,gen_pt[dir], tag[dir]);
#ifndef NO_LONG_LINKS
	restart_gather_field(src, sizeof(su3_vector), DIR3(dir), parity, 
				 gen_pt[DIR3(dir)], tag[DIR3(dir)]);
#endif
      }
  }
  
  /* Multiply by adjoint matrix at other sites */
  /* Use fat link for single link transport */
  FORSOMEPARITYDOMAIN_OMP( i, s, otherparity, private(fat4,long4) ){
    //NOPRE if( i < loopend-FETCH_UP ){
       //NOPRE fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
       //NOPRE prefetch_V(&(src[i+FETCH_UP]));
       //NOPRE prefetch_4MVVVV( 
		       //NOPRE fat4,
		       //NOPRE &(temp[0][i+FETCH_UP]),
		       //NOPRE &(temp[1][i+FETCH_UP]),
		       //NOPRE &(temp[2][i+FETCH_UP]),
		       //NOPRE &(temp[3][i+FETCH_UP]) );
#ifndef NO_LONG_LINKS
       //NOPRE long4 = &(t_longlink[4*(i+FETCH_UP)]);
       //NOPRE prefetch_4MVVVV( 
		       //NOPRE long4,
		       //NOPRE &(temp[4][i+FETCH_UP]),
		       //NOPRE &(temp[5][i+FETCH_UP]),
		       //NOPRE &(temp[6][i+FETCH_UP]),
		       //NOPRE &(temp[7][i+FETCH_UP]) );
#endif
    //NOPRE }

    fat4 = &(t_fatlink[4*i]);
    mult_adj_su3_mat_4vec( fat4, &(src[i]), &(temp[0][i]),
			   &(temp[1][i]), &(temp[2][i]), &(temp[3][i]) );
#ifndef NO_LONG_LINKS
    /* multiply by 3-link matrices too */
    long4 = &(t_longlink[4*i]);
    mult_adj_su3_mat_4vec( long4, &(src[i]),&(temp[4][i]),
			   &(temp[5][i]), &(temp[6][i]), &(temp[7][i]) );
#endif
  } END_LOOP_OMP
      
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
#ifndef NO_LONG_LINKS
    wait_gather(tag[DIR3(dir)]);
#endif
  }
  
  FORSOMEPARITYDOMAIN_OMP(i,s,parity, private(fat4,long4) ){
    //NOPRE if( i < loopend-FETCH_UP ){
      //NOPRE fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      //NOPRE prefetch_4MVVVV( 
		      //NOPRE fat4,
		      //NOPRE (su3_vector *)gen_pt[XUP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[YUP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
      //NOPRE prefetch_VVVV( 
		    //NOPRE (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
#ifndef NO_LONG_LINKS
      //NOPRE long4 = &(t_longlink[4*(i+FETCH_UP)]);
      //NOPRE prefetch_4MVVVV( 
		      //NOPRE long4,
		      //NOPRE (su3_vector *)gen_pt[X3UP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[Y3UP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[Z3UP][i+FETCH_UP],
		      //NOPRE (su3_vector *)gen_pt[T3UP][i+FETCH_UP] );
      //NOPRE prefetch_VVVV( 
		    //NOPRE (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
		    //NOPRE (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
#endif
    //NOPRE }
    
    fat4 = &(t_fatlink[4*i]);
    mult_su3_mat_vec_sum_4dir( fat4,
	       (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	       (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	       &(dest[i]) );
    
#ifndef NO_LONG_LINKS
    long4 = &(t_longlink[4*i]);
    mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    &(temp[8][i]));
#endif
  } END_LOOP_OMP
   
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[OPP_DIR(dir)]);
  }
  for(dir=X3UP; dir<=T3UP; dir++){
    wait_gather(tag[OPP_3_DIR(dir)]);
  }
  
  FORSOMEPARITYDOMAIN_OMP(i,s,parity, ){
    //NOPRE if( i < loopend-FETCH_UP ){
      //NOPRE prefetch_VVVVV( 
		     //NOPRE &(dest[i+FETCH_UP]),
		     //NOPRE (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
      //NOPRE prefetch_VVVVV( 
		     //NOPRE &(temp[8][i+FETCH_UP]), 
		     //NOPRE (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
		     //NOPRE (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
    //NOPRE }
    
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
  } END_LOOP_OMP
      
}

void 
dslash_fn_dir(su3_vector *src, su3_vector *dest, int parity,
	      fn_links_t *fn, int dir, int fb, 
	      Real wtfat, Real wtlong)
{
  register int i ;
  site *s;
  msg_tag *tag[2];
  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
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
  
    FORSOMEPARITYDOMAIN(i,s,parity)
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

    su3_vector *tvec1, *tvec2;

    tvec1 = create_v_field();
    if(do_long)tvec2 = create_v_field();

    FORSOMEPARITYDOMAIN(i,s,parity)
      {
	mult_adj_su3_mat_vec( fat+4*i+dir, src+i, tvec1+i ) ;
	if(do_long)
	  mult_adj_su3_mat_vec( lng+4*i+dir, src+i, tvec2+i ) ;
      } END_LOOP

    tag[0] = start_gather_field(tvec1, sizeof(su3_vector), OPP_DIR(dir), 
				parity, gen_pt[0] );
    if(do_long)
      tag[1] = start_gather_field(tvec2, sizeof(su3_vector), OPP_3_DIR(DIR3(dir)), 
				  parity, gen_pt[1] );
  
    wait_gather(tag[0]);
    if(do_long)wait_gather(tag[1]);

    /* Apply weights.  NOTE minus sign convention! */

    FORSOMEPARITYDOMAIN(i,s,parity)
      {
        scalar_mult_add_su3_vector( dest+i, (su3_vector *)gen_pt[0][i], -wtfat, dest+i ) ;    
        if(do_long)
	  scalar_mult_add_su3_vector( dest+i, (su3_vector *)gen_pt[1][i], -wtlong, dest+i ) ;    
      } END_LOOP

    if(do_long)destroy_v_field(tvec2);
    destroy_v_field(tvec1);

  }

  cleanup_gather(tag[0]);

  if(do_long)cleanup_gather(tag[1]);

}
	 
#if 0

/* d(D_slash)/d(u0) routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */
void ddslash_fn_du0_site( field_offset src, field_offset dest, int parity,
			  ferm_links_t *fn, ferm_links_t *fn_dmdu0) {
   register int i;
   register site *s;
   register int dir,otherparity=0;
   register su3_matrix *fat4;
#ifndef NO_LONG_LINKS
   register su3_matrix *long4;
   su3_matrix *t_longlink;
   su3_vector *templongvec, *templongv1;
#endif
   msg_tag *tag[16];
   su3_matrix *t_dfatlink_du0;
   su3_vector *tempvec;
   char myname[] = "ddslash_fn_site_site";

   if(!fn->fl.valid){
     printf("ddslash_fn_du0_site: invalid fn links!\n");
     terminate(1);
   }
#ifndef NO_LONG_LINKS
   t_longlink = fn->fl.lng;
#endif
   if(!fn_dmdu0->fl.valid){
     printf("ddslash_fn_du0_site: invalid fn_du0 links!\n");
     terminate(1);
   }
   t_dfatlink_du0 = fn_dmdu0->fl.fat;

   tempvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
   if(tempvec == NULL){
     printf("%s(%d)No room for temporary\n",myname, this_node);
     terminate(1);
   }
   
#ifndef NO_LONG_LINKS
   templongvec = (su3_vector *) malloc(sizeof(su3_vector)*4*sites_on_node);
   if(templongvec == NULL){
     printf("%s(%d)No room for temporary\n",myname, this_node);
     terminate(1);
   }
   
   templongv1 = create_v_field();
#endif
   
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
#ifndef NO_LONG_LINKS
	tag[DIR3(dir)] = start_gather_site( src, sizeof(su3_vector), DIR3(dir),
	    parity, gen_pt[DIR3(dir)] );
#endif
    }

    /* Multiply by adjoint matrix at other sites */
    /* Use fat link for single link transport */
    FORSOMEPARITYDOMAIN( i, s, otherparity ){
      if( i < loopend-FETCH_UP ){
	fat4 = &(t_dfatlink_du0[4*(i+FETCH_UP)]);
	prefetch_4MV4V( 
		       fat4,
		       (su3_vector *)F_PT(s+FETCH_UP,src),
		       tempvec+4*i+FETCH_UP);
#ifndef NO_LONG_LINKS
	long4 = &(t_longlink[4*(i+FETCH_UP)]);
	prefetch_4MV4V(
		       long4,
		       (su3_vector *)F_PT(s+FETCH_UP,src),
		       templongvec+4*i+FETCH_UP);
#endif
      }

      fat4 = &(t_dfatlink_du0[4*i]);
      mult_adj_su3_mat_vec_4dir( fat4,
	    (su3_vector *)F_PT(s,src), tempvec+4*i );
#ifndef NO_LONG_LINKS
      long4 = &(t_longlink[4*i]);
      /* multiply by 3-link matrices too */
      mult_adj_su3_mat_vec_4dir( long4,
	    (su3_vector *)F_PT(s,src), templongvec+4*i );
      for( dir=XUP; dir<=TUP; dir++ )
	scalar_mult_su3_vector( templongvec+4*i+dir, -2.0/u0,
				templongvec+4*i+dir );
#endif
    } END_LOOP

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
	tag[OPP_DIR(dir)] = 
	  declare_strided_gather( (char *)(tempvec+dir), 4*sizeof(su3_vector), 
				  sizeof(su3_vector), OPP_DIR( dir), parity, 
				  gen_pt[OPP_DIR(dir)] );
	prepare_gather(tag[OPP_DIR(dir)]);
	do_gather(tag[OPP_DIR(dir)]);
    }

#ifndef NO_LONG_LINKS
    /* Start 3-neighbour gathers from negative directions */
    for( dir=X3UP; dir <= T3UP; dir++){
	tag[OPP_3_DIR(dir)] = 
	  declare_strided_gather( (char *)(templongvec+INDEX_3RD(dir)), 
				  4*sizeof(su3_vector), 
				  sizeof(su3_vector), OPP_3_DIR(dir), 
				  parity, gen_pt[OPP_3_DIR(dir)] );
	prepare_gather(tag[OPP_3_DIR(dir)]);
	do_gather(tag[OPP_3_DIR(dir)]);
    }
#endif
    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
#ifndef NO_LONG_LINKS
	wait_gather(tag[DIR3(dir)]);
#endif
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
      fat4 = &(t_dfatlink_du0[4*i]);
      mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    (su3_vector *)F_PT(s,dest));

#ifndef NO_LONG_LINKS
      long4 = &(t_longlink[4*i]);
      mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
            templongv1+i);
      scalar_mult_su3_vector( templongv1+i, -2.0/u0,
				templongv1+i );
#endif

      if( i < loopend-FETCH_UP ){
	fat4 = &(t_dfatlink_du0[4*(i+FETCH_UP)]);
	prefetch_4MVVVV( 
              fat4,
	      (su3_vector *)gen_pt[XUP][i+FETCH_UP],
              (su3_vector *)gen_pt[YUP][i+FETCH_UP],
              (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
              (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
#ifndef NO_LONG_LINKS
	long4 = &(t_longlink[4*(i+FETCH_UP)]);
	prefetch_4MVVVV( 
              long4,
              (su3_vector *)gen_pt[X3UP][i+FETCH_UP],
              (su3_vector *)gen_pt[Y3UP][i+FETCH_UP],
              (su3_vector *)gen_pt[Z3UP][i+FETCH_UP],
              (su3_vector *)gen_pt[T3UP][i+FETCH_UP] );
        }
#endif
    } END_LOOP

    /* Wait gathers from negative directions, accumulate (negative) */
    /* and the same for the negative 3-rd neighbours */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    }
    for(dir=X3UP; dir<=T3UP; dir++){
	wait_gather(tag[OPP_3_DIR(dir)]);
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
        if( i < loopend-FETCH_UP ){
	  prefetch_VVVV( 
              (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
	  prefetch_VVVV( 
              (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
              (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
        }

        sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
#ifndef NO_LONG_LINKS
        sub_four_su3_vecs( templongv1+i,
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
        /* Now need to add these things together */
        add_su3_vector((su3_vector *)F_PT(s,dest), templongv1+i,
		       (su3_vector *)F_PT(s,dest));
#endif
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
			   ferm_links_t *fn, ferm_links_t *fn_dmdu0) {
   register int i;
   register site *s;
   register int dir,otherparity=0;
   msg_tag *tag[16];
   su3_vector *tempvec[4];
   register su3_matrix *fat4;
   su3_matrix *t_dfatlink_du0;
#ifndef NO_LONG_LINKS
   su3_vector *templongvec[4], *templongv1 ;
   register su3_matrix *long4;
   su3_matrix *t_longlink;
#endif
    
   for( dir=XUP; dir<=TUP; dir++ )
     {
       tempvec[dir]    =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
#ifndef NO_LONG_LINKS
       templongvec[dir]=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
#endif
     }
#ifndef NO_LONG_LINKS
   templongv1=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
#endif

   if(!fn->fl.valid){
     printf("ddslash_fn_du0_field: invalid fn links!\n");
     terminate(1);
   }
#ifndef NO_LONG_LINKS
   t_longlink = fn->fl.lng;
#endif
   if(!fn->fl.valid){
     printf("ddslash_fn_du0_field: invalid fn_du0 links!\n");
     terminate(1);
   }
   t_dfatlink_du0 = fn_dmdu0->fl.fat;

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
#ifndef NO_LONG_LINKS
     tag[DIR3(dir)] = start_gather_field( src, sizeof(su3_vector), 
					      DIR3(dir),parity, 
					      gen_pt[DIR3(dir)] );
#endif
   }

   /* Multiply by adjoint matrix at other sites */
   /* Use fat link for single link transport */
   FORSOMEPARITYDOMAIN( i, s, otherparity ){
     if( i < loopend-FETCH_UP ){
       fat4 = &(t_dfatlink_du0[4*(i+FETCH_UP)]);
       prefetch_V(&(src[i]));
       prefetch_4MVVVV(
		       fat4,
		       &(tempvec[0][i+FETCH_UP]),
		       &(tempvec[1][i+FETCH_UP]), 
		       &(tempvec[2][i+FETCH_UP]), 
		       &(tempvec[3][i+FETCH_UP])); 
       
#ifndef NO_LONG_LINKS
       long4 = &(t_longlink[4*(i+FETCH_UP)]);
       prefetch_4MVVVV(
		       long4,
		       &(templongvec[0][i+FETCH_UP]),
		       &(templongvec[1][i+FETCH_UP]), 
		       &(templongvec[2][i+FETCH_UP]), 
		       &(templongvec[3][i+FETCH_UP])); 
#endif
     }
     fat4 = &(t_dfatlink_du0[4*i]);
     mult_adj_su3_mat_4vec( fat4, &(src[i]), &(tempvec[0][i]),
			    &(tempvec[1][i]), &(tempvec[2][i]), 
			    &(tempvec[3][i]) );
#ifndef NO_LONG_LINKS
     /* multiply by 3-link matrices too */
     long4 = &(t_longlink[4*i]);
     mult_adj_su3_mat_4vec( long4, &(src[i]),&(templongvec[0][i]),
			    &(templongvec[1][i]), &(templongvec[2][i]), 
			    &(templongvec[3][i]) );
     for( dir=XUP; dir<=TUP; dir++ )
       scalar_mult_su3_vector( &(templongvec[dir][i]), -2.0/u0,
			       &(templongvec[dir][i]) );
#endif
   } END_LOOP

   /* Start gathers from negative directions */
   for( dir=XUP; dir <= TUP; dir++){
     tag[OPP_DIR(dir)] = start_gather_field( tempvec[dir],
	   sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
   }

#ifndef NO_LONG_LINKS
  /* Start 3-neighbour gathers from negative directions */
    for( dir=X3UP; dir <= T3UP; dir++){
      tag[OPP_3_DIR(dir)]=start_gather_field(templongvec[INDEX_3RD(dir)],
	sizeof(su3_vector), OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
    }
#endif
    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    /* wait for the 3-neighbours from positive directions, multiply */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
#ifndef NO_LONG_LINKS
	wait_gather(tag[DIR3(dir)]);
#endif
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
      if( i < loopend-FETCH_UP ){
	fat4 = &(t_dfatlink_du0[4*(i+FETCH_UP)]);
	prefetch_VV(
		   &(dest[i+FETCH_UP]),
		   &(templongv1[i+FETCH_UP]));
	prefetch_4MVVVV( 
		       fat4,
		       (su3_vector *)gen_pt[XUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[YUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
#ifndef NO_LONG_LINKS
	long4 = &(t_longlink[4*(i+FETCH_UP)]);
	prefetch_4MVVVV( 
		       long4,
		       (su3_vector *)gen_pt[X3UP][i+FETCH_UP],
		       (su3_vector *)gen_pt[Y3UP][i+FETCH_UP],
		       (su3_vector *)gen_pt[Z3UP][i+FETCH_UP],
		       (su3_vector *)gen_pt[T3UP][i+FETCH_UP] );
#endif       
      }
      fat4 = &(t_dfatlink_du0[4*i]);
      mult_su3_mat_vec_sum_4dir( fat4,
	    (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	    (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	    &(dest[i]) );

#ifndef NO_LONG_LINKS
      long4 = &(t_longlink[4*i]);
      mult_su3_mat_vec_sum_4dir( long4,
	    (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	    (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	    &(templongv1[i]));
      scalar_mult_su3_vector( &(templongv1[i]), -2.0/u0, &(templongv1[i]) );
#endif
    } END_LOOP
   
    /* Wait gathers from negative directions, accumulate (negative) */
    /* and the same for the negative 3-rd neighbours */
    for(dir=XUP; dir<=TUP; dir++){
      wait_gather(tag[OPP_DIR(dir)]);
    }
    for(dir=X3UP; dir<=T3UP; dir++){
      wait_gather(tag[OPP_3_DIR(dir)]);
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVV( 
		      (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
#ifndef NO_LONG_LINKS
	prefetch_VV(
		    &(dest[i+FETCH_UP]),
		    &(templongv1[i+FETCH_UP]));
	prefetch_VVVV( 
		      (su3_vector *)gen_pt[X3DOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[Y3DOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[Z3DOWN][i+FETCH_UP],
		      (su3_vector *)gen_pt[T3DOWN][i+FETCH_UP] );
#endif
      }
      sub_four_su3_vecs( &(dest[i]),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
#ifndef NO_LONG_LINKS
      sub_four_su3_vecs( &(templongv1[i]), 
	    (su3_vector *)(gen_pt[X3DOWN][i]),
	    (su3_vector *)(gen_pt[Y3DOWN][i]),
	    (su3_vector *)(gen_pt[Z3DOWN][i]),
	    (su3_vector *)(gen_pt[T3DOWN][i]) );
      /* Now need to add these things together */
      add_su3_vector(&(dest[i]), &(templongv1[i]),&(dest[i]));
#endif
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
#ifndef NO_LONG_LINKS
      free(templongvec[dir]);
#endif
    }
#ifndef NO_LONG_LINKS
    free(templongv1);
#endif
} /* end ddslash_fn_du0_field */

#endif /* DM_DU0 */
