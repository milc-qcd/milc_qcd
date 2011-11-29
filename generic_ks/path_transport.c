/****** path_transport.c  -- ******************/
/* MIMD version 7 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/********************************************************************/
/* parallel transport a vector to the current site along a path.
   For example, if the path is "XUP", bring in the vector from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   vectors rather than stuff in lattice.h */
/********************************************************************/
void path_transport_field( su3_vector *src, su3_vector *dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_dest = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_work = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = src[i];
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_vec( &(s->link[dir[j]]),
		    (su3_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_vec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        dest[i] = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        dest[i] = src[i];
    }
  }
} /* path_transport_field */


/********************************************************************/
/* transport a matrix to the current site along a path.
   This matrix has gauge transformation properties of a "connection" along the
   path --- the column index transforms at the far end, and the row index
   at the "current site".   Thus, at each step, only one matrix multiplication
   is required.  At the end, result is gauge covariant "at the current site".
   For example, if the path is "XUP", bring in the matrix from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   matrix rather than stuff in lattice.h */
/********************************************************************/
void path_transport_connection( su3_matrix *src, su3_matrix *dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_matrix *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_matrix *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tmp_dest = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tmp_work = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = src[i];
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src, sizeof(su3_matrix),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_nn( &(s->link[dir[j]]),
		    (su3_matrix *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_su3_an( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work, sizeof(su3_matrix),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_matrix *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        dest[i] = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        dest[i] = src[i];
    }
  }
} /* path_transport_connection */


/* THE FOLLOWING PROCEDURE SEEMS TO BE DISUSED - CD 7/11*/
void path_transport_connection_hisq( su3_matrix *src, su3_matrix **links, su3_matrix *dest,
    int parity, int *dir, int length ){
//TEST ME
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_matrix *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_matrix *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tmp_dest = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tmp_work = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = src[i];
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src, sizeof(su3_matrix),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_nn( &(links[dir[j]][i]), (su3_matrix *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_su3_an( &(links[OPP_DIR(dir[j])][i]), &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work, sizeof(su3_matrix),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_matrix *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        dest[i] = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        dest[i] = src[i];
    }
  }
} /* path_transport_connection_hisq */

// special case to transport a "connection" by one link, does both parities
void link_transport_connection( su3_matrix *src, su3_matrix *dest,
  su3_matrix *work, int dir ){
    register int i;
    register site *s;
    msg_tag *mtag0;

    if( GOES_FORWARDS(dir) ) {
	mtag0 = start_gather_field( src, sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir]), (su3_matrix *)(gen_pt[0][i]),
		&(dest[i]) );
	}
	cleanup_gather(mtag0);
    }

    else{ /* GOES_BACKWARDS(dir) */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[OPP_DIR(dir)]),
		&(src[i]), &(work[i]) );
	}
	mtag0 = start_gather_field( work, sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);
	FORALLSITES(i,s){
	    dest[i] = *(su3_matrix *)gen_pt[0][i];
	}
	cleanup_gather(mtag0);
    }
} /* link_transport_connection */
// special case to transport a "connection" by one link, does both parities
void link_transport_connection_hisq( su3_matrix *src, su3_matrix *links, su3_matrix *dest,
  su3_matrix *work, int dir ){
//TEST ME
    register int i;
    register site *s;
    msg_tag *mtag0;

    if( GOES_FORWARDS(dir) ) {
	mtag0 = start_gather_field( src, sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);
	FORALLSITES(i,s){
	  //	    mult_su3_nn( &(links[dir][i]), (su3_matrix *)(gen_pt[0][i]), &(dest[i]) );
	    mult_su3_nn( &(links[4*i+dir]), (su3_matrix *)(gen_pt[0][i]), &(dest[i]) );
	}
	cleanup_gather(mtag0);
    }

    else{ /* GOES_BACKWARDS(dir) */
	FORALLSITES(i,s){
	  //	    mult_su3_an( &(links[OPP_DIR(dir)][i]), &(src[i]), &(work[i]) );
	    mult_su3_an( &(links[4*i+OPP_DIR(dir)]), &(src[i]), &(work[i]) );
	}
	mtag0 = start_gather_field( work, sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);
	FORALLSITES(i,s){
	    dest[i] = *(su3_matrix *)gen_pt[0][i];
	}
	cleanup_gather(mtag0);
    }
} /* link_transport_connection_hisq */
// like link_transport, except doesn't multiply by link matrices.  use this, for example,
// when storing the intermediate HISQ force (a connection) at the lattice site
// associated with a link
void link_gather_connection_hisq( su3_matrix *src, su3_matrix *dest,
  su3_matrix *work, int dir ){
//TEST ME
    register int i;
    register site *s;
    msg_tag *mtag0;

    mtag0 = start_gather_field( src, sizeof(su3_matrix),
	dir, EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES(i,s){ dest[i] = *(su3_matrix *)(gen_pt[0][i]); }
    cleanup_gather(mtag0);
} /* link_gather_connection_hisq */




/********************************************************************/
/* Path transport a half_wilson_vector */
/********************************************************************/
void path_transport_hwv_field( half_wilson_vector *src, half_wilson_vector * dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    half_wilson_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    half_wilson_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (half_wilson_vector *)malloc(
      sites_on_node*sizeof(half_wilson_vector) );
    tmp_dest = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );
    tmp_work = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = src[i];
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_hwvec( &(s->link[dir[j]]),
		    (half_wilson_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_hwvec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(half_wilson_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        dest[i] = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        dest[i] = src[i];
    }
  }
} /* path_transport_hwv_field */

/********************************************************************/
/* parallel transport a vector to the current site along a path.
   For example, if the path is "XUP", bring in the vector from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   vectors rather than stuff in lattice.h */
/********************************************************************/
void path_transport( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );
    tmp_dest = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );
    tmp_work = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(su3_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_vec( &(s->link[dir[j]]),
		    (su3_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_vec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
    }
  }
} /* path_transport */

/********************************************************************/
/* Path transport a half_wilson_vector */
/********************************************************************/
void path_transport_hwv( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    half_wilson_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    half_wilson_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (half_wilson_vector *)special_alloc(
      sites_on_node*sizeof(half_wilson_vector) );
    tmp_dest = (half_wilson_vector *)special_alloc(
       sites_on_node*sizeof(half_wilson_vector) );
    tmp_work = (half_wilson_vector *)special_alloc(
       sites_on_node*sizeof(half_wilson_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(half_wilson_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_hwvec( &(s->link[dir[j]]),
		    (half_wilson_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_hwvec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(half_wilson_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) =
	  *(half_wilson_vector *)F_PT(s,src);
    }
  }
} /* path_transport_hwv_field */

/********************************************************************/
/* parallel transport a vector to the current site along a path.
   For example, if the path is "XUP", bring in the vector from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   vectors rather than stuff in lattice.h */
/********************************************************************/
void path_transport_site( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );
    tmp_dest = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );
    tmp_work = (su3_vector *)special_alloc( sites_on_node*sizeof(su3_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(su3_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_vec( &(s->link[dir[j]]),
		    (su3_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_vec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
    }
  }
} /* path_transport_site */

/********************************************************************/
/* Path transport a half_wilson_vector */
/********************************************************************/
void path_transport_hwv_site( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    half_wilson_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    half_wilson_vector *tmp_pt; /* scratch */
    int tmp_parity=0, tmp_otherparity=0; /* parity for this step */

  if( length > 0 ){
    tmp_src = (half_wilson_vector *)special_alloc(
      sites_on_node*sizeof(half_wilson_vector) );
    tmp_dest = (half_wilson_vector *)special_alloc(
       sites_on_node*sizeof(half_wilson_vector) );
    tmp_work = (half_wilson_vector *)special_alloc(
       sites_on_node*sizeof(half_wilson_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(half_wilson_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_field( tmp_src,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_hwvec( &(s->link[dir[j]]),
		    (half_wilson_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_hwvec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_field( tmp_work,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(half_wilson_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) =
	  *(half_wilson_vector *)F_PT(s,src);
    }
  }
} /* path_transport_hwv_site */

/* path_transport.c */
