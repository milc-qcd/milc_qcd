/********************** path_product.c ***************************/
/* MIMD version 6 */
/* Compute product of links along a specified path */
/* On return lattice[i].tempmat1 contains the product for the path
   ENDING at site i.  e.g., for the 1-link path XUP,
   lattice[i].tempmat1 = lattice[i-\hat x][XUP] */

#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/prefetch.h"
#define FETCH_UP 1

/* LOOPEND is required now -CD */
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

void path_product( const int *dir, const int length) {
    register int i;
    register site *s;
    msg_tag *mtag0;
    su3_matrix *tempmat2t, *tempmat3t;
    int j;
    /* a forward step leaves the answer in gen_pt[0], which points into
	link, tempmat1 or tempmat2, and backwards step in tempmat1 or tempmat2,
	After a forwards step, need to wait and clean a gather.
	  STEP	leaves answer in
	  even # forward	gen_pt[0]->tempmat1  (gen_pt[0]->link for step 0
	  even # backward	tempmat1
	  odd  # forward	gen_pt[0]->tempmat2
	  odd  # backward	tempmat2
	*/

    /* Trivial path case */
    if(length == 0){
      FORALLSITES(i,s){
	clear_su3mat(&(s->tempmat1));
	s->tempmat1.e[0][0].real = s->tempmat1.e[1][1].real 
	  = s->tempmat1.e[2][2].real = 1.;
      } END_LOOP
      return;
    }

    /* allocate temporary space */
    tempmat3t = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tempmat2t = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );

    /* j=0 */
    if( GOES_FORWARDS(dir[0]) )  {
	mtag0 = start_gather( F_OFFSET(link[dir[0]]), sizeof(su3_matrix),
	    OPP_DIR(dir[0]), EVENANDODD, gen_pt[0] );
    }
    else{  /* if GOES_BACKWARDS(dir[0]) */
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &((s+FETCH_UP)->tempmat1) );
	  }
	    su3_adjoint(&(s->link[OPP_DIR(dir[0])]),&(s->tempmat1) );
	} END_LOOP
    }

    for(j=1;j<length;j++) {
	if( j%2==1 ){
	    if( GOES_FORWARDS(dir[j]) ) {
	      if( GOES_FORWARDS(dir[j-1]) ){
	        wait_gather(mtag0);
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		    prefetch_M(  &(tempmat2t[i+FETCH_UP]) );
		  }
		  mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir[j]]),
		    &(tempmat2t[i]) );
	        } END_LOOP
	        cleanup_gather(mtag0);
	      }
	      else{ /* last link was backwards */
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( &((s+FETCH_UP)->tempmat1) );
		    prefetch_M(  &(tempmat2t[i+FETCH_UP]) );
		  }
		  mult_su3_nn( &(s->tempmat1),&(s->link[dir[j]]),
		    &(tempmat2t[i]) );
	        } END_LOOP
	      }
	      mtag0 = start_gather_from_temp( tempmat2t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }  /* for GOES_FORWARDS */

	    else{ /* GOES_BACKWARDS(dir[j]), which is an odd numbered step */
	      if( GOES_FORWARDS(dir[j-1]) ){
	        wait_gather(mtag0);
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  }
	          su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(tempmat3t[i]) );
	        } END_LOOP
	        cleanup_gather(mtag0);
	        mtag0 = start_gather_from_temp( tempmat3t, sizeof(su3_matrix),
		  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
	        mtag0 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
		  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	      }
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  prefetch_M(  &((s+FETCH_UP)->link[OPP_DIR(dir[j])]) );
		}
		  mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		    &(s->link[OPP_DIR(dir[j])]), &(tempmat2t[i]) );
	      } END_LOOP
	      cleanup_gather(mtag0);
	    } /* end for GOES_BACKWARDS */
	} /* end for j=odd */

	else{	/* j=even */
	  if( GOES_FORWARDS(dir[j]) ) {
	    if( GOES_FORWARDS(dir[j-1]) ){
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  prefetch_M(  &((s+FETCH_UP)->link[dir[j]]) );
		}
		mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir[j]]),
		    &(s->tempmat1) );
	      } END_LOOP
	      cleanup_gather(mtag0);
	    }
	    else{ /* last link goes backwards */
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( &(tempmat2t[i+FETCH_UP]) );
		  prefetch_M(  &((s+FETCH_UP)->link[dir[j]]) );
		}
		mult_su3_nn( &(tempmat2t[i]),&(s->link[dir[j]]),
		    &(s->tempmat1) );
	      } END_LOOP
	    }
	    mtag0 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	  }  /* for GOES_FORWARDS */

	  else{ /* GOES_BACKWARDS(dir[j]) */
	    if( GOES_FORWARDS(dir[j-1]) ){
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  prefetch_M(  &(tempmat3t[i+FETCH_UP]) );
		}
	        su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(tempmat3t[i]) ); 
	      } END_LOOP
	      cleanup_gather(mtag0);
	      mtag0 = start_gather_from_temp( tempmat3t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }
	    else{ /* last step was backwards */
	      mtag0 = start_gather_from_temp( tempmat2t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }
	    wait_gather(mtag0);
	    FORALLSITES(i,s){
	      if( i < loopend-FETCH_UP ){
		prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		prefetch_M(  &((s+FETCH_UP)->link[OPP_DIR(dir[j])]) );
	      }
	      mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		    &(s->link[OPP_DIR(dir[j])]), &(s->tempmat1) );
	    } END_LOOP
	    cleanup_gather(mtag0);
	  } /* for GOES_BACKWARDS */
	} /* for j=even */

    }  /* j=link in loop */

    /* Want to end in tempmat1 */
    if( length%2==0 ){  /* last step was odd */
      if( GOES_FORWARDS(dir[length-1]) ){
	wait_gather(mtag0);
	  FORALLSITES(i,s){
	    if( i < loopend-FETCH_UP ){
	      prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
	    }
	    su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(s->tempmat1) ); 
	} END_LOOP
	cleanup_gather(mtag0);
      }
      else{
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &(tempmat2t[i+FETCH_UP]) );
	  }
	  su3mat_copy(&(tempmat2t[i]),&(s->tempmat1) );
	} END_LOOP
      }
    }
    else{ /* odd length path */
      if( GOES_FORWARDS(dir[length-1]) ){
	wait_gather(mtag0);
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
	  }
	  su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(tempmat3t[i]) );
	} END_LOOP
	cleanup_gather(mtag0);
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &(tempmat3t[i+FETCH_UP]) );
	  }
	  su3mat_copy( &(tempmat3t[i]), &(s->tempmat1) );
	} END_LOOP
      }
      else{
      }
    }
    free(tempmat3t);
    free(tempmat2t);
} /* path */

#ifdef N_SUBL32
/* code from symanzik_sl32/dsdu_qhb.c ****************************/
/* U.M. Heller August 1997 */

/* This is a modification of "path_product" from gauge_stuff.c
   which works only on one sublattice. */
void path_prod_subl(const int *dir, const int length, const int subl)
{
register int i;
register site *s;
msg_tag *mtag0;
su3_matrix *tempmat2t, *tempmat3t;
int j, nsubl;

    /* A forward step leaves the answer in gen_pt[0], which points into
       link, tempmat1 or tempmat2, and backwards step in tempmat1 or tempmat2.
       After a forwards step, need to wait and clean a gather.
	STEP			leaves answer in
	even # forward		gen_pt[0]->tempmat1 (gen_pt[0]->link for step 0)
	even # backward		tempmat1
	odd  # forward		gen_pt[0]->tempmat2
	odd  # backward		tempmat2
    */

    /* allocate temporary space */
    tempmat3t = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
    tempmat2t = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );

    /* j=0 */
    if( GOES_FORWARDS(dir[0]) ) {
	nsubl = neighsubl[subl][dir[0]];
	mtag0 = start_gather( F_OFFSET(link[dir[0]]), sizeof(su3_matrix),
		OPP_DIR(dir[0]), nsubl, gen_pt[0] );
    }
    else{  /* if GOES_BACKWARDS(dir[0]) */
	nsubl = neighsubl[subl][dir[0]];
	FORSOMESUBLATTICE(i,s,nsubl){
	    su3_adjoint(&(s->link[OPP_DIR(dir[0])]), &(s->tempmat1) );
	}
    }

    for(j=1;j<length;j++) {
	if( j%2==1 ){
	    if( GOES_FORWARDS(dir[j]) ) {
	      if( GOES_FORWARDS(dir[j-1]) ){
		wait_gather(mtag0);
		FORSOMESUBLATTICE(i,s,nsubl){
		  mult_su3_nn( (su3_matrix *)(gen_pt[0][i]),
		    &(s->link[dir[j]]), &(tempmat2t[i]) );
		}
		cleanup_gather(mtag0);
	      }
	      else{ /* last link was backwards */
		FORSOMESUBLATTICE(i,s,nsubl){
		  mult_su3_nn( &(s->tempmat1),
		    &(s->link[dir[j]]), &(tempmat2t[i]) );
		}
	      }
	      nsubl = neighsubl[nsubl][dir[j]];
	      mtag0 = start_gather_from_temp( tempmat2t, sizeof(su3_matrix),
		      OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	    }  /* for GOES_FORWARDS */

	    else{ /* GOES_BACKWARDS(dir[j]), which is an odd numbered step */
	      if( GOES_FORWARDS(dir[j-1]) ){
		wait_gather(mtag0);
		FORSOMESUBLATTICE(i,s,nsubl){
		  su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat3t[i]) );
		}
		cleanup_gather(mtag0);
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather_from_temp( tempmat3t, sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      wait_gather(mtag0);
	      FORSOMESUBLATTICE(i,s,nsubl){
		mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		  &(s->link[OPP_DIR(dir[j])]), &(tempmat2t[i]) );
	      }
	      cleanup_gather(mtag0);
	    } /* end for GOES_BACKWARDS */
	} /* end for j=odd */

	else{   /* j=even */
	    if( GOES_FORWARDS(dir[j]) ) {
	      if( GOES_FORWARDS(dir[j-1]) ){
		wait_gather(mtag0);
		FORSOMESUBLATTICE(i,s,nsubl){
		  mult_su3_nn( (su3_matrix *)(gen_pt[0][i]),
		    &(s->link[dir[j]]), &(s->tempmat1) );
		}
		cleanup_gather(mtag0);
	      }
	      else{ /* last link was backwards */
		FORSOMESUBLATTICE(i,s,nsubl){
		  mult_su3_nn( &(tempmat2t[i]),
		    &(s->link[dir[j]]), &(s->tempmat1) );
		}
	      }
	      nsubl = neighsubl[nsubl][dir[j]];
	      mtag0 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
		      OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	    }  /* for GOES_FORWARDS */

	    else{ /* GOES_BACKWARDS(dir[j]), which is an even numbered step */
	      if( GOES_FORWARDS(dir[j-1]) ){
		wait_gather(mtag0);
		FORSOMESUBLATTICE(i,s,nsubl){
		  su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat3t[i]) );
		}
		cleanup_gather(mtag0);
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather_from_temp( tempmat3t, sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather_from_temp( tempmat2t, sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      wait_gather(mtag0);
	      FORSOMESUBLATTICE(i,s,nsubl){
		mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		  &(s->link[OPP_DIR(dir[j])]), &(s->tempmat1) );
	      }
	      cleanup_gather(mtag0);
	    } /* end for GOES_BACKWARDS */
	} /* end for j=even */

    }  /* j=link in loop */

    /* Want to end in tempmat1 */
    if( length%2==0 ){  /* last step was odd */
	if( GOES_FORWARDS(dir[length-1]) ){
	    wait_gather(mtag0);
	    FORSOMESUBLATTICE(i,s,nsubl){
	      su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat1) );
	    }
	    cleanup_gather(mtag0);
	}
	else{
	    FORSOMESUBLATTICE(i,s,nsubl){
	      su3mat_copy(&(tempmat2t[i]), &(s->tempmat1) );
	    }
	}
    }
    else{ /* odd length path: last step was even */
	if( GOES_FORWARDS(dir[length-1]) ){
	    wait_gather(mtag0);
	    FORSOMESUBLATTICE(i,s,nsubl){
	      su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat3t[i]) );
	    }
	    cleanup_gather(mtag0);
	    FORSOMESUBLATTICE(i,s,nsubl){
	      su3mat_copy(&(tempmat3t[i]), &(s->tempmat1) );
	    }
	}
    }

    free(tempmat3t);
    free(tempmat2t);

} /* path_prod_subl */

#endif /* N_SUBL32 */
