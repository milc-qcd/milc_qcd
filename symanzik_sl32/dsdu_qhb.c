/****** dsdu_qhb.c  -- compute the staple ******************/
/* This is a version for extended actions where 32 sublattices are
   needed to make the links independent. */
/* MIMD version 6 */
/* U.M. Heller August 1997 */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#include "symanzik_sl32_includes.h"
#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include "gauge_action.h"
#undef GAUGE_ACTION_PART1
extern	int loop_length[NLOOP];	/* lengths of various kinds of loops */
extern	int loop_num[NLOOP];	/* number of rotations/reflections  for each kind */
    /* table of directions, for each rotation and reflection of each kind of
       loop.  tabulated with "canonical" starting point and direction. */
extern	int loop_table[NLOOP][MAX_NUM][MAX_LENGTH];
    /* table of coefficients in action, for various "representations" (actually,
       powers of the trace) */
extern	Real loop_coeff[NLOOP][NREPS];

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void path_prod_subl(int *dir, int length, int subl);

void dsdu_qhb(int dir, int subl)
{
register site *st;
register int i;
int iloop, ln, k, j;
int dirs[MAX_LENGTH], length;
int path_dir[MAX_LENGTH], path_length;
su3_matrix tmat1;
int fsubl;


    FORSOMESUBLATTICE(i,st,subl) {
	clear_su3mat(&(st->staple));
    }

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	for(ln=0;ln<loop_num[iloop];ln++){
	    /* set up dirs.  we are looking at loop starting in "XUP"
	       direction, rotate so it starts in "dir" direction. */
	    for(k=0;k<length;k++){
		if( GOES_FORWARDS(loop_table[iloop][ln][k]) ){
		    dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
		}
		else {
		    dirs[k]=OPP_DIR(
			(dir+OPP_DIR(loop_table[iloop][ln][k]))%4 );
		}
	    }

	    path_length = length-1;	/* generalized "staple" */
	    /* The path starts at the forward end of the link */
	    fsubl = neighsubl[subl][dir];

	    /* check for links in direction of link to be updated.
	       Note the direction of the path - opposite the link. */
	    for(k=0;k<length;k++)if( dirs[k]==dir||dirs[k]==OPP_DIR(dir)) {
		if( GOES_FORWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[j] = dirs[(k+j+1)%length];
		}
		if( GOES_BACKWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[path_length-1-j] =
			OPP_DIR(dirs[(k+j+1)%length]);
		}
		path_prod_subl(path_dir, path_length, fsubl);

		/* We took the path in the other direction from our old
		   convention in order to get it to end up "at our site".
		   So now take adjoint */
		FORSOMESUBLATTICE(i,st,subl) {
		    su3_adjoint(&(st->tempmat1), &tmat1 );
		    scalar_mult_add_su3_matrix(&(st->staple), &tmat1,
			loop_coeff[iloop][0], &(st->staple) );
		}
	    } /* k (location in path) */
	} /* ln */
    } /* iloop */

    g_sync();

} /* dsdu_qhb */

/* This is a modification of "path_product" from gauge_stuff.c
   which works only on one sublattice. */
void path_prod_subl(int *dir, int length, int subl)
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
	mtag0 = start_gather_site( F_OFFSET(link[dir[0]]), sizeof(su3_matrix),
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
	      mtag0 = start_gather_field( tempmat2t, sizeof(su3_matrix),
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
		mtag0 = start_gather_field( tempmat3t, sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
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
	      mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
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
		mtag0 = start_gather_field( tempmat3t, sizeof(su3_matrix),
			OPP_DIR(dir[j]), nsubl, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
		nsubl = neighsubl[nsubl][dir[j]];
		mtag0 = start_gather_field( tempmat2t, sizeof(su3_matrix),
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

