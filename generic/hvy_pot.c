/****** hvy_pot.c  -- ******************/
/* Heavy quark potential
* MIMD version 7
* DT 6/10/97
* DT 7/17/97 modified to do all displacements
* CD 9/24/09 converted site variables to field variables
* Evaluate in different spatial directions, to check rotational
* invariance.  Gauge fix to Coulomb gauge to do spatial segments.
*
* Argument tells whether to use ordinary links or fat links.
* Gauge must be fixed to Coulomb gauge first.  DO THIS BEFORE CALLING!
*/
/* Suggest defining max_t = (nt/6)  maximum time value for loops */
/* and  max_x = (nx/2-2)	maximum spatial distance */

#include "generic_includes.h"	/* definitions files and prototypes */
void shiftmat( su3_matrix *src, su3_matrix *dest, int dir );

void hvy_pot( su3_matrix *links, int max_t, int max_x ) {
    register int i;
    int t_dist, x_dist, y_dist, z_dist;
    double wloop;
    msg_tag *mtag0;
    su3_matrix *tempmat1, *tempmat2, *staple, *oldmat, *newmat, *tt;
    char myname[] = "hvy_pot";

    node0_printf("hvy_pot(): MAX_T = %d, MAX_X = %d\n",max_t,max_x);

    tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }

    tempmat2 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat2 == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }

    staple = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(staple == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }

    /* Use tempmat1 to construct t-direction path from each point */
    for( t_dist=1; t_dist<= max_t; t_dist ++){
	if(t_dist==1 ){
	    FORALLFIELDSITES(i){
//		su3mat_copy( &(((su3_matrix *)(F_PT(s,links)))[TUP]),
//		    tempmat1+i );
	      su3mat_copy( links+4*i+TUP, tempmat1+i );
	    }
	}
	else{
	    mtag0 = start_gather_field( tempmat1, sizeof(su3_matrix),
	        TUP, EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);
	    FORALLFIELDSITES(i){
//	 	mult_su3_nn( &(((su3_matrix *)(F_PT(s,links)))[TUP]),
//		    (su3_matrix *)gen_pt[0][i], staple+i );
	      mult_su3_nn( links+4*i+TUP, (su3_matrix *)gen_pt[0][i], staple+i );
	    }
	    cleanup_gather( mtag0 );
	    FORALLFIELDSITES(i){
	 	su3mat_copy( staple+i, tempmat1+i );
	    }
	}
	/* now tempmat1 is path of length t_dist in TUP direction */
	oldmat = tempmat2;
	newmat = staple;	/* will switch these two */

	for( x_dist=0; x_dist<=max_x; x_dist++ ){
	    /**for( y_dist=0; y_dist<=x_dist; y_dist+= x_dist>0?x_dist:1 ){**/
	    for( y_dist=0; y_dist<=max_x; y_dist+= 1 ){
		/* now gather from spatial dirs, compute products of paths */
		FORALLFIELDSITES(i){
		    su3mat_copy( tempmat1+i, oldmat+i );
		}
		for(i=0;i<x_dist;i++){
		    shiftmat( oldmat, newmat, XUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
		for(i=0;i<y_dist;i++){
		    shiftmat( oldmat, newmat, YUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
		for( z_dist=0; z_dist<=max_x; z_dist++ ){
		    /* evaluate potential at this separation */
		    wloop = 0.0;
		    FORALLFIELDSITES(i){
			wloop += (double)realtrace_su3( tempmat1+i,
			    oldmat+i );
		    }
		    g_doublesum( &wloop );
		    node0_printf("POT_LOOP: %d %d %d %d \t%e\n",
			x_dist, y_dist, z_dist, t_dist, wloop/volume );

		    /* as we increment z, shift in z direction */
		    shiftmat( oldmat, newmat, ZUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		} /* z distance */
	    } /* y distance */
	} /*x distance */
    } /*t_dist*/

    free(staple);
    free(tempmat2);
    free(tempmat1);
} /* hvy_pot() */

/* shift, without parallel transport, a matrix from direction "dir" */
void shiftmat( su3_matrix *src, su3_matrix *dest, int dir ){
    register int i;
    msg_tag *mtag;
    mtag = start_gather_field( src, sizeof(su3_matrix),
        dir, EVENANDODD, gen_pt[0] );
    wait_gather(mtag);
    FORALLFIELDSITES(i){
        su3mat_copy( (su3_matrix *)gen_pt[0][i], dest+i );
    }
    cleanup_gather( mtag );
}
