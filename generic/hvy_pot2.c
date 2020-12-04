/****** hvy_pot.c  -- ******************/
/* Heavy quark potential
* MIMD version 7
* DT 6/10/97
* DT 7/17/97 modified to do all displacements
* CD 9/24/09 converted site variables to field variables
* AB 10/6/12 added Polyakov loop correlators
* Evaluate in different spatial directions, to check rotational
* invariance.  Gauge fix to Coulomb gauge to do spatial segments.
*
* Argument tells whether to use ordinary links or fat links.
* Gauge must be fixed to Coulomb gauge first.  DO THIS BEFORE CALLING!
*/
/* Suggest defining max_t = (nt/6)  maximum time value for loops */
/* and  max_x = (nx/2-2)	maximum spatial distance */

/* Two new versions implemented:
* 1) uses general gather to do x-y shifts all at once
* 2) assigns additional su3_mat valued field buffer that  
*    contains x-y shifted Wilson lines; 
*    x-shifts performed as general gathers (except x_dist==1),
*    y-shifts performed sequentially as single step gathers, 
*    after z-shifts: y-shifted field recovered from buffer
*
* Typical speedup is nearly 40% less execution time. JHW 2/24/15
*/

#define BUFMAT2

#include "generic_includes.h"	/* definitions files and prototypes */
void shiftmat( su3_matrix *src, su3_matrix *dest, int dir );

#ifdef BUFMAT0
#define NEWCODE
#endif
#ifdef BUFMAT1
#define NEWCODE
#define BUFMAT
#endif
#ifdef BUFMAT2
#define NEWCODE
#define BUFMAT
#endif

#ifdef NEWCODE
void general_shiftmat( su3_matrix *src, su3_matrix *dest, int *disp );
#endif

void hvy_pot( su3_matrix *links, int max_t, int max_x ) {
    register int i;
    int ngather=0;
    int t_dist, x_dist, y_dist, z_dist;
    double wloop;
    msg_tag *mtag0;
    su3_matrix *tempmat1, *tempmat2, *staple, *oldmat, *newmat, *tt;
    char myname[] = "hvy_pot";
#ifdef PLOOPCOR_MEAS
    double ploopcor;
    complex ctr1, ctr2, cprod;
#endif
#ifdef NEWCODE
    int disp[4]={0,0,0,0};
#endif
#ifdef BUFMAT
    su3_matrix *bufmat1=NULL,*buffer1=NULL;
#ifdef BUFMAT2
    su3_matrix *bufmat2=NULL,*buffer2=NULL;
#endif
#endif

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

#ifdef BUFMAT
    buffer1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
#ifdef BUFMAT2
    buffer2 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if (buffer1 == NULL || buffer2 == NULL )
#else
    if (buffer1 == NULL )
#endif
    {
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }
#endif
    /* Use tempmat1 to construct t-direction path from each point */
    for( t_dist=1; t_dist<= max_t; t_dist ++){
	if(t_dist==1 ){
          FORALLFIELDSITES(i){ su3mat_copy( links+4*i+TUP, tempmat1+i ); }
	}
	else{
	    mtag0 = start_gather_field( tempmat1, sizeof(su3_matrix),
	        TUP, EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);
	    FORALLFIELDSITES(i){
	      mult_su3_nn( links+4*i+TUP, (su3_matrix *)gen_pt[0][i], staple+i );
	    }
	    cleanup_gather( mtag0 ); ngather++;
	    FORALLFIELDSITES(i){ su3mat_copy( staple+i, tempmat1+i ); }
	}
	/* now tempmat1 is path of length t_dist in TUP direction */
	oldmat = tempmat2;
	newmat = staple;	/* will switch these two */
#ifdef BUFMAT
 	 bufmat1 = buffer1;
#ifdef BUFMAT2
 	 bufmat2 = buffer2;
#endif
#endif

	for( x_dist=0; x_dist<=max_x; x_dist++ ){
	    /**for( y_dist=0; y_dist<=x_dist; y_dist+= x_dist>0?x_dist:1 ){**/
	    for( y_dist=0; y_dist<=max_x; y_dist+= 1 ){
		/* now gather from spatial dirs, compute products of paths */
#ifndef NEWCODE
		FORALLFIELDSITES(i){
		    su3mat_copy( tempmat1+i, oldmat+i );
		}
		for(i=0;i<x_dist;i++){
		    shiftmat( oldmat, newmat, XUP ); ngather++;
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
		for(i=0;i<y_dist;i++){
		    shiftmat( oldmat, newmat, YUP ); ngather++;
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
#else
#ifdef BUFMAT0
                // uses general gather 
                FORALLFIELDSITES(i){ su3mat_copy( tempmat1+i, oldmat+i ); }
                disp[XUP]=x_dist; disp[YUP]=y_dist;
                general_shiftmat( oldmat, newmat, disp ); ngather++;
                tt=oldmat; oldmat=newmat; newmat=tt;
#endif //BUFMAT0
// minimal shifts, extra memory
#ifdef BUFMAT
                // shift matrix by x_dist steps only, if y_dist loop starts from zero, otherwise keep buffered
                if (y_dist==0) {
                  // put x-shifted wline into bufmat
                  // shift oldmat by x_dist steps in negative xdir
                  disp[XUP]=x_dist;
#ifdef BUFMAT2
                  if (x_dist>=1) { 
                    // shift bufmat2 by one further step in x_dist
                    shiftmat( bufmat2, bufmat1, XUP ); ngather++;
                    tt=bufmat2; bufmat2=bufmat1; bufmat1=tt;
                  } else {
                    // put wline into bufmat2 
                    FORALLFIELDSITES(i){ su3mat_copy( tempmat1+i, bufmat2+i ); }
                  } 
                  // copy shifted wline from bufmat2 to bufmat1
                  FORALLFIELDSITES(i){ su3mat_copy( bufmat2+i, bufmat1+i ); }
#else
                  // put wline into oldmat
                  FORALLFIELDSITES(i){ su3mat_copy( tempmat1+i, oldmat+i ); }
                  if (x_dist>=1) { 
                    if (x_dist>1) {
                      // shift oldmat by x_dist steps 
                      general_shiftmat( oldmat, newmat, disp ); ngather++;
                    } else { 
                      // shift bufmat2 by just one step in x_dist
                      shiftmat( oldmat, newmat, XUP ); ngather++;
                    }
                    tt=oldmat; oldmat=newmat; newmat=tt;
                  }
                  // copy shifted wline from oldmat to bufmat1
                  FORALLFIELDSITES(i){ su3mat_copy( oldmat+i, bufmat1+i ); }
#endif
                } else if (y_dist>0) {
                  // shift x-shifted buffer one further step in y_dist, keep xy-shifted matrix in buffer
                  shiftmat( bufmat1, newmat, YUP ); ngather++;
                  tt=bufmat1; bufmat1=newmat; newmat=tt;
                }
#endif //BUFMAT
#endif //ifndef NEWCODE
		for( z_dist=0; z_dist<=max_x; z_dist++ ){
		    /* evaluate potential at this separation */
		    wloop = 0.0;
#ifdef PLOOPCOR_MEAS
                    ploopcor = 0.0;
#endif
		    FORALLFIELDSITES(i){
#ifndef BUFMAT
			wloop += (double)realtrace_su3( tempmat1+i,
			    oldmat+i );
#else
			if (z_dist>0) { wloop += (double)realtrace_su3( tempmat1+i,oldmat+i ); }
			else { wloop += (double)realtrace_su3( tempmat1+i,bufmat1+i ); }
#endif //ifndef BUFMAT
#ifdef PLOOPCOR_MEAS
                        ctr1 = trace_su3( tempmat1+i );
#ifndef BUFMAT
                        ctr2 = trace_su3( oldmat+i );
#else
			if (z_dist>0) { ctr2 = trace_su3( oldmat+i ); }
			else {  ctr2 = trace_su3( bufmat1+i ); }
#endif //ifndef BUFMAT
                        CMUL_J( ctr1, ctr2, cprod );
                        ploopcor += cprod.real;
#endif
		    }
		    g_doublesum( &wloop );
		    node0_printf("POT_LOOP: %d %d %d %d \t%.6e\n",
			x_dist, y_dist, z_dist, t_dist, wloop/volume );
#ifdef PLOOPCOR_MEAS
                    g_doublesum( &ploopcor );
                    node0_printf("POL_CORR: %d %d %d %d \t%.6e\n",
                        x_dist, y_dist, z_dist, t_dist, ploopcor/volume );
#endif
		    /* as we increment z, shift in z direction */
#ifndef BUFMAT
		    shiftmat( oldmat, newmat, ZUP ); ngather++;
		    tt=oldmat; oldmat=newmat; newmat=tt;
#else
		    if (z_dist>0) { shiftmat( oldmat, newmat, ZUP ); ngather++; }
                    else  { shiftmat( bufmat1, newmat, ZUP ); ngather++;}
		    tt=oldmat; oldmat=newmat; newmat=tt;
#endif //ifndef BUFMAT
		} /* z distance */
	    } /* y distance */
	} /*x distance */
    } /*t_dist*/

#ifdef BUFMAT
    free(buffer1);
#ifdef BUFMAT2
    free(buffer2);
#endif
#endif
    free(staple);
    free(tempmat2);
    free(tempmat1);
#ifdef DEBUG
node0_printf("NUMBER OF GATHERS = %d\n",ngather);
#endif
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

#ifdef NEWCODE
/* shift, without parallel transport, a matrix from displacement "disp" */
void general_shiftmat( su3_matrix *src, su3_matrix *dest, int *disp ){
    register int i;
    msg_tag *mtag;
    mtag = start_general_gather_field( (void *)src, sizeof(su3_matrix),
        disp, EVENANDODD, gen_pt[0] );
    wait_general_gather(mtag);
    FORALLFIELDSITES(i){
        su3mat_copy( (su3_matrix *)gen_pt[0][i], dest+i );
    }
    cleanup_general_gather( mtag );
}
#endif
