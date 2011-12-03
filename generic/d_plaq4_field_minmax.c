/******************* d_plaq4_field_minimax.c **********************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */

/* AB 11/27/7 - measurement of average and min/max plaquette
   for the field */

/* Measure the average and min/max plaquette of the space-space and
   space-time plaquettes */

#include "generic_includes.h"

void d_plaquette_field_minmax(su3_matrix **U_field,
       double *ss_plaq,double *st_plaq,
       double *ss_plaq_min, double *st_plaq_min,
       double *ss_plaq_max, double *st_plaq_max) {
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register site *s;
register int first_pass_s,first_pass_t;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ss_sum,st_sum;
double rtrace_s, rtrace_t, ss_min, st_min, ss_max, st_max;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
    first_pass_s=1;
    first_pass_t=1;

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_field( U_field[dir2], sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_field( U_field[dir1], sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(U_field[dir1][i]);
		m4 = &(U_field[dir2][i]);
		mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP ) {
                  rtrace_t = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  st_sum += rtrace_t;
                }
		else {
                  rtrace_s = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  ss_sum += rtrace_s;
                }
//                printf("Plaq i=%d, dir1=%d, dir2=%d: %f %f\n",
//                  i,dir1,dir2,rtrace_s,rtrace_t);
                /* set min and max values on the first pass */
                if( dir1==TUP ) {
                  if( 1==first_pass_t ) {
                    st_min = rtrace_t;
                    st_max = rtrace_t;
                    first_pass_t = 0;
                  }
                  else {
                    if( rtrace_t < st_min ) st_min = rtrace_t;
                    if( rtrace_t > st_max ) st_max = rtrace_t;
                  }
                }
                else {
                  if( 1==first_pass_s ) {
                    ss_min = rtrace_s;
                    ss_max = rtrace_s;
                    first_pass_s = 0;
                  }
                  else {
                    if( rtrace_s < ss_min ) ss_min = rtrace_s;
                    if( rtrace_s > ss_max ) ss_max = rtrace_s;
                  }
                }
            }

            cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    /* DEBUGGING: output min/max values on each node before global operation */
    printf("MIN PLAQ (W links) on node %d:\t%f\t%f\n", mynode(), ss_min, st_min );
    printf("MAX PLAQ (W links) on node %d:\t%f\t%f\n", mynode(), ss_max, st_max );

    /* to find min: flip sign, find maximum, flip sign */
    ss_min = -ss_min;
    st_min = -st_min;
    g_doublemax( &ss_min );
    g_doublemax( &st_min );
    g_doublemax( &ss_max );
    g_doublemax( &st_max );

    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));

    *ss_plaq_min = -ss_min;
    *st_plaq_min = -st_min;
    *ss_plaq_max = ss_max;
    *st_plaq_max = st_max;

    free(su3mat);
} /* d_plaquette4 */

