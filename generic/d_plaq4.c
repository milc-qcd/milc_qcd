/************************** d_plaq4.c *******************************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */

/* Double precision version of "plaquette4.c" including optional
   Schroedinger functional - UMH - 1/27/00 */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

#include "generic_includes.h"

void d_plaquette(double *ss_plaq,double *st_plaq) {
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register site *s;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ss_sum,st_sum;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->link[dir1]);
		m4 = &(s->link[dir2]);
		mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
#ifdef SCHROED_FUN
		if(dir1==TUP ){
		    if(s->t==(nt-1)){
			mult_su3_nn( &su3mat[i],
			    &(s->boundary[dir2]), &mtmp);
		    }
		    else{
			mult_su3_nn( &su3mat[i],
			    (su3_matrix *)(gen_pt[0][i]), &mtmp);
		    }
		    st_sum +=
			realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
		}
		else if(s->t > 0){
		    mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
			&mtmp);
		    ss_sum +=
			realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
		}
#else
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP )st_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
		else          ss_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
#endif
	    }

	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
#ifdef SCHROED_FUN
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*(nt-1)));
#else
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
#endif
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));

    free(su3mat);
} /* d_plaquette4 */

