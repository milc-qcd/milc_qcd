/******************* d_plaq4_hist.c *******************************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */

/* AB 11/28/7 - measurement of (3-plaquette) histogram,
                histogram range is given for each power */

#include "generic_includes.h"

void d_plaquette_hist(int Npowers, int *Nhist, double **hist,
                      double **hist_bounds,
                      double *ss_plaq, double *st_plaq) {
/* measure (3-plaquette)^(n+1), where n=0,...,Npower-1,
   *Nhist - array with histogram sizes,
   **hist - first index is power "n",
   *hist_bounds - ranges for histogram */
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register int ipower,ihist;
register site *s;
register su3_matrix *m1,*m4;
double *plaq_power, *step_hist;
su3_matrix mtmp;
double ss_sum,st_sum;
double rtrace, rtrace3;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
#ifdef HISQ_DUMP_PLAQ_INTO_FILE
FILE *fp;
char plaq_file_name[300];
#endif /* HISQ_DUMP_PLAQ_INTO_FILE */

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    /* zero out the histogram */
    for(ipower=0;ipower<Npowers;ipower++) {
      for(ihist=0;ihist<Nhist[ipower];ihist++) {
        hist[ipower][ihist]=0.0;
      }
    }

    /* array with powers of (3-plaquette) */
    plaq_power=(double*)malloc(sizeof(double)*Npowers);

    /* array with step sizes */
    step_hist=(double*)malloc(sizeof(double)*Npowers);
    for(ipower=0;ipower<Npowers;ipower++) {
      step_hist[ipower]=
        (hist_bounds[ipower][1]-hist_bounds[ipower][0])/Nhist[ipower];
    }

#ifdef HISQ_DUMP_PLAQ_INTO_FILE
    sprintf( plaq_file_name, "plaq_fund_node%04d.dat", this_node );
    fp = fopen( plaq_file_name, "wt" );
#endif /* HISQ_DUMP_PLAQ_INTO_FILE */

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
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP ) {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  st_sum += rtrace;
                }
		else {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  ss_sum += rtrace;
                }
//                printf("Plaq i=%d, dir1=%d, dir2=%d: %f %f\n",
//                  i,dir1,dir2,rtrace_s,rtrace_t);
#ifdef HISQ_DUMP_PLAQ_INTO_FILE
                fprintf( fp, "%18.12g\n", rtrace );
#endif /* HISQ_DUMP_PLAQ_INTO_FILE */
                /* powers of (3-plaquette) */
                rtrace3=3.0-rtrace;
                plaq_power[0]=rtrace3;
                for(ipower=1;ipower<Npowers;ipower++) {
                  plaq_power[ipower]=plaq_power[ipower-1]*rtrace3;
                }
                /* find histogram entry */
                for(ipower=0;ipower<Npowers;ipower++) {
                  if( (plaq_power[ipower]>hist_bounds[ipower][0])
                   && (plaq_power[ipower]<hist_bounds[ipower][1]) ) {
                    ihist=(int)( 
                      (plaq_power[ipower]-hist_bounds[ipower][0])/
                      step_hist[ipower] );
                    hist[ipower][ihist]+=1.0;
                  }
                }
	    }

	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }

    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));

    /* global gather histograms */
    for(ipower=0;ipower<Npowers;ipower++) {
      plaq_power[ipower]=0;
      for(ihist=0;ihist<Nhist[ipower];ihist++) {
        g_doublesum( &(hist[ipower][ihist]) );
        plaq_power[ipower]+=hist[ipower][ihist];
      }
    }
    /* normalize histograms */
    for(ipower=0;ipower<Npowers;ipower++) {
      for(ihist=0;ihist<Nhist[ipower];ihist++) {
        hist[ipower][ihist]/=(plaq_power[ipower]*step_hist[ipower]);
      }
    }

#ifdef HISQ_DUMP_PLAQ_INTO_FILE
    fclose( fp );
#endif /* HISQ_DUMP_PLAQ_INTO_FILE */

    free(plaq_power);
    free(step_hist);
    free(su3mat);
} /* d_plaquette4 */

