/************************** coupling.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */

/* Measure the coupling in the Schroedinger functional, or more precisely
   k/g^2, with k=12*nx^2*{sin(\theta)+sin(2\theta)}, \theta=\pi/(3*nx*nt) */

#include "generic_schroed_includes.h"

void coupling(double *ds_deta, double *bd_plaq) {
register int i,dir;
register site *s;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ds_sum, plq_sum;
msg_tag *mtag0,*mtag1,*mtag2;
    ds_sum = 0.0;
    plq_sum = 0.0;
    for(dir=XUP;dir<TUP;dir++){

	mtag0 = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
	    TUP, EVENANDODD, gen_pt[0] );
	mtag1 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[1] );
	mtag2 = start_gather_site( F_OFFSET(boundary[dir]), sizeof(su3_matrix),
	    OPP_DIR(TUP), EVENANDODD, gen_pt[2] );

	/* Note: the derivative of the boundary fields at t=0 are in
	         boundary[dir] at t=0, while the derivative for the
	         boundary fields at t=nt are in boundary[dir] at t=nt-2! */

	FORALLSITES(i,s){
	    if(s->t == 0){
		m1 = &(s->boundary[dir]);
		m4 = &(s->link[TUP]);
		mult_su3_an(m4, m1, &(s->tempmat1));
		m1 = &(s->link[TUP]);
		m4 = &(s->link[dir]);
		mult_su3_an(m4, m1, &(s->staple));
	    }
	    else if(s->t == (nt-1)){
		m1 = &(s->link[TUP]);
		m4 = &(s->link[dir]);
		mult_su3_an(m4, m1, &(s->staple));
	    }
	}

	wait_gather(mtag0);
	wait_gather(mtag1);
	wait_gather(mtag2);

	FORALLSITES(i,s){
	    if(s->t == 0){
		mult_su3_nn( &(s->tempmat1),
		    (su3_matrix *)(gen_pt[1][i]), &mtmp);
		ds_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[0][i]), &mtmp);
		mult_su3_nn( &(s->staple),(su3_matrix *)(gen_pt[0][i]), &mtmp);
		plq_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
	    }
	    else if(s->t == (nt-1)){
		mult_su3_nn( &(s->staple),(su3_matrix *)(gen_pt[2][i]), &mtmp);
		ds_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
		mult_su3_nn( &(s->staple), &(s->boundary[dir]), &mtmp);
		plq_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
	    }
	}

	cleanup_gather(mtag0);
	cleanup_gather(mtag1);
	cleanup_gather(mtag2);
    }
    g_doublesum( &ds_sum );
    g_doublesum( &plq_sum );
    ds_sum *= (double)(beta/3.0);
    *ds_deta = -ds_sum;
    *bd_plaq = plq_sum /((double)(6*nx*ny*nz));
} /* coupling */

