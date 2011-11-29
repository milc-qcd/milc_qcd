/************************** gball_simp.c *******************************/
/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* version of 9/4/94 by UMH */
/* port to V5 CD */

/* Computes space-like plaquettes, summed over each time-slice,
   as simple glueball operators. */

#include "hvy_qpot_includes.h"

void gball_simp(int tot_smear) {
register int i,dir1,dir2,n,t;
int n_op,norm_vol,my_t;
register site *s;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
msg_tag *mtag0,*mtag1;
Real ftmp1,ftmp2;
complex *plaq_op, plaq;

    n_op = 3;
    norm_vol = 3*nx*ny*nz;

    plaq_op = (complex *)malloc(nt*n_op*sizeof(complex));
    for(n=0;n<n_op;n++) for(t=0;t<nt;t++){
	plaq_op[t+nt*n] = cmplx(0.0,0.0);
    }

    n = 0;
    for(dir1=XUP;dir1<=YUP;dir1++){
	for(dir2=dir1+1;dir2<=ZUP;dir2++){

	    mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->link[dir1]);
		m4 = &(s->link[dir2]);
		mult_su3_an(m4,m1,&(s->staple));
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_nn( &(s->staple), (su3_matrix *)(gen_pt[0][i]),
			     &mtmp);

		plaq = complextrace_su3( (su3_matrix *)(gen_pt[1][i]), &mtmp);

		my_t = s->t;
		plaq_op[my_t+n].real += plaq.real;
		plaq_op[my_t+n].imag += plaq.imag;
	    }

	    n += nt;

	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);

	}
    }

    /* Normalize and print the glueball operators */
    for(n=0;n<n_op;n++) for(t=0;t<nt;t++){
	ftmp1 = plaq_op[t+nt*n].real;
	g_floatsum( &ftmp1);
	ftmp1 /= (Real)(norm_vol);
	ftmp2 = plaq_op[t+nt*n].imag;
	g_floatsum( &ftmp2);
	ftmp2 /= (Real)(norm_vol);
	if(this_node == 0)printf("GBALL_OP1_%d  %d  %d  %e  %e\n",
	     tot_smear, n, t, (double)ftmp1, (double)ftmp2);
    }

    free( plaq_op);

} /* gball_simp */

