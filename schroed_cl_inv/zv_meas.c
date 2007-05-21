/************************** zv_meas.c ****************************/
/* MIMD version 7 */
/* UMH August 99 */

#include "schroed_cl_includes.h"

void zv_meas(
    field_offset src1,	/* src1 is type wilson_propagator (forward) */
    field_offset src2,	/* src2 is type wilson_propagator (backward) */
    Real *f_1, complex *f_V, Real Kappa)
{

register int i;
register site *s;

int my_t, src_s, src_c, si, ci;
Real ftmp;
complex cc;

half_wilson_vector hwv1, hwv2;
wilson_vector *K_prop[3][4], *K_tmp, wv1, wtmp;


    for(ci=0;ci<3;ci++) for(si=0;si<4;si++){
	K_prop[ci][si] = (wilson_vector *)malloc(sizeof(wilson_vector));
    }
    K_tmp = (wilson_vector *)malloc(sizeof(wilson_vector));

    /* Construct K, gamma_5*K and f_1 */
    ftmp = Kappa / (Real)(nx*ny*nz);
    *f_1 = 0.0;
    for(src_c=0;src_c<3;src_c++) for(src_s=0;src_s<4;src_s++){
	clear_wvec( K_tmp);
	FORALLSITES(i,s) if(s->t == (nt-1)){
	    wp_shrink(
		&(((wilson_propagator *)F_PT(s,src1))->c[src_c].d[src_s]),
		&hwv1, TUP, MINUS);
	    mult_adj_su3_mat_hwvec( &(s->link[TUP]), &hwv1, &hwv2);
	    wp_grow_add( &hwv2, K_tmp, TUP, MINUS);
	}
	g_wvectorsumfloat( K_tmp);
	scalar_mult_wvec( K_tmp, ftmp, K_tmp);
	*f_1 += magsq_wvec( K_tmp);
	mult_by_gamma( K_tmp, K_prop[src_c][src_s], GAMMAFIVE);
    }

    *f_1 /= 8.0*Kappa*Kappa;

    FORALLSITES(i,s) if(s->t > 0){

	my_t = s->t;
	for(src_c=0;src_c<3;src_c++) for(src_s=0;src_s<4;src_s++){

	    /* K_tmp = src2 * gamma_5 * K_prop */
	    clear_wvec( K_tmp);
	    for(ci=0;ci<3;ci++) for(si=0;si<4;si++){
		c_scalar_mult_add_wvec( K_tmp,
		    &(((wilson_propagator *)F_PT(s,src2))->c[ci].d[si]),
		    &K_prop[src_c][src_s]->d[si].c[ci], K_tmp);
/*    &((wilson_vector *)(K_prop[src_c][src_s])->d[si].c[ci]), K_tmp); */
		    /* K_prop[src_c][src_s]->d[si].c[ci], K_tmp); */
	    }

	    /* wv1 = gamma_5 * gamma_0 * src1 */
	    mult_by_gamma(
		&(((wilson_propagator *)F_PT(s,src1))->c[src_c].d[src_s]),
		&wtmp, TUP);
	    mult_by_gamma( &wtmp, &wv1, GAMMAFIVE);

	    /* Now construct f_V */
	    /* gamma_0 is negative of usual chiral basis */
	    cc = wvec2_dot( K_tmp, &wv1);
	    f_V[my_t].real -= cc.real;
	    f_V[my_t].imag -= cc.imag;
	}
    }

    /* Normalize f_V */
    ftmp = 4.0 * Kappa;
    for(my_t=0; my_t<nt; my_t++){
	f_V[my_t].real /= ftmp;
	f_V[my_t].imag /= ftmp;
    }

    for(ci=0;ci<3;ci++) for(si=0;si<4;si++){
	free(K_prop[ci][si]);
    }
    free(K_tmp);

} /* zv_meas */
