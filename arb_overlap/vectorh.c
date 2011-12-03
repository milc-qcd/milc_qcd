/************ vectorh.c *******************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

complex dotproduct(wilson_vector * w1, wilson_vector * w2)
{
    register int i;
    register site *s;
    complex ctmp;
    double_complex dctmp;

    dctmp = dcmplx(0., 0.);
    FORALLSITES(i, s) {
	ctmp = wvec_dot(&w1[i], &w2[i]);
	CSUM(dctmp, ctmp);
    }

    g_dcomplexsum(&dctmp);
    ctmp.real = dctmp.real;
    ctmp.imag = dctmp.imag;
    return ctmp;
}

Real vectornorm(wilson_vector * w1)
{
    register int i;
    register site *s;
    double dtmp;

    dtmp = 0.;
    FORALLSITES(i, s)
	dtmp += magsq_wvec(&w1[i]);
    g_doublesum(&dtmp);
    return (Real) sqrt(dtmp);
}

void project_out_eigen(wilson_vector* vec, wilson_vector** evec, int nvec)
{
    register int i, ivec;
    register site* s;
    double_complex dctmp,dc2;
    complex ctmp;



    for (ivec=0;ivec<nvec;ivec++)
    {
      dctmp=dcmplx(0.0,0.0);
	FORALLSITES(i,s)
	{
	    ctmp = wvec_dot(&(evec[ivec][i]), &vec[i]);
	    dc2.real=(double)ctmp.real; dc2.imag = (double)ctmp.imag;
	    CSUM(dctmp, dc2);
	}
        g_dcomplexsum(&dctmp);
/*
printf(" OVERLAP %d %e %e\n",ivec,dctmp.real,dctmp.imag);
*/
	/*
	CNEGATE(dctmp,ctmp);
	*/
	ctmp.real= -(Real)dctmp.real;
	ctmp.imag= -(Real)dctmp.imag;

	FORALLSITES(i, s)
	        c_scalar_mult_add_wvec(&(vec[i]), &(evec[ivec][i]), &ctmp,
	                            &(vec[i]));
    }
}
