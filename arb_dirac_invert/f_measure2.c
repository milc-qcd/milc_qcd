/**************** f_measure2.c ***************************************/

/* MIMD version 4 */
/* Wilson fermions */

/* Measure fermionic observables:

   Output is an F2MES line containing the expectation value of
   the real part of Tr(1/M),
   the imaginary part of  Tr(1/M),
   the real part of Tr(gamma_5/M),
   the imaginary part of Tr(gamma_5/M).

*/

#include "arb_dirac_inv_includes.h"
#include <string.h>

int f_measure2(int flag) {
/* local variables for accumulators */
register int i,j,k,dir;
register site *s;
msg_tag *tag0,*tag1;
register complex cc;
complex pbp;
complex pbg5p;
int iters;
wilson_vector wv0;
int MinCG,MaxCG;
Real size_r,RsdCG;
int x,y,z,t;


MaxCG=niter;
RsdCG=resid[0];

    /* gaussian random vector */
FORALLSITES(i,s){
        for(k=0;k<4;k++)for(j=0;j<3;j++){
            s->source.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->source.d[k].c[j].imag = gaussian_rand_no(&node_prn);
	}
    }

FORALLSITES(i,s){
s->chi=s->source;
}


    /* Invert, result in psi */
            /* compute the propagator.  Result in psi. */
#ifdef BI
            iters = bicgstab(F_OFFSET(chi),F_OFFSET(psi),
                                      MaxCG,RsdCG,&size_r,flag);
#else
iters = congrad_t(F_OFFSET(chi),F_OFFSET(psi),
                                      MaxCG,RsdCG,&size_r,flag);
#endif
        if(this_node==0)printf("size_r= %e, iters= %e\n",(double)size_r,
                (double)iters);


/*Temporary*/
/* Multiply by M and see if I get source back */
/* use dir as flag*/
/**
delta0( F_OFFSET(psi), F_OFFSET(mp), PLUS);
FORALLSITES(i,s){
    for(dir=0,j=0;j<4;j++)for(k=0;k<3;k++){
	if(s->source.d[j].c[k].real - s->mp.d[j].c[k].real > 2e-5 )dir=1;
	if(s->source.d[j].c[k].imag - s->mp.d[j].c[k].imag > 2e-5 )dir=1;
	if(dir)printf("%d %d %d  ( %.4e , %.4e )  ( %.4e , %.4e )\n",
	    i,j,k,s->source.d[j].c[k].real,s->source.d[j].c[k].imag,
	    s->mp.d[j].c[k].real,s->mp.d[j].c[k].imag);
    }
} *End temporary **/

    pbp = cmplx(0.0,0.0);
    pbg5p = cmplx(0.0,0.0);

    /* psi-bar-psi = source.psi */
    /* psi-bar-gamma-5 psi = source. gamma-5 psi */
    FORALLSITES(i,s){
        cc = wvec_dot( &(s->source), &(s->psi) );
	CSUM(pbp,cc);
        mult_by_gamma( &(s->psi),&wv0,  GAMMAFIVE);
        cc = wvec_dot( &(s->source), &wv0 );
	CSUM(pbg5p,cc);
    }

    g_complexsum( &pbp );
    g_complexsum( &pbg5p );

    CDIVREAL(pbp,volume,pbp);
    CDIVREAL(pbg5p,volume,pbg5p);
    if(this_node==0)printf("F2MES %e %e %e %e\n",
       (double)pbp.real, (double)pbp.imag,(double)pbg5p.real,
 (double)pbg5p.imag);
    fflush(stdout);
    return(iters);
}

