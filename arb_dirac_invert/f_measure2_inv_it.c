/**************** f_measure2_inv_it.c ***************************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */

/* Wilson fermions */

/* Measure fermionic observables:

   Output is an F2MES line containing the expectation value of
   the real part of Tr(1/M),
   the imaginary part of  Tr(1/M),
   the real part of Tr(gamma_5/M),
   the imaginary part of Tr(gamma_5/M).

And for eigenvalues, the volume averages of etabar eta and etabar gamma-5 eta.

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

/* code for direct tests */
complex num1,num2,den1,den2,eig;
int invert,flagg;

MaxCG=niter;
RsdCG=resid[0];

if(flag==0){

    /* gaussian random vector */
FORALLSITES(i,s){
        for(k=0;k<4;k++)for(j=0;j<3;j++){
            s->source.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->source.d[k].c[j].imag = gaussian_rand_no(&node_prn);
	}
    }

FORALLSITES(i,s){
s->psi=s->source;
}


} /* flag */
eig.real=1.0;

    /* new code for source overlaps */



   for(invert=0;invert<9;invert++){
if(this_node==0)printf("invert=%d\n",invert);



/* begin by setting chi_n = psi_{n-1}, psi_n=chi_n/eig */

FORALLSITES(i,s){
s->chi=s->psi;
}     
eig.real=1.0/eig.real;
   FORALLSITES(i,s) {
                scalar_mult_wvec(&(s->psi), eig.real,
                   &(s->psi) );
	}

flagg=flag;
if(invert>3)flagg=1; /* to iterate from a good solution */


printf("flagg=%d\n",flagg);

    /* Invert, result in psi */
#ifdef BI
            iters = bicgstab(F_OFFSET(chi),F_OFFSET(psi),
                                      MaxCG,RsdCG,&size_r,flagg);
#else

            iters = congrad_t(F_OFFSET(chi),F_OFFSET(psi),
                                      MaxCG,RsdCG,&size_r,flagg);
#endif
        if(this_node==0)printf("size_r= %e, iters= %e\n",(double)size_r,
                (double)iters);

    eig = cmplx(0.0,0.0);
    den1 = cmplx(0.0,0.0);

delta0( F_OFFSET(psi), F_OFFSET(chi), PLUS);

    /* den1 = psi.psi */
    FORALLSITES(i,s){
        cc = wvec_dot( &(s->psi), &(s->psi) );
	CSUM(den1,cc);
        cc = wvec_dot( &(s->psi), &(s->chi) );
	CSUM(eig,cc); 
    }
    g_complexsum( &den1 );
    g_complexsum( &eig );

    CDIVREAL(eig,den1.real,eig);

    if(this_node==0)printf("F2EIG %e %e\n",
       (double)eig.real, (double)eig.imag);

    den1.real=1.0/den1.real;
    den1.real=sqrt(den1.real);

        FORALLSITES(i,s) {
                scalar_mult_wvec(&(s->psi), den1.real,
                   &(s->psi) );
	}
    pbg5p = cmplx(0.0,0.0);

    /* psi is normalized. Check for chirality */

    FORALLSITES(i,s){
        cc = wvec_dot( &(s->psi), &(s->psi) );
        CSUM(pbp,cc);
        mult_by_gamma( &(s->psi),&wv0,  GAMMAFIVE);
        cc = wvec_dot( &(s->psi), &wv0 );
        CSUM(pbg5p,cc);
    }

   if(this_node==0)printf("F2CHI %e %e\n",
       (double)pbg5p.real, (double)pbg5p.imag);



   }

    return(iters);
}

