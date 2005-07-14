/******** meson_cont.c *************/
/* MIMD version 7 */
/* UMH April 96 */

#include "generic_wilson_includes.h"

void meson_cont(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop) 
/* src1 and src2 are type spin_wilson_vector */
{

register int i;
register site *s; 

int my_t;
int cf, sf, si;
int i_in,i_out;

complex g1,g2;

spin_wilson_vector localmat;       /* temporary storage for quark */
spin_wilson_vector quark;          /* temporary storage for quark */
spin_wilson_vector antiquark;      /* temporary storage for antiquark */



    FORALLSITES(i,s){
	my_t = s->t;

	/*first, dirac multiplication by the source gamma matrices (on left) */

	/*  antiquark = c.c. of quark propagator */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++){
	    CONJG(((spin_wilson_vector *)F_PT(s,src1))->d[si].d[sf].c[cf],
		antiquark.d[si].d[sf].c[cf]);
         }

	/* left multiply antiquark by source gamma matrices,
	   beginning with gamma_5 for quark -> antiquark */
	mult_swv_by_gamma_l( &antiquark, &localmat, GAMMAFIVE);

	/* right dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_swv_by_gamma_r( &localmat, &antiquark, GAMMAFIVE);


	/* left multiply by the particular source dirac matrices */
	for(i_in=0;i_in<n_in;i_in++){
	    mult_swv_by_gamma_l( &antiquark, &localmat, gamma_in[i_in]);
	    antiquark = localmat;
	}

	/* right dirac multiplication by the sink gamma matrices */
	for(i_out=0;i_out<n_out;i_out++) {
	    mult_swv_by_gamma_r( &antiquark, &localmat, gamma_out[i_out]);
	    antiquark = localmat;
	}

	/* copy into quark */
	/* 2/27/98 Simplify and avoid bug in Origin compiler UMH */
        quark = *(spin_wilson_vector *)F_PT(s,src2);

	/* trace over propagators */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++)
	{
	    g1 = antiquark.d[si].d[sf].c[cf];
	    g2 = quark.d[si].d[sf].c[cf];

	    prop[my_t].real += (g1.real*g2.real - g1.imag*g2.imag); 
	    prop[my_t].imag += (g1.real*g2.imag + g1.imag*g2.real); 
	}

    }

} /* meson_cont */

