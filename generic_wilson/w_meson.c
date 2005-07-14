/*********** w_meson.c *************/
/* MIMD version 7 */
/* UMH April 96 
 2/27/98 Simplify and avoid bug in Origin compiler UMH
 4/03/00 Split out meson_cont.c as a separate file CD 
 */

#include "generic_wilson_includes.h"

void w_meson(field_offset src1,field_offset src2,complex *prop[10])
/* src1 and src2 are type spin_wilson_vector */
{

int gamma_in[4],gamma_out[4];
int n_in,n_out;
void meson_cont(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);

/* gamma_in = source Dirac matrix, gamma_out = sink Dirac matrix */

/* note (gamma_5 gamma_mu)^\dagger = gamma_5 gamma_mu */

    /* PION */
    n_in=1;n_out=1;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= GAMMAFIVE;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[0]);


    /* PS505 */
    n_in=1;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= TUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[1]);


    /* PS055 */
    n_in=2;n_out=1;
    gamma_in[0]= TUP;
    gamma_in[1]= GAMMAFIVE;
    gamma_out[0]= GAMMAFIVE;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[2]);


    /* PS0505 */
    n_in=2;n_out=2;
    gamma_in[0]= TUP;
    gamma_in[1]= GAMMAFIVE;
    gamma_out[0]= TUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[3]);


    /* RHO33 (we use also the other polarizations) */
    n_in=1;n_out=1;
    gamma_in[0]= ZUP;
    gamma_out[0]= ZUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[0]= XUP;
    gamma_out[0]= XUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[0]= YUP;
    gamma_out[0]= YUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);


    /* RHO0303 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= TUP;
    gamma_in[1]= ZUP;
    gamma_out[0]= TUP;
    gamma_out[1]= ZUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);

    gamma_in[1]= XUP;
    gamma_out[1]= XUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);

    gamma_in[1]= YUP;
    gamma_out[1]= YUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);


    /* SCALAR */
    n_in=0;n_out=0;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[6]);


    /* SCALA0 */
    n_in=1;n_out=1;
    gamma_in[0]= TUP;
    gamma_out[0]= TUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[7]);


    /* PV35 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_in[1]= ZUP;
    gamma_out[0]= ZUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);

    gamma_in[1]= XUP;
    gamma_out[0]= XUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);

    gamma_in[1]= YUP;
    gamma_out[0]= YUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);


    /* B12 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= XUP;
    gamma_in[1]= YUP;
    gamma_out[0]= YUP;
    gamma_out[1]= XUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

    gamma_in[0]= YUP;
    gamma_in[1]= ZUP;
    gamma_out[0]= ZUP;
    gamma_out[1]= YUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

    gamma_in[0]= ZUP;
    gamma_in[1]= XUP;
    gamma_out[0]= XUP;
    gamma_out[1]= ZUP;
    meson_cont(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

} /* w_meson */



