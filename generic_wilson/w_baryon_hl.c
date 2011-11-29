/********** w_baryon_hl.c *************/
/* MIMD version 7 */

/* CD  April 07 switched from field_offset to wilson_prop_field  */
/* UMH April 96 */

/* Heavy-light baryon spectrum for Wilson hadrons */

/* Construct baryon propagators for the Sigma^+, the Sigma^{*+} and
   the Lambda with non-degenerate "u" and "s" quarks. For the Sigma^+,
   which becomes the Proton for degenerate quarks, we take

   |S, s_z=1/2> = (s C gamma_5 u) "u_up"

   for the Sigma^{*+}, which becomes the Delta^+ for degenerate quarks,

   |S*, s_z=3/2> = 2*(s C gamma_- u) "u_up" + (u C gamma_- u) "s_up"

   and for the Lambda, for which we consider the "u" and "d" quark
   as degenerate!

   |L, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
                  + (u C gamma_5 s) "d_up".

   We have put "q_up" in quotes, since this is meant in the Dirac basis,
   not in the 'DeGrand-Rossi' chiral basis used in the program!
   A spin-up quark in the Dirac basis corresponds to
   1/sqrt(2) * ( - q_1 - q_3 ) in this chiral basis. We shall neglect
   the sign and the 1/sqrt(2) here and will use spin 1 (q_1) only. In
   Dirac basis this corresponds to an equal mixure of forward propagating
   baryon and backward propagating anti-baryon.

   For all baryons we compute a 'B0' that differs from the 'B' above
   by insertion of a gamma_0 between C and the gamma_{5,-}.
*/

#include "generic_wilson_includes.h"

#define Nc 3
#define Ns 4

void w_baryon_hl(wilson_prop_field *src1, wilson_prop_field *src2,
		 wilson_prop_field *src3, complex *prop[6]) 
{

  int i,j,k,t;
  
  int eps[3][3][3], chi_b[4][4];
  complex *prop_tmp, *prop_l;

  if(src1->nc != 3 || src2->nc != 3 || src3->nc != 3){
    node0_printf("w_baryon_hl: requires 3 colors for each field\n");
    terminate(1);
  }

    prop_tmp = (complex *)malloc(nt*sizeof(complex));
    prop_l = (complex *)malloc(nt*sizeof(complex));

    /* initialize epsilon factors */
    for( i=0; i< 3; i++ )
    for( j=0; j< 3; j++ )
    for( k=0; k< 3; k++ ){
	eps[i][j][k]=0;
    }
    eps[0][1][2]=1;eps[1][2][0]=1;eps[2][0][1]=1;
    eps[0][2][1]= -1;eps[1][0][2]= -1;eps[2][1][0]= -1;

    /* PROTON (used also for LAMBDA) */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][0] = 1;
    chi_b[0][1] =  -1;
    chi_b[3][2] = 1;
    chi_b[2][3] =  -1; 

    for(t=0;t<nt;t++){
	prop_l[t] = cmplx(0.0,0.0);
	prop_tmp[t] = cmplx(0.0,0.0);
    }

    baryon_cont1(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);

    for(t=0;t<nt;t++){
	prop[0][t].real += prop_l[t].real + prop_tmp[t].real;
	prop[0][t].imag += prop_l[t].imag + prop_tmp[t].imag;
	prop_l[t].real -= prop_tmp[t].real;
	prop_l[t].imag -= prop_tmp[t].imag;
    }

    /* LAMBDA */
    baryon_cont1(src2, src3, src1, chi_b, eps, prop[4]);
    baryon_cont2(src2, src3, src1, chi_b, eps, prop[4]);
    baryon_cont2(src2, src1, src3, chi_b, eps, prop[4]);

    for(t=0;t<nt;t++){
	prop[4][t].real += 0.5*prop_l[t].real;
	prop[4][t].imag += 0.5*prop_l[t].imag;
    }

    /* PROTON0 (used also for LAMBDA0) */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[0][3] = 1;
    chi_b[1][2] =  -1;
    chi_b[2][1] = 1;
    chi_b[3][0] =  -1; 

    for(t=0;t<nt;t++){
	prop_l[t] = cmplx(0.0,0.0);
	prop_tmp[t] = cmplx(0.0,0.0);
    }

    baryon_cont1(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);

    for(t=0;t<nt;t++){
	prop[1][t].real += prop_l[t].real + prop_tmp[t].real;
	prop[1][t].imag += prop_l[t].imag + prop_tmp[t].imag;
	prop_l[t].real -= prop_tmp[t].real;
	prop_l[t].imag -= prop_tmp[t].imag;
    }

    /* LAMBDA0 */
    baryon_cont1(src2, src3, src1, chi_b, eps, prop[5]);
    baryon_cont2(src2, src3, src1, chi_b, eps, prop[5]);
    baryon_cont2(src2, src1, src3, chi_b, eps, prop[5]);

    for(t=0;t<nt;t++){
	prop[5][t].real += 0.5*prop_l[t].real;
	prop[5][t].imag += 0.5*prop_l[t].imag;
    }

    /* DELTA */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][3] =  -1; 
    chi_b[3][1] =  -1;

    for(t=0;t<nt;t++)prop_l[t] = cmplx(0.0,0.0);
    baryon_cont1(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src2, src3, src1, chi_b, eps, prop_l);
    baryon_cont2(src2, src1, src3, chi_b, eps, prop_l);

    for(t=0;t<nt;t++)prop_tmp[t] = cmplx(0.0,0.0);
    baryon_cont1(src2, src3, src1, chi_b, eps, prop_tmp);

    for(t=0;t<nt;t++){
	prop[2][t].real += (2.0*prop_l[t].real+prop_tmp[t].real)/3.0;
	prop[2][t].imag += (2.0*prop_l[t].imag+prop_tmp[t].imag)/3.0;
    }

    /* DELTA0 */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][1] = 1;
    chi_b[3][3] = 1;

    for(t=0;t<nt;t++)prop_l[t] = cmplx(0.0,0.0);
    baryon_cont1(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_l);
    baryon_cont2(src2, src3, src1, chi_b, eps, prop_l);
    baryon_cont2(src2, src1, src3, chi_b, eps, prop_l);

    for(t=0;t<nt;t++)prop_tmp[t] = cmplx(0.0,0.0);
    baryon_cont1(src2, src3, src1, chi_b, eps, prop_tmp);

    for(t=0;t<nt;t++){
	prop[3][t].real += (2.0*prop_l[t].real+prop_tmp[t].real)/3.0;
	prop[3][t].imag += (2.0*prop_l[t].imag+prop_tmp[t].imag)/3.0;
    }

    free(prop_l);
    free(prop_tmp);

}  /* w_baryon_hl */

#undef Nc
#undef Ns
