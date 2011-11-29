/********** w_baryon.c *************/
/* MIMD version 7 */

/* CD  April 07 switched from field_offset to wilson_propagator * */
/* UMH April 96 */

/* Baryon spectrum for Wilson hadrons */

/* Construct baryon propagators for the Proton and the Delta^+ with
   degenerate "u" and "d" quarks. For the Proton we take

   |P, s_z=1/2> = (d C gamma_5 u) "u_up"

   and for the Delta^+

   |D, s_z=3/2> = 2*(d C gamma_- u) "u_up" + (u C gamma_- u) "d_up".

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

void w_baryon(wilson_prop_field *src1, wilson_prop_field *src2,
	      wilson_prop_field *src3, complex *prop[4]) 
{

  int i,j,k,t;
  
  int eps[3][3][3], chi_b[4][4];
  complex *prop_tmp;
  
  if(src1->nc != 3 || src2->nc != 3 || src3->nc != 3){
    node0_printf("w_baryon: requires 3 colors for each field\n");
    terminate(1);
  }

    prop_tmp = (complex *)malloc(nt*sizeof(complex));

    /* initialize epsilon factors */
    for( i=0; i< 3; i++ )
    for( j=0; j< 3; j++ )
    for( k=0; k< 3; k++ ){
	eps[i][j][k]=0;
    }
    eps[0][1][2]=1;eps[1][2][0]=1;eps[2][0][1]=1;
    eps[0][2][1]= -1;eps[1][0][2]= -1;eps[2][1][0]= -1;

    /* PROTON */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][0] = 1;
    chi_b[0][1] =  -1;
    chi_b[3][2] = 1;
    chi_b[2][3] =  -1; 

    baryon_cont1(src1, src2, src3, chi_b, eps, prop[0]);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop[0]);

    /* PROTON0 */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[0][3] = 1;
    chi_b[1][2] =  -1;
    chi_b[2][1] = 1;
    chi_b[3][0] =  -1; 

    baryon_cont1(src1, src2, src3, chi_b, eps, prop[1]);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop[1]);

    /* DELTA */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][3] =  -1; 
    chi_b[3][1] =  -1;

    baryon_cont1(src1, src2, src3, chi_b, eps, prop[2]);

    for(t=0;t<nt;t++)prop_tmp[t] = cmplx(0.0,0.0);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);
    for(t=0;t<nt;t++){
	prop[2][t].real += 2.0*prop_tmp[t].real;
	prop[2][t].imag += 2.0*prop_tmp[t].imag;
    }

    /* DELTA0 */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][1] = 1;
    chi_b[3][3] = 1;

    baryon_cont1(src1, src2, src3, chi_b, eps, prop[3]);

    for(t=0;t<nt;t++)prop_tmp[t] = cmplx(0.0,0.0);
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);
    for(t=0;t<nt;t++){
	prop[3][t].real += 2.0*prop_tmp[t].real;
	prop[3][t].imag += 2.0*prop_tmp[t].imag;
    }

    free(prop_tmp);

}  /* w_baryon */

