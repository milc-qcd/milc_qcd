/********** w_baryon_p.c *************/
/* MIMD version 6 */
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

#include "arb_dirac_inv_includes.h"

#define Nc 3
#define Ns 4

void baryon_cont1(field_offset src1, field_offset src2, field_offset src3, 
		  int chi_b[4][4], int eps[3][3][3], Real prop[MAX_P][MAX_NT]);
void baryon_cont2(field_offset src1, field_offset src2, field_offset src3, 
		  int chi_b[4][4], int eps[3][3][3], Real prop[MAX_P][MAX_NT]);

void w_baryon_p(field_offset src1,field_offset src2,field_offset src3,
	      Real prop[4][MAX_P][MAX_NT]) 
/* srci's are type wilson_propagator */
{

int i,j,k,t;

int eps[3][3][3], chi_b[4][4];
Real prop_tmp[MAX_P][MAX_NT];


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

    for(i=0;i<MAX_P;i++)for(t=0;t<nt;t++)prop_tmp[i][t] = 0.0;
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);
    for(i=0;i<MAX_P;i++)for(t=0;t<nt;t++){
	prop[2][i][t] += 2.0*prop_tmp[i][t];
    }

    /* DELTA0 */
    for( i=0; i< 4; i++ )
    for( j=0; j< 4; j++ ){
	chi_b[i][j] = 0;
    }
    chi_b[1][1] = 1;
    chi_b[3][3] = 1;

    baryon_cont1(src1, src2, src3, chi_b, eps, prop[3]);

    for(i=0;i<MAX_P;i++)for(t=0;t<nt;t++)prop_tmp[i][t] = 0.0;
    baryon_cont2(src1, src2, src3, chi_b, eps, prop_tmp);
    for(i=0;i<MAX_P;i++)for(t=0;t<nt;t++){
	prop[3][i][t] += 2.0*prop_tmp[i][t];
    }

}  /* w_baryon */

/******** baryon_cont1.c *************/
/* MIMD version 6 */
/* UMH April 96 */

/* Construct baryon propagator contraction 1, where the first two quark
   propagators form a loop with contracted spin indices.
   In all baryons the colour components are contracted with the totally
   antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   The 'wave functions' , i.e. the C*Gamma, are encoded in
   chi_b(i,j), where i and j label the spin indices of the
   first and second quark propagator.
*/

void baryon_cont1(field_offset src1, field_offset src2, field_offset src3, 
		  int chi_b[4][4], int eps[3][3][3], Real prop[MAX_P][MAX_NT])
/* src1-3 are type wilson_propagator */
{

register int i;
register site *s;

int my_t;

int ci_1, ci_2, ci_3, si_1, si_2, si_3;
int cf_1, cf_2, cf_3, sf_1, sf_2, sf_3;
int  chi_i, chi_f, eps_f, eps_i,j;
Real factor;
complex diquark, diquark_temp;

/* for nonzero momentum */
Real cx,cy,cz,cxy,cyz,cxz,c111;

    FORALLSITES(i,s){

	my_t = s->t;

	/* Sum over source and sink colors of quark 3 */
	for(ci_3=0;ci_3<Nc;ci_3++)for(cf_3=0;cf_3<Nc;cf_3++){

	  diquark = cmplx(0.0,0.0);

	  /* Sum over source spins of quarks 1 and 2 */
	  /* They will form the "di_quark" */
	  for(si_1=0;si_1<Ns;si_1++)for(si_2=0;si_2<Ns;si_2++){

	    chi_i = chi_b[si_1][si_2];
	    if( chi_i != 0 ){

	      /* Sum over sink spins of quarks 1 and 2 */
	      for(sf_1=0;sf_1<Ns;sf_1++)for(sf_2=0;sf_2<Ns;sf_2++){

		chi_f = chi_b[sf_1][sf_2];
		if( chi_f != 0 ){

		  /* Sum over source colors of quarks 1 and 2 */
		  for(ci_1=0;ci_1<Nc;ci_1++)
		  if(ci_1 != ci_3) for(ci_2=0;ci_2<Nc;ci_2++){

		    eps_i = eps[ci_1][ci_2][ci_3];
		    if( eps_i != 0 ){

		      /* Sum over sink colors of quarks 1 and 2 */
		      for(cf_1=0;cf_1<Nc;cf_1++)
		      if(cf_1 != cf_3) for(cf_2=0;cf_2<Nc;cf_2++){

			eps_f = eps[cf_1][cf_2][cf_3];
			if( eps_f != 0 ){

			  factor = (Real)(eps_f*eps_i*chi_i*chi_f);
			  CMUL(
((wilson_propagator *)F_PT(s,src1))->c[cf_1].d[sf_1].d[si_1].c[ci_1],
((wilson_propagator *)F_PT(s,src2))->c[cf_2].d[sf_2].d[si_2].c[ci_2],
				diquark_temp);
			  diquark.real += factor*diquark_temp.real;
			  diquark.imag += factor*diquark_temp.imag;

			}  /* eps_f */
		      }  /* sum cf_1, cf_2 */
		    }  /* eps_i */
		  }  /* sum ci_1, ci_2 */
		}  /* chi_f */
	      }  /* sum sf_1, sf_2 */
	    }  /* chi_i */
	  } /* sum si_1, si_2 */

	  /* Sum over source and sink spin of uncontracted quark 3 */
	  /* Actually just use spin 1 */
	  si_3 = sf_3 = 1;

	    CMUL(diquark,
((wilson_propagator *)F_PT(s,src3))->c[cf_3].d[sf_3].d[si_3].c[ci_3],
		diquark_temp);


                for(j=0;j<3;j++){
                cz=cos(2.0*PI/(Real)nz*(Real)(s->z)*(Real)j);
                cx=cos(2.0*PI/(Real)nx*(Real)(s->x)*(Real)j);
                cy=cos(2.0*PI/(Real)ny*(Real)(s->y)*(Real)j);

 
                prop[j][my_t] += diquark_temp.real*(cx+cy+cz)/3.0;
                }
                cxy=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y));
                cxz=cos(2.0*PI/(Real)nz*(Real)(s->x +s->z));
                cyz=cos(2.0*PI/(Real)nz*(Real)(s->y +s->z));
                c111=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y + s->z));
                prop[3][my_t] += diquark_temp.real*(cxy+cyz+cxz)/3.0;
                prop[4][my_t] += diquark_temp.real*c111;



	  /* }  */ /* sum sf_3, si_3 */
	}  /* sum cf_3, ci_3 */

    }  /* FORALLSITES */

}  /* baryon_cont1 */

/******** baryon_cont2.c *************/
/* MIMD version 6 */
/* UMH April 96 */

/* Construct baryon propagator contraction 2, with no contracted loop,
   where at the source quark 2 and 3 are interchanged compared to the
   contraction where the first two propagators form a loop with
   contracted Dirac indices.
   In all baryons the colour components are contracted with the totally
   antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   The 'wave functions' , i.e. the C*Gamma, are encoded in
   chi_b(i,j), where i and j label the spin indices of the
   first and second quark propagator.
*/

void baryon_cont2(field_offset src1, field_offset src2, field_offset src3, 
		  int chi_b[4][4], int eps[3][3][3], Real prop[MAX_P][MAX_NT])
/* src1-3 are type wilson_propagator */
{

register int i;
register site *s;

int my_t;

int ci_1, ci_2, ci_3, si_1, si_2, si_3;
int cf_1, cf_2, cf_3, sf_1, sf_2, sf_3;
int  chi_i, chi_f, eps_f, eps_i,j;
Real factor;
complex diquark, diquark_temp;

/* for nonzero momentum */
Real cx,cy,cz,cxy,cyz,cxz,c111;

    FORALLSITES(i,s){

	my_t = s->t;

	/* Sum over source and sink colors of quark 3 */
	for(ci_3=0;ci_3<Nc;ci_3++)for(cf_3=0;cf_3<Nc;cf_3++){

	  /* Sum over source spin of quark 3, which is connected there! */
	  for(si_3=0;si_3<Ns;si_3++){

	    diquark = cmplx(0.0,0.0);

	    /* Sum over source spin of connected quark 1 */
	    /* Quark 1 and 2, connected at the sink, will form the "di_quark" */
	    for(si_1=0;si_1<Ns;si_1++){

	      chi_i = chi_b[si_1][si_3];
	      if( chi_i != 0 ){

		/* Sum over sink spins of quarks 1 and 2 */
		for(sf_1=0;sf_1<Ns;sf_1++)for(sf_2=0;sf_2<Ns;sf_2++){

		  chi_f = chi_b[sf_1][sf_2];
		  if( chi_f != 0 ){

		    /* Sum over source colors of quarks 1 and 2 */
		    for(ci_1=0;ci_1<Nc;ci_1++)
		    if(ci_1 != ci_3) for(ci_2=0;ci_2<Nc;ci_2++){

		      eps_i = eps[ci_1][ci_2][ci_3];
		      if( eps_i != 0 ){

			/* Sum over sink colors of quarks 1 and 2 */
			for(cf_1=0;cf_1<Nc;cf_1++)
			if(cf_1 != cf_3) for(cf_2=0;cf_2<Nc;cf_2++){

			  eps_f = eps[cf_1][cf_2][cf_3];
			  if( eps_f != 0 ){

			    /* Sum over source spin of quark 2 */
			    /* Actually just use spin 1 */
			    si_2 = 1;

			      factor = (Real)(eps_f*eps_i*chi_i*chi_f);
			      CMUL(
((wilson_propagator *)F_PT(s,src1))->c[cf_1].d[sf_1].d[si_1].c[ci_1],
((wilson_propagator *)F_PT(s,src2))->c[cf_2].d[sf_2].d[si_2].c[ci_2],
				  diquark_temp);
			      diquark.real += factor*diquark_temp.real;
			      diquark.imag += factor*diquark_temp.imag;

			    /* }  */ /* sum si_2 */
			  }  /* eps_f */
			}  /* sum cf_1, cf_2 */
		      }  /* eps_i */
		    }  /* sum ci_1, ci_2 */
		  }  /* chi_f */
		}  /* sum sf_1, sf_2 */
	      }  /* chi_i */
	    }  /* sum si_1 */

	    /* Sum over sink spin of uncontracted quark 3 */
	    /* Actually just use spin 1 */
	    sf_3 = 1;

	      CMUL(diquark,
((wilson_propagator *)F_PT(s,src3))->c[cf_3].d[sf_3].d[si_3].c[ci_3],
		  diquark_temp);
                for(j=0;j<3;j++){
                cz=cos(2.0*PI/(Real)nz*(Real)(s->z)*(Real)j);
                cx=cos(2.0*PI/(Real)nx*(Real)(s->x)*(Real)j);
                cy=cos(2.0*PI/(Real)ny*(Real)(s->y)*(Real)j);

 
                prop[j][my_t] += diquark_temp.real*(cx+cy+cz)/3.0;
                }
                cxy=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y));
                cxz=cos(2.0*PI/(Real)nz*(Real)(s->x +s->z));
                cyz=cos(2.0*PI/(Real)nz*(Real)(s->y +s->z));
                c111=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y + s->z));
                prop[3][my_t] += diquark_temp.real*(cxy+cyz+cxz)/3.0;
                prop[4][my_t] += diquark_temp.real*c111;

	    /* }  */ /* sum sf_3 */
	  }  /* sum si_3 */
	}  /* sum cf_3, ci_3 */

    }  /* FORALLSITES */

}  /* baryon_cont2 */

#undef Nc
#undef Ns
