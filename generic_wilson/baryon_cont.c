/******** baryon_cont.c *************/
/* MIMD version 7 */
/* UMH April 96 */

/* Construct baryon propagator contraction 1, where the first two quark
   propagators form a loop with contracted spin indices.
   In all baryons the colour components are contracted with the totally
   antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   The 'wave functions' , i.e. the C*Gamma, are encoded in
   chi_b(i,j), where i and j label the spin indices of the
   first and second quark propagator.
*/

#include "generic_wilson_includes.h"

#define Nc 3
#define Ns 4

void baryon_cont1(wilson_prop_field *src1, wilson_prop_field *src2, 
		  wilson_prop_field *src3, 
		  int chi_b[4][4], int eps[3][3][3], complex *prop)
{

  register int i;
  register site *s;
  
  int my_t;
  
  int ci_1, ci_2, ci_3, si_1, si_2, si_3;
  int cf_1, cf_2, cf_3, sf_1, sf_2, sf_3;
  int  chi_i, chi_f, eps_f, eps_i;
  Real factor;
  complex diquark, diquark_temp;

  if(src1->nc != 3 || src2->nc != 3 || src3->nc != 3){
    node0_printf("baryon_cont1: requires 3 colors for each field\n");
    terminate(1);
  }

  
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
			       src1->swv[cf_1][i].d[sf_1].d[si_1].c[ci_1],
			       src2->swv[cf_2][i].d[sf_2].d[si_2].c[ci_2],
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
		 src3->swv[cf_3][i].d[sf_3].d[si_3].c[ci_3],
		 diquark_temp);
	    prop[my_t].real += diquark_temp.real;
	    prop[my_t].imag += diquark_temp.imag;

	  /* }  */ /* sum sf_3, si_3 */
	}  /* sum cf_3, ci_3 */

    }  /* FORALLSITES */

}  /* baryon_cont1 */

/******** baryon_cont2.c *************/
/* MIMD version 7 */
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

void baryon_cont2(wilson_prop_field *src1, wilson_prop_field *src2, 
		  wilson_prop_field *src3, 
		  int chi_b[4][4], int eps[3][3][3], complex *prop)
{

  register int i;
  register site *s;
  
  int my_t;
  
  int ci_1, ci_2, ci_3, si_1, si_2, si_3;
  int cf_1, cf_2, cf_3, sf_1, sf_2, sf_3;
  int  chi_i, chi_f, eps_f, eps_i;
  Real factor;
  complex diquark, diquark_temp;

  if(src1->nc != 3 || src2->nc != 3 || src3->nc != 3){
    node0_printf("baryon_cont2: requires 3 colors for each field\n");
    terminate(1);
  }



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
				   src1->swv[cf_1][i].d[sf_1].d[si_1].c[ci_1],
				   src2->swv[cf_2][i].d[sf_2].d[si_2].c[ci_2],
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
		   src3->swv[cf_3][i].d[sf_3].d[si_3].c[ci_3],
		   diquark_temp);
	      prop[my_t].real += diquark_temp.real;
	      prop[my_t].imag += diquark_temp.imag;
	      
	    /* }  */ /* sum sf_3 */
	  }  /* sum si_3 */
	}  /* sum cf_3, ci_3 */

    }  /* FORALLSITES */

}  /* baryon_cont2 */

#undef Nc
#undef Ns
