/************* hoverlap.c ******************************/
/* MIMD version 7 */

#include "arb_ov_includes.h"




void hoverlap(field_offset src,
	      field_offset dest)
{

register int i;
register site *s;

wilson_vector wtmp;
Real term1,term2;

/*
term1=R0*(1.0-m0);
term2=m0+R0*(1.0-m0);
*/
/* massless overlap */
term1=term2=R0;

     /* compute the step function term */

      step(src,dest); 

      /* and add in the ''gamma-5'' term, multiply sum by R0  */


    FORALLSITES(i,s){ 
              mult_by_gamma((wilson_vector *)F_PT(s,src),&wtmp,
			    GAMMAFIVE);

              scalar_mult_wvec(&(wtmp),term2,&wtmp);

              scalar_mult_wvec((wilson_vector *)F_PT(s,dest),term1,
			      (wilson_vector *)F_PT(s,dest) );

	      add_wilson_vector((wilson_vector *)F_PT(s,dest),
			       &wtmp,(wilson_vector *)F_PT(s,dest));
    }

}
void hoverlap_field(wilson_vector* src,
	            wilson_vector* dest)
{

register int i;
register site *s;

wilson_vector wtmp;
/* massless overlap */

/* compute the step function term */

    step_field(src,dest); 

/* and add in the ''gamma-5'' term, multiply sum by R0  */


    FORALLSITES(i,s){ 
              mult_by_gamma(&src[i],&wtmp,
			    GAMMAFIVE);
	      add_wilson_vector(&dest[i], &wtmp,&dest[i]);
              scalar_mult_wvec(&(dest[i]),R0,&(dest[i]));
    }

}
