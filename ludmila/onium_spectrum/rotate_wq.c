#include "generic_wilson_includes.h"

//the rotated quark ends up in the quark_propagator_copy!

void rotate_w_quark(field_offset src, field_offset dest, float d1)
{
  int spin,color,spin1, color1,i;
  site *s;
  float fconst;
  wilson_propagator *qsrc;
  wilson_propagator *qdest;
  wilson_vector qtemp;
  wilson_propagator qtemp1;

  fconst = d1 * 0.25;
   FORALLSITES(i,s)
    {
      qsrc = (wilson_propagator *)F_PT(s,src);
      for(color=0;color<3;color++){
      mult_swv_by_gamma_l( &(qsrc->c[color]), &(qtemp1.c[color]), GAMMAFIVE);
      mult_swv_by_gamma_r( &(qtemp1.c[color]), &(qsrc->c[color]), GAMMAFIVE);
      }

    }  

  for(spin=0;spin<4;spin++)for(color=0;color<3;color++){
    FORALLSITES(i,s)
      {
	qsrc = (wilson_propagator *)F_PT(s,src);
	s->psi = qsrc->c[color].d[spin];
      }
    
    /* Do Wilson Dslash3 on the psi field */
    dslash3(F_OFFSET(psi), F_OFFSET(mp), PLUS,  EVENANDODD); 
    //Do Wilson Dslash3 on the psi field
    dslash3(F_OFFSET(psi), F_OFFSET(tmp), MINUS, EVENANDODD);
 
    /* From subtraction we get 2*Dslash3 */
    
    FORALLSITES(i,s){
      qsrc = (wilson_propagator *)F_PT(s,src);
      qdest = (wilson_propagator *)F_PT(s,dest);
      sub_wilson_vector(&(s->mp), &(s->tmp), &qtemp);
        
      scalar_mult_wvec(&qtemp, fconst, &qtemp);
      
      add_wilson_vector(&qtemp, &(qsrc->c[color].d[spin]),
			&(qdest->c[color].d[spin]));
		
    }
  }


}


