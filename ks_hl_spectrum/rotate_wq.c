#include "ks_hl_spectrum_includes.h"

void rotate_w_quark(field_offset src, field_offset dest, Real d1)
{
  int spin,color,spin1, color1,i;
  site *s;
  Real fconst;
  wilson_propagator *qsrc;
  wilson_propagator *qdest;
  wilson_vector qtemp;
  wilson_propagator qtemp1;

  fconst = d1 * 0.25;

  /* gamma_5 multiplication, needed for the milc-fermi formats compatibility of the results*/

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
    
    /* Do Wilson Dslash_w_3D on the psi field */
    dslash_w_3D_site(F_OFFSET(psi), F_OFFSET(mp), PLUS,  EVENANDODD);
 
    dslash_w_3D_site(F_OFFSET(psi), F_OFFSET(tmp), MINUS, EVENANDODD);
   
    /* From subtraction we get 2*Dslash_W_3D */
    
    FORALLSITES(i,s){
      qsrc = (wilson_propagator *)F_PT(s,src);
      qdest = (wilson_propagator *)F_PT(s,dest);
 
      sub_wilson_vector(&(s->mp), &(s->tmp), &qtemp);
      scalar_mult_wvec(&qtemp, fconst, &qtemp);      
      add_wilson_vector(&qtemp, &(qsrc->c[color].d[spin]),
			&(qdest->c[color].d[spin]));
		
    }
  }
  /* Transposing the propagator */
  FORALLSITES(i,s)
    {
      qdest = (wilson_propagator *)F_PT(s,dest);
      qtemp1 = *(wilson_propagator *)F_PT(s,dest);
	for(spin=0;spin<4;spin++)for(color=0;color<3;color++)for(spin1=0;spin1<4;spin1++)for(color1=0;color1<3;color1++)
	  {
	    qdest->c[color].d[spin1].d[spin].c[color1].real = qtemp1.c[color].d[spin].d[spin1].c[color1].real;
	    qdest->c[color].d[spin1].d[spin].c[color1].imag = qtemp1.c[color].d[spin].d[spin1].c[color1].imag;
	  }
    }
  cleanup_dslash_w_3D_temps();
}


