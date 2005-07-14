/******* cppacs_help5.c -vector routine for CPPACS */
/* scalar_mult_add_latvec in main call in congrad */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* unroll to color loops */
/* MIMD version 7 */
/* Kogut-Susskind fermions */

#include "ks_dyn_includes.h"

void cppacs_help5( int parity, Real b ){
    register int is,c,stop;
    register site *s;
    switch(parity){
	case EVEN: is=0; stop=even_sites_on_node; s=&(lattice[0]); break;
	case ODD: is=even_sites_on_node; stop=sites_on_node; 
	    s=&(lattice[even_sites_on_node]); break;
	case EVENANDODD: is=0; stop=sites_on_node; s=&(lattice[0]); break;
    }
    /*voption vec*/
    for( ; is<stop; is++,s++ ){
	s->cg_p.c[0].real = s->resid.c[0].real + b * s->cg_p.c[0].real;
	s->cg_p.c[0].imag = s->resid.c[0].imag + b * s->cg_p.c[0].imag;
	s->cg_p.c[1].real = s->resid.c[1].real + b * s->cg_p.c[1].real;
	s->cg_p.c[1].imag = s->resid.c[1].imag + b * s->cg_p.c[1].imag;
	s->cg_p.c[2].real = s->resid.c[2].real + b * s->cg_p.c[2].real;
	s->cg_p.c[2].imag = s->resid.c[2].imag + b * s->cg_p.c[2].imag;
    }
}
