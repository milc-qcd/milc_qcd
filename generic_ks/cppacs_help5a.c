/******* cppacs_help5a.c -vector routine for CPPACS */
/* MIMD version 7 */
/* scalar_mult_add_latvec in main call in congrad */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

#include "ks_dyn_includes.h"

void cppacs_help5( int parity, Real b ){
    register int is,c;
    switch(parity){
	case EVEN:  for( is=0; is<even_sites_on_node; is++ ){
	  for(c=0;c<3;c++){
	    lattice[is].cg_p.c[c].real = lattice[is].resid.c[c].real
		+ b *  lattice[is].cg_p.c[c].real;
	    lattice[is].cg_p.c[c].imag = lattice[is].resid.c[c].imag
		+ b *  lattice[is].cg_p.c[c].imag;
	  }
/**scalar_mult_add_su3_vector( &lattice[is].resid,
&lattice[is].cg_p, b,  &lattice[is].cg_p);**/
	} break;
	case ODD:  for( is=even_sites_on_node; is<sites_on_node; is++ ){
	  for(c=0;c<3;c++){
	    lattice[is].cg_p.c[c].real = lattice[is].resid.c[c].real
		+ b *  lattice[is].cg_p.c[c].real;
	    lattice[is].cg_p.c[c].imag = lattice[is].resid.c[c].imag
		+ b *  lattice[is].cg_p.c[c].imag;
	  }
	} break;
	case EVENANDODD:  for( is=0; is<sites_on_node; is++ ){
	  for(c=0;c<3;c++){
	    lattice[is].cg_p.c[c].real = lattice[is].resid.c[c].real
		+ b *  lattice[is].cg_p.c[c].real;
	    lattice[is].cg_p.c[c].imag = lattice[is].resid.c[c].imag
		+ b *  lattice[is].cg_p.c[c].imag;
	  }
	} break;
    }
}
