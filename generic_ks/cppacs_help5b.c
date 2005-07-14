/******* cppacs_help5b.c -vector routine for CPPACS */
/* MIMD version 7 */
/* scalar_mult_add_latvec in main call in congrad */
/* unroll to color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

#include "ks_dyn_includes.h"

void cppacs_help5( int parity, Real b ){
    register int is,c;
    switch(parity){
	case EVEN:
		/*voption vec*/
	    for( is=0; is<even_sites_on_node; is++ ){
	    lattice[is].cg_p.c[0].real = lattice[is].resid.c[0].real
		+ b *  lattice[is].cg_p.c[0].real;
	    lattice[is].cg_p.c[0].imag = lattice[is].resid.c[0].imag
		+ b *  lattice[is].cg_p.c[0].imag;
	    lattice[is].cg_p.c[1].real = lattice[is].resid.c[1].real
		+ b *  lattice[is].cg_p.c[1].real;
	    lattice[is].cg_p.c[1].imag = lattice[is].resid.c[1].imag
		+ b *  lattice[is].cg_p.c[1].imag;
	    lattice[is].cg_p.c[2].real = lattice[is].resid.c[2].real
		+ b *  lattice[is].cg_p.c[2].real;
	    lattice[is].cg_p.c[2].imag = lattice[is].resid.c[2].imag
		+ b *  lattice[is].cg_p.c[2].imag;
	} break;
	case ODD:
		/*voption vec*/
 	    for( is=even_sites_on_node; is<sites_on_node; is++ ){
	    lattice[is].cg_p.c[0].real = lattice[is].resid.c[0].real
		+ b *  lattice[is].cg_p.c[0].real;
	    lattice[is].cg_p.c[0].imag = lattice[is].resid.c[0].imag
		+ b *  lattice[is].cg_p.c[0].imag;
	    lattice[is].cg_p.c[1].real = lattice[is].resid.c[1].real
		+ b *  lattice[is].cg_p.c[1].real;
	    lattice[is].cg_p.c[1].imag = lattice[is].resid.c[1].imag
		+ b *  lattice[is].cg_p.c[1].imag;
	    lattice[is].cg_p.c[2].real = lattice[is].resid.c[2].real
		+ b *  lattice[is].cg_p.c[2].real;
	    lattice[is].cg_p.c[2].imag = lattice[is].resid.c[2].imag
		+ b *  lattice[is].cg_p.c[2].imag;
	} break;
	case EVENANDODD:
		/*voption vec*/
	    for( is=0; is<sites_on_node; is++ ){
	    lattice[is].cg_p.c[0].real = lattice[is].resid.c[0].real
		+ b *  lattice[is].cg_p.c[0].real;
	    lattice[is].cg_p.c[0].imag = lattice[is].resid.c[0].imag
		+ b *  lattice[is].cg_p.c[0].imag;
	    lattice[is].cg_p.c[1].real = lattice[is].resid.c[1].real
		+ b *  lattice[is].cg_p.c[1].real;
	    lattice[is].cg_p.c[1].imag = lattice[is].resid.c[1].imag
		+ b *  lattice[is].cg_p.c[1].imag;
	    lattice[is].cg_p.c[2].real = lattice[is].resid.c[2].real
		+ b *  lattice[is].cg_p.c[2].real;
	    lattice[is].cg_p.c[2].imag = lattice[is].resid.c[2].imag
		+ b *  lattice[is].cg_p.c[2].imag;
	} break;
    }
}
