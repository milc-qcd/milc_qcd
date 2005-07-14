/******* cppacs_help3.c -vector routine for CPPACS */
/* loop number 3 in main call in congrad */
/* unroll to color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

/*
FORSOMEPARITY{
}
            scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4,
                &(s->ttt) );
            pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
*/

#include "ks_dyn_includes.h"

void cppacs_help3( int parity, Real scalar, double *pkp ){
    register int is,start,stop;
    register site *s,*s0;
    register Real pkptmp;
    switch(parity){
	case EVEN: start=0; stop=even_sites_on_node; s0=&(lattice[0]); break;
	case ODD: start=even_sites_on_node; stop=sites_on_node; 
	    s0=&(lattice[even_sites_on_node]); break;
	case EVENANDODD: start=0; stop=sites_on_node; s0=&(lattice[0]); break;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
      /* scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4, &(s->ttt) );*/
	s->ttt.c[0].real = s->ttt.c[0].real + scalar * s->cg_p.c[0].real;
	s->ttt.c[0].imag = s->ttt.c[0].imag + scalar * s->cg_p.c[0].imag;
	s->ttt.c[1].real = s->ttt.c[1].real + scalar * s->cg_p.c[1].real;
	s->ttt.c[1].imag = s->ttt.c[1].imag + scalar * s->cg_p.c[1].imag;
	s->ttt.c[2].real = s->ttt.c[2].real + scalar * s->cg_p.c[2].real;
	s->ttt.c[2].imag = s->ttt.c[2].imag + scalar * s->cg_p.c[2].imag;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
      /* pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );*/
        pkptmp=0.0;
	pkptmp += s->cg_p.c[0].real * s->ttt.c[0].real;
	pkptmp += s->cg_p.c[0].imag * s->ttt.c[0].imag;
	pkptmp += s->cg_p.c[1].real * s->ttt.c[1].real;
	pkptmp += s->cg_p.c[1].imag * s->ttt.c[1].imag;
	pkptmp += s->cg_p.c[2].real * s->ttt.c[2].real;
	pkptmp += s->cg_p.c[2].imag * s->ttt.c[2].imag;
    *pkp += (double)pkptmp;
    }
}
