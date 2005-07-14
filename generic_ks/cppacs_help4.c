/******* cppacs_help4.c -vector routine for CPPACS */
/* loop number 4 in main call in congrad */
/* unroll to color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

/*
FORSOMEPARITY{
    scalar_mult_add_su3_vector( xxx_pt, cg_p_pt, a, xxx_pt );
    scalar_mult_add_su3_vector( resid_pt, ttt_pt, a, resid_pt );
    rsq += (double)magsq_su3vec( resid_pt );
    xxx_pt = (su3_vector *)( (char *)xxx_pt + stride );
    cg_p_pt = (su3_vector *)( (char *)cg_p_pt + stride );
    ttt_pt = (su3_vector *)( (char *)ttt_pt + stride );
    resid_pt = (su3_vector *)( (char *)resid_pt + stride );
}
*/

#include "ks_dyn_includes.h"

void cppacs_help4( int parity, Real a, double *rsq ){
    register int is,start,stop;
    register site *s,*s0;
    register Real rsqtmp;
    switch(parity){
	case EVEN: start=0; stop=even_sites_on_node; s0=&(lattice[0]); break;
	case ODD: start=even_sites_on_node; stop=sites_on_node; 
	    s0=&(lattice[even_sites_on_node]); break;
	case EVENANDODD: start=0; stop=sites_on_node; s0=&(lattice[0]); break;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
      /* scalar_mult_add_su3_vector( &(s->xxx), &(s->cg_p), a, &(s->xxx) );*/
	s->xxx.c[0].real = s->xxx.c[0].real + a * s->cg_p.c[0].real;
	s->xxx.c[0].imag = s->xxx.c[0].imag + a * s->cg_p.c[0].imag;
	s->xxx.c[1].real = s->xxx.c[1].real + a * s->cg_p.c[1].real;
	s->xxx.c[1].imag = s->xxx.c[1].imag + a * s->cg_p.c[1].imag;
	s->xxx.c[2].real = s->xxx.c[2].real + a * s->cg_p.c[2].real;
	s->xxx.c[2].imag = s->xxx.c[2].imag + a * s->cg_p.c[2].imag;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
      /* scalar_mult_add_su3_vector( &(s->resid), &(s->ttt), a, &(s->resid));*/
	s->resid.c[0].real = s->resid.c[0].real + a * s->ttt.c[0].real;
	s->resid.c[0].imag = s->resid.c[0].imag + a * s->ttt.c[0].imag;
	s->resid.c[1].real = s->resid.c[1].real + a * s->ttt.c[1].real;
	s->resid.c[1].imag = s->resid.c[1].imag + a * s->ttt.c[1].imag;
	s->resid.c[2].real = s->resid.c[2].real + a * s->ttt.c[2].real;
	s->resid.c[2].imag = s->resid.c[2].imag + a * s->ttt.c[2].imag;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
      /* rsq += (double)magsq_su3vec( &(s->resid) );*/
        rsqtmp=0.0;
	rsqtmp += s->resid.c[0].real * s->resid.c[0].real;
	rsqtmp += s->resid.c[0].imag * s->resid.c[0].imag;
	rsqtmp += s->resid.c[1].real * s->resid.c[1].real;
	rsqtmp += s->resid.c[1].imag * s->resid.c[1].imag;
	rsqtmp += s->resid.c[2].real * s->resid.c[2].real;
	rsqtmp += s->resid.c[2].imag * s->resid.c[2].imag;
    *rsq += (double)rsqtmp;
    }
}
