/******* cppacs_help2.c -vector routine for CPPACS */
/* sub4vecs loop in dslash_special */
/* unroll color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

/*
FORSOMEPARITY{
sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
            (su3_vector *)(gen_pt[XDOWN][i]),
            (su3_vector *)(gen_pt[YDOWN][i]),
            (su3_vector *)(gen_pt[ZDOWN][i]),
            (su3_vector *)(gen_pt[TDOWN][i]) );
}
*/

#include "ks_dyn_includes.h"

void cppacs_help2( int parity, field_offset dest ){
    register int is,start,stop;
    register site *s,*s0;
    register su3_vector *a,*bx,*by,*bz,*bt;
    switch(parity){
	case EVEN: start=0; stop=even_sites_on_node; s0=&(lattice[0]); break;
	case ODD: start=even_sites_on_node; stop=sites_on_node; 
	    s0=&(lattice[even_sites_on_node]); break;
	case EVENANDODD: start=0; stop=sites_on_node; s0=&(lattice[0]); break;
    }
    /*voption vec*/
    for( is=start,s=s0 ; is<stop; is++,s++ ){
	a = (su3_vector *)F_PT(s,dest);
	bx =  (su3_vector *)(gen_pt[XDOWN][is]);
	by =  (su3_vector *)(gen_pt[YDOWN][is]);
	bz =  (su3_vector *)(gen_pt[ZDOWN][is]);
	bt =  (su3_vector *)(gen_pt[TDOWN][is]);
	a->c[0].real -= ( bx->c[0].real + by->c[0].real
		+ bz->c[0].real + bt->c[0].real );
	a->c[0].imag -= ( bx->c[0].imag + by->c[0].imag
		+ bz->c[0].imag + bt->c[0].imag );
	a->c[1].real -= ( bx->c[1].real + by->c[1].real
		+ bz->c[1].real + bt->c[1].real );
	a->c[1].imag -= ( bx->c[1].imag + by->c[1].imag
		+ bz->c[1].imag + bt->c[1].imag );
	a->c[2].real -= ( bx->c[2].real + by->c[2].real
		+ bz->c[2].real + bt->c[2].real );
	a->c[2].imag -= ( bx->c[2].imag + by->c[2].imag
		+ bz->c[2].imag + bt->c[2].imag );
    }
}
