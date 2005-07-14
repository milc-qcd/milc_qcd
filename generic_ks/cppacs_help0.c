/******* cppacs_help0.c -vector routine for CPPACS */
/*  mult_adj_su3_mat_vec_4dir loop in dslash_special */
/* unroll color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

/*
FORSOMEPARITY{
 mult_adj_su3_mat_vec_4dir( s->link, (su3_vector *)F_PT(s,src), s->tempvec );
}
*/

#include "ks_dyn_includes.h"

void cppacs_help0( int parity, field_offset src ){
    register int is,start,stop;
    register site *s,*s0;
    register int dir,color;
    register su3_matrix *mat;
    register su3_vector *vec;
    register complex tdest;
    switch(parity){
	case EVEN: start=0; stop=even_sites_on_node; s0=&(lattice[0]); break;
	case ODD: start=even_sites_on_node; stop=sites_on_node; 
	    s0=&(lattice[even_sites_on_node]); break;
	case EVENANDODD: start=0; stop=sites_on_node; s0=&(lattice[0]); break;
    }
    for( dir=XUP; dir<=TUP; dir++)for(color=0;color<3;color++){
        /*voption vec*/
        for( is=start,s=s0 ; is<stop; is++,s++ ){
	    mat = &(s->link[dir]);
	    vec = (su3_vector *)F_PT(s,src );

	    /*mult_adj_su3_mat_vec( mat, vec, dest );*/
	    tdest.real = mat->e[0][color].real * vec->c[0].real
				+ mat->e[0][color].imag * vec->c[0].imag
				+ mat->e[1][color].real * vec->c[1].real
				+ mat->e[1][color].imag * vec->c[1].imag
				+ mat->e[2][color].real * vec->c[2].real
				+ mat->e[2][color].imag * vec->c[2].imag ;
	    tdest.imag = mat->e[0][color].real * vec->c[0].imag
				- mat->e[0][color].imag * vec->c[0].real
				+ mat->e[1][color].real * vec->c[1].imag
				- mat->e[1][color].imag * vec->c[1].real
				+ mat->e[2][color].real * vec->c[2].imag
				- mat->e[2][color].imag * vec->c[2].real ;
	    s->tempvec[dir].c[color].real = tdest.real;
	    s->tempvec[dir].c[color].imag = tdest.imag;
        }
    }
}
