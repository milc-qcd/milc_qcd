/******* cppacs_help1.c -vector routine for CPPACS */
/*  mult_su3_mat_vec_sum_4dir loop in dslash_special */
/* unroll color loops */
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions */

/*
FORSOMEPARITY{
 mult_su3_mat_vec_sum_4dir( s->link,
            (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
            (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
            (su3_vector *)F_PT(s,dest));
}
*/

#include "ks_dyn_includes.h"

void cppacs_help1( int parity, field_offset dest ){
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
    dir=XUP;  /* First one is not summed */
      for(color=0;color<3;color++){
        /*voption vec*/
        for( is=start,s=s0 ; is<stop; is++,s++ ){
	    mat = &(s->link[dir]);
	    vec = (su3_vector *)gen_pt[dir][is];

	    /*mult_su3_mat_vec( mat, vec, dest );*/
	    tdest.real = mat->e[color][0].real * vec->c[0].real
			- mat->e[color][0].imag * vec->c[0].imag
			+ mat->e[color][1].real * vec->c[1].real
			- mat->e[color][1].imag * vec->c[1].imag
			+ mat->e[color][2].real * vec->c[2].real
			- mat->e[color][2].imag * vec->c[2].imag ;
	    tdest.imag = mat->e[color][0].real * vec->c[0].imag
			+ mat->e[color][0].imag * vec->c[0].real
			+ mat->e[color][1].real * vec->c[1].imag
			+ mat->e[color][1].imag * vec->c[1].real
			+ mat->e[color][2].real * vec->c[2].imag
			+ mat->e[color][2].imag * vec->c[2].real ;
	    ((su3_vector *)(F_PT(s,dest)))->c[color].real = tdest.real;
	    ((su3_vector *)(F_PT(s,dest)))->c[color].imag = tdest.imag;
        }
      }
    for( dir=YUP; dir<=TUP; dir++){
      for(color=0;color<3;color++){
        /*voption vec*/
        for( is=start,s=s0 ; is<stop; is++,s++ ){
	    mat = &(s->link[dir]);
	    vec = (su3_vector *)gen_pt[dir][is];

	    /*mult_su3_mat_vec_sum( mat, vec, dest );*/
	    tdest.real = mat->e[color][0].real * vec->c[0].real
			- mat->e[color][0].imag * vec->c[0].imag
			+ mat->e[color][1].real * vec->c[1].real
			- mat->e[color][1].imag * vec->c[1].imag
			+ mat->e[color][2].real * vec->c[2].real
			- mat->e[color][2].imag * vec->c[2].imag ;
	    tdest.imag = mat->e[color][0].real * vec->c[0].imag
			+ mat->e[color][0].imag * vec->c[0].real
			+ mat->e[color][1].real * vec->c[1].imag
			+ mat->e[color][1].imag * vec->c[1].real
			+ mat->e[color][2].real * vec->c[2].imag
			+ mat->e[color][2].imag * vec->c[2].real ;
	    ((su3_vector *)(F_PT(s,dest)))->c[color].real += tdest.real;
	    ((su3_vector *)(F_PT(s,dest)))->c[color].imag += tdest.imag;
        }
      }
    }
}
