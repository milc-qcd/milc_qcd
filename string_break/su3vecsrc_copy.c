/*****************  su3vecsrc_copy.c ***************************/

/* Contains su3vecsrc_copy, su3_vec_to_src, su3_src_to_vec and
   su3vecsrc_outer_prod.
*/

#include "../include/complex.h"
#include "../include/su3.h"

#include <lattice.h>

/* Copy a su3 vector:  b <- a   */
void su3vecsrc_copy( su3_vector_src *a, su3_vector_src *b, int num_src ){
register int i,ks;
    for(i=0;i<3;i++)for(ks=0;ks<num_src;ks++){
	b->s[ks].c[i].real = a->s[ks].c[i].real;
	b->s[ks].c[i].imag = a->s[ks].c[i].imag;
    }
}


/* Copy a su3 vector:  b <- a   */
void su3_vec_to_src( su3_vector *a, su3_vector_src *b, int num_src ){
register int i,ks;
    for(i=0;i<3;i++)for(ks=0;ks<num_src;ks++){
	b->s[ks].c[i].real = a[ks].c[i].real;
	b->s[ks].c[i].imag = a[ks].c[i].imag;
    }
}


/* Copy a su3 vector:  b <- a   */
void su3_src_to_vec( su3_vector_src *a, su3_vector *b, int j ){
register int i;
    for(i=0;i<3;i++){
	b->c[i].real = a->s[j].c[i].real;
	b->c[i].imag = a->s[j].c[i].imag;
    }
}


void su3vecsrc_outer_prod( su3_vector *a, su3_vector_src *b,
			   su3_matrix *c, int num_src){
register int i, j, ks;
register complex x, y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	x.real = x.imag = 0.0;
	for(ks=0;ks<num_src;ks++){
	    CMUL_J( a[ks].c[i], b->s[ks].c[j], y);
	    CSUM(x, y);
	}
	c->e[i][j] = x;
    }
}

/* su3vecsrc_copy */
