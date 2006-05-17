/******************* mult_ldu1.c ***************************************/
/* MIMD version 6 */

/* MIMD version 4 */
/* version of 12/29/94 by UMH */

/* Multiply a Wilson vector (spin,color) with a block-diagonal hermition
   matrix stored as a complex lower triangular matrix (without diagonal)
   and a real diagonal. The blocks are spin indices 0,1 and 2,3. */

#include "arb_dirac_eig_includes.h"
/* To simplify the code with the above block structure, introduce
   the somewhat dirty structure, equivalent to a wilson_vector: */
typedef struct { complex b[2][6]; } wilson_block_vector;

void mult_ldu1_site(field_offset src,field_offset dest,
field_offset triang,field_offset diag, int parity)
 {
register int i;
register site *s;
register int b,j,k,jk,kj;
complex ctmp;



    FORSOMEPARITY(i,s,parity){
	for(b=0;b<2;b++)for(j=0;j<6;j++){
	    /* diagonal part */
	    CMULREAL(((wilson_block_vector *)F_PT(s,src))->b[b][j],
		((diagonal *)F_PT(s,diag))->di[b][j],
		((wilson_block_vector *)F_PT(s,dest))->b[b][j]);

	    /* lower triangular part */
	    jk=j*(j-1)/2;
	    for(k=0;k<j;k++){
		CMUL(((triangular *)F_PT(s,triang))->tr[b][jk],
		    ((wilson_block_vector *)F_PT(s,src))->b[b][k],ctmp);
		CADD(ctmp,((wilson_block_vector *)F_PT(s,dest))->b[b][j],
		    ((wilson_block_vector *)F_PT(s,dest))->b[b][j]);
		jk++;
	    }
	}

	for(b=0;b<2;b++)for(k=0;k<6;k++){

	    /* upper triangular part */
	    kj=k*(k-1)/2;
	    for(j=0;j<k;j++){
		CMULJ_(((triangular *)F_PT(s,triang))->tr[b][kj],
		    ((wilson_block_vector *)F_PT(s,src))->b[b][k],ctmp);
		CADD(ctmp,((wilson_block_vector *)F_PT(s,dest))->b[b][j],
		    ((wilson_block_vector *)F_PT(s,dest))->b[b][j]);
		kj++;
	    }
	}
    }

} /* mult_ldu1 */

