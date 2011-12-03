/******************* mult_ldu2.c ***************************************/

/* MIMD version 7 */
/* version of 12/29/94 by UMH */

/* Multiply a Wilson vector (spin,color) with a block-diagonal hermition
   matrix stored as a complex lower triangular matrix (without diagonal)
   and a real diagonal. The blocks are spin indices 0,1 and 2,3. */

#include "arb_ov_includes.h"



/* To simplify the code with the above block structure, introduce
   the somewhat dirty structure, equivalent to a wilson_vector: */
typedef struct { complex b[2][6]; } wilson_block_vector;

void mult_ldu1_field(wilson_vector *src,wilson_vector *dest,
triangular *triang, diagonal *diag, int parity)
 {
register int i;
register site *s;
register int b,j,k,jk,kj;

wilson_block_vector *src1,*dest1;
wilson_block_vector s1, d1;
triangular tt;

src1=(wilson_block_vector*)src;
dest1=(wilson_block_vector*)dest;

    FORALLSITES(i,s){

	s1=src1[i];
	tt=triang[i];

        for(b=0;b<2;b++)for(j=0;j<6;j++){

             (d1.b[b][j]).real = (diag[i].di[b][j]) * (s1.b[b][j]).real; 
	     (d1.b[b][j]).imag = (diag[i].di[b][j]) * (s1.b[b][j]).imag; 




            jk=j*(j-1)/2;
            for(k=0;k<j;k++){
                 (d1.b[b][j]).real += (tt.tr[b][jk]).real*(s1.b[b][k]).real - 
		        (tt.tr[b][jk]).imag*(src1[i].b[b][k]).imag; 
		 (d1.b[b][j]).imag += (tt.tr[b][jk]).real*(s1.b[b][k]).imag + 
			(tt.tr[b][jk]).imag*(src1[i].b[b][k]).real;

                jk++;
            }
        }

        for(b=0;b<2;b++)for(k=0;k<6;k++){


            kj=k*(k-1)/2;
            for(j=0;j<k;j++){
                 (d1.b[b][j]).real += (tt.tr[b][kj]).real*(s1.b[b][k]).real + 
		                           (tt.tr[b][kj]).imag*(s1.b[b][k]).imag; 
		 (d1.b[b][j]).imag += (tt.tr[b][kj]).real*(s1.b[b][k]).imag - 
		                           (tt.tr[b][kj]).imag*(s1.b[b][k]).real; 

                kj++;
            }
        }
	
	(dest1[i])=d1;

    }

} /* mult_ldu */

void mult_ldu1(field_offset src,field_offset dest,
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

} /* mult_ldu */

