/******************* make_clov1.c **************************************/

/* MIMD version 4 */
/* version of 1/4/95 by UMH */

/* Prepare the "clover" term */


#include "arb_dirac_eig_includes.h"
#include <string.h>



void make_clov1()
 {
register int i;
register site *s;

register int b,j,k,jk,jk2;
register complex ctmp;

Real Clov_c;

field_offset f_mn;
f_mn = F_OFFSET(tmp);

/* The clover term
*                                         ( X | 0 )
* sum_{mu<nu} i sigma_{mu,nu} F_{mu,nu} = ( ----- )
*                                         ( 0 | Y )
*
* is block diagonal in the Weyl basis used, and hermitian.

This version does not have the 1 term, in contrast to the ``usual''
MILC code!


* we will store it in a lower complex triangular (block) matrix (without
* diagonal) and a real diagonal. The blocks go over Dirac indices 0,1 and 2,3.
*
* Here X = sigma_1 (F_23 - F_14) + sigma_2 (F_13 + F_24) + sigma_3 (F_12 - F_34)
* and  Y = sigma_1 (F_23 + F_14) + sigma_2 (F_13 - F_24) + sigma_3 (F_12 + F_34)
*
* This has to be multiplied with the "clover coefficient and subtracted from 1|.
*
* Note that F_mn = f_mn/i = -i f_mn!
*
* (The indices above start from 1, in the program they start from 0,
* i.e. subtract -1 from the above) */

Clov_c=clover_term;
if(this_node==0)printf("entered clover setup with term=%e\n",Clov_c);

    f_mu_nu_site(f_mn, 0, 1);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    jk=0;
    for(j=0;j<3;j++){
	FORALLSITES(i,s){
	    s->clov_diag.di[0][j] = 
		-((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[0][j+3] = 
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[1][j] = 
		-((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[1][j+3] = 
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	}
	jk2=(j+3)*(j+2)/2+3;
	for(k=0;k<j;k++){
	    FORALLSITES(i,s){
		TIMESPLUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[0][jk]);
		TIMESMINUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[0][jk2]);
		TIMESPLUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[1][jk]);
		TIMESMINUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[1][jk2]);
	    }
	    jk++;
	    jk2++;
	}
    }

    f_mu_nu_site(f_mn, 2, 3);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    jk=0;
    for(j=0;j<3;j++){
	FORALLSITES(i,s){
	    s->clov_diag.di[0][j] +=
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[0][j+3] -=
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[1][j] -=
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	    s->clov_diag.di[1][j+3] +=
		((su3_matrix *)F_PT(s,f_mn))->e[j][j].imag;
	}
	jk2=(j+3)*(j+2)/2+3;
	for(k=0;k<j;k++){
	    FORALLSITES(i,s){
		TIMESMINUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k], ctmp);
		CADD( s->clov.tr[0][jk], ctmp, s->clov.tr[0][jk]);
		CSUB( s->clov.tr[0][jk2], ctmp, s->clov.tr[0][jk2]);
		CSUB( s->clov.tr[1][jk], ctmp, s->clov.tr[1][jk]);
		CADD( s->clov.tr[1][jk2], ctmp, s->clov.tr[1][jk2]);
	    }
	    jk++;
	    jk2++;
	}
    }

    f_mu_nu_site(f_mn, 1, 2);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    for(j=0;j<3;j++){
	jk=(j+3)*(j+2)/2;
	for(k=0;k<3;k++){
	    FORALLSITES(i,s){
		TIMESPLUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[0][jk]);
		TIMESPLUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k],
		    s->clov.tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu_site(f_mn, 0, 3);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    for(j=0;j<3;j++){
	jk=(j+3)*(j+2)/2;
	for(k=0;k<3;k++){
	    FORALLSITES(i,s){
		TIMESMINUSI( ((su3_matrix *)F_PT(s,f_mn))->e[j][k], ctmp);
		CADD( s->clov.tr[0][jk], ctmp, s->clov.tr[0][jk]);
		CSUB( s->clov.tr[1][jk], ctmp, s->clov.tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu_site(f_mn, 0, 2);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    for(j=0;j<3;j++){
	jk=(j+3)*(j+2)/2;
	for(k=0;k<3;k++){
	    FORALLSITES(i,s){
		CSUB( s->clov.tr[0][jk],
		    ((su3_matrix *)F_PT(s,f_mn))->e[j][k], s->clov.tr[0][jk]);
		CSUB( s->clov.tr[1][jk],
		    ((su3_matrix *)F_PT(s,f_mn))->e[j][k], s->clov.tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu_site(f_mn, 1, 3);
    FORALLSITES(i,s){
        scalar_mult_su3_matrix( ((su3_matrix *)F_PT(s,f_mn)), Clov_c,
	    ((su3_matrix *)F_PT(s,f_mn)) );
    }
    for(j=0;j<3;j++){
	jk=(j+3)*(j+2)/2;
	for(k=0;k<3;k++){
	    FORALLSITES(i,s){
		CSUB( s->clov.tr[0][jk],
		    ((su3_matrix *)F_PT(s,f_mn))->e[j][k], s->clov.tr[0][jk]);
		CADD( s->clov.tr[1][jk],
		    ((su3_matrix *)F_PT(s,f_mn))->e[j][k], s->clov.tr[1][jk]);
	    }
	    jk++;
	}
    }



} /* make_clov */
