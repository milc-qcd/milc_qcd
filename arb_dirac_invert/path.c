/************************** path.c *******************************/
/* MIMD version 3 */
/* An arbitrary path walker */
/*  uses tempmat1.e[0][0] to accumulate, and ordinary gathers*/

#include "arb_dirac_inv_includes.h"

void path(int *dir,int *sign,int length)
{
register int i;
register site *s;
msg_tag *mtag0, *mtag1;
int j;


/* j=0 */
	if(sign[0]>0)  {
	    mtag0 = start_gather_site( F_OFFSET(link[dir[0]]), sizeof(su3_matrix),
		OPP_DIR(dir[0]), EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);

	      FORALLSITES(i,s){
	      su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(s->tempmat1) );
	      }

	    cleanup_gather(mtag0);
	}

	if(sign[0]<0) { 
	      FORALLSITES(i,s){
	      su3_adjoint(&(s->link[dir[0]]),&(s->tempmat1) );
	      }
	}


	for(j=1;j<length;j++) {
	if(sign[j] > 0) {

	      FORALLSITES(i,s){
		mult_su3_nn( &(s->tempmat1),&(s->link[dir[j]]),
		    &(s->tempmat2) );
	      }
	    mtag0 = start_gather_site( F_OFFSET(tempmat2), sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);

	      FORALLSITES(i,s){
	      su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(s->tempmat1) ); 
	      }
	    cleanup_gather(mtag0);
	}

	if(sign[j] < 0) {
	    mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
		dir[j], EVENANDODD, gen_pt[1] );
	    wait_gather(mtag0);

	    FORALLSITES(i,s){
		mult_su3_na((su3_matrix *)(gen_pt[1][i]),
			&(s->link[dir[j]]), &(s->tempmat2) );
	    }

	    FORALLSITES(i,s){
	      su3mat_copy(&(s->tempmat2),&(s->tempmat1) );
	    }
	    cleanup_gather(mtag0);
	}

      }



} /* path */

