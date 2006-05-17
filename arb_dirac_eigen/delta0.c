/******* delta0.c -   FP fermions ****/
/* MIMD version 6 */

/*
dest= delta0*src


isign=1 for delta0, -1 for delta0^\dagger.
*/



#include "arb_dirac_eig_includes.h"
#include <string.h>




void delta0(field_offset src,field_offset dest,  int isign)
{
	register int i;
	register site *s;
	msg_tag *tag[2];
	wilson_vector wvec1,wvec2;
	int ipath,n[4],k,mu;
	Real inr[4];




	/* begin with the origin */
	FORALLSITES(i,s) {
		scalar_mult_wvec((wilson_vector *)F_PT(s,src), lambda[0],
		    (wilson_vector *)F_PT(s,dest));
	}


/* multiply the clover term  */

	if(clover_term != 0.0){
        mult_ldu1_site(src,F_OFFSET(htmp[1]),
                F_OFFSET(clov), F_OFFSET(clov_diag), EVENANDODD);
	 FORALLSITES(i,s) {
			scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,dest)),
			    &(s->htmp[1]), 1.0, 
			    ((wilson_vector *)F_PT(s,dest)) );
	 }
	}


	/* now loop over the offsets and collect them */
	for(ipath=0;ipath<off_max;ipath++){

		k=label[ipath];
/*  printf("ipath=%d k=%d\n",ipath,k); */

		if(rho[k] != 0.0)if(lambda[k]!= 0.0) {
/*
printf("ipath=%d k=%d  %e   %e \n",ipath,k,rho[k],lambda[k]); 
*/


			for(mu=0;mu<4;mu++) n[mu] = -offset[ipath][mu];
			for(mu=0;mu<4;mu++) inr[mu]=isign*n[mu]*rho[k];


        /* Take  src displaced in up direction, gather it to "our site" */
 /*printf("ahead %d %d %d %d\n",offset[ipath][0], offset[ipath][1],
offset[ipath][2],offset[ipath][3]);  */

			tag[0] = start_general_gather_site( src,
			    sizeof(wilson_vector), offset[ipath],
				 EVENANDODD, gen_pt[0] );


        /* Take  src displaced in up direction, gathered,
                multiply it by link matrix,  and add to dest */
			wait_general_gather(tag[0]);
/* printf("mult 1\n"); */

			FORALLSITES(i,s) {
			mult_mat_wilson_vec( &(s->blocked_link[ipath]), 
				(wilson_vector * )(gen_pt[0][i]), &wvec1 ); 
			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec1,lambda[k], 
			    (wilson_vector *)F_PT(s,dest));

		/* multiply by the four Dirac matrices ahead*/
			for(mu=0;mu<4;mu++)if(n[mu]!=0){
			mult_by_gamma(&wvec1,&wvec2,mu);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec2, -inr[mu], 
			    (wilson_vector *)F_PT(s,dest));
			} /* dirac */

				}
		cleanup_general_gather(tag[0]);

        /* Take  src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */


			FORALLSITES(i,s){
			  mult_adj_mat_wilson_vec(&(s->blocked_link[ipath]),
			  (wilson_vector *)F_PT(s,src), &(s->htmp[1]));
			}

/* printf("behind %d %d %d %d\n",n[0],n[1],n[2],n[3]); */
			tag[1] = start_general_gather_site( F_OFFSET(htmp[1]),
			    sizeof(wilson_vector), n, EVENANDODD, gen_pt[1] );
			wait_general_gather(tag[1]);
/* printf("mult 2\n"); */

	/* multiply by the scalar term behind*/

			FORALLSITES(i,s) {
			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    (wilson_vector *)gen_pt[1][i],lambda[k], 
			    (wilson_vector *)F_PT(s,dest));

		/* multiply by the four Dirac matrices behind*/
			for(mu=0;mu<4;mu++)if(n[mu]!=0){
			mult_by_gamma((wilson_vector *)gen_pt[1][i],&wvec1,mu);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec1,inr[mu], 
			    (wilson_vector *)F_PT(s,dest));
			}/* dirac */
			} 

		cleanup_general_gather(tag[1]);
		} /* if (rho and lambda) */


	} /* offsets */


}
/* delta0.c */
