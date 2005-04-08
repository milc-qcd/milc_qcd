/******* pauli.c -   FP fermions ****/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */

/*
dest= pauli*src
``pauli'' is a full-sized Pauli term done by lookup tables.
We run over ipath=0 to off_max_p-1, with an offset_p[ipath]
labelling the offset, pauli_link[ipath] the appropriate gauge field,
 sigma[ipath] the Dirac matrix, and k the label of the path
(to set the overall scale of the term).  The oppositely-oriented terms
contribute with a negative sign.
*/



#include "arb_dirac_inv_includes.h"

void pauli(field_offset src,field_offset dest)
{
	register int i;
	register site *s;
	msg_tag *tag[2];
	wilson_vector wvec1,wvec2;
	int ipath,n[4],k,mu;

	void mult_by_sigma(wilson_vector *src, wilson_vector *dest, int dir);


/* zero dest */

        FORALLSITES(i,s) {
            clear_wvec((wilson_vector *)F_PT(s,dest) );
        }

	/*  loop over the offsets and collect them */
	for(ipath=0;ipath<off_max_p;ipath++){

		k=label_p[ipath];

			for(mu=0;mu<4;mu++) n[mu] = -offset_p[ipath][mu];


        /* Take  src displaced in up direction, gather it to "our site" */

/*
printf("ahead %d %d %d %d\n",offset_p[ipath][0], offset_p[ipath][1],
offset_p[ipath][2],offset_p[ipath][3]); 
printf("pauli_term %e  sigma_p %d\n",pauli_term[k],sigma_p[ipath]);
*/

			tag[0] = start_general_gather_site( src,
			    sizeof(wilson_vector), offset_p[ipath],
				 EVENANDODD, gen_pt[0] );

        /* Take  src displaced in up direction, gathered,
                multiply it by link matrix,  and add to dest */

			wait_general_gather(tag[0]);


		/* multiply by the sigma matrices ahead*/
			  FORALLSITES(i,s) {
			mult_mat_wvec( &(s->pauli_link[ipath]), 
				(wilson_vector *)(gen_pt[0][i]), &wvec1 ); 
			mult_by_sigma(&wvec1,&wvec2,sigma_p[ipath]);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec2, pauli_term[k], 
			    (wilson_vector *)F_PT(s,dest));
					}


		cleanup_general_gather(tag[0]);

        /* Take  src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */


			FORALLSITES(i,s){
			  mult_adj_mat_wvec(&(s->pauli_link[ipath]),
			  F_PT(s,src), &(s->htmp[0]));
			}
/*
 printf("behind %d %d %d %d\n",n[0],n[1],n[2],n[3]); 
*/
			tag[1] = start_general_gather_site( F_OFFSET(htmp[0]),
			    sizeof(wilson_vector), n, EVENANDODD, gen_pt[1] );


        /* Take  src displaced in down direction,
        process it, and add to dest */
			wait_general_gather(tag[1]);


		/* multiply by the sigma matrices behind*/
			  FORALLSITES(i,s) {
			mult_by_sigma((wilson_vector *)gen_pt[1][i],&wvec1,
			sigma_p[ipath]);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec1,-pauli_term[k], 
			    (wilson_vector *)F_PT(s,dest));
					}


		cleanup_general_gather(tag[1]);


	} /* offsets */


}
/* pauli.c */

void mult_by_sigma(wilson_vector *src, wilson_vector *dest, int dir)
{
/* for the time being, assume that dir is one of 12,13,14,23,24,34,
mod, and just do two Dirac mults (wasteful) */
register int i; /* color */

  switch(dir){
    case 12:
        for(i=0;i<3;i++){
            TIMESMINUSI(  src->d[0].c[i], dest->d[0].c[i] );
            TIMESPLUSI(  src->d[1].c[i], dest->d[1].c[i] );
            TIMESMINUSI( src->d[2].c[i], dest->d[2].c[i] );
            TIMESPLUSI( src->d[3].c[i], dest->d[3].c[i] );
        }
        break;

    case 13:
        for(i=0;i<3;i++){
            TIMESMINUSONE( src->d[1].c[i], dest->d[0].c[i] );
            TIMESPLUSONE(  src->d[0].c[i], dest->d[1].c[i] );
            TIMESMINUSONE(  src->d[3].c[i], dest->d[2].c[i] );
            TIMESPLUSONE( src->d[2].c[i], dest->d[3].c[i] );
        }
        break;

    case 14:
        for(i=0;i<3;i++){
            TIMESPLUSI(  src->d[1].c[i], dest->d[0].c[i] );
            TIMESPLUSI(  src->d[0].c[i], dest->d[1].c[i] );
            TIMESMINUSI( src->d[3].c[i], dest->d[2].c[i] );
            TIMESMINUSI( src->d[2].c[i], dest->d[3].c[i] );
        }
        break;

     case 23:
        for(i=0;i<3;i++){
            TIMESMINUSI(  src->d[1].c[i], dest->d[0].c[i] );
            TIMESMINUSI(  src->d[0].c[i], dest->d[1].c[i] );
            TIMESMINUSI( src->d[3].c[i], dest->d[2].c[i] );
            TIMESMINUSI( src->d[2].c[i], dest->d[3].c[i] );
        }
        break;

    case 24:
        for(i=0;i<3;i++){
            TIMESMINUSONE( src->d[1].c[i], dest->d[0].c[i] );
            TIMESPLUSONE(  src->d[0].c[i], dest->d[1].c[i] );
            TIMESPLUSONE(  src->d[3].c[i], dest->d[2].c[i] );
            TIMESMINUSONE( src->d[2].c[i], dest->d[3].c[i] );
        }
        break;

    case 34:
        for(i=0;i<3;i++){
            TIMESPLUSI(  src->d[0].c[i], dest->d[0].c[i] );
            TIMESMINUSI(  src->d[1].c[i], dest->d[1].c[i] );
            TIMESMINUSI( src->d[2].c[i], dest->d[2].c[i] );
            TIMESPLUSI( src->d[3].c[i], dest->d[3].c[i] );
        }
        break;

   default:
        printf("BAD CALL TO MULT_BY_SIGMA()\n");
  }
}
