/************************** wl_1l_1corr_new.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 12/28/98  PL  */

/* Computes the operator consisting of a U-shaped product time-like links
   and a light propagator on gauge configurations in the axial gauge.
   Note: here the time-like gauge fields of the last time slice is copied
   in all other time slices as well, instead of the unit matrix. 

   This version: first collect all random propagators
                 into single vector (su3_vec_src) before communicating


         ____
        |    |
        ^    V      where  ~~~~  = light prop
        |    |
         ~~~~
*/


#include "string_break_includes.h"

void wl_1l_1corr(int tot_smear, int step) {

register int i,dir,r,t;
int nth,nxh, j;
register site *s;
su3_matrix tmat1,tmat2;
su3_vector tvec;
msg_tag *mtag[4];
complex *wl1_mes, *wl1_mes_t, cc;


    if( nx != ny || nx != nz){
	if(this_node == 0)printf("wl_lprop gives wrong results for nx!=ny!=nz");
	return;
    }

    nth = nt/2;  nxh = nx/2;
    wl1_mes = (complex *)malloc(num_src*nth*nxh*sizeof(complex));
    wl1_mes_t = (complex *)malloc(nth*nxh*sizeof(complex));

    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wl1_mes_t[r+nxh*t] = cmplx(0.0,0.0);
    }

    for(t=0; t<nth; t++) for(r=0;r<nxh;r++)for(j=0;j<num_src;j++){
	wl1_mes[r+nxh*t+nxh*nth*j] = cmplx(0.0,0.0);
    }

    for(dir=XUP;dir<=ZUP;dir++){

	FORALLSITES(i,s){
	    su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]), num_src);
	}

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[dir]), &(s->s_link));
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Start gather of forward time-like links */
	mtag[0] = start_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );

	/* gather light propagators in direction dir    */
	mtag[1] = start_gather_site( F_OFFSET(dtmpvecs[1].n[0]),
	    sizeof(su3_vector_src), dir, EVENANDODD, gen_pt[1] );

	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r>0 ){
		wait_gather( mtag[2]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    mult_su3_nn( &(s->link[dir]), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    wait_gather( mtag[1]);
	    FORALLSITES(i,s){
		su3vecsrc_copy((su3_vector_src *)(gen_pt[1][i]),
		    &(s->resid_src), num_src);
	    }
	    FORALLSITES(i,s){
		su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[0]),
		    num_src);
	    }

	    /* Inbetween gather space-links for next r, if still needed. */
	    if( r==0 ){
		restart_gather_site( F_OFFSET(dtmpvecs[1].n[0]),
		    sizeof(su3_vector_src),
		    dir, EVENANDODD, gen_pt[1], mtag[1] );

		mtag[2] = start_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[2] );
	    }
	    else if( r<(nxh-1) ){
		restart_gather_site( F_OFFSET(dtmpvecs[1].n[0]),
		    sizeof(su3_vector_src),
		    dir, EVENANDODD, gen_pt[1], mtag[1] );

		restart_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[2], mtag[2] );
	    }
	    else{
		cleanup_gather( mtag[1]);
		cleanup_gather( mtag[2]);
	    }

	    /* Collect forward time-like links. */
	    wait_gather( mtag[0]);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }

	    /* Recursively compute the "Wilson loops" of different time extent*/
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), &(s->s_link_f));
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
			TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

               
		/*   average over sources                                 */

		for(j=0; j<num_src; j++){

		    /* first rewrite current prop and source in standard format  */

		    FORALLSITES(i,s){
			su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
			    &(s->tempvec[2]), j);
		    }

		    /* construct operator by multiplying U-shape Wilson loop
		       with light propagator and then by random source  */
		    FORALLSITES(i,s){
			if( ((s->t)+t+1)>=nt ){
			    mult_su3_nn( &(s->link[TUP]), &(s->s_link_f),
				&tmat1);
			    mult_su3_na( &tmat1, &(s->t_link_f), &tmat2);
			    mult_su3_mat_vec(&tmat2, &(s->tempvec[2]), &tvec);
			}
			else{
			    mult_su3_mat_vec(&(s->s_link_f), &(s->tempvec[2]),
				&tvec);
			}

			cc = su3_dot(&(s->g_rand[j]), &tvec);

			wl1_mes[r+nxh*t+nxh*nth*j].real += cc.real;
			wl1_mes[r+nxh*t+nxh*nth*j].imag += cc.imag;

		    }

		} /* end loop over j */

	    } /* end loop over t */

	    /* Start gather of forward time-like links for next r. */
	    if( r<(nxh-1) ){
		restart_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[0], mtag[0] );
	    }
	    else{
		cleanup_gather( mtag[0]);
	    }

	} /* end loop over r */

    } /* end loop over dir */

    /* Normalize and print the Wilson loops */
    for(j=0;j<num_src;j++) for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	g_floatsum( &wl1_mes[r+nxh*t+nxh*nth*j].real );
	wl1_mes[r+nxh*t+nxh*nth*j].real /= (Real)(3*volume);
	g_floatsum( &wl1_mes[r+nxh*t+nxh*nth*j].imag );
	wl1_mes[r+nxh*t+nxh*nth*j].imag /= (Real)(3*volume);

     /*
	if(this_node == 0)
	    printf("WL_1LC1_%d_OP_%d %d %d %d %e %e\n", step, tot_smear, j, r, t,
		(double)wl1_mes[r+nxh*t+nxh*nth*j].real,
		(double)wl1_mes[r+nxh*t+nxh*nth*j].imag);

     */

	wl1_mes_t[r+nxh*t].real += wl1_mes[r+nxh*t+nxh*nth*j].real;
	wl1_mes_t[r+nxh*t].imag += wl1_mes[r+nxh*t+nxh*nth*j].imag;
    }



    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wl1_mes_t[r+nxh*t].real /= (Real)num_src;
	wl1_mes_t[r+nxh*t].imag /= (Real)num_src;
	if(this_node == 0)
	    printf("WL_1LC1_%d_OP_T%d %d %d %e %e\n", step, tot_smear, r, t,
		(double)wl1_mes_t[r+nxh*t].real,
		(double)wl1_mes_t[r+nxh*t].imag);
    }


    free( wl1_mes);
    free( wl1_mes_t);

} /*  wl_1l_1corr.c  */

