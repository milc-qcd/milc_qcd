/************************** wl_2l_1corr_new.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 08/14/99  PL  */

/* Computes the combination of the product of 2 light quarks
   with 2 legs of time-like links to form a kind of "Wilson loop".
   Note: here the time-like gauge fields of the last time slice is copied
   in all other time slices as well, instead of the unit matrix. */
/*  
         ~~~~   
        |    |
        ^    V      where  ~~~~  = light prop 
        |    | 
         ~~~~ 
*/


#include "string_break_includes.h"

void wl_2l_1corr(int tot_smear, int step) {

register int i,dir,r,t;
int nth,nxh,j;
register site *s;
su3_vector tvec1, tvec2;
su3_matrix tmat1, tmat2;
dble_su3_vec_src *dpt;
msg_tag *mtag[5];
complex  cc, cc1, cc2;
complex *wl2_mes, *wl2_mes_d;
Real xnsrc;

    if( nx != ny || nx != nz){
	if(this_node == 0)printf("wl_2lprop gives wrong results for nx!=ny!=nz");
	return;
    }

    xnsrc=num_src*(num_src-1); 
    nth = nt/2;  nxh = nx/2;
    wl2_mes = (complex *)malloc(nth*nxh*sizeof(complex));
    wl2_mes_d = (complex *)malloc(nth*nxh*sizeof(complex));

    /*  initialize                       */


     for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wl2_mes[r+nxh*t] = cmplx(0.0,0.0);
	wl2_mes_d[r+nxh*t] = cmplx(0.0,0.0);
     }

    for(dir=XUP;dir<=ZUP;dir++){

	FORALLSITES(i,s){
	    su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
	    su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[0].n[1]), num_src);
	}

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Start gather of forward time-like links */
	mtag[4] = start_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	   dir, EVENANDODD, gen_pt[4] );

	/* gather g_rand and light prop in direction dir  */
	mtag[0] = start_gather_site( F_OFFSET(dtmpvecs[0]),
	    sizeof(dble_su3_vec_src), dir, EVENANDODD, gen_pt[0]);


	for(r=0;r<nxh;r++){

	    /* Start gather for light prop for next time slice */
	    /* Do actual gathering below together with g_rand */
	    FORALLSITES(i,s){
		su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]),
		    num_src);
	    }

	    /* Collect forward random numbers g_rand in direction dir */
	    /* and make SU(3) matrix from qprop(0) * g_rand(R)^\dagger
	       in s_link_f */
	    wait_gather( mtag[0]);
	    FORALLSITES(i,s){
               dpt = (dble_su3_vec_src*)(gen_pt[0][i]);
                su3vecsrc_copy( &(dpt->n[0]),
                       &(s->resid_src), num_src);
	    }
	    FORALLSITES(i,s){
		su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[0]),
		    num_src);
		su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[1].n[1]),
		    num_src);
		su3vecsrc_outer_prod( &(s->qprop[0]), &(s->dtmpvecs[0].n[0]),
		    &(s->s_link_f), num_src);
	    }

	    /* Start gather of g_rand and light prop in ==TUP== direction */
	    mtag[1] = start_gather_site( F_OFFSET(dtmpvecs[1]),
		sizeof(dble_su3_vec_src), TUP, EVENANDODD, gen_pt[1]);

	    /* Start gather s_link_f for next time slice */
	    mtag[3] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[3]);

	    /* Collect light props in direction dir */
	    /* and make SU(3) matrix from g_rand(0) * qprop(R)^\dagger
	       in s_link */
	    FORALLSITES(i,s){
               dpt = (dble_su3_vec_src*)(gen_pt[0][i]);
                su3vecsrc_copy( &(dpt->n[1]),
                       &(s->resid_src), num_src);
	    }
	    FORALLSITES(i,s){
		su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[1]),
		    num_src);
		su3vecsrc_outer_prod( &(s->g_rand[0]), &(s->dtmpvecs[0].n[1]),
		    &(s->s_link), num_src);
	    }

	    /* Collect forward time-like links                 */
	    wait_gather( mtag[4]);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }


	    /* gather next light propagators, random numbers and links
	       in direction dir if needed  */
	    if( r<(nxh-1) ){
		restart_gather_site( F_OFFSET(dtmpvecs[0]), sizeof(dble_su3_vec_src),
		    dir, EVENANDODD, gen_pt[0], mtag[0] );

		restart_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[4], mtag[4] );
	    }
	    else{
		cleanup_gather( mtag[0]);
		cleanup_gather( mtag[4]);
	    }


	    /* Recursively compute the mixed "Wilson loops" 
	       of different time extent */

	    for(t=0;t<nth;t++){

		/* gather qprop, g_rand and s_link_f in TUP dir */
		wait_gather(mtag[1]);
		FORALLSITES(i,s){
                  dpt = (dble_su3_vec_src*)(gen_pt[1][i]);
                  su3vecsrc_copy( &(dpt->n[0]),
                         &(s->resid_src), num_src);
		}
		FORALLSITES(i,s){
		    su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[1].n[0]),
			num_src);
		}
		FORALLSITES(i,s){
                  dpt = (dble_su3_vec_src*)(gen_pt[1][i]);
                  su3vecsrc_copy( &(dpt->n[1]),
                         &(s->resid_src), num_src);
		}
		FORALLSITES(i,s){
		    su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[1].n[1]),
			num_src);
		}

		wait_gather(mtag[3]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[3][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), &(s->s_link_f));
		}


		/* start gather for next t if needed   */
		if( t<(nth-1) ){
		    restart_gather_site( F_OFFSET(dtmpvecs[1]),
			sizeof(dble_su3_vec_src), TUP, EVENANDODD,
			gen_pt[1], mtag[1] );

		    restart_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
			TUP, EVENANDODD,
			gen_pt[3], mtag[3] );
		}
		else{
		    cleanup_gather( mtag[1]);
		    cleanup_gather( mtag[3]);
		}


		/* construct operator by multiplying product of time-like
		links with light propagators and random sources */ 
		FORALLSITES(i,s){

		    /* "Diagonal" contribution */
		    for(j=0; j<num_src; j++){
			su3_src_to_vec(&(s->dtmpvecs[0].n[1]),
			    &(s->tempvec[1]), j);
			su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
			    &(s->tempvec[2]), j);
			su3_src_to_vec(&(s->dtmpvecs[1].n[1]),
			    &(s->tempvec[3]), j);
			if( ((s->t)+t+1)>=nt ){
			    mult_su3_mat_vec(&(s->link[TUP]),
				&(s->tempvec[2]), &tvec1);
			    cc1 = su3_dot(&(s->g_rand[j]), &tvec1);
			    mult_adj_su3_mat_vec(&(s->t_link_f),
				&(s->tempvec[1]), &tvec2);
			    cc2 = su3_dot(&(s->tempvec[3]), &tvec2);
			}
			else{
			    cc1 = su3_dot(&(s->g_rand[j]), &(s->tempvec[2]));
			    cc2 = su3_dot(&(s->tempvec[3]), &(s->tempvec[1]));
			}
			cc=cmul(&cc1,&cc2);
			wl2_mes_d[r+nxh*t].real += cc.real;
			wl2_mes_d[r+nxh*t].imag += cc.imag;
		    }

		    /* "num_src^2" contribution */
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), &(s->s_link_f), &tmat1);
			mult_su3_na( &tmat1, &(s->t_link_f), &tmat2);
			cc = complextrace_su3( &(s->s_link), &tmat2);
		    }
		    else{
			cc = complextrace_su3( &(s->s_link), &(s->s_link_f));
		    }
		    wl2_mes[r+nxh*t].real += cc.real;
		    wl2_mes[r+nxh*t].imag += cc.imag;
		}

	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over dir */

    /* Subtract diagonal contribution and normalize */
    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wl2_mes[r+nxh*t].real -= wl2_mes_d[r+nxh*t].real;
	wl2_mes[r+nxh*t].imag -= wl2_mes_d[r+nxh*t].imag;
	g_floatsum( &wl2_mes[r+nxh*t].real );
	wl2_mes[r+nxh*t].real /= (Real)(3*volume);
	g_floatsum( &wl2_mes[r+nxh*t].imag );
	wl2_mes[r+nxh*t].imag /= (Real)(3*volume);
	wl2_mes[r+nxh*t].real /= xnsrc;
	wl2_mes[r+nxh*t].imag /= xnsrc;
	if(this_node == 0)
	    printf("WL_2LC1_%d_OP_T%d %d %d %e %e\n", step, tot_smear, r, t,
		(double)wl2_mes[r+nxh*t].real, (double)wl2_mes[r+nxh*t].imag);
    }

 
    free( wl2_mes);
    free( wl2_mes_d);

} /* wl_2l_1corr.c */


