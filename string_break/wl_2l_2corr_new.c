/************************** wl_2l_2corr_new.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 01/03/99  PL  */

/* Computes the combination of the product of 2 light quarks
   with 2 legs of time-like links to form a kind of "Wilson loop".
   Note: here the time-like gauge fields of the last time slice is copied
   in all other time slices as well, instead of the unit matrix. */
/*  
         ~     ~ 
        | ~   ~ |
        ^ ~   ~ V      where  ~~~~  = light prop 
        | ~   ~ | 
         ~     ~
      
        0       R  
*/


#include "string_break_includes.h"

void wl_2l_2corr(int tot_smear, int step) {

register int i,dir,r,t;
int nth,nxh,j,k;
register site *s;
su3_vector tvec1,tvec2;
msg_tag *mtag[3];
dble_su3_vec_src *dpt;
complex  cc, cc1, cc2;
complex *wl2_mes, *wl2_mes_d;
Real xnsrc;
field_offset meson1, meson2, meson_f, mes_t;

meson1 = F_OFFSET(xxx.c[0]);
meson2 = F_OFFSET(xxx.c[1]);
meson_f = F_OFFSET(ttt.c[0]);
mes_t = F_OFFSET(ttt.c[1]);

    if( nx != ny || nx != nz){
	if(this_node == 0)printf("wl_2lprop gives wrong results for nx!=ny!=nz");
	return;
    }

    xnsrc=num_src*(num_src-1); 
    nth = nt/2;  nxh = nx/2;
    wl2_mes = (complex *)malloc(nth*nxh*sizeof(complex));
    wl2_mes_d = (complex *)malloc(nth*nxh*sizeof(complex));

    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wl2_mes[r+nxh*t] = cmplx(0.0,0.0);
	wl2_mes_d[r+nxh*t] = cmplx(0.0,0.0);
    }


    FORALLSITES(i,s){
	su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[0].n[0]), num_src);
	su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[1]), num_src);
    }

    mtag[0] = start_gather_site( F_OFFSET(dtmpvecs[0]), sizeof(dble_su3_vec_src),
	TUP, EVENANDODD, gen_pt[0]);

    for(t=0;t<nth;t++){

	/* Construct the static-light meson and anti-meson */
	wait_gather( mtag[0]);
	FORALLSITES(i,s){
           dpt = (dble_su3_vec_src*)(gen_pt[0][i]);
           su3vecsrc_copy( &(dpt->n[0]),
                  &(s->resid_src), num_src);
	}
	FORALLSITES(i,s){
	    su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[0]), num_src);
	}
	FORALLSITES(i,s){
           dpt = (dble_su3_vec_src*)(gen_pt[0][i]);
           su3vecsrc_copy( &(dpt->n[1]),
                  &(s->resid_src), num_src);
	}
	FORALLSITES(i,s){
	    su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[1]), num_src);
	}

	if( t<(nth-1) ){
	    restart_gather_site( F_OFFSET(dtmpvecs[0]), sizeof(dble_su3_vec_src),
		TUP, EVENANDODD, gen_pt[0], mtag[0] );
	}
	else{
	    cleanup_gather( mtag[0]);
	}

	/* Construct static-light meson and anti-meson */
	FORALLSITES(i,s){
	    (*(complex *)F_PT(s,meson1)) = cmplx(0.0,0.0);
	    (*(complex *)F_PT(s,meson2)) = cmplx(0.0,0.0);
	    for(j=0; j<num_src; j++){ 
		su3_src_to_vec(&(s->dtmpvecs[0].n[0]), &(s->tempvec[0]), j);
		su3_src_to_vec(&(s->dtmpvecs[0].n[1]), &(s->tempvec[1]), j);
		if( ((s->t)+t+1)>=nt ){
		    mult_su3_mat_vec(&(s->link[TUP]), &(s->tempvec[0]), &tvec1);
		    cc1 = su3_dot(&(s->g_rand[j]), &tvec1);
		    mult_adj_su3_mat_vec(&(s->link[TUP]), &(s->qprop[j]),
			&tvec2);
		    cc2 = su3_dot(&(s->tempvec[1]), &tvec2);
		}
		else
		{
		    cc1 = su3_dot(&(s->g_rand[j]), &(s->tempvec[0]));
		    cc2 = su3_dot(&(s->tempvec[1]), &(s->qprop[j]));
		}
		CSUM((*(complex *)F_PT(s,meson1)), cc1);
		CSUM((*(complex *)F_PT(s,meson2)), cc2);
		s->dtmpvecs[1].n[0].s[j].c[0] = cc1;
		s->dtmpvecs[1].n[0].s[j].c[1] = cc2;
	    }
	}

	for(dir=XUP;dir<=ZUP;dir++){

	    FORALLSITES(i,s){
		(*(complex *)F_PT(s,meson_f)) = (*(complex *)F_PT(s,meson2));
		for(j=0; j<num_src; j++){ 
		    s->dtmpvecs[1].n[1].s[j].c[0] =
			s->dtmpvecs[1].n[0].s[j].c[1];
		}
	    }

	    /* Start gather of forward meson */
	    mtag[1] = start_gather_site( F_OFFSET(dtmpvecs[1].n[1]),
		sizeof(su3_vector_src), dir, EVENANDODD, gen_pt[1] );
   	    mtag[2] = start_gather_site( meson_f, sizeof(complex),
		dir, EVENANDODD, gen_pt[2] );

	    for(r=0;r<nxh;r++){
		/* Collect forward meson  and compute meson -
		   anti-meson correlation function */
		wait_gather( mtag[1]);
		wait_gather( mtag[2]);
		FORALLSITES(i,s){
		    (*(complex *)F_PT(s,mes_t)) =
			(*(complex *)(gen_pt[2][i]));
		    for(j=0; j<num_src; j++){ 
			s->resid_src.s[j].c[2] =
			    (*(su3_vector_src *)(gen_pt[1][i])).s[j].c[0];
		    }
		}
		FORALLSITES(i,s){
		    (*(complex *)F_PT(s,meson_f)) =
			(*(complex *)F_PT(s,mes_t));

		    cc = cmul(((complex *)F_PT(s,meson1)),
			((complex *)F_PT(s,meson_f)) );
		    wl2_mes[r+nxh*t].real += cc.real;
		    wl2_mes[r+nxh*t].imag += cc.imag;
		    for(j=0; j<num_src; j++){ 
			s->dtmpvecs[1].n[1].s[j].c[0] = s->resid_src.s[j].c[2];
			cc = cmul(&(s->dtmpvecs[1].n[0].s[j].c[0]),
			    &(s->dtmpvecs[1].n[1].s[j].c[0]));
			wl2_mes_d[r+nxh*t].real += cc.real;
			wl2_mes_d[r+nxh*t].imag += cc.imag;
		    }
		}

		if( r<(nxh-1) ){
		    restart_gather_site( F_OFFSET(dtmpvecs[1].n[1]),
			sizeof(su3_vector_src), dir, EVENANDODD,
			gen_pt[1], mtag[1] );

		    restart_gather_site( meson_f, sizeof(complex),
			dir, EVENANDODD, gen_pt[2], mtag[2] );
		}
		else{
		    cleanup_gather( mtag[1]);
		    cleanup_gather( mtag[2]);
		}

	    } /* end loop over r */

	} /* end loop over dir */

    } /* end loop over t */


    /* Normalize and print the Wilson-light operator  */
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
	    printf("WL_2LC2_%d_OP_T%d %d %d %e %e\n", step, tot_smear, r, t,
		(double)wl2_mes[r+nxh*t].real, (double)wl2_mes[r+nxh*t].imag);
    }

    free( wl2_mes);
    free( wl2_mes_d);

} /* wl_2l_2corr.c */

