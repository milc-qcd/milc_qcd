/************************** wl_2l_2corr_offaxn_red1_new.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 05/24/99  PL  */

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

     05/24/99:  NOTE:  calculate ONLY sqrt(5) and sqrt(8) correlations
                Also:  nrmax = 2 + 2


     08/20/99:   Introduce new vector "su3_vec_src"
                 (and "dble_su3_vec_src") to send and receive
                 all components of g_rand and qprop at the
                 same time
                
                 Here: use O(N) algorithm to calculate 
                 sum over sources j and k (sum_k_j)

*/


#include "string_break_includes.h"

void wl_2l_2corr_offax(int tot_smear, int step) {

register int i,dir1,dir2,dir3,r,t,r_off;
int nth,nxh,nrmax,j,k,num_lp,rmax;
register site *s;
dble_su3_vec_src *dpt;
su3_vector tvec1,tvec2;
msg_tag *mtag[8],*gmtag;
complex  cc,cc1,cc2;
complex *wl2_mes, *wl2_mes_d;
Real ftmpr, ftmpi, fns, xnsrc;
field_offset meson1, meson2, meson_f, mes_t;
int disp[4];    /* displacement vector for general gather */


meson1 = F_OFFSET(xxx.c[0]);
meson2 = F_OFFSET(xxx.c[1]);
meson_f = F_OFFSET(ttt.c[0]);
mes_t = F_OFFSET(ttt.c[1]);

    if( nx != ny || nx != nz){
     if(this_node == 0)printf("wl_2lprop gives wrong results for nx!=ny!=nz");
	return;
    }

    xnsrc = num_src*(num_src-1);
    nth = nt/2;  nxh = nx/2;  nrmax = 4;
    wl2_mes = (complex *)malloc(nth*nrmax*sizeof(complex));
    wl2_mes_d = (complex *)malloc(nth*nrmax*sizeof(complex));

    for(t=0;t<nth;t++) for(r=0;r<nrmax;r++){
        wl2_mes[r+nrmax*t] = cmplx(0.0,0.0);
        wl2_mes_d[r+nrmax*t] = cmplx(0.0,0.0);
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
            su3vecsrc_copy( &(dpt->n[0]),&(s->resid_src), num_src);
        }
        FORALLSITES(i,s){
            su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[0]), num_src);
        }
        FORALLSITES(i,s){
            dpt = (dble_su3_vec_src*)(gen_pt[0][i]);
            su3vecsrc_copy( &(dpt->n[1]),&(s->resid_src), num_src);
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


	/* Loop over the different loops and orientations */

	for(num_lp=0; num_lp<30; num_lp++) {

	for(i=XUP;i<=TUP;i++)disp[i]=0;

	if (num_lp < 24) {
	  /* off-axis "sqrt(5)" loop */
	  r_off = 0;
	  rmax = 2;
          switch (num_lp){
	    case 0:
	      disp[XUP] = 2;
	      disp[YUP] = 1;
	      break;
	    case 1:
	      disp[XUP] = 2;
	      disp[ZUP] = 1;
	      break;
	    case 2:
	      disp[YUP] = 2;
	      disp[XUP] = 1;
	      break;
	    case 3:
	      disp[YUP] = 2;
	      disp[ZUP] = 1;
	      break;
	    case 4:
	      disp[ZUP] = 2;
	      disp[XUP] = 1;
	      break;
	    case 5:
	      disp[ZUP] = 2;
	      disp[YUP] = 1;
	      break;
	    case 6:
	      disp[XUP] = 2;
	      disp[YUP] = -1;
	      break;
	    case 7:
	      disp[XUP] = 2;
	      disp[ZUP] = -1;
	      break;
	    case 8:
	      disp[YUP] = 2;
	      disp[XUP] = -1;
	      break;
	    case 9:
	      disp[YUP] = 2;
	      disp[ZUP] = -1;
	      break;
	    case 10:
	      disp[ZUP] = 2;
	      disp[XUP] = -1;
	      break;
	    case 11:
	      disp[ZUP] = 2;
	      disp[YUP] = -1;
	      break;
	    case 12:
	      disp[XUP] = 1;
	      disp[YUP] = 2;
	      break;
	    case 13:
	      disp[XUP] = 1;
	      disp[ZUP] = 2;
	      break;
	    case 14:
	      disp[YUP] = 1;
	      disp[XUP] = 2;
	      break;
	    case 15:
	      disp[YUP] = 1;
	      disp[ZUP] = 2;
	      break;
	    case 16:
	      disp[ZUP] = 1;
	      disp[XUP] = 2;
	      break;
	    case 17:
	      disp[ZUP] = 1;
	      disp[YUP] = 2;
	      break;
	    case 18:
	      disp[XUP] = 1;
	      disp[YUP] = -2;
	      break;
	    case 19:
	      disp[XUP] = 1;
	      disp[ZUP] = -2;
	      break;
	    case 20:
	      disp[YUP] = 1;
	      disp[XUP] = -2;
	      break;
	    case 21:
	      disp[YUP] = 1;
	      disp[ZUP] = -2;
	      break;
	    case 22:
	      disp[ZUP] = 1;
	      disp[XUP] = -2;
	      break;
	    case 23:
	      disp[ZUP] = 1;
	      disp[YUP] = -2;
	      break;
	  }
	}
	else {
	  /* off-axis "sqrt(8)" loop */
	  r_off = 2;
	  rmax = 2;
          switch (num_lp){
            case 24:
              disp[XUP] = 2;
              disp[YUP] = 2;
              break;
            case 25:
              disp[XUP] = 2;
              disp[ZUP] = 2;
              break;
            case 26:
              disp[YUP] = 2;
              disp[ZUP] = 2;
              break;
            case 27:
              disp[XUP] = 2;
              disp[YUP] = -2;
              break;
            case 28:
              disp[XUP] = 2;
              disp[ZUP] = -2;
              break;
            case 29:
              disp[YUP] = 2;
              disp[ZUP] = -2;
              break;
	  }
	}

            FORALLSITES(i,s){
                (*(complex *)F_PT(s,meson_f)) = (*(complex *)F_PT(s,meson2));
                for(j=0; j<num_src; j++){
                    s->dtmpvecs[1].n[1].s[j].c[0] =
                        s->dtmpvecs[1].n[0].s[j].c[1];
                }
            }


	    /* Start gather of forward meson */

	    gmtag = start_general_gather_site( meson_f, sizeof(complex),
		disp, EVENANDODD, gen_pt[4] );

   	    wait_general_gather( gmtag);
	     FORALLSITES(i,s){
	       ((complex *)F_PT(s,mes_t))->real =
	 	  ((complex *)(gen_pt[4][i]))->real;
  	      ((complex *)F_PT(s,mes_t))->imag =
	 	  ((complex *)(gen_pt[4][i]))->imag;
		}
	    cleanup_general_gather( gmtag);

	    gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[1].n[1]),
                 sizeof(su3_vector_src),disp,EVENANDODD, gen_pt[1] );


	    for(r=0;r<rmax;r++){

		/* Collect forward meson  and compute meson -
		   anti-meson correlation function */

		wait_general_gather( gmtag);

                FORALLSITES(i,s){
                    for(j=0; j<num_src; j++){
                        s->resid_src.s[j].c[2] =
                            (*(su3_vector_src *)(gen_pt[1][i])).s[j].c[0];
                    }
                }
	        cleanup_general_gather( gmtag);

                FORALLSITES(i,s){
                    (*(complex *)F_PT(s,meson_f)) =
                        (*(complex *)F_PT(s,mes_t));

                    cc = cmul(((complex *)F_PT(s,meson1)),
                        ((complex *)F_PT(s,meson_f)) );
                    wl2_mes[r+r_off+nrmax*t].real += cc.real;
                    wl2_mes[r+r_off+nrmax*t].imag += cc.imag;
                    for(j=0; j<num_src; j++){
                        s->dtmpvecs[1].n[1].s[j].c[0] = s->resid_src.s[j].c[2];
                        cc = cmul(&(s->dtmpvecs[1].n[0].s[j].c[0]),
                            &(s->dtmpvecs[1].n[1].s[j].c[0]));
                        wl2_mes_d[r+r_off+nrmax*t].real += cc.real;
                        wl2_mes_d[r+r_off+nrmax*t].imag += cc.imag;
                    }

                }

                if( r<(rmax-1) ){

	         gmtag = start_general_gather_site( meson_f, sizeof(complex),
		   disp, EVENANDODD, gen_pt[4] );

   	         wait_general_gather( gmtag);
	         FORALLSITES(i,s){
	            ((complex *)F_PT(s,mes_t))->real =
	 	      ((complex *)(gen_pt[4][i]))->real;
  	            ((complex *)F_PT(s,mes_t))->imag =
	 	      ((complex *)(gen_pt[4][i]))->imag;
		 }
	         cleanup_general_gather( gmtag);


	        gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[1].n[1]),
                   sizeof(su3_vector_src),disp,EVENANDODD, gen_pt[1] );

                }

	    } /* end loop over r */

	} /* end loop over num_lp */

    } /* end loop over t */

    /* Normalize and print the Wilson-light operator  */
    for(t=0;t<nth;t++) {

	r_off = 0;
	for(r=0;r<2;r++){
        wl2_mes[r+r_off+nrmax*t].real -= wl2_mes_d[r+r_off+nrmax*t].real;
        wl2_mes[r+r_off+nrmax*t].imag -= wl2_mes_d[r+r_off+nrmax*t].imag;
        g_floatsum( &wl2_mes[r+r_off+nrmax*t].real );
        wl2_mes[r+r_off+nrmax*t].real /= (Real)(24*volume);
        g_floatsum( &wl2_mes[r+r_off+nrmax*t].imag );
        wl2_mes[r+r_off+nrmax*t].imag /= (Real)(24*volume);
        wl2_mes[r+r_off+nrmax*t].real /= xnsrc;
        wl2_mes[r+r_off+nrmax*t].imag /= xnsrc;
        if(this_node == 0)
            printf("WL_2LC2_OA_T_%d_OP_%d %d %d %e %e\n", 
                 step, tot_smear, r, t,
                (double)wl2_mes[r+r_off+nrmax*t].real, 
                (double)wl2_mes[r+r_off+nrmax*t].imag);
        }


	r_off = 2;
	for(r=0;r<2;r++){
        wl2_mes[r+r_off+nrmax*t].real -= wl2_mes_d[r+r_off+nrmax*t].real;
        wl2_mes[r+r_off+nrmax*t].imag -= wl2_mes_d[r+r_off+nrmax*t].imag;
        g_floatsum( &wl2_mes[r+r_off+nrmax*t].real );
        wl2_mes[r+r_off+nrmax*t].real /= (Real)(6*volume);
        g_floatsum( &wl2_mes[r+r_off+nrmax*t].imag );
        wl2_mes[r+r_off+nrmax*t].imag /= (Real)(6*volume);
        wl2_mes[r+r_off+nrmax*t].real /= xnsrc;
        wl2_mes[r+r_off+nrmax*t].imag /= xnsrc;
        if(this_node == 0)
            printf("WL_2LC2_OA_T_%d_OP_%d %d %d %e %e\n", step, tot_smear,
                 r+r_off, t,
                (double)wl2_mes[r+r_off+nrmax*t].real, 
                (double)wl2_mes[r+r_off+nrmax*t].imag);
        }

        }


      free( wl2_mes);
      free( wl2_mes_d);


} /* wl_2l_2corr_offax */


