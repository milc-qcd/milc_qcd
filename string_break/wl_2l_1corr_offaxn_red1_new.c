/******************** wl_2l_1corr_offaxn_red1_new.c ********************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 08/25/99  PL  */

/* Computes the combination of the product of 2 light quarks
   with 2 legs of time-like links to form a kind of "Wilson loop".
   Note: here the time-like gauge fields of the last time slice is copied
   in all other time slices as well, instead of the unit matrix.

   ==  NOTE: here off-axis correlations calculated (ala UMH) == 

         ~~~~
        |    |
        ^    V      where  ~~~~  = light prop
        |    |
         ~~~~

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

void wl_2l_1corr_offax(int tot_smear, int step) {

register int i,dir1,dir2,dir3,r,t,r_off;
int nth,nxh,nrmax,j,k,num_lp,rmax;
register site *s;
su3_vector tvec1,tvec2;
su3_matrix tmat1,tmat2;
dble_su3_vec_src *dpt;
msg_tag *mtag[8],*gmtag;
complex  cc, cc1, cc2;
complex *wl2_mes, *wl2_mes_d;
Real ftmpr, ftmpi, xnsrc;
int disp[4];    /* displacement vector for general gather */

    if( nx != ny || nx != nz){
	if(this_node == 0)printf("w_loop2 gives wrong results for nx!=ny!=nz");
        return;
    }

    xnsrc=num_src*(num_src-1);
    nth = nt/2;  nxh = nx/2; nrmax=4;
    wl2_mes = (complex *)malloc(nth*nrmax*sizeof(complex));
    wl2_mes_d = (complex *)malloc(nth*nrmax*sizeof(complex));

    for(t=0;t<nth;t++) for(r=0;r<nrmax;r++){
        wl2_mes[r+nrmax*t] = cmplx(0.0,0.0);
        wl2_mes_d[r+nrmax*t] = cmplx(0.0,0.0);
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
            su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
        }


        FORALLSITES(i,s){
          su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[0].n[1]), num_src);
          su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
        }

        gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[0]),
           sizeof(dble_su3_vec_src),disp, EVENANDODD, gen_pt[0]);


        /* Recursively construct the space-like segments and compute
           the Wilson loops with that segment */

        for(r=0;r<rmax;r++){


	    /* Start gather for light prop for next time slice */
	    /* Do actual gathering below together with g_rand */
	    FORALLSITES(i,s){
		su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]),
		    num_src);
	    }

	    /* Collect forward random numbers g_rand in direction dir */
	    /* and make SU(3) matrix from qprop(0) * g_rand(R)^\dagger
	       in s_link_f */
            wait_general_gather( gmtag);
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

            cleanup_general_gather( gmtag);

            /* Gather time-like links across the diagonal. */
            gmtag = start_general_gather_site( F_OFFSET(t_link_f),
                sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

	    FORALLSITES(i,s){
		su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[0].n[1]),
		    num_src);
		su3vecsrc_outer_prod( &(s->g_rand[0]), &(s->dtmpvecs[0].n[1]),
		    &(s->s_link), num_src);
	    }

            /* Collect time-like links */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
                su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
            }
            cleanup_general_gather( gmtag);


	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }

	    /* gather next light propagators, random numbers and links
	       in direction dir if needed  */
	    if( r<(rmax-1) ){
		gmtag = start_general_gather_site(F_OFFSET(dtmpvecs[0]),
                    sizeof(dble_su3_vec_src),
		    disp, EVENANDODD, gen_pt[0] );

             }

	    /* Recursively compute the mixed "Wilson loops" 
	       of different time extent */

	    for(t=0;t<nth;t++){

		/* gather qprop, g_rand and s_link_f in TUP dir */
		wait_gather(mtag[1]);
		FORALLSITES(i,s){
		    su3vecsrc_copy( (su3_vector_src *)(gen_pt[1][i]),
			&(s->resid_src), num_src);
		}
		FORALLSITES(i,s){
		    su3vecsrc_copy( &(s->resid_src), &(s->dtmpvecs[1].n[0]),
			num_src);
		}
		FORALLSITES(i,s){
		    su3vecsrc_copy(
			(su3_vector_src *)(gen_pt[1][i]+sizeof(su3_vector_src)),
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
			wl2_mes_d[r+r_off+nrmax*t].real += cc.real;
			wl2_mes_d[r+r_off+nrmax*t].imag += cc.imag;
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
		    wl2_mes[r+r_off+nrmax*t].real += cc.real;
		    wl2_mes[r+r_off+nrmax*t].imag += cc.imag;
		}

	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over lp  */


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
            printf("WL_2LC1_OA_T_%d_OP_%d %d %d %e %e\n", 
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
            printf("WL_2LC1_OA_T_%d_OP_%d %d %d %e %e\n", step, tot_smear,
                 r+r_off, t,
                (double)wl2_mes[r+r_off+nrmax*t].real,
                (double)wl2_mes[r+r_off+nrmax*t].imag);
        }

        }


      free( wl2_mes);
      free( wl2_mes_d);


} /* wl_2l_1corr_offax */



