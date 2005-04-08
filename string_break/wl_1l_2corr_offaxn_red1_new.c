/******************* wl_1l_2corr_offaxn_red1_new.c  *********************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* v. 08/20/99  PL  */

/* Computes the operator consisting of a U-shaped product time-like links
   and a light propagator on gauge configurations in the axial gauge.
   Note: here the time-like gauge fields of the last time slice is copied
   in all other time slices as well, instead of the unit matrix. 

   ==  NOTE: here off-axis correlations calculated (ala UMH) ==

         ~~~~
        |    |
        ^    V      where  ~~~~  = light prop
        |    |
         -<-- 

     05/24/99:   NOTE: this version calculates ONLY sqrt(5)
                    and sqrt(8) correlations
                 Also:  nrmax=2+2

     08/20/99:   Introduce new vector "su3_vec_src"            
                 (and "dble_su3_vec_src") to send and receive
                 all components of g_rand and qprop at the
                 same time 


     03/26/00:  correction sqrt(5) loops  (dir1,-dir2)


*/

#include "string_break_includes.h"

void wl_1l_2corr_offax(int tot_smear, int step) {

register int i,dir1,dir2,dir3,r,t,r_off;
int nth,nxh,nrmax,j;
register site *s;
su3_vector tvec;
su3_matrix tmat1,tmat2;
dble_su3_vec_src *dpt;
msg_tag *mtag[8],*gmtag;
complex *wils_loop2, *wils_loop2_t, cc;
Real ftmp;
int disp[4];    /* displacement vector for general gather */

    if( nx != ny || nx != nz){
	if(this_node == 0)printf("w_loop2 gives wrong results for nx!=ny!=nz");
        return;
    }

    nth = nt/2;  
    nxh = 2;   /*   note   */
    nrmax = 4;
    wils_loop2 = (complex *)malloc(num_src*nth*nrmax*sizeof(complex));
    wils_loop2_t = (complex *)malloc(nth*nrmax*sizeof(complex));

    for(t=0;t<nth;t++) for(r=0;r<nrmax;r++){
	wils_loop2_t[r+nrmax*t] = cmplx(0.0,0.0);
    }

    for(t=0; t<nth; t++) for(r=0;r<nrmax;r++)for(j=0;j<num_src;j++){
        wils_loop2[r+nrmax*t+nrmax*nth*j] = cmplx(0.0,0.0);
    }


    /*  "sqrt(5)"  Correlations               */
    /*  ------------------------------------  */

    r_off = 0;

    for(dir1=XUP;dir1<=ZUP;dir1++) for(dir2=XUP;dir2<=ZUP;dir2++)
    if( dir1 != dir2){

	/* Do off-axis "sqrt(5)" loops in (dir1,dir2)-plane */
        /* ------------------------------------------------ */

        FORALLSITES(i,s){
            su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
        }

	/* First construct the "diagonal" link in (2*dir1,dir2) direction */

	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	/* Start gather of dir2-link from "2*dir1" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Gather the double links from dir2 direction */
	mtag[dir2] = start_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );

	/* Make first corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Make second corner and add to first */
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = 1;
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

        /* Recursively construct the space-like segments and compute
	  the Wilson loops with that segment */

	for(r=0;r<2;r++){

           FORALLSITES(i,s){
               su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]), num_src);
           }

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
                wait_general_gather( gmtag);
                FORALLSITES(i,s){
                    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
                }
                cleanup_general_gather( gmtag);

                FORALLSITES(i,s){
                    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
                }

                /* Inbetween gather time-like links across the diagonal. */
                gmtag = start_general_gather_site( F_OFFSET(t_link_f),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


           /* Collect forward time-like links */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
                su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
            }
            FORALLSITES(i,s){
                su3mat_copy( &(s->staple), &(s->t_link_f));
            }
            cleanup_general_gather( gmtag);


           /*  gather random numbers g_rand[j]  */

            gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[0].n[0]),
              sizeof(su3_vector_src),disp, EVENANDODD, gen_pt[5] );

            /* Collect forward light props       */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
              dpt = (dble_su3_vec_src*)(gen_pt[5][i]);
              su3vecsrc_copy( &(dpt->n[0]),&(s->resid_src), num_src);
            }
            FORALLSITES(i,s){
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[0].n[0]),
                    num_src);
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                    num_src);
            }
            cleanup_general_gather( gmtag);
  

            /* Start gather of g_rand and qprop in ==TUP== direction */
             mtag[7] = start_gather_site( F_OFFSET(dtmpvecs[1]), 
                sizeof(dble_su3_vec_src),TUP, EVENANDODD, gen_pt[7]);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

                /* Collect random numbers and light propagators in dir TUP    */
                wait_gather( mtag[7]);
                FORALLSITES(i,s){
                    su3vecsrc_copy((su3_vector_src *)(gen_pt[7][i]),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[0]),
                        num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(
                        (su3_vector_src *)(gen_pt[7][i]+sizeof(su3_vector_src)),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                        num_src);
                }


                /* Start gather for next t, if still needed. */
                if( t<(nth-1) ){
                    restart_gather_site( F_OFFSET(dtmpvecs[1]),
                        sizeof(dble_su3_vec_src),
                        TUP, EVENANDODD, gen_pt[7], mtag[7] );
                }
                else{
                    cleanup_gather( mtag[7]);
                }

                /*   average over sources                      */


             for(j=0; j<num_src; j++){

                  /* first rewrite current prop in standard format */

                    FORALLSITES(i,s){
                        su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
                            &(s->tempvec[3]), j);
                        su3_src_to_vec(&(s->dtmpvecs[1].n[1]),
                            &(s->tempvec[2]), j);
                    }

	       /* Finally, compute the Wilson loops. */

               /* construct operator by multiplying U-shape Wilson loop
                 with light propagator and then by random source  */
                  FORALLSITES(i,s){
                   if( ((s->t)+t+1)>=nt ){
                    mult_su3_an( &(s->s_link), &(s->link[TUP]),&tmat1);
                    mult_su3_an( &(s->t_link_f), &tmat1, &tmat2);
                    mult_su3_mat_vec(&tmat2,&(s->tempvec[3]),&tvec);
                   }
                   else
                   {
                    mult_adj_su3_mat_vec(&(s->s_link), &(s->tempvec[3]),&tvec);
                   }
                    cc = su3_dot(&(s->tempvec[2]),&tvec);
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real += cc.real;
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag += cc.imag;
                }
              }

	    } /* end loop over t */

	} /* end loop over r */


	/* Do off-axis "sqrt(5)" loops in (dir1,-dir2)-plane */
        /* ------------------------------------------------ */

        FORALLSITES(i,s){
            su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
        }

	/* First construct the "diagonal" link in the (2*dir1,-dir2) dir */

	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	/* Gather dir2-link from across the diagonal */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Make one corner and then gather it */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->s_link), &(s->staple));
	}
	mtag[dir2] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[dir2] );

	/* Make second corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Collect first corner and add to second */
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of time-like links across the diagonal. */
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<2;r++){

           FORALLSITES(i,s){
               su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]), num_src);
           }


	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
                wait_general_gather( gmtag);
                FORALLSITES(i,s){
                    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
                }
                cleanup_general_gather( gmtag);

                FORALLSITES(i,s){
                    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
                }

                /* Inbetween gather time-like links across the diagonal. */
                gmtag = start_general_gather_site( F_OFFSET(t_link_f),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


           /* Collect forward time-like links */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
                su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
            }
            FORALLSITES(i,s){
                su3mat_copy( &(s->staple), &(s->t_link_f));
            }
            cleanup_general_gather( gmtag);


           /*  gather random numbers g_rand[j]  */

            gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[0].n[0]),
              sizeof(su3_vector_src),disp, EVENANDODD, gen_pt[5] );

            /* Collect forward light props       */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
              dpt = (dble_su3_vec_src*)(gen_pt[5][i]);
              su3vecsrc_copy( &(dpt->n[0]),&(s->resid_src), num_src);
            }
            FORALLSITES(i,s){
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[0].n[0]),
                    num_src);
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                    num_src);
            }
            cleanup_general_gather( gmtag);


            /* Start gather of g_rand and qprop in ==TUP== direction */
             mtag[7] = start_gather_site( F_OFFSET(dtmpvecs[1]),
                sizeof(dble_su3_vec_src),TUP, EVENANDODD, gen_pt[7]);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }

	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

                /* Collect random numbers and light propagators in dir TUP    */
                wait_gather( mtag[7]);
                FORALLSITES(i,s){
                    su3vecsrc_copy((su3_vector_src *)(gen_pt[7][i]),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[0]),
                        num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(
                        (su3_vector_src *)(gen_pt[7][i]+sizeof(su3_vector_src)),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                        num_src);
                }


                /* Start gather for next t, if still needed. */
                if( t<(nth-1) ){
                    restart_gather_site( F_OFFSET(dtmpvecs[1]),
                        sizeof(dble_su3_vec_src),
                        TUP, EVENANDODD, gen_pt[7], mtag[7] );
                }
                else{
                    cleanup_gather( mtag[7]);
                }

                /*   average over sources                      */


             for(j=0; j<num_src; j++){

                  /* first rewrite current prop in standard format */

                    FORALLSITES(i,s){
                        su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
                            &(s->tempvec[3]), j);
                        su3_src_to_vec(&(s->dtmpvecs[1].n[1]),
                            &(s->tempvec[2]), j);
                    }

	       /* Finally, compute the Wilson loops. */

               /* construct operator by multiplying U-shape Wilson loop
                 with light propagator and then by random source  */
                  FORALLSITES(i,s){
                   if( ((s->t)+t+1)>=nt ){
                    mult_su3_an( &(s->s_link), &(s->link[TUP]),&tmat1);
                    mult_su3_an( &(s->t_link_f), &tmat1, &tmat2);
                    mult_su3_mat_vec(&tmat2,&(s->tempvec[3]),&tvec);
                   }
                   else
                   {
                    mult_adj_su3_mat_vec(&(s->s_link), &(s->tempvec[3]),&tvec);
                   }
                    cc = su3_dot(&(s->tempvec[2]),&tvec);
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real -= cc.real;
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag -= cc.imag;
                }
              }


	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over dir1 != dir2 */



    /*  "sqrt(8)"  Correlations               */
    /*  ------------------------------------  */

    /* Off-set for next bunch of Wilson loops */
    r_off = 2;

    for(dir1=XUP;dir1<=YUP;dir1++) for(dir2=dir1+1;dir2<=ZUP;dir2++){

	/* Do off-axis "sqrt(8)" loops in (dir1,dir2)-plane */
        /* ------------------------------------------------ */

        FORALLSITES(i,s){
            su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
        }

	/* First construct the "diagonal" link in (2*dir1,2*dir2) direction */
	/* First construct the "double links" */
	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );
	mtag[dir2] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );

	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Start gather of dir1-double-link from site "2*dir2" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = 2;
	gmtag = start_general_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	wait_gather( mtag[dir2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->s_link_f));
	}
	cleanup_gather( mtag[dir2]);

	/* Make first corner */
	wait_general_gather( gmtag);
        FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link_f), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
        }
	cleanup_general_gather( gmtag);

	/* Start gather of dir2-double-link from site "2*dir1" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	gmtag = start_general_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make second corner and add to first */
	ftmp = 0.5;
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_general_gather( gmtag);


	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = 2;
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );


        /* Recursively construct the space-like segments and compute
	  the Wilson loops with that segment */

	for(r=0;r<2;r++){

           FORALLSITES(i,s){
               su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]), num_src);
           }


	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
                wait_general_gather( gmtag);
                FORALLSITES(i,s){
                    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
                }
                cleanup_general_gather( gmtag);

                FORALLSITES(i,s){
                    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
                }

                /* Inbetween gather time-like links across the diagonal. */
                gmtag = start_general_gather_site( F_OFFSET(t_link_f),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


           /* Collect forward time-like links */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
                su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
            }
            FORALLSITES(i,s){
                su3mat_copy( &(s->staple), &(s->t_link_f));
            }
            cleanup_general_gather( gmtag);


            gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[0].n[0]),
              sizeof(su3_vector_src),disp, EVENANDODD, gen_pt[5] );

            /* Collect forward light props       */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
              dpt = (dble_su3_vec_src*)(gen_pt[5][i]);
              su3vecsrc_copy( &(dpt->n[0]),&(s->resid_src), num_src);
            }
            FORALLSITES(i,s){
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[0].n[0]),
                    num_src);
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                    num_src);
            }
            cleanup_general_gather( gmtag);

            /* Start gather of g_rand and qprop in ==TUP== direction */
             mtag[7] = start_gather_site( F_OFFSET(dtmpvecs[1]),
                sizeof(dble_su3_vec_src),TUP, EVENANDODD, gen_pt[7]);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

                /* Collect random numbers and light propagators in dir TUP    */
                wait_gather( mtag[7]);
                FORALLSITES(i,s){
                    su3vecsrc_copy((su3_vector_src *)(gen_pt[7][i]),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[0]),
                        num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(
                        (su3_vector_src *)(gen_pt[7][i]+sizeof(su3_vector_src)),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                        num_src);
                }


                /* Start gather for next t, if still needed. */
                if( t<(nth-1) ){
                    restart_gather_site( F_OFFSET(dtmpvecs[1]),
                        sizeof(dble_su3_vec_src),
                        TUP, EVENANDODD, gen_pt[7], mtag[7] );
                }
                else{
                    cleanup_gather( mtag[7]);
                }

                /*   average over sources                      */


             for(j=0; j<num_src; j++){

                  /* first rewrite current prop in standard format */

                    FORALLSITES(i,s){
                        su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
                            &(s->tempvec[3]), j);
                        su3_src_to_vec(&(s->dtmpvecs[1].n[1]),
                            &(s->tempvec[2]), j);
                    }

	       /* Finally, compute the Wilson loops. */

               /* construct operator by multiplying U-shape Wilson loop
                 with light propagator and then by random source  */
                  FORALLSITES(i,s){
                   if( ((s->t)+t+1)>=nt ){
                    mult_su3_an( &(s->s_link), &(s->link[TUP]),&tmat1);
                    mult_su3_an( &(s->t_link_f), &tmat1, &tmat2);
                    mult_su3_mat_vec(&tmat2,&(s->tempvec[3]),&tvec);
                   }
                   else
                   {
                    mult_adj_su3_mat_vec(&(s->s_link), &(s->tempvec[3]),&tvec);
                   }
                    cc = su3_dot(&(s->tempvec[2]),&tvec);
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real += cc.real;
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag += cc.imag;
                }
              }


	    } /* end loop over t */

	} /* end loop over r */


	/* Do off-axis "sqrt(8)" loops in (dir1,-dir2)-plane */
        /* ------------------------------------------------ */

        FORALLSITES(i,s){
            su3_vec_to_src( &(s->g_rand[0]), &(s->dtmpvecs[0].n[0]), num_src);
        }

	/* First construct the "diagonal" link in the (2*dir1,-2*dir2) dir */
	/* First construct the "double links" */
	mtag[dir2] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );
	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	wait_gather( mtag[dir2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->s_link_f));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of dir2-double-link from site "(2*dir1,-2*dir2)" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Make first corner */
        FORALLSITES(i,s){
	    mult_su3_an( &(s->s_link_f), &(s->s_link), &(s->staple));
        }

	/* Make second corner */
	wait_general_gather( gmtag);
        FORALLSITES(i,s){
	    mult_su3_na( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
        }
	cleanup_general_gather( gmtag);

	/* Start gather first corner from site "-2*dir2" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Collect first corner and add to second */
	ftmp = 0.5;
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[4][i]),
		 &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_general_gather( gmtag);


	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] =  2;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

        /* Recursively construct the space-like segments and compute
	  the Wilson loops with that segment */

	for(r=0;r<2;r++){

           FORALLSITES(i,s){
               su3_vec_to_src( &(s->qprop[0]), &(s->dtmpvecs[1].n[0]), num_src);
           }


	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
                wait_general_gather( gmtag);
                FORALLSITES(i,s){
                    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
                }
                cleanup_general_gather( gmtag);

                FORALLSITES(i,s){
                    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
                }

                /* Inbetween gather time-like links across the diagonal. */
                gmtag = start_general_gather_site( F_OFFSET(t_link_f),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


           /* Collect forward time-like links */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
                su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
            }
            FORALLSITES(i,s){
                su3mat_copy( &(s->staple), &(s->t_link_f));
            }
            cleanup_general_gather( gmtag);


           /*  gather random numbers g_rand[j]  */

            gmtag = start_general_gather_site( F_OFFSET(dtmpvecs[0].n[0]),
              sizeof(su3_vector_src),disp, EVENANDODD, gen_pt[5] );

            /* Collect forward light props       */
            wait_general_gather( gmtag);
            FORALLSITES(i,s){
              dpt = (dble_su3_vec_src*)(gen_pt[5][i]);
              su3vecsrc_copy( &(dpt->n[0]),&(s->resid_src), num_src);
            }
            FORALLSITES(i,s){
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[0].n[0]),
                    num_src);
                su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                    num_src);
            }
            cleanup_general_gather( gmtag);
 

            /* Start gather of g_rand and qprop in ==TUP== direction */
             mtag[7] = start_gather_site( F_OFFSET(dtmpvecs[1]),
                sizeof(dble_su3_vec_src),TUP, EVENANDODD, gen_pt[7]);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

                /* Collect random numbers and light propagators in dir TUP    */
                wait_gather( mtag[7]);
                FORALLSITES(i,s){
                    su3vecsrc_copy((su3_vector_src *)(gen_pt[7][i]),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[0]),
                        num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(
                        (su3_vector_src *)(gen_pt[7][i]+sizeof(su3_vector_src)),
                        &(s->resid_src), num_src);
                }
                FORALLSITES(i,s){
                    su3vecsrc_copy(&(s->resid_src), &(s->dtmpvecs[1].n[1]),
                        num_src);
                }


                /* Start gather for next t, if still needed. */
                if( t<(nth-1) ){
                    restart_gather_site( F_OFFSET(dtmpvecs[1]),
                        sizeof(dble_su3_vec_src),
                        TUP, EVENANDODD, gen_pt[7], mtag[7] );
                }
                else{
                    cleanup_gather( mtag[7]);
                }

                /*   average over sources                      */


             for(j=0; j<num_src; j++){

                  /* first rewrite current prop in standard format */

                    FORALLSITES(i,s){
                        su3_src_to_vec(&(s->dtmpvecs[1].n[0]),
                            &(s->tempvec[3]), j);
                        su3_src_to_vec(&(s->dtmpvecs[1].n[1]),
                            &(s->tempvec[2]), j);
                    }

	       /* Finally, compute the Wilson loops. */

               /* construct operator by multiplying U-shape Wilson loop
                 with light propagator and then by random source  */
                  FORALLSITES(i,s){
                   if( ((s->t)+t+1)>=nt ){
                    mult_su3_an( &(s->s_link), &(s->link[TUP]),&tmat1);
                    mult_su3_an( &(s->t_link_f), &tmat1, &tmat2);
                    mult_su3_mat_vec(&tmat2,&(s->tempvec[3]),&tvec);
                   }
                   else
                   {
                    mult_adj_su3_mat_vec(&(s->s_link), &(s->tempvec[3]),&tvec);
                   }
                    cc = su3_dot(&(s->tempvec[2]),&tvec);
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real += cc.real;
                     wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag += cc.imag;
                }
              }


	    } /* end loop over t */

	} /* end loop over r */


    } /* end loop over dir1 != dir2 */


    /* Normalize and print the Wilson loops */
         for(j=0;j<num_src;j++) for(t=0;t<nth;t++){

	for(r=0;r<2;r++){
	    g_floatsum(&wils_loop2[r+nrmax*t+nrmax*nth*j].real);
            wils_loop2[r+nrmax*t+nrmax*nth*j].real /= (Real)(12*volume);
	    g_floatsum(&wils_loop2[r+nrmax*t+nrmax*nth*j].imag);
            wils_loop2[r+nrmax*t+nrmax*nth*j].imag /= (Real)(12*volume);
/*
            if(this_node == 0)
            printf("WL_1LC2_OA_%d_OP_%d %d %d %d %e %e\n", step,tot_smear,j,r,
              t,(double)wils_loop2[r+nrmax*t+nrmax*nth*j].real,
              (double)wils_loop2[r+nrmax*t+nrmax*nth*j].imag);
*/
             wils_loop2_t[r+nrmax*t].real 
               += wils_loop2[r+nrmax*t+nrmax*nth*j].real;
             wils_loop2_t[r+nrmax*t].imag 
                += wils_loop2[r+nrmax*t+nrmax*nth*j].imag;
        }

	r_off = 2;
	for(r=0;r<2;r++){
	    g_floatsum(&wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real);
            wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real /= (Real)(6*volume);
	    g_floatsum(&wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag);
            wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag /= (Real)(6*volume);
/*
            if(this_node == 0)
            printf("WL_1LC2_OA_%d_OP_%d %d %d %d %e %e\n", 
              step,tot_smear,j,r+r_off,
              t,(double)wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real,
              (double)wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag);
*/
             wils_loop2_t[r+r_off+nrmax*t].real 
               += wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].real;
             wils_loop2_t[r+r_off+nrmax*t].imag 
                += wils_loop2[r+r_off+nrmax*t+nrmax*nth*j].imag;
        }

     }


    /* Normalize and print the Wilson loops */

     for(t=0;t<nth;t++){

	for(r=0;r<nrmax;r++){
           wils_loop2_t[r+nrmax*t].real /= (Real)num_src;
           wils_loop2_t[r+nrmax*t].imag /= (Real)num_src;
           if(this_node == 0)
            printf("WL_1LC2_OA_T_%d_OP_%d %d %d %e %e\n", step,tot_smear,r,t,
              (double)wils_loop2_t[r+nrmax*t].real,
              (double)wils_loop2_t[r+nrmax*t].imag);
        }
     }
	   

    free( wils_loop2);
    free( wils_loop2_t);

} /* wl_1l_2corr_offax  */


