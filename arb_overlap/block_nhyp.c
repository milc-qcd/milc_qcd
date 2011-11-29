/********** block_nhyp.c ******************************/
/* MIMD version 7       */
/* adapted YS Sept 2008 */

/* Reference:
* Hypercubic smeared links for dynamical fermions.
* By Anna Hasenfratz, Roland Hoffmann, Stefan Schaefer.
* JHEP 0705:029,2007. [hep-lat/0702028]
*/

#include "arb_ov_includes.h"

/* parts of this code are specific to SU(3)
   valid for any NCOL: calculation of Omega, Staple, and Q
   valid for SU(3) only: calculation of Q^{-1/2}, including compute_fhb()
*/
void scalar_add_diag_su3(su3_matrix *a, Real s);


void staple_nhyp( int dir1, int dir2, su3_matrix *lnk1, su3_matrix *lnk2,
                  su3_matrix *stp );

void block_nhyp1();
void block_nhyp2();
void block_nhyp3();


/* do n smearing levels
   where n<=3 is the value of SMEAR_LEVEL (set in defines.h)
*/
void block_nhyp()
{
#if (SMEAR_LEVEL==3)
    block_nhyp1();
#endif
#if (SMEAR_LEVEL>1)
    block_nhyp2();
#endif
    block_nhyp3();
} /* block_nhyp */


void block_nhyp3()
{
    register int dir, dir2, i;
    register site *st;
    Real f[3];   /* related code is specific to SU(3) */
    Real ftmp1,ftmp2;
    su3_matrix tmat, Omega, eQ,  Q, Q2;

    ftmp1=alpha_smear[0]/(6.*(1.-alpha_smear[0]));
    ftmp2=1.-alpha_smear[0];

    for(dir=XUP;dir<=TUP;dir++){

	/* compute the staple */
	FORALLDYNLINKS(i,st,dir)  clear_su3mat(&Staple3[dir][i]);
	for(dir2=XUP;dir2<=TUP;dir2++) if(dir2!=dir){
#if (SMEAR_LEVEL>1)
	    staple_nhyp(dir,dir2,hyplink2[dir2][dir],
		        hyplink2[dir][dir2],Staple3[dir]);
#else /* one-level only */
	    staple_nhyp(dir,dir2,gauge_field_thin[dir],
                        gauge_field_thin[dir2],Staple3[dir]);
#endif
	}

	FORALLDYNLINKS(i,st,dir){
	    /* make Omega  */
	    scalar_mult_add_su3_matrix(gauge_field_thin[dir]+i,
                                         Staple3[dir]+i,ftmp1 ,&Q);
	    scalar_mult_su3_matrix(&Q,ftmp2,&Omega);
	    Staple3[dir][i]=Omega;
	    mult_su3_an(&Omega,&Omega,&Q);
            /* IR regulator, see clover_xxx/defines.h               */
            scalar_add_diag_su3(&Q,IR_STAB);
#ifndef NHYP_DEBUG
	    compute_fhb(&Q,f,NULL, 0);
#else
            compute_fhb(&Omega,&Q,f,NULL, 0);
#endif

	    /* make Q**2 */
	    mult_su3_nn(&Q,&Q,&Q2);

	    /* compute Q^(-1/2) via Eq. 19  */
	    scalar_mult_su3_matrix(&Q,f[1],&tmat);
	    scalar_mult_add_su3_matrix(&tmat,&Q2,f[2],&eQ);
	    scalar_add_diag_su3(&eQ,f[0]);

	    /* multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)  */
	    mult_su3_nn(&Omega,&eQ,gauge_field[dir]+i);
	}

    } /* dir */
} /* block_nhyp3 */


#if (SMEAR_LEVEL>1)
void block_nhyp2()
{
    register int dir, dir2, dir3, dir4, i;
    register site *st;
    Real f[3], ftmp1, ftmp2;
    su3_matrix tmat, Omega, eQ,  Q, Q2;

    ftmp1=alpha_smear[1]/(4.*(1.-alpha_smear[1]));
    ftmp2=(1.-alpha_smear[1]);

    for(dir=XUP;dir<=TUP;dir++){
	for(dir2=XUP;dir2<=TUP;dir2++) if(dir2!=dir){

	    /* compute the staple */
	    FORALLDYNLINKS(i,st,dir) clear_su3mat(Staple2[dir2][dir]+i);

	    for(dir3=0;dir3<4;dir3++) if(dir3!=dir && dir3!=dir2){
	        for(dir4=XUP;dir4<=TUP;dir4++){
                    if(dir4!=dir && dir4!=dir2 && dir4 !=dir3) break;
                }
#if (SMEAR_LEVEL==3)
                staple_nhyp(dir,dir3, hyplink1[dir4][dir],
                            hyplink1[dir4][dir3],Staple2[dir2][dir]);
#else /* SMEAR_LEVEL==2 */
                staple_nhyp(dir,dir3, gauge_field_thin[dir],
                            gauge_field_thin[dir3],Staple2[dir2][dir]);
#endif
	    }

	    FORALLDYNLINKS(i,st,dir){

		/* make Omega  */
		scalar_mult_add_su3_matrix(gauge_field_thin[dir]+i,
                                             Staple2[dir2][dir]+i,ftmp1, &Q);
		scalar_mult_su3_matrix(&Q,ftmp2,&Omega);
		Staple2[dir2][dir][i]=Omega;

		mult_su3_an(&Omega,&Omega,&Q);
                scalar_add_diag_su3(&Q,IR_STAB);
#ifndef NHYP_DEBUG
	        compute_fhb(&Q,f,NULL, 0);
#else
                compute_fhb(&Omega,&Q,f,NULL, 0);
#endif

		/* make Q**2 */
		mult_su3_nn(&Q,&Q,&Q2);

		/* compute Q^(-1/2) via Eq. 19  */
		scalar_mult_su3_matrix(&Q,f[1],&tmat);
		scalar_mult_add_su3_matrix(&tmat,&Q2,f[2],&eQ);
		scalar_add_diag_su3(&eQ,f[0]);

		/* multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)  */
		mult_su3_nn(&Omega,&eQ,hyplink2[dir2][dir]+i);

	    }

	} /* dir2 */
    } /* dir */
} /* block_nhyp2 */
#endif /* SMEAR_LEVEL>1 */

#if (SMEAR_LEVEL==3)
void block_nhyp1()
{
    register int dir1, dir2, i;
    register site *st;
    Real f[3], ftmp1, ftmp2;
    su3_matrix tmat, Omega, eQ,  Q, Q2;

    ftmp1=alpha_smear[2]/(2.*(1.-alpha_smear[2]));
    ftmp2=(1.-alpha_smear[2]);

    /* dir1 is the direction of the original link
       dir2 is the other direction that defines the staple        */

    for(dir1=XUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<=TUP;dir2++) if(dir1!=dir2){
	    FORALLDYNLINKS(i,st,dir1) clear_su3mat(Staple1[dir2][dir1]+i);

            /* compute the staple */
	    staple_nhyp(dir1,dir2,gauge_field_thin[dir1],
                        gauge_field_thin[dir2],Staple1[dir2][dir1]);

            FORALLDYNLINKS(i,st,dir1){
                /* make Omega  */
                scalar_mult_add_su3_matrix(gauge_field_thin[dir1]+i,
                                       Staple1[dir2][dir1]+i,ftmp1 ,&Q);
                scalar_mult_su3_matrix(&Q,ftmp2,&Omega);
                Staple1[dir2][dir1][i]=Omega;

                mult_su3_an(&Omega,&Omega,&Q);
                scalar_add_diag_su3(&Q,IR_STAB);
#ifndef NHYP_DEBUG
	        compute_fhb(&Q,f,NULL, 0);
#else
                compute_fhb(&Omega,&Q,f,NULL, 0);
#endif

                /* make Q**2 */
                mult_su3_nn(&Q,&Q,&Q2);

                /* compute Q^(-1/2) via Eq. 19  */
                scalar_mult_su3_matrix(&Q,f[1],&tmat);
                scalar_mult_add_su3_matrix(&tmat,&Q2,f[2],&eQ);
                scalar_add_diag_su3(&eQ,f[0]);

                /* multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)  */
                mult_su3_nn(&Omega,&eQ,hyplink1[dir2][dir1]+i);
            }
	} /* dir2 */
    } /* dir1 */
} /* block_nhyp1 */
#endif /* SMEAR_LEVEL=3 */

void staple_nhyp( int dir1, int dir2, su3_matrix *lnk1, su3_matrix *lnk2,
                  su3_matrix *stp )
{
    register int i;
    register site *st;
    msg_tag *tag0,*tag1,*tag2;
    su3_matrix tmat1,tmat2;

    /* dir1 is the direction of the original link
       dir2 is the other direction that defines the staple        */

    /* get blocked_link[dir2] from direction dir1 */
    tag0 = start_gather_field( lnk2, sizeof(su3_matrix), dir1,
                               EVENANDODD, gen_pt[0] );

    /* get blocked_link[dir1] from direction dir2 */
    tag1 = start_gather_field( lnk1, sizeof(su3_matrix), dir2,
                               EVENANDODD, gen_pt[1] );

    /* start working on the lower staple while we wait for the gathers.
       the lower staple is prepared at x-dir2 and stored in tempmat_nhyp1,
       then gathered to x.
    */

    FORALLSITES(i,st){
	mult_su3_an( lnk2+i, lnk1+i, tempmat_nhyp1+i );
    }

    wait_gather(tag0);
    wait_gather(tag1);

    /* fix boundary values for SF */
    FORALLSITES(i,st) {
        gen_pt[0][i] = CHOOSE_NBR(i,st,dir1,linkf_bndr_up[dir2],0);
        gen_pt[1][i] = CHOOSE_NBR(i,st,dir2,linkf_bndr_up[dir1],1);
    }

    /* finish lower staple */
    FORALLSITES(i,st){
	mult_su3_nn( tempmat_nhyp1+i, (su3_matrix *)gen_pt[0][i], &tmat1 );
        su3mat_copy( &tmat1, tempmat_nhyp1+i );
    }

    /* gather staple from direction -dir2 to "home" site */
    tag2 = start_gather_field( tempmat_nhyp1, sizeof(su3_matrix),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[2] );

    /* calculate upper staple, add it */
    FORALLDYNLINKS(i,st,dir1){
	mult_su3_nn( lnk2+i, (su3_matrix *)gen_pt[1][i], &tmat1 );
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
	add_su3_matrix( stp+i, &tmat2, stp+i );
    }

    /* finally add the lower staple.
       for SF, content of tempmat_nhyp1 is bogus if t=0, dir1=0,1,2,
       but we don't use it.
    */

    wait_gather(tag2);

    FORALLDYNLINKS(i,st,dir1){
	add_su3_matrix( stp+i, (su3_matrix *)gen_pt[2][i], stp+i );
    }

    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
}
void scalar_add_diag_su3(su3_matrix *a, Real s){
  register int i;

  for(i=0;i<3;i++){
    a->e[i][i].real += s;
  }
}
