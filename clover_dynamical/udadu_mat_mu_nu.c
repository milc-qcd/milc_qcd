/******************* udadu_mat_mu_nu.c ****************************************/

/* MIMD version 7 */
/* version of May 1995 by MBW */
/* calculates udadu_mu_nu and multiplies by an su3_matrix
   for the Tr_{dirac} \sigma_mu_nu A_odd^{-1} U dA_odd/dU term of
   the fermion force.
*/

#include "cl_dyn_includes.h"

void udadu_mat_mu_nu( field_offset matsrc, field_offset matdest, 
		     int mu, int nu )
{
int disp[4];
register int i;
register site *s;
msg_tag *tag[8];
su3_matrix tmp1, utmp, uAcntr;		/* for middle steps */

/**********  First work on upper leaf  **********/
    /* get link[nu] from direction +mu */
    tag[0] = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        mu, EVENANDODD, gen_pt[0] );

    /* get link[mu] from direction +nu */
    tag[1] = start_gather_site( F_OFFSET(link[mu]), sizeof(su3_matrix),
        nu, EVENANDODD, gen_pt[1] );

 /*****  Upper leaf - same parity *****/
    /* get matsrc from +mu +nu */
    for(i=XUP;i<=TUP;i++) disp[i]=0;
    disp[mu] = 1;  disp[nu] = 1;
    tag[2] = start_general_gather_site( matsrc, sizeof(su3_matrix),
        disp, ODD, gen_pt[2] );

    wait_gather(tag[0]);
    wait_gather(tag[1]);

    /* utmp = plaquette */
    FORODDSITES(i,s) {
	mult_su3_nn( &(s->link[mu]), (su3_matrix *)(gen_pt[0][i]),
	    &utmp );
	mult_su3_na( &utmp, (su3_matrix *)(gen_pt[1][i]), &tmp1 );
	mult_su3_na( &tmp1, &(s->link[nu]), &utmp );

	mult_su3_nn( &utmp, ((su3_matrix *)F_PT(s,matsrc)),
		((su3_matrix *)F_PT(s,matdest)) );
    }

    wait_general_gather(tag[2]);

    /* utmp = upperleftcorner, uAcntr = lowerrightcorner*matsrc */
    FORODDSITES(i,s) {
	mult_su3_nn( &(s->link[mu]), (su3_matrix *)(gen_pt[0][i]), &tmp1 );
	mult_su3_nn( &tmp1, (su3_matrix *)(gen_pt[2][i]), &uAcntr );

	mult_su3_nn( &(s->link[nu]), (su3_matrix *)(gen_pt[1][i]), &utmp );

	mult_su3_na( &uAcntr, &utmp, &tmp1);
	add_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

 /***** upper leaf - different parity *****/

    /* get matsrc from +nu */
    tag[3] = start_gather_site( matsrc, sizeof(su3_matrix), nu, EVEN, gen_pt[3] );
    wait_gather(tag[3]);

    /* uAcntr = staple*matsrc, utmp = link[nu]^dagger */
    FOREVENSITES(i,s) {
	mult_su3_nn( &(s->link[mu]), (su3_matrix *)(gen_pt[0][i]),
	    &tmp1 );
	mult_su3_na( &tmp1, (su3_matrix *)(gen_pt[1][i]), &utmp );
	mult_su3_nn( &utmp, (su3_matrix *)(gen_pt[3][i]), &uAcntr );
	mult_su3_na( &uAcntr,&(s->link[nu]),((su3_matrix *)F_PT(s,matdest)));
    }

    cleanup_gather(tag[3]);
    /* get matsrc from +mu */
    tag[3] = start_gather_site( matsrc, sizeof(su3_matrix), mu, EVEN, gen_pt[3] );
    wait_gather(tag[3]);

    /* uAcntr = link[mu]*matsrc, utmp = staple */
    FOREVENSITES(i,s) {
	mult_su3_nn( &(s->link[mu]), (su3_matrix *)(gen_pt[3][i]), &uAcntr );

	mult_su3_nn( &(s->link[nu]), (su3_matrix *)(gen_pt[1][i]), &tmp1 );
	mult_su3_na( &tmp1, (su3_matrix *)(gen_pt[0][i]), &utmp );

	mult_su3_na( &uAcntr, &utmp, &tmp1);
	add_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

/**********  The Lower leaf  **********/

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    /* get link[nu] from direction -nu */
    tag[0] = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        OPP_DIR(nu), EVENANDODD, gen_pt[0] );

    /* get link[mu] from direction -nu */
    tag[1] = start_gather_site( F_OFFSET(link[mu]), sizeof(su3_matrix),
        OPP_DIR(nu), EVENANDODD, gen_pt[1] );

    /* get link[nu] from direction -nu +mu */
    /* disp[mu] = 1; already */
    disp[nu] = -1;
    tag[4] = start_general_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        disp, EVENANDODD, gen_pt[4] );

    wait_gather(tag[0]);
    wait_gather(tag[1]);
    wait_general_gather(tag[4]);

 /***** lower leaf, different parity *****/

    /* uAcntr = staple*matsrc, utmp = link[mu]^dagger */
    FOREVENSITES(i,s) {
	mult_su3_an( (su3_matrix *)(gen_pt[0][i]),
	    (su3_matrix *)(gen_pt[1][i]), &uAcntr );
	mult_su3_nn( &uAcntr, (su3_matrix *)(gen_pt[4][i]), &tmp1 );
	mult_su3_nn( &tmp1, (su3_matrix *)(gen_pt[3][i]), &uAcntr );

	mult_su3_na( &uAcntr, &(s->link[mu]), &tmp1 );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

    cleanup_gather(tag[3]);
    /* get matsrc from -nu */
    tag[3] = start_gather_site( matsrc, sizeof(su3_matrix), OPP_DIR(nu),
	EVEN, gen_pt[3] );
    wait_gather(tag[3]);

    /* uAcntr = link[nu]^dagger*matsrc, utmp=staple */
    FOREVENSITES(i,s) {
	mult_su3_an( (su3_matrix *)(gen_pt[0][i]),
	    (su3_matrix *)(gen_pt[3][i]), &uAcntr );

	mult_su3_nn( (su3_matrix *)(gen_pt[1][i]),
	    (su3_matrix *)(gen_pt[4][i]), &tmp1 );
	mult_su3_na( &tmp1, &(s->link[mu]), &utmp );

	mult_su3_nn( &uAcntr, &utmp, &tmp1 );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

 /***** lower leaf, same parity *****/

    cleanup_general_gather(tag[2]);
    /* get matsrc from +mu -nu */
    /* disp[mu] = 1;  disp[nu] = -1; already */
    tag[2] = start_general_gather_site( matsrc, sizeof(su3_matrix), disp,
	ODD, gen_pt[2] );

    /* utmp= plaq */
    FORODDSITES(i,s) {
	mult_su3_an( (su3_matrix *)(gen_pt[0][i]),
	    (su3_matrix *)(gen_pt[1][i]), &utmp );
	mult_su3_nn( &utmp, (su3_matrix *)(gen_pt[4][i]), &tmp1 );
	mult_su3_na( &tmp1, &(s->link[mu]), &utmp );

	mult_su3_nn( ((su3_matrix *)F_PT(s,matsrc)), &utmp, &tmp1 );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

    wait_general_gather(tag[2]);

    /* uAcntr = lowerleftcorner*matsrc, utmp=upperrightcorner */
    FORODDSITES(i,s) {
	mult_su3_an( (su3_matrix *)(gen_pt[0][i]),
	    (su3_matrix *)(gen_pt[1][i]), &tmp1 );
	mult_su3_nn( &tmp1, (su3_matrix *)(gen_pt[2][i]), &uAcntr );

	mult_su3_na( &(s->link[mu]), (su3_matrix *)(gen_pt[4][i]), &utmp );

	mult_su3_na( &uAcntr, &utmp, &tmp1 );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,matdest)), &tmp1,
	    ((su3_matrix *)F_PT(s,matdest)) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    cleanup_general_gather(tag[2]);
    cleanup_gather(tag[3]);
    cleanup_general_gather(tag[4]);

} /* end udadu_mat_mu_nu */
