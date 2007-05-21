/******************* ape_smear_SF.c ****************************************/

/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* For Schroedinger functional - special b.c. */
/* version of 8/24/99 by UMH */

/* Perform one iteration of "APE smearing" on all links.
   For the projection back to SU(3) use multiples of 3 hits over the
   3 different SU(2) subgroups until the procedure converges.
   Note: only this complete projection makes the smearing truly
   gauge covariant. */

#include "schroed_cl_includes.h"
#define TOL 1e-5
#define MAXCOUNT 50

void ape_smear_SF() {
register int i,dir,diro;
register site *s;
msg_tag *tag0,*tag1,*tag2,*tag4;
su3_matrix tmat1,tmat2;
int disp[4];	/* displacement vector for general gather */
int index,ind1,ind2;
Real a0,a1,a2,a3,asq;
Real smear_fac, alpha6;
su2_matrix u;
void left_su2_hit_n();
int count;
double old_tr, new_tr;
Real conver;
field_offset tmpmat;

    tmpmat = F_OFFSET(psi);

    smear_fac = 1.0 - alpha;
    alpha6 = alpha / 6.0;

    /* Loop over the directions and APE-smear the links.
       The results will temporarily be stored in link_tmp
       and to link[dir] at the end. */

    for(dir=XUP;dir<=TUP;dir++){

	/* Direct link with weight smear_fac */
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix( &(s->link[dir]), smear_fac,
		&(s->link_tmp[dir]));
	}

	/* Loop over other directions, computing contribution from
	   plaquettes in the (dir,diro)-plane */
	for(diro=XUP;diro<=TUP;diro++) if(diro != dir){

	    /* get link[diro] from displacement +dir-diro */
	    for(i=XUP;i<=TUP;i++)disp[i]=0;
	    disp[dir] = 1;
	    disp[diro] = -1;
	    tag4 = start_general_gather_site( F_OFFSET(link[diro]),
		sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

	    /* get link[diro] from direction dir */
	    tag0 = start_gather_site( F_OFFSET(link[diro]), sizeof(su3_matrix),
		dir, EVENANDODD, gen_pt[0] );

	    /* get link[dir] from direction diro */
	    tag1 = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
		diro, EVENANDODD, gen_pt[1] );

	    /* Make one corner of lower staple and gather "up" */
	    FORALLSITES(i,s){
		if(dir==TUP == s->t==(nt-1)){
		    mult_su3_an( &(s->link[diro]), &(s->link[dir]), &tmat1 );
		    mult_su3_nn( &tmat1, &(s->boundary[diro]),
			((su3_matrix *)F_PT(s,tmpmat)) );
		}
		else{
		    mult_su3_an( &(s->link[diro]), &(s->link[dir]),
			((su3_matrix *)F_PT(s,tmpmat)) );
		}
	    }
	    tag2 = start_gather_site( tmpmat, sizeof(su3_matrix),
		OPP_DIR(diro), EVENANDODD, gen_pt[2] );

	    /* Make upper staple */
	    wait_gather(tag0);
	    wait_gather(tag1);
	    FORALLSITES(i,s){
		if(diro==TUP && s->t==(nt-1)){
		    mult_su3_nn( &(s->link[diro]), &(s->boundary[dir]),
			&tmat1 );
		}
		else{
		    mult_su3_nn( &(s->link[diro]), (su3_matrix *)(gen_pt[1][i]),
			&tmat1 );
		}
		if(dir==TUP && s->t==(nt-1)){
		    mult_su3_na( &tmat1, &(s->boundary[diro]), &tmat2 );
		}
		else{
		    mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[0][i]), &tmat2 );
		}
		scalar_mult_add_su3_matrix( &(s->link_tmp[dir]), &tmat2,
		    alpha6, &(s->link_tmp[dir]));
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);

	    /* Finish making lower staple */
	    wait_gather(tag2);
	    wait_general_gather(tag4);
	    FORALLSITES(i,s){
		if(dir==TUP && s->t==(nt-1)){
		    scalar_mult_add_su3_matrix( &(s->link_tmp[dir]),
			(su3_matrix *)(gen_pt[2][i]),
			alpha6, &(s->link_tmp[dir]));
		}
		else{
		    mult_su3_nn( (su3_matrix *)(gen_pt[2][i]),
			(su3_matrix *)(gen_pt[4][i]), &tmat1 );
		    scalar_mult_add_su3_matrix( &(s->link_tmp[dir]), &tmat1,
			alpha6, &(s->link_tmp[dir]));
		}
	    }
	    cleanup_gather(tag2);
	    cleanup_general_gather(tag4);
	} /* end loop over diro */

	/* Now do hits in the SU(2) subgroup to "normalize" link_tmp[dir] */
	/* Start with the "old" link and put result into tmpmat */
	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[dir]), ((su3_matrix *)F_PT(s,tmpmat)) );
	}

	old_tr = (double)0.0;
	FORALLSITES(i,s){
	    old_tr += (double)realtrace_su3(((su3_matrix *)F_PT(s,tmpmat)),
		&(s->link_tmp[dir]) );
	}
	g_doublesum( &old_tr );
	old_tr /= (double)(3*volume);
	conver = 1.0;

	for(count=0; count<MAXCOUNT && conver>TOL; count++){

	    for(index=0;index<3;index++){

		/* Pick the SU(2) subgroup */
		ind2 = (index+1) % 3;
		if( ind2 > index ){
		    ind1 = index;
		}
		else{
		    ind1 = ind2;
		    ind2 = index;
		}

		/* FORALLSITES(i,s){ */
		FORALLSITES(i,s) if(dir==TUP || s->t>0){
	            mult_su3_na( ((su3_matrix *)F_PT(s,tmpmat)),
			&(s->link_tmp[dir]), &tmat1 );

		    /* Extract SU(2) subgroup in Pauli matrix representation,
		       a0 + i * sum_j a_j sigma_j, from the SU(3) matrix tmat1 */
		    a0 = tmat1.e[ind1][ind1].real + tmat1.e[ind2][ind2].real;
		    a1 = tmat1.e[ind1][ind2].imag + tmat1.e[ind2][ind1].imag;
		    a2 = tmat1.e[ind1][ind2].real - tmat1.e[ind2][ind1].real;
		    a3 = tmat1.e[ind1][ind1].imag - tmat1.e[ind2][ind2].imag;

		    /* Normalize and put complex conjugate into u */
		    asq = a0*a0 + a1*a1 + a2*a2 + a3*a3;
		    asq = sqrt((double)asq);
		    a0 = a0/asq; a1 = a1/asq; a2 = a2/asq; a3 = a3/asq;
		    u.e[0][0] = cmplx( a0,-a3);
		    u.e[0][1] = cmplx(-a2,-a1);
		    u.e[1][0] = cmplx( a2,-a1);
		    u.e[1][1] = cmplx( a0, a3);

		    /* Do the SU(2) hit */
		    left_su2_hit_n( &u, ind1, ind2,
			((su3_matrix *)F_PT(s,tmpmat)) );

		} /* end loop over sites */

	    } /* end loop over index */

	    new_tr = (double)0.0;
	    FORALLSITES(i,s){
		new_tr += (double)realtrace_su3(((su3_matrix *)F_PT(s,tmpmat)),
		    &(s->link_tmp[dir]) );
	    }
	    new_tr /= (double)(3*volume);
	    g_doublesum( &new_tr );
	    conver = (Real)((new_tr-old_tr)/old_tr);
	    old_tr = new_tr;

	} /* end loop over projection convergence */

	if( conver>TOL && this_node==0 ){
	    printf("No convergence in SMEAR for dir =%d: conver = %e\n",
		dir, (double)conver);
	}

	/* And temporarily store the new link */
	/* FORALLSITES(i,s){ */
	FORALLSITES(i,s) if(dir==TUP || s->t>0){
	    su3mat_copy( ((su3_matrix *)F_PT(s,tmpmat)), &(s->link_tmp[dir]) );
	}

    } /* end loop over dir */

    /* Now copy the smeared links, overwriting the previous ones */
    for(dir=XUP;dir<=TUP;dir++){
	/* FORALLSITES(i,s){ */
	FORALLSITES(i,s) if(dir==TUP || s->t>0){
	    su3mat_copy( &(s->link_tmp[dir]), &(s->link[dir]) );
	}
    }

} /* smearing */
