/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 6 */
/* Wilson fermions */

#include "wi_dyn_includes.h"
void gauge_force(Real eps);
void fermion_force(Real eps);

void update_h( Real eps ) {
    /* gauge field force */
    gauge_force( eps );
    /* fermionic force */
    fermion_force( eps );
} /* update_h */

/* update the momenta with the gauge force */
void gauge_force(Real eps) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
int start;
su3_matrix tmat1,tmat2;
register Real eb3;

/**double dtime,dclock();
dtime = -dclock();**/

    eb3 = eps*beta/3.0;
    /* Loop over directions, update mom[dir1] */
    for(dir1=XUP; dir1<=TUP; dir1++){
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	start=1; /* indicates staple sum not initialized */
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1){

	    /* get link[dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    wait_gather(tag0);
	    FORALLSITES(i,st){
	        mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1 );
	        mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
		    &(st->tempmat1) );
	    }

	    /* Gather this partial result "up to home site" */
	    tag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	    /* begin the computation of the "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
	    if(start){	/* this is the first contribution to staple */
	        FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i],
		        &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->staple) );
		}
		start=0;
	    }
	    else{
	        FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]), (su3_matrix *)gen_pt[2][i],
			&tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
		    add_su3_matrix( &(st->staple),&tmat2,&(st->staple));
	        }
	    }

	    wait_gather(tag1);
	    FORALLSITES(i,st){
		add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[1][i],
		    &(st->staple));
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}
	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	    mult_su3_na( &(st->link[dir1]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir1]), &tmat2 );
	    scalar_mult_sub_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );
	}
    }
/**dtime += dclock();
if(this_node==0)printf("G_FORCE: time = %e mflops = %e\n",
dtime, (double)(10848.0*volume/(1.0e6*dtime*numnodes())) );**/
}

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in psi */
/* pseudofermion vector is chi */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
 
  without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
  with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe

   With LU we will set
	psi(odd) = Dslash_oe * psi(even)
   and
	p(odd) = Dslash_adjoint * p(even)
   With this convention, the formula for the force on the even and
   odd sites is the same, and except for the "ferm_epsilon" is the
   same as the expression without the LU preconditioned matrix,
   so we have only one fermion force code.
*/
#ifdef LU
#define FORMYSITES FOREVENSITES
#else
#define FORMYSITES FORALLSITES
#endif

void fermion_force(Real eps) {
register int i,dir;
register site *st;
msg_tag *tag0,*tag1;
wilson_vector tvec1,tvec2;
half_wilson_vector hvec;
su3_matrix temp1,temp2;
Real ferm_epsilon;

/**double dtime,dclock();
dtime = -dclock();**/
 
    /* First compute M*psi in temporary vector p */
    /* set ferm_epsilon */
#ifdef LU
    ferm_epsilon = nflavors*eps*kappa*kappa;
    dslash_w_site( F_OFFSET(psi), F_OFFSET(psi), PLUS, ODD );
    dslash_w_site( F_OFFSET(psi), F_OFFSET(p)  , PLUS, EVEN);
    FOREVENSITES(i,st)
        scalar_mult_add_wvec( &(st->psi), &(st->p), -kappa*kappa, &(st->p) );

    /* psi(odd) = Dslash * psi(even) already, set
         p(odd) = Dslash_adjoint * p(even) . */
    dslash_w_site( F_OFFSET(p  ), F_OFFSET(p  ), MINUS, ODD );
#else
    ferm_epsilon = nflavors*eps*kappa;
    dslash_w_site( F_OFFSET(psi), F_OFFSET(p), PLUS, EVENANDODD );
    FORALLSITES(i,st)
	scalar_mult_add_wvec( &(st->psi), &(st->p), -kappa, &(st->p) );
#endif

    for(dir=XUP;dir<=TUP;dir++){

        FORALLSITES(i,st){
	    wp_shrink( &(st->psi), &(st->htmp[0]), dir, PLUS);
	    wp_shrink( &(st->p), &(st->htmp[1]), dir, MINUS);
	}
	tag0 = start_gather_site( F_OFFSET(htmp[0]), sizeof(half_wilson_vector), dir,
	    EVENANDODD, gen_pt[0] );
	tag1 = start_gather_site( F_OFFSET(htmp[1]), sizeof(half_wilson_vector), dir,
	    EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,st){
	    /* psi and p parallel transported in from positive directions */
	    mult_su3_mat_hwvec( &(st->link[dir]),
		(half_wilson_vector *)gen_pt[0][i], &hvec);
	    wp_grow( &hvec, &tvec1, dir, PLUS);
	    /* if LU is defined, tvec1 = 1+gamma_mu)*U*DSLASH*psi,
		otherwise it is (1+gamma_mu)*U*psi(x+mu_hat) */

	    mult_su3_mat_hwvec( &(st->link[dir]),
		(half_wilson_vector *)gen_pt[1][i], &hvec);
	    wp_grow( &hvec, &tvec2, dir, MINUS);
	    /* if LU is defined, tvec2 = (1-gamma_mu)*U*Dslash_adjoint*p,
		otherwise it is (1-gamma_mu)*U*p(x+mu_hat) */

	    su3_projector_w( &tvec1, &(st->p), &temp1 );
	    su3_projector_w( &tvec2, &(st->psi), &temp2 );
	    add_su3_matrix( &temp1, &temp2, &temp1);

	    uncompress_anti_hermitian( &(st->mom[dir]), &temp2 );
	    scalar_mult_add_su3_matrix( &temp2, &temp1, ferm_epsilon, &temp1 );
	    make_anti_hermitian( &temp1, &(st->mom[dir]) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
    }
/**dtime += dclock();
if(this_node==0)printf("F_FORCE: time = %e mflops = %e\n",
dtime, (double)(5584.0*volume/(1.0e6*dtime*numnodes())) );**/
}
