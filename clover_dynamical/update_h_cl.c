/****** update_h_cl.c  -- update the momentum matrices ******************/
/* MIMD version 7.  April 1995 MBW */
/* Clover fermions */

#include "cl_dyn_includes.h"
void gauge_force(Real eps);
void fermion_force(Real eps);

void update_h(Real eps){
    /* gauge field force */
    gauge_force(eps);
    /* fermionic force */
    fermion_force(eps);
} /* update_h */

/* update the momenta with the gauge force */
void gauge_force(Real eps){
register int i,dir;
register site *st;
su3_matrix tmat1,tmat2;
register Real eb3;

int j,k;
int dirs[10],sign[10],length;
int path_dir[10],path_sign[10],path_length;

msg_tag *mtag0;
int ln,iloop;

int ncount;

    eb3 = eps*beta/3.0;


    /* Loop over directions, update mom[dir] */
    for(dir=XUP; dir<=TUP; dir++){

	FORALLSITES(i,st)for(j=0;j<3;j++)for(k=0;k<3;k++){
			st->staple.e[j][k]=cmplx(0.0,0.0);
	}

ncount=0;
   for(iloop=0;iloop<nloop;iloop++){
      length=loop_length[iloop];
      /* Note: changed loop_ch[iloop][ln] =1 to 
	 loop_ch[iloop][ln] ==1 CD 6/6/98 */
      for(ln=0;ln<loop_num[iloop];ln++)if(loop_ch[iloop][ln] ==1){
 /* set up dirs and sign  */
       for(k=0;k<length;k++){
                if(loop_table[iloop][ln][k] < 4 ){sign[k]=1;
                       dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;}
                else {sign[k]=-1;
                       dirs[k]=(7+dir-loop_table[iloop][ln][k] )% 4;}
                             }

	path_length= length-1;
	for(k=0;k<length;k++)if(dirs[k]==dir) {

		if(sign[k]== -1)
		for(j=0;j<path_length;j++) {
		path_dir[j] = dirs[(k+j+1)%length];
		path_sign[j] = sign[(k+j+1)%length];
		}

		if(sign[k]== 1)
		for(j=0;j<path_length;j++) {
		path_dir[path_length-1-j] = dirs[(k+j+1)%length];
		path_sign[path_length-1-j]= -sign[(k+j+1)%length];
		}
	path(path_dir,path_sign,path_length);

/* now we must gather the path from direction dir, due to 
	the way path.c works*/
	mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
                dir, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);

	FORALLSITES(i,st){
              su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(st->tempmat2) ); 
	}
	cleanup_gather(mtag0);


	single_action(dir,loop_term[ncount]);
	ncount++;

/*
printf("ncount=%d\n",ncount);
*/

	} /* k (location in path) */
	}} /* ln and iloop */



	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	    mult_su3_na( &(st->link[dir]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir]), &tmat2 );
	    scalar_mult_sub_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir]) );
	}
    }
} /* gauge_force.c */

void single_action(int dir,Real *coeff)
{
site *st;
Real action,act2,new_term;
int i,j;


	  FORALLSITES(i,st){
/* first we compute the fundamental term */

	new_term= coeff[0];


/* now we add in the higher representations */
	if(nreps > 1){

		act2=1.0;

		action = 3.0 - realtrace_su3(&(st->link[dir]),&(st->tempmat2)); 


			for(j=1;j<nreps;j++){
			act2 *= action;
			new_term += coeff[j]*act2*(Real)(j+1);
			}
	    } 

	scalar_mult_add_su3_matrix( &(st->staple), &(st->tempmat2),
		new_term, &(st->staple) );

	} /* sites */

}

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in psi */
/* pseudofermion vector is chi */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
 
  without LU:
   	M = A - kappa Dslash
  with LU:
  	M = A_e - kappa^2 * Dslash_eo (A_o)^{-1}* Dslash_oe

*/

void fermion_force(Real eps){
register int i, mu, nu;
register site *st;
msg_tag *tag0,*tag1;
wilson_vector tvec1,tvec2;
half_wilson_vector hvec;
su3_matrix temp1,temp2;
Real ferm_epsilon, KAP, CKU0=-clov_c*kappa/(8.0*u0*u0*u0);

/* printf("CKU0 = %g\n",CKU0); */
/**double dtime,dclock();
dtime = -dclock();**/

ferm_epsilon = nflavors*eps;

#ifdef LU
    KAP = -kappa*kappa;
    dslash_w_site( F_OFFSET(psi), F_OFFSET(tmp), PLUS, ODD );
    mult_this_ldu_site( gen_clov, F_OFFSET(tmp), F_OFFSET(psi), ODD );
    dslash_w_site( F_OFFSET(psi), F_OFFSET(p), PLUS, EVEN );
    mult_this_ldu_site( gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVEN );
    FOREVENSITES(i,st)
        scalar_mult_add_wvec( &(st->tmp), &(st->p), KAP, &(st->p) );
    dslash_w_site( F_OFFSET(p), F_OFFSET(tmp), MINUS, ODD );
    cleanup_dslash_wtemps();
    mult_this_ldu_site( gen_clov, F_OFFSET(tmp), F_OFFSET(p), ODD );

    /* M = A_even - kappa^2 * Dslash * A_odd^{-1} * Dslash
       psi(even) = psi(even)
       psi(odd)  = A_odd^{-1} * Dslash * psi(even) 
         p(even) = M * psi(even)
         p(odd)  = A_odd^{-1} * Dslash_adjoint * M * psi(even). */
#else
    KAP = -kappa;
    dslash_w_site( F_OFFSET(psi), F_OFFSET(p), PLUS, EVENANDODD );
    cleanup_dslash_wtemps();
    mult_this_ldu_site( gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVENANDODD );
    FORALLSITES(i,st) {
        scalar_mult_add_wvec( &(st->tmp), &(st->p), KAP, &(st->p) );
}
    /* M = A - kappa * Dslash
       psi = psi
       p = M * psi          */
#endif /*LU*/

    for(mu=XUP;mu<=TUP;mu++) {
/* printf("mu = %d\n",mu); */

        FORALLSITES(i,st){
	    wp_shrink( &(st->psi), &(st->htmp[0]), mu, PLUS);
	    wp_shrink( &(st->p), &(st->htmp[1]), mu, MINUS);
	}
	tag0 = start_gather_site( F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
	     mu, EVENANDODD, gen_pt[0] );
	tag1 = start_gather_site( F_OFFSET(htmp[1]), sizeof(half_wilson_vector),
	     mu, EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);

	/* First the U d(dslash)/dU terms (do for only one value of nu) */
	FORALLSITES(i,st) {
	    /* psi and p parallel transported in from positive directions */
	    mult_su3_mat_hwvec( &(st->link[mu]), 
				(half_wilson_vector *)gen_pt[0][i], &hvec);
	    wp_grow( &hvec, &tvec1, mu, PLUS);
	    /* i even => tvec1 = (1+gamma_mu)*U*Aodd^(-1)*D*psi,
	   i odd  => tvec1 = (1+gamma_mu)*U*psi */

	    mult_su3_mat_hwvec( &(st->link[mu]), 
				(half_wilson_vector *)gen_pt[1][i], &hvec);
	    wp_grow( &hvec, &tvec2, mu, MINUS);
	    /* i even => tvec2 = (1-gamma_mu)*U*Aodd^(-1)*D_adj*M*psi,
	       i odd  => tvec2 = (1-gamma_mu)*U*M*psi */

	    su3_projector_w( &tvec1, &(st->p), &temp1 );
	    su3_projector_w( &tvec2, &(st->psi), &temp2 );
	    add_su3_matrix( &temp1, &temp2, &(st->tempmat1) );
	    scalar_mult_su3_matrix( &(st->tempmat1), KAP, &(st->tempmat1) );
	} /* end loop over sites */
	cleanup_gather(tag0);
	cleanup_gather(tag1);

	/* Now the U dA/dU terms */
        for(nu=XUP;nu<=TUP;nu++) if(nu!=mu) {

#ifdef LU
	    /* U dA_odd/dU from U dM/dU */
	    udadu_mu_nu( F_OFFSET(p), F_OFFSET(psi), F_OFFSET(tempmat2),
		mu, nu, ODD );
	    FORALLSITES(i,st) {
		scalar_mult_add_su3_matrix( &(st->tempmat1), 
		    &(st->tempmat2), -KAP*CKU0, &(st->tempmat1) );
	    }
	    /* U dA_odd/dU from U dM^dagger/dU */
	    udadu_mu_nu( F_OFFSET(psi), F_OFFSET(p), F_OFFSET(tempmat2),
		mu, nu, ODD );
	    FORALLSITES(i,st) {
		scalar_mult_add_su3_matrix( &(st->tempmat1),
		    &(st->tempmat2), -KAP*CKU0, &(st->tempmat1) );
	    }
            /* U dA_even/dU from U dM/dU */
	    udadu_mu_nu( F_OFFSET(p), F_OFFSET(psi), F_OFFSET(tempmat2),
		mu, nu, EVEN );
	    FORALLSITES(i,st) {
                scalar_mult_add_su3_matrix( &(st->tempmat1), 
		    &(st->tempmat2), CKU0, &(st->tempmat1) );
	    }
	    /* U dA_even/dU from U dM^dagger/dU */
	    udadu_mu_nu( F_OFFSET(psi), F_OFFSET(p), F_OFFSET(tempmat2),
		mu, nu, EVEN );
	    FORALLSITES(i,st) {
		scalar_mult_add_su3_matrix( &(st->tempmat1),
		    &(st->tempmat2), CKU0, &(st->tempmat1) );
	    }
#else /*if not LU*/
	/* Now the U dA/dU terms */
	    /* U dA/dU from U dM/dU */
	    udadu_mu_nu( F_OFFSET(p), F_OFFSET(psi), F_OFFSET(tempmat2),
		mu, nu, EVENANDODD );
	    FORALLSITES(i,st) {
                scalar_mult_add_su3_matrix( &(st->tempmat1), 
		    &(st->tempmat2), CKU0, &(st->tempmat1) );
	    }
	    /* U dA/dU from U dM^dagger/dU */
	    udadu_mu_nu( F_OFFSET(psi), F_OFFSET(p), F_OFFSET(tempmat2),
		mu, nu, EVENANDODD );
	    FORALLSITES(i,st) {
		scalar_mult_add_su3_matrix( &(st->tempmat1),
		    &(st->tempmat2), CKU0, &(st->tempmat1) );
	    }
#endif /*LU/else*/

/* Add the Tr_{dirac} \sigma_mu_nu A_odd^{-1} U dA_odd/dU term */
#ifdef LU
	/* staple is used as a temporary su3_matrix */
        tr_sigma_this_ldu_mu_nu_site( gen_clov, F_OFFSET(staple), mu, nu );
	udadu_mat_mu_nu( F_OFFSET(staple), F_OFFSET(tempmat2), mu, nu );
	FORALLSITES(i,st) {
	    scalar_mult_add_su3_matrix( &(st->tempmat1),
		&(st->tempmat2), 2.0*CKU0, &(st->tempmat1) );
	}
#endif /*LU*/

        } /* end loop over nu & endif( nu != mu )*/

	FORALLSITES(i,st) {
	    uncompress_anti_hermitian( &(st->mom[mu]), &temp2 );
	    scalar_mult_sub_su3_matrix( &temp2,
		&(st->tempmat1), ferm_epsilon, &temp1 );
	    make_anti_hermitian( &temp1, &(st->mom[mu]) );
	}

    } /* end loop over mu */
/**dtime += dclock();
if(this_node==0)printf("F_FORCE: time = %e mflops = %e\n",
dtime, (double)(5584.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* end fermion_force */
