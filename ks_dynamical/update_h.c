/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 6 */

#include "ks_dyn_includes.h"
void gauge_force( Real eps );
void fermion_force( Real eps );

void update_h( Real eps ) {
    /* gauge field force */
    gauge_force(eps);
    /* fermionic force */
    /* First compute M*xxx in temporary vector ttt */
	/* The diagonal term in M doesn't matter */
    load_ferm_links(&fn_links);
    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links );
    fermion_force(eps);
} /* update_h */

/* update the momenta with the gauge force */
void gauge_force( Real eps ) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
int start;
su3_matrix tmat1,tmat2;
register Real eb3;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
gf_av=gf_max=0.0;
**END TEMP**/

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
		    mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->staple) );
		}
		start=0;
	    }
	    else{
	        FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
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
	    scalar_mult_add_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );
/** FIND AVERAGE AND MAXIMUM GAUGE FORCE **/
/**TEMP **
for(gf_x=0,gf_i=0;gf_i<3;gf_i++)for(gf_j=0;gf_j<3;gf_j++)
gf_x+=cabs_sq(&(tmat1.e[gf_i][gf_j]));
gf_x *= beta/3.0;
gf_av += gf_x;
if(gf_x > gf_max) gf_max = gf_x;
** END TEMP **/
	}
    }
/**TEMP**
g_floatsum( &gf_av );
g_floatmax( &gf_max );
gf_av /= (4*volume);
if(this_node==0)printf("GF: %e %e\n",gf_av,gf_max);
**END TEMP **/
/**dtime += dclock();
if(this_node==0)printf("G_FORCE: time = %e mflops = %e\n",
dtime, (double)(10848.0*volume/(1.0e6*dtime*numnodes())) );**/
}

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   xxx, and dslash_site(xxx,ttt) has been run. */
void fermion_force( Real eps ) {
register int i,dir;
register site *st;
msg_tag *tag0,*tag1;
su3_vector tvec;
su3_matrix temp1,temp2,temp3;
Real ferm_epsilon;
/**TEMP**
Real ff_x,ff_av,ff_max;
int ff_i,ff_j;
ff_av=ff_max=0.0;
**END TEMP**/

/**double dtime,dclock();
dtime = -dclock();**/

    ferm_epsilon = (nflavors/2.0)*eps;
    /* For even sites, gather ttt  get first one befor entering loop */
    tag0 = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector), XUP, EVEN,
	gen_pt[0] );

    for(dir=XUP;dir<=TUP;dir++){

	/* For odd sites, gather xxx */
	tag1 = start_gather_site( F_OFFSET(xxx), sizeof(su3_vector), dir, ODD,
	    gen_pt[1] );
	wait_gather(tag0);
	FOREVENSITES(i,st){
	    mult_su3_mat_vec( &(st->link[dir]), (su3_vector *)gen_pt[0][i],
		&tvec);
	    su3_projector( &tvec, &(st->xxx), &temp1 );
	    uncompress_anti_hermitian( &(st->mom[dir]), &temp2 );
	    scalar_mult_add_su3_matrix( &temp2, &temp1, ferm_epsilon, &temp3 );
	    make_anti_hermitian( &temp3, &(st->mom[dir]) );
/** FIND AVERAGE AND MAXIMUM FERMION FORCE **/
/**TEMP **
for(ff_x=0,ff_i=0;ff_i<3;ff_i++)for(ff_j=0;ff_j<3;ff_j++)
ff_x+=cabs_sq(&(temp1.e[ff_i][ff_j]));
ff_x *= (nflavors/2.0);
ff_av += ff_x;
if(ff_x > ff_max) ff_max = ff_x;
** END TEMP **/
	}
	cleanup_gather(tag0);

	/* For even sites, gather ttt */
	if(dir<TUP){
	    tag0 = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector),
		dir+1, EVEN, gen_pt[0] );
	}

	wait_gather(tag1);
	FORODDSITES(i,st){
	    mult_su3_mat_vec( &(st->link[dir]), (su3_vector *)gen_pt[1][i],
		&tvec);
	    su3_projector( &(st->ttt), &tvec, &temp1 );
	    uncompress_anti_hermitian( &(st->mom[dir]), &temp2 );
	    scalar_mult_add_su3_matrix( &temp2, &temp1, ferm_epsilon, &temp3 );
	    make_anti_hermitian( &temp3, &(st->mom[dir]) );
/** FIND AVERAGE AND MAXIMUM FERMION FORCE **/
/**TEMP **
for(ff_x=0,ff_i=0;ff_i<3;ff_i++)for(ff_j=0;ff_j<3;ff_j++)
ff_x+=cabs_sq(&(temp1.e[ff_i][ff_j]));
ff_x *= (nflavors/2.0);
ff_av += ff_x;
if(ff_x > ff_max) ff_max = ff_x;
** END TEMP **/
	}
	cleanup_gather(tag1);
    }
/**TEMP**
g_floatsum( &ff_av );
g_floatmax( &ff_max );
ff_av /= (4*volume);
if(this_node==0)printf("FF: %e %e\n",ff_av,ff_max);
**END TEMP **/
/**dtime += dclock();
if(this_node==0)printf("F_FORCE: time = %e mflops = %e\n",
dtime, (double)(672.0*volume/(1.0e6*dtime*numnodes())) );**/
}
