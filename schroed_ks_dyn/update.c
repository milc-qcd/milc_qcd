/********** update.c ****************************************************/
/* MIMD version 7 */
/* NOTE: NEEDS UPGRADING TO ASQTAD */
/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "schroed_ks_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>    /* For "finite" */
#endif

#ifdef PHI_ALGORITHM
static void predict_next_xxx(Real *oldtime, Real *newtime, Real *nexttime);
#endif

void dump_v_site(char *id, field_offset f){
  int i,j;
  site *s;
  su3_vector *v;

  printf("Dump of %s\n",id);
  FORALLSITES(i,s){
    printf("%d %d %d %d ", s->x, s->y, s->z, s->t);
    v = (su3_vector *)F_PT(s,f);
    for(j = 0; j < 3; j++)
      printf("( %e, %e) ", v->c[j].real, v->c[j].imag);
    printf("\n");
  }
}

int update()  {
int step, iters=0;
Real final_rsq;
Real cg_time;	/* simulation time for last two CG's */
#ifdef PHI_ALGORITHM
Real old_cg_time,next_cg_time;	/* simulation time for last two CG's */
#endif
#ifdef HMC_ALGORITHM
double startaction,endaction,change;
Real xrandom;
#endif
  imp_ferm_links_t** fn;

    /* refresh the momenta */
    ranmom();

    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){

#ifdef PHI_ALGORITHM
	/* generate a pseudofermion configuration only at start*/
	if(step==1){
	  restore_fermion_links_from_site(fn_links, PRECISION);
	  fn = get_fm_links(fn_links);
	  clear_latvec( F_OFFSET(phi), EVENANDODD );
	  grsource_imp( F_OFFSET(phi), mass, EVEN, fn[0]);
	  clear_latvec( F_OFFSET(xxx), EVENANDODD );
	  clear_latvec( F_OFFSET(ttt), EVENANDODD );
	  //	  grsource(EVEN); 
	  old_cg_time = cg_time = -1.0e6;
	}

#ifdef HMC_ALGORITHM
	/* find action */
	/* do conjugate gradient to get (Madj M)inverse * phi */
	if(step==1){
	    /* do conjugate gradient to get (Madj M)inverse * phi */
	  restore_fermion_links_from_site(fn_links, PRECISION);
	  fn = get_fm_links(fn_links);
	  iters += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			      niter, nrestart, rsqmin, PRECISION, EVEN,
			      &final_rsq, fn[0]);
	    cg_time = 0.0;
	    startaction=d_action();
	    /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
}
#endif

	/* update U's to middle of interval */
	update_u(0.5*epsilon);

	/* save conjugate gradient solution, predict next one */
	next_cg_time = ((Real)step-0.5)*epsilon;
	predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);

#else /* "R" algorithm */
	/* first update the U's to special time interval */
	update_u(epsilon*(0.5-nflavors/8.0));

	/* generate a pseudofermion configuration */
	restore_fermion_links_from_site(fn_links, PRECISION);
	fn = get_fm_links(fn_links);
	clear_latvec( F_OFFSET(phi), EVENANDODD );
     	grsource_imp( F_OFFSET(phi), mass, EVEN, fn[0]);
	clear_latvec( F_OFFSET(xxx), EVENANDODD );
	clear_latvec( F_OFFSET(ttt), EVENANDODD );
	//	grsource(EVEN); 
	cg_time = -1.0e6;

	/* update U's to middle of interval */
	update_u(epsilon*nflavors/8.0);
#endif

	/* do conjugate gradient to get (Madj M)inverse * phi */
	restore_fermion_links_from_site(fn_links, PRECISION);
	fn = get_fm_links(fn_links);
	iters += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, nrestart, rsqmin, PRECISION, EVEN, 
			    &final_rsq, fn[0]);
	cg_time = ((Real)step - 0.5)*epsilon;


	/* now update H by full time interval */
	update_h(epsilon);

	/* update U's by half time step to get to even time */
	update_u(epsilon*0.5);

	/* reunitarize the gauge field */
	rephase_sf( OFF );
	reunitarize();
	rephase_sf( ON );

    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    next_cg_time = steps*epsilon;
    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
    restore_fermion_links_from_site(fn_links, PRECISION);
    fn = get_fm_links(fn_links);
    iters += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			niter, nrestart, rsqmin, PRECISION, EVEN, 
			&final_rsq, fn[0]);
    cg_time = steps*epsilon;
    endaction=d_action();
    change = endaction-startaction;

#ifndef HAVE_IEEEFP_H
    /* Floating overflow gives +- NaN for change - need to reject if -NaN */
    if(change < -1.0e20 ){
	if(this_node==0)printf(
	    "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
	change = 1.0e20;
    }
#else
    /* Carleton's change for detecting NaN */
    /* Reject configurations giving overflows */
    if(!finite((double)change)){
	if(this_node==0)printf(
	    "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
	change = 1.0e20;
    }
#endif

    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( -change ) < (double)xrandom ){
	if(steps > 0)
	    gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
	invalidate_fermion_links(fn_links);
	if(this_node==0)printf("REJECT: delta S = %e\n", change);
    }
    else {
	if(this_node==0)printf("ACCEPT: delta S = %e\n", change);
    }
#endif

    if(steps > 0)return (iters/steps);
    else return(-99);
}

#ifdef PHI_ALGORITHM
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
static void predict_next_xxx(Real *oldtime, Real *newtime, Real *nexttime) {
register int i;
register site *s;
register Real x;
su3_vector tvec;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
        FOREVENSITES(i,s) if(s->t > 0){
	    s->old_xxx = s->xxx;
        }
    }
    else  {
        FOREVENSITES(i,s) if(s->t > 0){
            sub_su3_vector( &(s->xxx), &(s->old_xxx), &tvec);
	    s->old_xxx = s->xxx;
            scalar_mult_add_su3_vector( &(s->xxx), &tvec,x, &(s->xxx) );
        }
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}
#endif
