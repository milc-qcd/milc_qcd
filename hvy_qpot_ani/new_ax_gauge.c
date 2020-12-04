/************************** new_ax_gauge.c *******************************/
/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* version of 3/13/94 by UMH */
/* 2/19/98 Version 5 port CD */
/* 3/25/20 OMP support JHW*/

/* Gauge transform, in place, to axial gauge.
   Then store time-like gauge fields of last time slice in all
   other time slices as well, instead of the unit matrix! */

#include "hvy_qpot_includes.h"
#include "../include/openmp_defs.h"

void ax_gauge() {
register int i,dir,t,j,k;
register site *s;
su3_matrix tmat;
msg_tag *mtag[4];

    mtag[TUP] = start_gather_site( F_OFFSET(link[xc[TUP]]), sizeof(su3_matrix),
	xc[TDOWN], EVENANDODD, gen_pt[TUP] );

    /* Put gauge transformation into staple; it is unity for t=0 */
    FORALLSITES_OMP(i,s, private(j,k) ){ if( (site_coord(s,xc[TUP]))==0 ){
	for(j=0; j<3; j++)  {
	    for(k=0; k<3; k++)  {
		if(j != k)  {
		    s->staple.e[j][k] = cmplx(0.0,0.0);
		}
		else  {
		    s->staple.e[j][k].real = 1.0;
		    s->staple.e[j][k].imag = 0.0;
		}
	    }
	}
    } } END_LOOP_OMP;

    /* recursively determine the gauge transformations needed */
    for(t=1;t<nc[TUP];t++){

	wait_gather( mtag[TUP]);

	FORALLSITES_OMP(i,s, private(tmat) ){ if( (site_coord(s,xc[TUP]))==t ){
	    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
	    su3mat_copy( &(s->link[xc[TUP]]), &tmat);
	    mult_su3_nn( &(s->staple), &tmat, &(s->link[xc[TUP]]));
	} } END_LOOP_OMP;

	if(t<(nc[TUP]-1)){
	    restart_gather_site( F_OFFSET(link[xc[TUP]]), sizeof(su3_matrix),
		xc[TDOWN], EVENANDODD, gen_pt[TUP], mtag[TUP] );
	}
	else{
	    cleanup_gather( mtag[TUP]);
	    mtag[TUP] = start_gather_site( F_OFFSET(link[xc[TUP]]),
		sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[TUP] );
	}
    }

    /* Now do the gauge transformation of the space-like links */
    for(dir=XUP;dir<=ZUP;dir++){
	mtag[dir] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    xc[dir], EVENANDODD, gen_pt[dir] );
    }

    for(dir=XUP;dir<=ZUP;dir++){
	FORALLSITES_OMP(i,s, private(tmat) ){
	    su3mat_copy( &(s->link[xc[dir]]), &tmat);
	    mult_su3_nn( &(s->staple), &tmat, &(s->link[xc[dir]]));
	} END_LOOP_OMP;

	wait_gather( mtag[dir]);

	FORALLSITES_OMP(i,s, private(tmat) ){
	    su3mat_copy( &(s->link[xc[dir]]), &tmat);
	    mult_su3_na( &tmat, (su3_matrix *)(gen_pt[dir][i]),
		&(s->link[xc[dir]]));
	} END_LOOP_OMP;
    }

    for(dir=XUP;dir<=ZUP;dir++){
	cleanup_gather( mtag[dir]);
    }

    /* recursively copy the last time-like links to all time-slices */
    for(t=nc[TUP]-2;t>=0;t--){

	wait_gather( mtag[TUP]);

	FORALLSITES_OMP(i,s, default(shared) ){ if( (site_coord(s,xc[TUP]))==t ){
             su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]),&(s->link[xc[TUP]]));
	} } END_LOOP_OMP;

	if(t>0) restart_gather_site( F_OFFSET(link[xc[TUP]]), sizeof(su3_matrix),
		    xc[TUP], EVENANDODD, gen_pt[TUP], mtag[TUP] );
    }

    cleanup_gather( mtag[TUP]);

} /* new_ax_gauge */

