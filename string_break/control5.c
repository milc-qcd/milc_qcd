/************************ control5.c ******************************/
/* MIMD version 6 */
/* Main procedure for pure gauge SU3 */

/* original code by UMH */
/* 5/25/98 Version 5 port CD */

/* This version does some required number of smearing iterations,
   then transforms to axial gauge and computes (time-like) Wilson loops
   for the computation of the heavy quark potential. */
/* This version:  
   only step=1 implemented (i.e. random sources on ALL
   sites, and not step=2 (even sites only). Reason: for wl2_*        
   correlators, many entries identically equal to zero due to
   construction 
   Also: reinstate KS phases before calculating correlators
   */ 


#define CONTROL
#include "string_break_includes.h"

int main(int argc, char *argv[])  {
int todo, sm_lev;
int prompt;
double dtime;

int m_iters, avm_iters, meascount;
int i, j, k, t0, step;
int dir;
site *s;
Real rsq;

initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

        dtime = -dclock();
 
	/* fix to axial gauge */
	if( startflag != CONTINUE){
#ifdef RAN_GAUGE
	    rand_gauge(F_OFFSET(rgt));
#endif
	    ax_gauge();
	    tot_smear = 0;
	    /* Save these for inversion */
	    /* First put in KS phases */
	    rephase( ON );
	    for(dir=XUP; dir<TUP; dir++) FORALLSITES(i,s){
		su3mat_copy( &(s->link[dir]), &(s->link_save[dir]));
	    }
	    FORALLSITES(i,s){
		if(s->t == nt-1){
		    su3mat_copy( &(s->link[TUP]), &(s->link_save[TUP]));
		}
		else{
		    for(j=0; j<3; j++) for(k=0; k<3; k++){
			if(j != k){
			    s->link_save[TUP].e[j][k] = cmplx(0.0,0.0);
			}
			else{
			    s->link_save[TUP].e[j][k].real = 1.0;
			    s->link_save[TUP].e[j][k].imag = 0.0;
			}
		    }
		}
	    }
	    /* Take KS phases out again */
	    /* The unsmeared links, in "link_save" have the KS phases. */
	    rephase( OFF );
	}

	meascount = avm_iters = 0;
	if(this_node==0)printf("END OF HEADER\n");

	/* Loop over the different smearing levels */
	for(sm_lev=0; sm_lev < no_smear_level; sm_lev++ ){

	    /* Do the smearing iterations */
	    for(todo=smear_num[sm_lev]; todo > 0; --todo ){
		smearing();
	    }
	    if(this_node==0)printf("SMEARING COMPLETED\n");  fflush(stdout);
	    tot_smear += smear_num[sm_lev];

	    rephase( ON );

 	    /*   use ONLY step=1 here - i.e all sites used 
                (step=2 for even sites only)  */
	    for(step=1; step<2; step++){
	    /* Create "fuzzed" random Gaussian sources in qprop */
	    t0 = nt;
	    k = off_axis_flag;
	    for(j=0; j<num_src; j++)
	    {
	     off_axis_flag = j;
		fuz_source( F_OFFSET(g_rand[j]), F_OFFSET(qprop[j]), r0,
		    t0, step);
	    }
	    off_axis_flag = k;

	    /* Compute the light quark propagators */
	    /* Permute smeared and unsmeared (axial gauge) gauge fields */
	    for(dir=XUP; dir<=TUP; dir++) FORALLSITES(i,s){
		su3mat_copy( &(s->link_save[dir]), &(s->staple));
		su3mat_copy( &(s->link[dir]), &(s->link_save[dir]));
		su3mat_copy( &(s->staple), &(s->link[dir]));
	    }

	    /* NOTE: NEED TO CHANGE THIS TO BUILD ks_act_paths AND
	       CALL mat_invert_uml INSTEAD OF ks_congrad with EVENANDODD */
	    load_ferm_links(&fn_links);
	    for(j=0; j<num_src; j++){
		/* Put source, m^\dag qprop, in phi, result in xxx */
		dslash_fn_site( F_OFFSET(qprop[j]), F_OFFSET(phi), EVENANDODD,
				fn_links);
		scalar_mult_latvec( F_OFFSET(phi), -1.0, F_OFFSET(phi),
		    EVENANDODD);
		scalar_mult_add_latvec( F_OFFSET(phi), F_OFFSET(qprop[j]),
		    2.0*mass, F_OFFSET(phi), EVENANDODD);
		clear_latvec( F_OFFSET(xxx), EVENANDODD);
		
		m_iters = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				   niter, rsqmin, PRECISION, EVENANDODD, &rsq,
				   fn_links);
		avm_iters += m_iters;
		++meascount;
		copy_latvec( F_OFFSET(xxx), F_OFFSET(qprop[j]), EVENANDODD);
	    }


	    /* Permute smeared and unsmeared (axial gauge) gauge fields */
	    for(dir=XUP; dir<=TUP; dir++) FORALLSITES(i,s){
		su3mat_copy( &(s->link_save[dir]), &(s->staple));
		su3mat_copy( &(s->link[dir]), &(s->link_save[dir]));
		su3mat_copy( &(s->staple), &(s->link[dir]));
	    }

	    for(j=0; j<num_src; j++)
		fuz_prop( F_OFFSET(qprop[j]), r0);

	    stat_li_mesons(tot_smear, step);


	    /* Compute on-axis Wilson-(1)light (WL) loops */
	    wl_1l_1corr(tot_smear, step);
	    fflush(stdout);

	    /* Compute on-axis Wilson-(1)light (WL) loops */
	    wl_1l_2corr(tot_smear, step);
	    fflush(stdout);

	    /* Compute on-axis Wilson-(2)light (WLL) loops */
	    wl_2l_1corr(tot_smear, step);
	    fflush(stdout);

	    /* Compute on-axis Wilson-(2)light (WLL) loops */
	    wl_2l_2corr(tot_smear, step);
	    fflush(stdout);

	    /* Compute off-axis Wilson-(1)light (WL) loops */
	    wl_1l_1corr_offax(tot_smear, step);
	    fflush(stdout);

	    /* Compute off-axis Wilson-(1)light (WL) loops */
	    wl_1l_2corr_offax(tot_smear, step);
	    fflush(stdout);

	    /* Compute off-axis Wilson-(1)light (WL) loops */
	    wl_2l_1corr_offax(tot_smear, step);
	    fflush(stdout);

	    /* Compute off-axis Wilson-(1)light (WL) loops */
	    wl_2l_2corr_offax(tot_smear, step);
	    fflush(stdout);

}

/**	    t0 = 0;
	    k = off_axis_flag;
	    for(j=0; j<num_src; j++)
	    {
	     off_axis_flag = j;
	     fuz_source( F_OFFSET(g_rand[j]), F_OFFSET(qprop[j]), r0, t0, 2);
	    }
	    off_axis_flag = k; **/

	    /* Compute the light quark propagators */
	    /* Permute smeared and unsmeared (axial gauge) gauge fields */
/**	    for(dir=XUP; dir<=TUP; dir++) FORALLSITES(i,s){
		su3mat_copy( &(s->link_save[dir]), &(s->staple));
		su3mat_copy( &(s->link[dir]), &(s->link_save[dir]));
		su3mat_copy( &(s->staple), &(s->link[dir]));
	    }

	    for(j=0; j<num_src; j++){ **/
		/* Put source, m^\dag qprop, in phi, result in xxx */
/**		dslash_site( F_OFFSET(qprop[j]), F_OFFSET(phi), EVENANDODD);
		scalar_mult_latvec( F_OFFSET(phi), -1.0, F_OFFSET(phi),
		    EVENANDODD);
		scalar_mult_add_latvec( F_OFFSET(phi), F_OFFSET(qprop[j]),
		    2.0*mass, F_OFFSET(phi), EVENANDODD);
		clear_latvec( F_OFFSET(xxx), EVENANDODD);
		m_iters = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				   niter, rsqmin, PRECISION, EVENANDODD, &rsq);
		avm_iters += m_iters;
		++meascount;
		copy_latvec( F_OFFSET(xxx), F_OFFSET(qprop[j]), EVENANDODD);
	    } **/

	    /* Permute smeared and unsmeared (axial gauge) gauge fields */
/**	    for(dir=XUP; dir<=TUP; dir++) FORALLSITES(i,s){
		su3mat_copy( &(s->link_save[dir]), &(s->staple));
		su3mat_copy( &(s->link[dir]), &(s->link_save[dir]));
		su3mat_copy( &(s->staple), &(s->link[dir]));
		} **/

	    /* Compute point-sink light-light mesons */
/**	    li_li_mesons(tot_smear, t0, "POINT"); **/

	    /* Compute smeared-sink light-light mesons */
/**	    for(j=0; j<num_src; j++)
		fuz_prop( F_OFFSET(qprop[j]), r0);

	    li_li_mesons(tot_smear, t0, "SMEAR"); **/
	    /* Compute on-axis time-like Wilson loops */
	    w_loop1(tot_smear);
  	    w_loop2(tot_smear);

	    rephase( OFF );
	}
        if(this_node==0)printf("RUNNING COMPLETED\n");
	if(meascount>0)
	    if(this_node==0)printf("average cg iters for measurement= %e\n",
		(double)avm_iters/meascount);

        dtime += dclock();
        if(this_node==0){
            printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
        }
        fflush(stdout);
	dtime = -dclock();

	/* save lattice if requested */
	if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
	}
    }
    return 0;
}
