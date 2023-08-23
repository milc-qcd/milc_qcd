/****** imp_gauge_action.c  -- ******************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* E.W. modified 7/22, split off gauge action */

/**#define GATIME**/ /* For timing gauge action calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

#ifdef ANISOTROPY
    /* for each rotation/reflection, an integer indicating if the path
       is spatial (=0) or temporal (=1) */
extern int **loop_st;
#endif

double imp_gauge_action_cpu() {
    register int i;
    int rep;
    register site *s;
    complex trace;
    double g_action;
    double action,act2,total_action;
    su3_matrix *tempmat1;
    su3_matrix *links;
    int length;

    /* these are for loop_table  */
    int ln,iloop;

    /* get loop variables from functions */
    const int max_length = get_max_length();
    const int nloop = get_nloop();
    const int nreps = get_nreps();
    const int *loop_length = get_loop_length();
    const int *loop_num = get_loop_num();
    int ***loop_table = get_loop_table();
    Real **loop_coeff = get_loop_coeff();

    g_action=0.0;


    tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
      printf("imp_gauge_action: Can't malloc temporary\n");
      terminate(1);
    }

    links = create_G_from_site();
    
    /* gauge action */
    for(iloop=0;iloop<nloop;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product_fields(links, loop_table[iloop][ln] , length, tempmat1 );

	    FORALLSITES_OMP(i,s,private(trace,action,total_action,act2,rep) reduction(+:g_action)){
		trace=trace_su3( &tempmat1[i] );
		action =  3.0 - (double)trace.real;
		/* need the "3 -" for higher characters */
#ifndef ANISOTROPY
        	total_action= (double)loop_coeff[iloop][0]*action;
#else
		/* NOTE: for anisotropic case every loop is multiplied by
		   the corresponding spatial (beta[0]) or temporal (beta[1])
		   coupling, while in the isotropic case all loops are
		   added together and then multiplied by beta outside
                   of this function */
        	total_action= (double)loop_coeff[iloop][0]*action
                              *beta[loop_st[iloop][ln]];
		/* loop_st[iloop][ln] is either 0 or 1 */
#endif
        	act2=action;
		for(rep=1;rep<nreps;rep++){
		    act2 *= action;
		    total_action += (double)loop_coeff[iloop][rep]*act2;
		}

        	g_action  += total_action;

	    } END_LOOP_OMP; /* sites */

	} /* ln */
    } /* iloop */

    g_doublesum( &g_action );
    destroy_G(links);
    special_free(tempmat1);
    return( g_action );
} /* imp_gauge_action_cpu */

/* for backwards compatibility */
double imp_gauge_action() {
    return imp_gauge_action_cpu();
} /* imp_gauge_action */


