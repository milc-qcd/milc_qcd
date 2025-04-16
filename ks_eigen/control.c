/******************** control.c *****************************/
/* MIMD version 7 */
/* Main procedure for SU3 eigenvalues with improved dynamical fermions */

#define CONTROL
#include "ks_eig_includes.h"	/* definitions files and prototypes */

#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#endif
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"
#endif

#ifdef HAVE_GRID
#include "../include/generic_grid.h"
#endif

int main( int argc, char **argv ){
  register site *s;
  int i,si;
  int prompt;
  double dtime;
  double tottime;
  //  su3_vector **eigVec ;
  //  double *eigVal ;
  int total_R_iters ;
  double *resid = NULL;
  double chirality, chir_ev, chir_od ;
  eigVec = NULL;

  initialize_machine(&argc,&argv);
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  /* loop over input sets */
  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;

    tottime = -dclock();
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif
#ifdef HYPISQ_SVD_COUNTER
    hypisq_svd_counter = 0;
#endif

    /* Fix the gauge */
    
    if( param.fixflag == COULOMB_GAUGE_FIX)
      {
	if(this_node == 0) 
	  printf("Fixing to Coulomb gauge\n");
	
	rephase( OFF );
	STARTTIME;
	gaugefix(TUP,(Real)1.8,500,GAUGE_FIX_TOL);
	ENDTIME("gauge fix");
      }
    else
      if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");

    /* save lattice if requested */
    if( param.saveflag != FORGET ){
      rephase( OFF );
      savelat_p = save_lattice( param.saveflag, param.savefile, 
				param.stringLFN );
      rephase( ON );
    }

#ifdef U1_FIELD
    if( param.save_u1flag != FORGET ){
      save_u1_lattice( param.save_u1flag, param.save_u1file );
    }

    /* Insert U(1) phases into the gauge links in the site structure */
    u1phase_on(param.charge, u1_A);
    invalidate_fermion_links(fn_links);
#endif
    
    /* Recompute FN links from the links in the site structure */
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    
    /* Eigenpair calculation */
    STARTTIME;
    
    imp_ferm_links_t *fn = get_fm_links(fn_links, 0);

    /* Move KS phases and apply time boundary condition, based on the
       coordinate origin and time_bc */
    Real bdry_phase[4] = {0.,0.,0.,param.time_bc};
    /* Set values in the structure fn */
    set_boundary_twist_fn(fn, bdry_phase, param.coord_origin);
    /* Apply the operation */
    boundary_twist_fn(fn, ON);

    /* Calculate eigenpairs on even sites */
    param.eigen_param.parity = EVEN;
    total_R_iters=ks_eigensolve(eigVec, eigVal, &param.eigen_param, 1);
    fflush(stdout);

    /* Construct eigenpairs on odd sites */
    construct_eigen_other_parity(eigVec, eigVal, &param.eigen_param, fn);
    
    /* Calculate and print the residues and norms of the eigenvectors */
    resid = (double *)malloc(param.eigen_param.Nvecs*sizeof(double));
    node0_printf("Even site residuals\n");
    check_eigres( resid, eigVec, eigVal, param.eigen_param.Nvecs, EVEN, fn );
    node0_printf("Odd site residuals\n");
    check_eigres( resid, eigVec, eigVal, param.eigen_param.Nvecs, ODD, fn );
    fflush(stdout);

  /* Unapply twisted boundary conditions on the fermion links and
     restore conventional KS phases and antiperiodic BC, if
     changed. */
    boundary_twist_fn(fn, OFF);

    /* save eigenvectors if requested */
    int status = save_ks_eigen(param.ks_eigen_saveflag, param.ks_eigen_savefile,
      param.eigen_param.Nvecs, eigVal, eigVec, resid, 1);

    if(status != 0){
      node0_printf("ERROR writing eigenvectors\n"); fflush(stdout);
    }

    /* print eigenvalues of iDslash */
    node0_printf("The above were eigenvalues of -Dslash^2 in MILC normalization\n");
    node0_printf("Here we also list eigenvalues of iDslash in continuum normalization\n");
    for(i=0;i<param.eigen_param.Nvecs;i++)
      { 
	if ( eigVal[i] > 0.0 ){ chirality = sqrt(eigVal[i]) / 2.0; }
	else { chirality = 0.0; }
	node0_printf("eigenval(%i): %10g\n",i,chirality) ;
      } 

#ifdef U1_FIELD
    /* Unapply the U(1) field phases */
    u1phase_off();
    invalidate_fermion_links(fn_links);
#endif
    ENDTIME("calculate Dirac eigenpairs");
    
#ifdef CHIRALITY
    measure_chirality(eigVec[i], &chir_ev, EVEN);
    measure_chirality(eigVec[i], &chir_od, ODD);
    chirality = (chir_ev + chir_od) / 2.0;
    node0_printf("Chirality(%i) -- even, odd, total: %10g, %10g, %10g\n",
		 i,chir_ev,chir_od,chirality) ;
#endif

    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    tottime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",tottime);
      printf("total Rayleigh iters = %d\n",total_R_iters);
#ifdef HISQ_SVD_COUNTER
      printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
      printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
    }
    fflush(stdout);

    destroy_fn_links(fn);
  
  } /* readin(prompt) */

#ifdef EO
  cleanup_dslash_temps();
#endif
  /**
     for(i=0;i<Nvecs;i++)
     {
     sprintf(label,"DENSITY(%i)",i) ;
     print_densities(eigVec[i], label, ny/2,nz/2,nt/2, EVEN) ;
     }
  **/

  /* Clean up eigen storage */
  if(eigVec != NULL){
    for(i = 0; i < param.eigen_param.Nvecs; i++) free(eigVec[i]);
    free(eigVal); free(eigVec); free(resid);
  }
  invalidate_fermion_links(fn_links);
  
  /* Destroy fermion links (created in readin() */
  
#if FERM_ACTION == HISQ
  destroy_fermion_links_hisq(fn_links);
#else
  destroy_fermion_links(fn_links);
#endif
  fn_links = NULL;

  free_lattice();
  
#ifdef HAVE_QUDA
  finalize_quda();
#endif

#ifdef HAVE_GRID
  finalize_grid();
#endif
  
  normal_exit(0);
  return 0;
}

