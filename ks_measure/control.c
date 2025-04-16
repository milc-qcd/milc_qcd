/***************** control.c *****************************************/

/* Main procedure for measuring the staggered fermion chiral
   condensate and related quantities */

/* MIMD version 7 */

/* The chiral condensate is computed on a supplied background field */

/* Modifications ... */

// CD 2/5/16 Merged Hiroshi Ohno's support for deflated and eigcg solves */
   
//  $Log: control.c,v $
//  Revision 1.3  2012/05/08 20:40:16  detar
//  Call qudaFinalize to allow writing optimization file
//
//  Revision 1.2  2012/02/16 16:54:33  detar
//  Print hisq_svd_counter diagnostics only from node 0
//
//  Revision 1.1  2011/12/02 04:38:14  detar
//  Add
//
   

#define CONTROL
#include "ks_measure_includes.h"
#include <string.h>
#ifdef HAVE_QUDA
#include "../include/generic_quda.h"
#include <quda_milc_interface.h>
#endif
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#endif

extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

int main(int argc, char *argv[])
{
  int prompt;
  int k;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif

  int Nvecs_curr = 0;
  double *resid = NULL;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("set up");

#if 0
  gethostname(hostname, 128);
  
#endif

  /* loop over input sets */

  while( 1 ){

    starttime = dclock();
    if(readin(prompt) != 0)break;
    
    /* Skip calculations if just proofreading input parameters */
    if(prompt == 2)continue;

    total_iters=0;
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif


    /**************************************************************/
    /* Compute Dirac eigenpairs           */
    if(param.eigen_param.Nvecs > 0){

#if EIGMODE != EIGCG

      STARTTIME;
      
      /* First set of fn links is always charge 0 and Naik epsilon 0 */
      imp_ferm_links_t *fn = get_fm_links(fn_links, 0);

      /* Move KS phases and apply time boundary condition, based on the
	 coordinate origin and time_bc */
      Real bdry_phase[4] = {0.,0.,0.,param.time_bc};
      /* Set values in the structure fn */
      set_boundary_twist_fn(fn, bdry_phase, param.coord_origin);
      /* Apply the operation */
      boundary_twist_fn(fn, ON);
      
      Nvecs_curr = Nvecs_tot = param.eigen_param.Nvecs;
      
      /* compute eigenpairs if requested */
      if(param.ks_eigen_startflag == FRESH){
	int total_R_iters;
	total_R_iters=ks_eigensolve(eigVec, eigVal, &param.eigen_param, 1);
	node0_printf("total Rayleigh iters = %d\n", total_R_iters); fflush(stdout);
	construct_eigen_other_parity(eigVec, eigVal, &param.eigen_param, fn);
      }

      /* Check the eigenvectors */

      /* Calculate and print the residues and norms of the eigenvectors */
      resid = (double *)malloc(Nvecs_curr*sizeof(double));
      construct_eigen_other_parity(eigVec, eigVal, &param.eigen_param, fn);
      node0_printf("Even site residuals\n");
      check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn );
      node0_printf("Odd site residuals\n");
      check_eigres( resid, eigVec, eigVal, Nvecs_curr, ODD, fn );
      
      /* Unapply twisted boundary conditions on the fermion links and
	 restore conventional KS phases and antiperiodic BC, if
	 changed. */
      boundary_twist_fn(fn, OFF);
      
      /* print eigenvalues of iDslash */
      node0_printf("The above were eigenvalues of -Dslash^2 in MILC normalization\n");
      node0_printf("Here we also list eigenvalues of iDslash in continuum normalization\n");
      for(int i=0;i<Nvecs_curr;i++){ 
	if ( eigVal[i] > 0.0 ){
	  node0_printf("eigenval(%i): %10g\n", i, 0.5*sqrt(eigVal[i]));
	}
	else{
	  eigVal[i] = 0.0;
	  node0_printf("eigenval(%i): %10g\n", i, 0.0);
	}
      }

      destroy_fn_links(fn);

      ENDTIME("calculate/reload Dirac eigenpairs"); fflush(stdout);
      
#endif
    }
    
    /**************************************************************/
    /* Compute chiral condensate and other observables            */
    
    STARTTIME;
    
    for(k = 0; k < param.num_set; k++){

#if 0 /* If needed for debugging */
      /* (The Kalkreuter routine uses the random number generator to
	 initialize the eigenvector search, so, if you want to compare
	 results with and without deflation, you need to re-initialize
	 here.) */
      initialize_site_prn_from_seed(iseed);
#endif

      int num_pbp_masses = param.num_pbp_masses[k];
      Real masses[num_pbp_masses];
      Real charges[num_pbp_masses];
      imp_ferm_links_t *fn_mass[num_pbp_masses];
      int i0 = param.begin_pbp_masses[k];
      Real bdry_phase[4] = {0.,0.,0.,param.time_bc};
      
      /* Make table of FN links and masses and set boundary phases if
	 requested */
      imp_ferm_links_t *fn_pt[MAX_NAIK];
      for(int naik_index = 0; naik_index < MAX_NAIK; naik_index++)
	fn_pt[naik_index] = NULL;

      for(int j = 0; j < num_pbp_masses; j++){
	int naik_index = param.ksp_pbp[i0+j].naik_term_epsilon_index;
	int charge_index = param.ksp_pbp[i0+j].charge_index;
	if(fn_pt[naik_index] == NULL)
	  fn_pt[naik_index] = get_fm_links(fn_links_charge[charge_index], naik_index);
	fn_mass[j] = fn_pt[naik_index];
	masses[j] = param.ksp_pbp[i0+j].mass;
	charges[j] = param.ksp_pbp[i0+j].charge;

	/* Move KS phases and apply time boundary condition, based on the
	   coordinate origin and time_bc */
	/* Set values in the structure fn */
	set_boundary_twist_fn(fn_mass[j], bdry_phase, param.coord_origin);
	/* Apply the operation if not already done */
	if(twist_status(fn_mass[j]) == OFF)boundary_twist_fn(fn_mass[j], ON);
      }
      
#ifdef CURRENT_DISC
      if(param.truncate_diff[k])
	f_meas_current_diff( num_pbp_masses, param.npbp_reps[k],
			     param.thinning[k],
			     &param.qic_pbp[i0], &param.qic_pbp_sloppy[i0],
			     masses, charges, fn_mass, u1_A,
			     &param.pbp_filenames[i0] );
      else
	f_meas_current( num_pbp_masses, param.npbp_reps[k],
			param.thinning[k],
			&param.qic_pbp[i0], masses, charges,
			fn_mass, u1_A, &param.pbp_filenames[i0] );
#else
      // THESE NEED FIXING NOW
      if(num_pbp_masses == 1)
	f_meas_imp_field( param.npbp_reps[k], &param.qic_pbp[i0],
			  param.ksp_pbp[i0].mass, 
			  param.ksp_pbp[i0].naik_term_epsilon_index,
			  fn_links);
      else
	f_meas_imp_multi( param.num_pbp_masses[k], param.npbp_reps[k],
			  &param.qic_pbp[i0], &param.ksp_pbp[i0], fn_links);
      
#ifdef D_CHEM_POT
      if(num_pbp_masses == 1)
	Deriv_O6_field(param.npbp_reps[k], &param.qic_pbp[i0],
		       param.ksp_pbp[i0].mass,
		       fn_links, param.ksp_pbp[i0].naik_term_epsilon_index, 
		       param.ksp_pbp[i0].naik_term_epsilon);
      else
	Deriv_O6_multi( param.num_pbp_masses[k], param.npbp_reps[k],
			&param.qic_pbp[i0], &param.ksp_pbp[i0], fn_links);
#endif
#endif

      
      /* Unapply twisted boundary conditions on the fermion links and
	 restore conventional KS phases and antiperiodic BC, if
	 changed. */
#if 0 /* Not needed because we destroy their parents in the next lines */
      for(int j = 0; j < num_pbp_masses; j++)
	if(twist_status(fn_mass[j]) == ON)
	   boundary_twist_fn(fn_mass[j], OFF);
#endif

      for(int naik_index = 0; naik_index < MAX_NAIK; naik_index++)
	if(fn_pt[naik_index] != NULL)
	  destroy_fn_links(fn_pt[naik_index]);
          
    } /* k num_set */
 
#ifdef HISQ_SVD_COUNTER
    node0_printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
    
    ENDTIME("calculate observables");
    
    /* save lattice if requested */
    if( param.saveflag != FORGET ){
      STARTTIME;
      rephase( OFF );
      save_lattice( param.saveflag, param.savefile, param.stringLFN );
      rephase( ON );
      ENDTIME("save lattice");
    }

#ifdef U1_FIELD
    if( param.save_u1flag != FORGET ){
      save_u1_lattice( param.save_u1flag, param.save_u1file );
    }
#endif

#if EIGMODE == EIGCG

    STARTTIME;

    Nvecs_curr = param.eigcgp.Nvecs_curr;

    imp_ferm_links_t *fn = get_fm_links(fn_links, 0);
    resid = (double *)malloc(Nvecs_curr*sizeof(double));

    if(param.ks_eigen_startflag == FRESH)
      calc_eigenpairs(eigVal, eigVec, &param.eigcgp, EVEN);

    check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn );

    destroy_fn_links(fn);
    
    if(param.eigcgp.H != NULL) free(param.eigcgp.H);

    ENDTIME("compute eigenvectors");
#endif

    if(param.eigen_param.Nvecs > 0){
      /* save eigenvectors if requested */
      int status = save_ks_eigen(param.ks_eigen_saveflag, param.ks_eigen_savefile,
				 Nvecs_curr, eigVal, eigVec, resid, 1);
      if(status != 0){
	node0_printf("ERROR writing eigenvectors\n");
      }
      
      /* Clean up eigen storage */
      for(int i = 0; i < Nvecs_tot; i++) free(eigVec[i]);
      free(eigVal); free(eigVec); free(resid);
    }
    
    node0_printf("RUNNING COMPLETED\n");
    endtime = dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

    destroy_ape_links_4D(ape_links);

    /* Destroy fermion links (created in readin() */

    for(int j = 0; j < n_charges; j++){
#if FERM_ACTION == HISQ
      destroy_fermion_links_hisq(fn_links_charge[j]);
#else
      destroy_fermion_links(fn_links_charge[j]);
#endif
      fn_links_charge[j] = NULL;
    }
  } /* readin(prompt) */

#ifdef HAVE_QUDA
  finalize_quda();
#endif

  normal_exit(0);
  return 0;
}
