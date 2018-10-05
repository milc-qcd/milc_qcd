/***************** control.c *****************************************/

/* Main procedure for measuring the staggered fermion chiral
   condensate and related quantities 

   Control Flow

   if DEFLATION: 
      compute eigenpairs, if not provided, by 
        1) if HAVE_PRIMME, using cheb_PRIMME
	2) otherwise, using Kalkreuter method (current default)
   if CURR_DISC:
      compute current for disconnected correlations either for single or multi-mass current
      * uses deflation method if eigenpairs are avaiable.
   else:
      compute chiral condensate (psi_bar psi)
   if CHEM_POT:
      compute chemical potential using previously computed psi_bar psi
   if EIGMODE == EIGCG:
      compute eigenpairs by conjugate gradient method

   Note:
   i)   If HAVEPRIMME, cheb_PRIMME will be invoked
   ii)  Otherwise, Kalkreuter will be invoked
   iii) Kalkreuter is equated by a macro to
          Kalkreuter_Ritz from eign_stuff_Ritz.c if HAVE_QOP is true
          Kalkreuter from eigen_stuff_qdp.c if HAVE_QOP is false but HAVE_QDP is true
	  Kalkreuter_Ritz from eign_stuff_Ritz.c if HAVE_QOP and HAVE_QDP are both false

*/

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
#include <quda_milc_interface.h>
#endif

#ifndef PARITY
#define PARITY EVEN /* Default value */
#endif 
#ifndef ISCORR
#define ISCORR 0
#endif
#ifndef IsColorThin
#define IsColorThin 0
#endif

int main(int argc, char *argv[])
{
  int prompt;
  int k;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif

#if EIGMODE == EIGCG || EIGMODE == DEFLATION
  int Nvecs_curr;
  double *resid = NULL;
  imp_ferm_links_t **fn;
#else
  int Nvecs_curr = 0;
#endif
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

    
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("set up");

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


#if EIGMODE == DEFLATION
    /**************************************************************/
    /* Compute Dirac eigenpairs           */

    STARTTIME;

    active_parity = EVEN;  /* EVEN is required */
    fn = get_fm_links(fn_links);
    Nvecs_curr = Nvecs_tot = param.Nvecs;

    /* compute eigenpairs if requested */
    if(param.ks_eigen_startflag == FRESH){
#ifdef PRIMME
      int mxv = cheb_PRIMME(eigVec, eigVal, param.eigenval_tol, Nvecs_curr, param.MaxIter, param.chebInfo, 0);
#else // default eigensolver
      int total_R_iters=Kalkreuter(eigVec, eigVal, param.eigenval_tol, param.error_decr,
			       Nvecs_curr, param.MaxIter, param.Restart, param.Kiters, 
			       param.ks_eigen_startflag == FRESH);
      node0_printf("total Rayleigh iters = %d\n", total_R_iters); fflush(stdout);
#endif //END: PRIMME
      construct_eigen_odd(eigVec, eigVal, Nvecs_curr, fn[0]);
    }

#if 1 /* If needed for debugging */
      /* (The Kalkreuter routine uses the random number generator to
	 initialize the eigenvector search, so, if you want to compare
	 first results with and without deflation, you need to
	 re-initialize here.) */
      initialize_site_prn_from_seed(iseed);
#endif
    
    /* Calculate and print the residues and norms of the eigenvectors */
    resid = (double *)malloc(Nvecs_curr*sizeof(double));
    node0_printf("Even site residuals\n");
    check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn[0] );
    construct_eigen_odd(eigVec, eigVal, Nvecs_curr, fn[0]);
    node0_printf("Odd site residuals\n");
    check_eigres( resid, eigVec, eigVal, Nvecs_curr, ODD, fn[0] );

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

    ENDTIME("calculate Dirac eigenpairs"); fflush(stdout);

#endif // END: EIGMODE == DEFLATION
    
    /**************************************************************/
    /* Compute chiral condensate and other observables            */

    STARTTIME;

    for(k = 0; k < param.num_set; k++){
      int num_pbp_masses = param.num_pbp_masses[k];
      int i0 = param.begin_pbp_masses[k];

      restore_fermion_links_from_site(fn_links, param.qic_pbp[i0].prec);

#ifdef CURRENT_DISC
	
	// initialize_site_prn_from_seed(iseed); /* Use the same random number sequence for all sets */

#ifdef PT2PT /* Exact computation using pt. to pt. sources */
      f_meas_current_pt2pt( param.num_pbp_masses[k], &param.qic_pbp[i0], &param.ksp_pbp[i0], fn_links,
			    &param.pbp_filenames[i0]);
#else /* Stochastic estimation */
      
      /* Compute the bias correction current */
      if(param.truncate_diff[k]){
	if(strstr(param.repeated_sloppy_cg_mode[k], "none")!=NULL);
	else 
	  f_meas_current_diff( param.num_pbp_masses[k], param.npbp_reps[k], param.nwrite[k], param.thinning[k], 
			       &param.qic_pbp[i0], &param.qic_pbp_sloppy[i0],
			       eigVec, eigVal, Nvecs_curr,
			       &param.ksp_pbp[i0], fn_links, 
			       param.slp_prec[k], param.num_pbp_slp_prec_sets[k],
			       &param.pbp_filenames[i0], 
			       PARITY, ISCORR, param.numMassPair[k], IsColorThin);
      }
      else {
	/* Compute the sloppy current */
	if(strstr(param.repeated_sloppy_cg_mode[k], "none")!=NULL);
	else 
	  f_meas_current( param.num_pbp_masses[k], param.npbp_reps[k], param.nwrite[k],	param.thinning[k], &param.qic_pbp[i0],
			  eigVec, eigVal, Nvecs_curr,
			  &param.ksp_pbp[i0], fn_links, 
			  param.slp_prec[k], param.num_pbp_slp_prec_sets[k],
			  &param.pbp_filenames[i0],
			  PARITY, param.numMassPair[k], IsColorThin);
      }
#endif
#else // Not CURRENT_DISK
      f_meas_imp_multi( param.num_pbp_masses[k], param.npbp_reps[k], &param.qic_pbp[i0], 
			&param.ksp_pbp[i0], fn_links);
#endif //END: CURRENT_DISC
#ifdef D_CHEM_POT
      Deriv_O6_multi( param.num_pbp_masses[k], param.npbp_reps[k], &param.qic_pbp[i0],
		      &param.ksp_pbp[i0], fn_links);
#endif
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

#if EIGMODE == EIGCG

    STARTTIME;

    active_parity = EVEN;
    Nvecs_curr = param.eigcgp.Nvecs_curr;

    fn = get_fm_links(fn_links);
    resid = (double *)malloc(Nvecs_curr*sizeof(double));

    if(param.ks_eigen_startflag == FRESH)
      calc_eigenpairs(eigVal, eigVec, &param.eigcgp, active_parity);

    check_eigres( resid, eigVec, eigVal, Nvecs_curr, active_parity, fn[0] );

    if(param.eigcgp.H != NULL) free(param.eigcgp.H);

    ENDTIME("compute eigenvectors");
#endif

    /**************************************************************/
    /* clean up */

#if EIGMODE == EIGCG || EIGMODE == DEFLATION

    /* save eigenvectors if requested */
    int status = save_ks_eigen(param.ks_eigen_saveflag, param.ks_eigen_savefile,
			       Nvecs_curr, eigVal, eigVec, resid, 1);
    if(status != 0){
      node0_printf("ERROR writing eigenvectors\n");
    }

    /* Clean up eigen storage */
    for(int i = 0; i < Nvecs_tot; i++) free(eigVec[i]);
    free(eigVal); free(eigVec); free(resid);

#endif

    node0_printf("RUNNING COMPLETED\n");
    endtime = dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_matrix_mult = %d\n",total_iters);
    fflush(stdout);

    destroy_ape_links_4D(ape_links);

    /* Destroy fermion links (created in readin() */

#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;

  } /* readin(prompt) */

#ifdef HAVE_QUDA
  qudaFinalize();
#endif

  normal_exit(0);
  return 0;
}
