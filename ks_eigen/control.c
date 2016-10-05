/******************** control.c *****************************/
/* MIMD version 7 */
/* Main procedure for SU3 eigenvalues with improved dynamical fermions */

#define CONTROL
#include "ks_eig_includes.h"	/* definitions files and prototypes */

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
  register site *s;
  int i,si;
  int prompt;
  double dtime;
  su3_vector **eigVec ;
  su3_vector *tmp ;
  double *eigVal ;
  int total_R_iters ;
  double chirality, chir_ev, chir_od ;
  imp_ferm_links_t **fn;

  initialize_machine(&argc,&argv);
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  /* set up */
  prompt = setup();

  /* loop over input sets */
  while( readin(prompt) == 0){
    
    dtime = -dclock();
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif
#ifdef HYPISQ_SVD_COUNTER
    hypisq_svd_counter = 0;
#endif
    restore_fermion_links_from_site(fn_links, PRECISION);
    /* call fermion_variable measuring routines */
    /* results are printed in output file */
    f_meas_imp( 1, PRECISION, F_OFFSET(phi), F_OFFSET(xxx), mass,
		0, fn_links);

    eigVal = (double *)malloc(Nvecs*sizeof(double));
    eigVec = (su3_vector **)malloc(Nvecs*sizeof(su3_vector*));
    for(i=0;i<Nvecs;i++)
      eigVec[i]=
	(su3_vector*)malloc(sites_on_node*sizeof(su3_vector));
    
    fn = get_fm_links(fn_links);
    active_parity = EVEN;
    total_R_iters=Kalkreuter(eigVec, eigVal, eigenval_tol, 
			     error_decr, Nvecs, MaxIter, Restart, 
			     Kiters, 1);

    node0_printf("The above where eigenvalues of -Dslash^2 in MILC normalization\n");
    node0_printf("Here we also list eigenvalues of iDslash in continuum normalization\n");
    for(i=0;i<Nvecs;i++)
      { 
	if ( eigVal[i] > 0.0 ){ chirality = sqrt(eigVal[i]) / 2.0; }
	else { chirality = 0.0; }
	node0_printf("eigenval(%i): %10g\n",i,chirality) ;
      } 

    tmp = (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));
    for(i=0;i<Nvecs;i++)
      { 
//	if ( eigVal[i] > 10.0*eigenval_tol )
//	  {
	    /* Construct to odd part of the vector.                 *
	     * Note that the true odd part of the eigenvector is    *
	     *  i/sqrt(eigVal) Dslash Psi. But since I only compute *
	     * the chirality the i factor is irrelevant (-i)*i=1!!  */
	    dslash_field(eigVec[i], tmp, ODD, fn[0]) ;
	    FORSOMEPARITY(si,s,ODD){ 
	      scalar_mult_su3_vector( &(tmp[si]),
				      1.0/sqrt(eigVal[i]), 
				      &(eigVec[i][si]) ) ;
	    }
	
//	    measure_chirality(eigVec[i], &chirality, EVENANDODD);
	    /* Here I divide by 2 since the EVEN vector is normalized to
	     * 1. The EVENANDODD vector is normalized to 2. I could have
	     * normalized the EVENANDODD vector to 1 and then not devide
	     * by to.  The measure_chirality routine assumes vectors
	     * normalized to 1.  */
//	    node0_printf("Chirality(%i): %g\n",i,chirality/2) ;
	    measure_chirality(eigVec[i], &chir_ev, EVEN);
	    measure_chirality(eigVec[i], &chir_od, ODD);
	    chirality = (chir_ev + chir_od) / 2.0;
	    node0_printf("Chirality(%i) -- even, odd, total: %10g, %10g, %10g\n",
			 i,chir_ev,chir_od,chirality) ;
//	  }
//	else
//	  {
	    /* This is considered a "zero mode", and treated as such */
//	    measure_chirality(eigVec[i], &chirality, EVEN);
//	    node0_printf("Chirality(%i): %g\n",i,chirality) ;
//	  }
      }
#ifdef EO
    cleanup_dslash_temps();
#endif
    free(tmp);
    /**
       for(i=0;i<Nvecs;i++)
       {
       sprintf(label,"DENSITY(%i)",i) ;
       print_densities(eigVec[i], label, ny/2,nz/2,nt/2, EVEN) ;
       }
    **/
    for(i=0;i<Nvecs;i++)
      free(eigVec[i]) ;
    free(eigVec) ;
    free(eigVal) ;
    invalidate_fermion_links(fn_links);
    fflush(stdout);
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
      printf("total Rayleigh iters = %d\n",total_R_iters);
#ifdef HISQ_SVD_COUNTER
      printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
      printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
    }
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

