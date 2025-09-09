/********************** check_fermion_force.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the fermion force calculation */

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#else
BOMB THE COMPILE
#endif


void check_fermion_force( char phifile[MAX_MASS][MAXFILENAME], int phiflag,
			  char *ansfile, int ansflag, int n_naik, ks_param *ksp)
{
  Real diff, maxdiff, norm, maxnorm, reldiff;
  int i, dir;
  site *s;
  Real eps = 1.;
  int nflavors = 4;
  su3_matrix tmat, diffmat;
  char *filexml;
  char recxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test fermion force field</title>";
#if (MILC_PRECISION == 1)
  Real tol = 1e-3;
#else
  Real tol = 1e-8;
#endif
  int ff_prec = MILC_PRECISION;  /* Just use prevailing precision for now */
  /* Supports only asqtad at the moment */
  Real *residues = (Real *)malloc(n_naik*sizeof(Real));;
  if(residues == NULL){
    node0_printf("No room for residues\n");
    terminate(1);
  }

  su3_matrix *ansmom = (su3_matrix *)malloc(4*sites_on_node*sizeof(su3_matrix));
  if(ansmom == NULL){
    node0_printf("No room for ansmom\n");
    terminate(1);
  }

  /* Get test rational function parameters from file */
  params_rhmc *rf;
  char filename[] = "rat.2flavor";
  int n_pseudo = 1;
  rf = load_rhmc_params(filename, n_pseudo);
  if(rf == NULL)terminate(1);
			 
  /* Determine the maximum rational fcn order */
  max_rat_order = 0;
  for(i = 0; i < n_pseudo; i++){
    if(rf[i].MD.order > max_rat_order)max_rat_order = rf[i].MD.order;
    if(rf[i].GR.order > max_rat_order)max_rat_order = rf[i].GR.order;
    if(rf[i].FA.order > max_rat_order)max_rat_order = rf[i].FA.order;
  }
  if(mynode()==0)printf("Maximum rational func order is %d\n",max_rat_order);
  fflush(stdout);

  /* Determine the number of different Naik masses
     and fill in n_orders_naik and n_pseudo_naik        */
  Real current_naik_epsilon = rf[0].naik_term_epsilon;
  int tmporder = 0;
  int n_naiks = 0;
  double eps_naik[MAX_NAIK];
  imp_ferm_links_t *fn;

  n_order_naik_total = 0;
  for( i=0; i<n_pseudo; i++ ) {
    if( rf[i].naik_term_epsilon != current_naik_epsilon ) {
      if( tmporder > 0 ) {
        n_orders_naik[n_naiks] = tmporder;
	eps_naik[n_naiks] = current_naik_epsilon;
        current_naik_epsilon = rf[i].naik_term_epsilon;
        n_naiks++;
        n_order_naik_total += tmporder;
        tmporder = 0;
      }
    }
    tmporder += rf[i].MD.order;
    n_pseudo_naik[n_naiks]++;
  }
  if( tmporder > 0 ) {
    n_orders_naik[n_naiks] = tmporder;
    eps_naik[n_naiks] = current_naik_epsilon;
    n_order_naik_total += tmporder;
    n_naiks++;
  }

  su3_vector **phi = (su3_vector **)malloc(MAX_N_PSEUDO*sizeof(su3_vector *));
  for(i = 0; i < MAX_N_PSEUDO; i++)
    phi[i] = create_v_field();

  int n_multi_x = max_rat_order;
  su3_vector **multi_x = (su3_vector **)malloc(n_multi_x*sizeof(su3_vector *));
  for(i=0;i<n_multi_x;i++)
    multi_x[i] = create_v_field();

  Real *allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  su3_vector *sumvec = create_v_field();
  
  /* Make a random source in phi if we don't reload it */

  if(phiflag == RELOAD_SERIAL){
    for(int inaik = 0; inaik < n_naik; inaik++){
      restore_ks_vector_scidac_to_field(phifile[inaik], QIO_SERIAL, phi[inaik], 1);
      fflush(stdout);
    }
  }  else {
    
    int iphi = 0;
    for(int inaik = 0; inaik < n_naik; inaik++){
      
      /* For each pseudofermion belonging to this Naik epsilon
	 Generate g_rand random. Compute  phi =  (GR rat func) * g_rand
	 Compute multi_x = (Mdag M + beta_jphi)^{-1} g_rand */
      
      Real rsqmin_gr = param.qic[inaik].resid * param.qic[inaik].resid;
      int  niter_gr  = param.qic[inaik].max;
      int  prec_gr   = MILC_PRECISION;
      fn = get_fm_links(fn_links, inaik);
      
      for( int jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
	node0_printf("Generating random sources\n");
	grsource_imp_rhmc_field( phi[iphi], &rf[iphi].GR, EVEN,
				 multi_x, sumvec, rsqmin_gr, niter_gr,
				 prec_gr, fn, inaik, 
				 rf[iphi].naik_term_epsilon);
	iphi++;
      }
      destroy_fn_links(fn);
    }

    /* Construct  */
    
    tmporder = 0;
    iphi = 0;
    n_naik = fermion_links_get_n_naiks(fn_links);
    

    for( int inaik=0; inaik < n_naik; inaik++ ) {
      fn  = get_fm_links(fn_links, inaik);
      for( int jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
	
	// Add the current pseudofermion to the current set
	int order      = rf[iphi].MD.order;
	Real *residues = rf[iphi].MD.res;
	Real *roots    = rf[iphi].MD.pole;
	Real rsqmin_md = param.qic[inaik].resid * param.qic[inaik].resid;
	int  niter_md  = param.qic[inaik].max;
	
	/* Compute multi_x = M^\dagger M)^{-1} phi on even sites */
	Real final_rsq;
	int iters = 0;
	iters += ks_ratinv_field( phi[iphi], multi_x+tmporder, roots, residues,
				  order, niter_md, rsqmin_md, MILC_PRECISION, EVEN,
				  &final_rsq, fn, 
				  inaik, rf[iphi].naik_term_epsilon );
	
	/*Then compute multi_x = M * multi_x  = Dslash * multi_x on odd sites */
	for(int j=0;j<order;j++){
	  dslash_fn_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD,
			   fn);
	  allresidues[tmporder+j] = residues[j+1];
	  // remember that residues[0] is constant, no force contribution.
	}
	tmporder += order;
	iphi++;
      } /* jphi */
      destroy_fn_links(fn);
    } /* inaik */
  } /* phiflag != RELOAD_SERIAL */

  /* Compute the fermion force */

  node0_printf("Computing the fermion force\n"); fflush(stdout);
  
  /* Just to be safe, clear the answer */
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      s->mom[dir].m00im = 0.0;
      s->mom[dir].m11im = 0.0;
      s->mom[dir].m22im = 0.0;
      s->mom[dir].m01.real = 0.0;
      s->mom[dir].m01.imag = 0.0;
      s->mom[dir].m02.real = 0.0;
      s->mom[dir].m02.imag = 0.0;
      s->mom[dir].m12.real = 0.0;
      s->mom[dir].m12.imag = 0.0;
    }
  }

  eo_fermion_force_multi( 0.02, allresidues, multi_x,
			  n_order_naik_total, ff_prec, fn_links );

  /* If the answer file is given, read it for comparison */
  if(ansflag == RELOAD_SERIAL){
    restore_color_matrix_scidac_to_field(ansfile, ansmom, 4, MILC_PRECISION);
    node0_printf("Checking the answer\n"); fflush(stdout);
  }

  /* Clean up */

  free(residues);

  for(i=0;i<n_multi_x;i++)
    destroy_v_field(multi_x[i]);
  free(multi_x);

  /* Unpack the answer and compare if possible */
  maxdiff = 0;
  maxnorm = 0;
  norm = 0;
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat );
      /* If we have loaded an answer file, do the comparison */
      if(ansflag == RELOAD_SERIAL){
	sub_su3_matrix( ansmom + 4*i + dir, &tmat, &diffmat);
	diff = sqrt(realtrace_su3( &diffmat, &diffmat ));
	norm = sqrt(realtrace_su3( &tmat, &tmat));
	// printf("DIFF %g %g\n",norm,diff);
	if(diff > tol * norm){
	  printf("Intolerable relative difference %e node %d site %d\n",
		 diff/norm,this_node,i);
	  dumpmat(ansmom + 4*i + dir);
	  dumpmat(&tmat);
	}
	if(maxdiff < diff)maxdiff = diff;
	if(maxnorm < norm)maxnorm = norm;
      }
      /* In any case, copy the new result to the answer matrix */
      ansmom[4*i + dir] = tmat;
    }
  }

  if(ansflag == RELOAD_SERIAL){
    g_floatmax(&maxdiff);
    g_floatmax(&maxnorm);
    if(maxnorm > 0){
      reldiff = maxdiff/maxnorm;
      node0_printf("Relative difference %e\n",reldiff);
    }
    else
      node0_printf("Absolute difference %e but norm is 0???\n",maxdiff);
  }      

  /* Save source and answer if requested */
  if(phiflag == SAVE_SERIAL || phiflag == SAVE_PARTFILE_SCIDAC)
    for(i = 0; i < n_naik; i++){
#ifdef HAVE_QIO
      if(phiflag == SAVE_SERIAL)
	save_ks_vector_scidac_from_field(phifile[i], "check fermion force",
					"source color vector field", 
					QIO_SINGLEFILE, QIO_SERIAL, phi[i], 1, MILC_PRECISION);
      else if(phiflag == SAVE_PARTFILE_SCIDAC)
	save_ks_vector_scidac_from_field(phifile[i], "check fermion force",
					"source color vector field",
					QIO_PARTFILE, QIO_SERIAL, phi[i], 1, MILC_PRECISION);
    }

  if(ansflag == SAVE_SERIAL){
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
	recxml, QIO_SINGLEFILE, ansmom, 4, MILC_PRECISION, NULL);
    free_QCDML(filexml);
  }
  else if(ansflag == SAVE_PARTFILE_SCIDAC){
    node0_printf("Saving the momentum matrix\n");
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
	recxml, QIO_PARTFILE, ansmom, 4, MILC_PRECISION, NULL);
    free_QCDML(filexml);
  }
#endif

  free(ansmom);
  for(i = 0; i < MAX_N_PSEUDO; i++)
    destroy_v_field(phi[i]);
  free(phi);
}      
