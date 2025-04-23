/***************** control.c *****************************************/

/* Main procedure for staggered fermion spectroscopy     	     */
/* MIMD version 7 */

/* This version computes propagators for staggered fermions on a supplied
   background field config and ties them together according to the
   input parameters. */

/* Modifications ... */
   
//  $Log: control.c,v $
//  Revision 1.8  2013/12/28 20:57:27  detar
//  Fix improper type for creating ksprop quark field
//
//  Revision 1.7  2013/12/24 05:39:00  detar
//  Add combo operation
//
//  Revision 1.6  2012/11/24 05:14:20  detar
//  Add support for U(1) fields and for future HYPISQ action
//
//  Revision 1.5  2012/05/08 20:39:54  detar
//  Call qudaFinalize to allow writing optimization file.
//
//  Revision 1.4  2012/04/25 03:23:21  detar
//  Fix rephase flag
//
//  Revision 1.3  2012/01/21 21:34:36  detar
//  Move start time to beginning.  Remake APE links after gauge fixing.
//
//  Revision 1.2  2011/12/03 03:43:39  detar
//  Cosmetic
//
//  Revision 1.1  2011/11/30 22:11:38  detar
//  Add
//
   

#define CONTROL
#include "ks_spectrum_includes.h"
#include <string.h>
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#endif
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef HAVE_GRID
#include "../include/generic_grid.h"
#endif

int main(int argc, char *argv[])
{
  char myname[] = "main";
  int prompt;
  int i, is, j, k, iq0, iq1, iq2;
#ifdef KS_LEAN
  int oldiq0, oldiq1, oldiq2, oldip0;
  int quark_nc[MAX_QK];
#endif
  int naik_index, naik_index0, naik_index1;
  double mass;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  ks_prop_field *source[MAX_SOURCE];
  ks_prop_field *prop[MAX_PROP];
  ks_prop_field *quark[MAX_QK];
  int prop_nc[MAX_PROP];
  int Nvecs_curr;
  double *resid = NULL;
  
  initialize_machine(&argc,&argv);

  for(i = 0; i < MAX_SOURCE; i++)source[i] = NULL;
  for(i = 0; i < MAX_PROP; i++)prop[i] = NULL;
  for(i = 0; i < MAX_QK; i++)quark[i] = NULL;

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

  starttime=dclock();
    
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;

    total_iters=0;
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif

    STARTTIME;
    
    /**************************************************************/
    /* Do whatever is needed to get lattice */
    if( param.startflag == CONTINUE ){
      rephase( OFF );
    }
    if( param.startflag != CONTINUE ){
      startlat_p = reload_lattice( param.startflag, param.startfile );
    }
    /* if a lattice was read in, put in KS phases and AP boundary condition */
    phases_in = OFF;
    rephase( ON );

#ifdef U1_FIELD
    /* Read the U(1) gauge field, if wanted */
    start_u1lat_p = reload_u1_lattice( param.start_u1flag, param.start_u1file);
#endif

    ENDTIME("read lattice");

    /**************************************************************/
    /* Fix the gauge, but not if we are "continuing"              */
    
    if( param.fixflag == COULOMB_GAUGE_FIX && ! (param.startflag == CONTINUE) )
      {
	if(this_node == 0) 
	  printf("Fixing to Coulomb gauge\n");

	rephase( OFF );

	STARTTIME;
	gaugefix(TUP,(Real)1.8,500,GAUGE_FIX_TOL);
	//gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
	ENDTIME("gauge fix");

#if 0
	/* (Re)construct APE smeared links after gauge fixing.  
	   No KS phases here! */
	destroy_ape_links_4D(ape_links);
	ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
	if(param.time_bc == 0)apply_apbc( ape_links, param.coord_origin[3] );
	refresh_ape_links = 1;  // To signal refreshing of any cached links
	ape_links_ks_phases = OFF;  
	/* By default, the phases are ON */
	rephase_field_offset( ape_links, ON, &ape_links_ks_phases, param.coord_origin );
#endif
	
	rephase( ON );
	invalidate_fermion_links(fn_links);

      }
    else
      if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
    
    /**************************************************************/
    /* Construct APE smeared links without KS phases, but with
       conventional antiperiodic bc.  This is the same initial setup as
       the gauge field itself.  Later the phases are adjusted according
       to boundary phases and momentum twists.  If we are reading from a
       file, we assume it was saved with the same conventions.
    */
#ifdef APE_LINKS_FILE
    
    if(param.start_ape_flag == FRESH){

      /* Do APE smearing */
      rephase( OFF );

      ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
      if(param.time_bc == 0)apply_apbc( ape_links, param.coord_origin[3] );
      refresh_ape_links = 1;
      ape_links_ks_phases = OFF;
      /* By default, the KS phases in the APE links are ON */
      rephase_field_offset( ape_links, ON, &ape_links_ks_phases, param.coord_origin );
      
      rephase( ON );

    } else {

      /* Reload APE links from a file */
      ape_links = create_G();
      reload_apelinks( param.start_ape_flag, ape_links, param.start_ape_file );
      ape_links_ks_phases = ON;  /* Because we save them with phases on */

    }
    
#else

    rephase( OFF );

    /* Do APE smearing */
    ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
    if(param.time_bc == 0)apply_apbc( ape_links, param.coord_origin[3] );
    refresh_ape_links = 1;
    ape_links_ks_phases = OFF;
    /* By default, the phases are ON */
    rephase_field_offset( ape_links, ON, &ape_links_ks_phases, param.coord_origin );
    
    rephase( ON );

#endif
  
    /**************************************************************/
    /* save lattice if requested */

    STARTTIME;
    
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
#endif

#ifdef APE_LINKS_FILE

    /* Save the APE links to a file if requested */
    if(param.save_ape_flag != FORGET){
      if(ape_links == NULL){
	node0_printf("ERROR: main: Requested saving empty APE links\n");
	terminate(1);
      }
      save_apelinks( param.save_ape_flag, ape_links, param.save_ape_file );
    }

#endif

    ENDTIME("save lattice");

    /**************************************************************/
    /* Set up fermion links */
    
    STARTTIME;
    
#ifdef DBLSTORE_FN
    /* We want to double-store the links for optimization */
    fermion_links_want_back(1);
#endif
    
    /* Don't need to save HISQ auxiliary links */
    fermion_links_want_aux(0);
    
#if FERM_ACTION == HISQ
    
#ifdef DM_DEPS
    fermion_links_want_deps(1);
#endif
    
    fn_links = create_fermion_links_from_site(MILC_PRECISION, param.n_naiks, param.eps_naik);
    
#else
    
#ifdef DM_DU0
    fermion_links_want_du0(1);
#endif
    
    fn_links = create_fermion_links_from_site(MILC_PRECISION, 0, NULL);
    
#endif
    
    ENDTIME("create fermion links");

    imp_ferm_links_t *fn = get_fm_links(fn_links, 0);
    
    /**************************************************************/
    /* Set up eigenpairs, if requested */

    STARTTIME;
      
#if EIGMODE == EIGCG
    int Nvecs_max = param.eigcgp.Nvecs_max;
    if(param.ks_eigen_startflag == FRESH)
      Nvecs_tot = ((Nvecs_max - 1)/param.eigcgp.Nvecs)*param.eigcgp.Nvecs
	+ param.eigcgp.m;
    else
      Nvecs_tot = Nvecs_max;
    
    Nvecs_alloc = Nvecs_tot;
    eigVal = (double *)malloc(Nvecs_alloc*sizeof(double));
    eigVec = (su3_vector **)malloc(Nvecs_alloc*sizeof(su3_vector *));
    for(i = 0; i < Nvecs_alloc; i++)
      eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    
    /* Do whatever is needed to get eigenpairs */
    int status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile,
				 &Nvecs_tot, eigVal, eigVec, fn, 1);
    
    if(param.fixflag != NO_GAUGE_FIX){
      node0_printf("WARNING: Gauge fixing does not readjust the eigenvectors\n");
    }
    if(status != 0) normal_exit(0);
    
    if(param.ks_eigen_startflag != FRESH){
      param.eigcgp.Nvecs = 0;
      param.eigcgp.Nvecs_curr = Nvecs_tot;
      param.eigcgp.H = (double_complex *)malloc(Nvecs_max*Nvecs_max
						*sizeof(double_complex));
      for(i = 0; i < Nvecs_max; i++){
	for(k = 0; k < i; k++)
	  param.eigcgp.H[k + Nvecs_max*i] = dcmplx((double)0.0, (double)0.0);
	param.eigcgp.H[(Nvecs_max+1)*i] = dcmplx(eigVal[i], (double)0.0);
      }
    }
#endif
    
#if EIGMODE != EIGCG
    /* If using QUDA for deflation, then eigenvectors are loaded directly by QUDA and not MILC */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) && defined(USE_EIG_GPU) )
    if(param.eigen_param.Nvecs > 0){
      /* malloc for eigenpairs */
      eigVal = (Real *)malloc(param.eigen_param.Nvecs*sizeof(double));
      eigVec = (su3_vector **)malloc(param.eigen_param.Nvecs*sizeof(su3_vector *));
      for(i=0; i < param.eigen_param.Nvecs; i++){
	eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
	if(eigVec[i] == NULL){
	  printf("No room for eigenvector\n");
	  terminate(1);
	}
      }
      
      /* Do whatever is needed to get eigenpairs */
      int status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
				   &param.eigen_param.Nvecs, eigVal, eigVec, fn, 1);
      if(param.fixflag != NO_GAUGE_FIX){
	node0_printf("WARNING: Gauge fixing does not readjust the eigenvectors");
      }
    }
#endif
#endif
    
    /**************************************************************/
    /* Compute Dirac eigenpairs           */
    
    Nvecs_curr = Nvecs_tot = param.eigen_param.Nvecs;
      
    if(param.eigen_param.Nvecs > 0){
      
#if EIGMODE != EIGCG
      
      param.eigen_param.parity = EVEN;  /* Required */

      /* Move KS phases and apply time boundary condition, based on the
	 coordinate origin and time_bc */
      Real bdry_phase[4] = {0.,0.,0.,(Real)param.time_bc};
      /* Set values in the structure fn */
      set_boundary_twist_fn(fn, bdry_phase, param.coord_origin);
      /* Apply the operation */
      boundary_twist_fn(fn, ON);

      // If using QUDA deflated CG + asking for QUDA to do the eigensolve, then
      // the eigensolver is called from within QUDA's CG solver...not from MILC
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) && defined(USE_EIG_GPU) ) 
      /* compute eigenpairs if requested */
      if(param.ks_eigen_startflag == FRESH){
	int total_R_iters;
	total_R_iters=ks_eigensolve(eigVec, eigVal, &param.eigen_param, 1);
	construct_eigen_other_parity(eigVec, eigVal, &param.eigen_param, fn);
	node0_printf("total Rayleigh iters = %d\n", total_R_iters);
	
#if 0 /* If needed for debugging */
      /* (The ks_eigensolve routine uses the random number generator to
	 initialize the eigenvector search, so, if you want to compare
	 first results with and without deflation, you need to
	 re-initialize here.) */
	initialize_site_prn_from_seed(iseed);
#endif
      }
#endif
      /* Check the eigenvectors */

      /* If using QUDA for deflation, then eigenvectors are loaded directly by QUDA and not checked by MILC */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) && defined(USE_EIG_GPU) )
      /* Calculate and print the residues and norms of the eigenvectors */
      resid = (double *)malloc(Nvecs_curr*sizeof(double));
      node0_printf("Even site residuals\n");
      check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn );
      construct_eigen_other_parity(eigVec, eigVal, &param.eigen_param, fn);
      node0_printf("Odd site residuals\n");
      check_eigres( resid, eigVec, eigVal, Nvecs_curr, ODD, fn );
#endif
      /* Unapply twisted boundary conditions on the fermion links and
	 restore conventional KS phases and antiperiodic BC, if
	 changed. */
      boundary_twist_fn(fn, OFF);
     
      /* If using QUDA for deflation, then eigenvalues are printed by QUDA */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) ) 
      /* print eigenvalues of iDslash */
      node0_printf("The above were eigenvalues of -Dslash^2 in MILC normalization\n");
      node0_printf("Here we also list eigenvalues of iDslash in continuum normalization\n");
      for(i=0;i<Nvecs_curr;i++){ 
	if ( eigVal[i] > 0.0 ){
	  node0_printf("eigenval(%i): %10g\n", i, 0.5*sqrt(eigVal[i]));
	}
	else{
	  eigVal[i] = 0.0;
	  node0_printf("eigenval(%i): %10g\n", i, 0.0);
	}
      }
#endif
#endif
    }
    
    ENDTIME("calculate Dirac eigenpairs");

    /**************************************************************/
    /* Compute chiral condensate and related quantities           */
    
    STARTTIME;
    
    /* Make fermion links if not already done */
    
    for(i = 0; i < param.num_pbp_masses; i++){
#ifdef U1_FIELD
      u1phase_on(param.charge_pbp[i], u1_A);
      invalidate_fermion_links(fn_links);
      restore_fermion_links_from_site(fn_links, param.qic_pbp[i].prec);
#endif

      naik_index = param.ksp_pbp[i].naik_term_epsilon_index;
      mass = param.ksp_pbp[i].mass;

      f_meas_imp_field( param.npbp_reps, &param.qic_pbp[i], mass,
      			naik_index, fn_links);
#ifdef D_CHEM_POT
      
      Deriv_O6_field( param.npbp_reps, &param.qic_pbp[i], mass,
      		      fn_links, naik_index, 
		      param.ksp_pbp[i].naik_term_epsilon);
#endif
#ifdef U1_FIELD
      /* Unapply the U(1) field phases */
      u1phase_off();
      invalidate_fermion_links(fn_links);
      restore_fermion_links_from_site(fn_links, param.qic_pbp[i].prec);
      fn = get_fm_links(fn_links, 0);
#endif
    }

    ENDTIME("calculate pbp, etc");

    if(this_node==0)printf("END OF HEADER\n");
    

    /**************************************************************/

    /* Create sources */

    STARTTIME;

    /* Base sources */

    for(k=0; k<param.num_base_source; k++){
      quark_source *qs = &param.src_qs[k];
      source[k] = create_ksp_field(qs->ncolor);

      /* Open file with metadata */
      if(qs->saveflag != FORGET){
	char *fileinfo = create_ks_XML();
	w_source_open_ks(qs, fileinfo);
	free(fileinfo);
      }
      
      for(int color = 0; color < qs->ncolor; color++){

	/* Create a base source */
	if(v_source_field(source[k]->v[color], qs)){
	  printf("%s(%d): error getting source\n",myname,this_node);
	  terminate(1);
	}

	/* Write the source, if requested */
	if(qs->saveflag != FORGET){
	  qs->color = color;
	  if(w_source_ks( source[k]->v[color], qs ) != 0)
	    node0_printf("Error writing source\n");
	}
      } /* color */
  
      if(qs->saveflag != FORGET) w_source_close(qs);
    }


    /* Modified sources */
    for(is=param.num_base_source; is<param.num_base_source+param.num_modified_source; is++){

      quark_source *qs = &param.src_qs[is];
      int p = param.parent_source[is];
      source[is] = create_ksp_field(source[p]->nc);
      
      if(qs->saveflag != FORGET){
	char *fileinfo = create_ks_XML();
	w_source_open_ks(qs, fileinfo);
	free(fileinfo);
      }
      
      /* Copy parent source */
      copy_ksp_field(source[is],  source[p]);

      for(int color = 0; color < qs->ncolor; color++){

	node0_printf("Creating modified source %d for color %d from parent source %d\n",
		     is, color, p);
	/* Apply operator*/
        v_field_op(source[is]->v[color], &(param.src_qs_op[is]), qs->subset, qs->t0);

	/* Write the source, if requested */
	if(qs->saveflag != FORGET){
	  qs->color = color;
	  if(w_source_ks( source[is]->v[color], qs ) != 0)
	    node0_printf("Error writing source\n");
	}
      } /* color */

      if(qs->saveflag != FORGET) w_source_close(qs);

    } /* is */


#if defined(HAVE_QUDA) && defined(USE_GSMEAR_GPU)
    // delete 2-link current used for smearing
    gauss_smear_delete_2link_QUDA();
#endif

    ENDTIME("create sources");

    /**************************************************************/


    /* Loop over sets of propagators */

    STARTTIME;

    /* Temporary lists */
    ks_prop_field *tmp_source[MAX_PROP];
    quark_source *tmp_src_qs[MAX_PROP];
    
    for(k=0; k<param.num_set; k++){
      int num_prop = param.end_prop[k] - param.begin_prop[k] + 1;
      if(num_prop <= 0)continue;  /* Ignore set if zero */
      int i0 = param.begin_prop[k];

      for(i=param.begin_prop[k]; i <= param.end_prop[k]; i++){

	/**************************************************************/
	/* Read and/or generate quark propagator */
	
	is = param.source[i];  /* source index for this propagator */
	//	prop_nc[i] = param.src_qs[is].ncolor;
	/* Get number of colors from the parent */
	int p = param.parent_source[is];
	if(p == BASE_SOURCE_PARENT)
	  prop_nc[i] = param.src_qs[is].ncolor;
	else
	  prop_nc[i] = source[p]->nc;
	
	/* Allocate propagator */
	  
	prop[i] = create_ksp_field(prop_nc[i]);
	if(prop[i] == NULL){
	  printf("main(%d): No room for prop\n",this_node);
	  terminate(1);
	}
	tmp_source[i] = source[is];  /* Pointer copy */
	tmp_src_qs[i] = &param.src_qs[is]; 

	node0_printf("Mass= %g source %s ",
		     (double)param.ksp[i].mass,
		     param.src_qs[is].descrp);
#ifdef U1_FIELD
	node0_printf("Q %g ",param.charge[k]);
#endif
	node0_printf("residue= %g rel= %g\n",
		     (double)param.qic[i].resid,
		     (double)param.qic[i].relresid);
	fflush(stdout);
      }
      
      /* We pass the beginning addresses of the set data */
      
      total_iters += solve_ksprop(param.set_type[k],
				  param.inv_type[k],
				  num_prop,
				  param.startflag_ks + i0,
				  param.startfile_ks + i0,
				  param.saveflag_ks + i0,
				  param.savefile_ks + i0,
				  prop + i0,
				  tmp_source + i0,
				  tmp_src_qs + i0,
				  param.qic + i0, 
				  param.ksp + i0,
				  param.charge[k],
				  param.bdry_phase[i0],
				  param.coord_origin,
				  param.check[i0]);
      
    } /* sets */
    ENDTIME("compute propagators");
    
    
    /*****************************************************************/
    /* Complete the quark propagators by applying the sink operators
       to either the raw propagator or by building on an existing quark
       propagator */
    
#ifdef KS_LEAN
    oldip0 = -1;
    oldiq0 = -1;
    oldiq1 = -1;
#endif
    for(j=0; j<param.num_qk; j++){
      STARTTIME;
      i = param.prop_for_qk[j];
      
      if(param.parent_type[j] == PROP_TYPE){
#ifdef KS_LEAN
	/* Restore clover prop[i] from file. */
	/* But first destroy the old one, unless we still need it */
	
	/* If we saved the old prop to a file, we can safely free it */
	if(oldip0 >= 0 && oldip0 != i)
	  if(param.saveflag_ks[oldip0] != FORGET){
	    destroy_ksp_field(prop[oldip0]);  prop[oldip0] = NULL;
	    node0_printf("destroy prop[%d]\n",oldip0);
	  }
	
	/* If we saved the old quark 0 to a file, we can safely free it */
	if(oldiq0 >= 0)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }
	
	/* If we saved the old quark 1 to a file, we can safely free it */
	if(oldiq1 >= 0)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }
	
	/* Fetch the prop we currently need from a file */
	if(prop[i] == NULL)
	  prop[i] = reread_ksprop_to_ksp_field(param.saveflag_ks[i], 
					       param.savefile_ks[i],
					       prop_nc[i]);
#endif
	/* Apply sink operator quark[j] <- Op[j] prop[i] */
	quark[j] = create_ksp_field_copy(prop[i]);
	ksp_sink_op(&param.snk_qs_op[j], quark[j]);
#ifdef KS_LEAN
	quark_nc[j] = quark[j]->nc;
	oldip0 = i;
	oldiq0 = -1;
#endif
	
	/* Can we delete any props now? */
	/* For each prop, scan ahead to see if it is no longer needed. */
	for(i = 0; i <= param.end_prop[param.num_set-1]; i++)
	  if(prop[i] != NULL){
	    int wont_need = 1;
	    for(k = j + 1; k < param.num_qk; k++)
	      if(param.parent_type[k] == PROP_TYPE &&
		 param.prop_for_qk[k] == i)
		wont_need = 0;  /* Still need this one */
	    if(wont_need){
	      destroy_ksp_field(prop[i]); prop[i] = NULL;
	      node0_printf("destroy prop[%d]\n",i);
	    }
	  }
      }
      else if(param.parent_type[j] == QUARK_TYPE) { /* QUARK_TYPE */
#ifdef KS_LEAN
	/* Restore quark[i] from file */
	/* But first destroy the old ones, unless we still need one of them */
	
	/* In this case we won't need the old prop */
	if(oldip0 >= 0)
	  if(param.saveflag_ks[oldip0] != FORGET){
	    destroy_ksp_field(prop[oldip0]); prop[oldip0] = NULL;
	    node0_printf("destroy prop[%d]\n",oldip0);
	  }
	
	if(oldiq0 >= 0 && oldiq0 != i)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }
	
	if(oldiq1 >= 0 && oldiq1 != i)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }
	
	if(quark[i] == NULL)
	  quark[i] = reread_ksprop_to_ksp_field(param.saveflag_q[i], 
						param.savefile_q[i],
						quark_nc[i]);
	
#endif
	/* Apply sink operator quark[j] <- Op[j] quark[i] */
	quark[j] = create_ksp_field_copy(quark[i]);
	ksp_sink_op(&param.snk_qs_op[j], quark[j]);
#ifdef KS_LEAN
	quark_nc[j] = quark[j]->nc;
	oldip0 = -1;
	oldiq0 = i;
#endif
      } else { /* COMBO_TYPE */
	int k;
	int nc = quark[param.combo_qk_index[j][0]]->nc;
	/* Create a zero field */
	quark[j] = create_ksp_field(nc);
	/* Compute the requested linear combination */
	for(k = 0; k < param.num_combo[j]; k++){
	  ks_prop_field *q = quark[param.combo_qk_index[j][k]];
	  if(nc != q->nc){
	    printf("Error: Attempting to combine an inconsistent number of colors: %d != %d\n",nc, q->nc);
	    terminate(1);
	  }
	  scalar_mult_add_ksprop_field(quark[j], q, param.combo_coeff[j][k], quark[j]);
	}
      }      
      /* Save the resulting quark[j] if requested */
      dump_ksprop_from_ksp_field( param.saveflag_q[j], 
				  param.savefile_q[j], quark[j]);
#ifdef KS_LEAN
      oldiq1 = j;
#endif
      ENDTIME("generate sink operator");
    }
#ifdef KS_LEAN
    /* Free remaining memory */
    if(oldip0 >= 0)
      if(param.saveflag_ks[oldip0] != FORGET){
	destroy_ksp_field(prop[oldip0]); prop[oldip0] = NULL;
	node0_printf("destroy prop[%d]\n",oldip0);
      }
    
    if(oldiq0 >= 0)
      if(param.saveflag_q[oldiq0] != FORGET){
	destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    
    if(oldiq1 >= 0)
      if(param.saveflag_q[oldiq1] != FORGET){
	destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif

#if defined(HAVE_QUDA) && defined(USE_GSMEAR_GPU)
    // delete 2-link current used for smearing
    gauss_smear_delete_2link_QUDA();
#endif
    
    /* Now destroy all remaining propagator fields */
    
    if(param.num_set > 0){
      for(i = 0; i <= param.end_prop[param.num_set-1]; i++){
	if(prop[i] != NULL){
	  node0_printf("destroy prop[%d]\n",i);
	  destroy_ksp_field(prop[i]);
	  prop[i] = NULL;
	}
      }
    }

    // also destroy the source fields here
    for(is=0; is<param.num_base_source+param.num_modified_source; is++){
      if(source[is] != NULL)node0_printf("destroy source[%d]\n",is);
      destroy_ksp_field(source[is]); source[is] = NULL;
    }
    
    
    
    /****************************************************************/
    /* Compute the meson propagators */
    
    STARTTIME;
    for(i = 0; i < param.num_pair; i++){
      
      /* Index for the quarks making up this meson */
      iq0 = param.qkpair[i][0];
      iq1 = param.qkpair[i][1];
      
      /* Naik indices for the quarks */
      naik_index0 = param.naik_index[iq0];
      naik_index1 = param.naik_index[iq1];
      
      node0_printf("Mesons for quarks %d and %d\n",iq0,iq1);
      
#ifdef KS_LEAN
      /* Restore quarks from file and free old memory */
      /* We try to reuse props that are already in memory, so we don't
         destroy them immediately, but wait to see if we need
         them again for the next pair. */
      if(i > 0 && oldiq0 != iq0 && oldiq0 != iq1)
	if(param.saveflag_q[oldiq0] != FORGET){
	  destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq0);
	}
      
      if(i > 0 && oldiq1 != iq0 && oldiq1 != iq1)
	if(param.saveflag_q[oldiq1] != FORGET){
	  destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq1);
	}
      
      if(quark[iq0] == NULL)
	quark[iq0] = 
	  reread_ksprop_to_ksp_field(param.saveflag_q[iq0], 
				     param.savefile_q[iq0],
				     quark_nc[iq0]);
      
      if(quark[iq1] == NULL){
	quark[iq1] = 
	  reread_ksprop_to_ksp_field(param.saveflag_q[iq1], 
				     param.savefile_q[iq1],
				     quark_nc[iq1]);
      }
#endif
      
      /* Tie together to generate hadron spectrum */
      spectrum_ks(quark[iq0], naik_index0, quark[iq1], naik_index1, i);
      
#ifdef KS_LEAN
      /* Remember, in case we need to free memory */
      oldiq0 = iq0;
      oldiq1 = iq1;
#endif
    }
#ifdef KS_LEAN
    /* Free any remaining quark prop memory */
    if(quark[oldiq0] != NULL)
      if(param.saveflag_q[oldiq0] != FORGET){
	destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    if(quark[oldiq1] != NULL)
      if(param.saveflag_q[oldiq1] != FORGET){
	destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif
    ENDTIME("tie meson correlators");
    
    /****************************************************************/
    /* Compute the baryon propagators */
    
    STARTTIME;
    for(i = 0; i < param.num_triplet; i++){
      
      /* Index for the quarks making up this meson */
      iq0 = param.qktriplet[i][0];
      iq1 = param.qktriplet[i][1];
      iq2 = param.qktriplet[i][2];
      
      node0_printf("Baryons for quarks %d, %d, and %d\n",iq0,iq1,iq2);
      
#ifdef KS_LEAN
      /* Restore quarks from file and free old memory */
      /* We try to reuse props that are already in memory, so we don't
         destroy them immediately, but wait to see if we need
         them again for the next pair. */
      if(i > 0 && oldiq0 != iq0 && oldiq0 != iq1 && oldiq0 != iq2)
	if(param.saveflag_q[oldiq0] != FORGET){
	  destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq0);
	}
      
      if(i > 0 && oldiq1 != iq0 && oldiq1 != iq1 && oldiq1 != iq2)
	if(param.saveflag_q[oldiq1] != FORGET){
	  destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq1);
	}
      
      if(i > 0 && oldiq2 != iq0 && oldiq2 != iq1 && oldiq2 != iq2)
	if(param.saveflag_q[oldiq2] != FORGET){
	  destroy_ksp_field(quark[oldiq2]); quark[oldiq2] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq2);
	}
      
      if(quark[iq0] == NULL)
	quark[iq0] = 
	  reread_ksprop_to_ksp_field(param.saveflag_q[iq0], 
				     param.savefile_q[iq0],
				     quark_nc[iq0]);
      
      if(quark[iq1] == NULL)
	quark[iq1] = 
	  reread_ksprop_to_ksp_field(param.saveflag_q[iq1], 
				     param.savefile_q[iq1],
				     quark_nc[iq1]);
      if(quark[iq2] == NULL)
	quark[iq2] = 
	  reread_ksprop_to_ksp_field(param.saveflag_q[iq2], 
				     param.savefile_q[iq2],
				     quark_nc[iq2]);
#endif
      
      /* Tie together to generate hadron spectrum */
      spectrum_ks_baryon(quark[iq0], quark[iq1], quark[iq2], i);
      
#ifdef KS_LEAN
      /* Remember, in case we need to free memory */
      oldiq0 = iq0;
      oldiq1 = iq1;
      oldiq2 = iq2;
#endif
    }
#ifdef KS_LEAN
    /* Free any remaining quark prop memory */
    if(quark[oldiq0] != NULL)
      if(param.saveflag_q[oldiq0] != FORGET){
	destroy_ksp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    if(quark[oldiq1] != NULL)
      if(param.saveflag_q[oldiq1] != FORGET){
	destroy_ksp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
    if(quark[oldiq2] != NULL)
      if(param.saveflag_q[oldiq2] != FORGET){
	destroy_ksp_field(quark[oldiq2]); quark[oldiq2] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq2);
      }
#endif
    ENDTIME("tie baryon correlators");
    
#if EIGMODE == EIGCG

    Nvecs_curr = param.eigcgp.Nvecs_curr;
      
    if(param.eigcgp.Nvecs_max > 0){
      STARTTIME;
      
      resid = (double *)malloc(Nvecs_curr*sizeof(double));
      
      if(param.ks_eigen_startflag == FRESH)
	calc_eigenpairs(eigVal, eigVec, &param.eigcgp, EVEN);
      
      check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn );
      
      if(param.eigcgp.H != NULL) free(param.eigcgp.H);
      
      ENDTIME("compute eigenvectors");
    }
#endif

    if(param.eigen_param.Nvecs > 0){

      /* If using QUDA for deflation, then eigenvectors are loaded and saved directly by QUDA and not MILC */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
      STARTTIME;
      
      /* save eigenvectors if requested */
      int status = save_ks_eigen(param.ks_eigen_saveflag, param.ks_eigen_savefile,
				 Nvecs_curr, eigVal, eigVec, resid, 1);
      if(status != 0){
	node0_printf("ERROR writing eigenvectors\n");
      }
      
      /* Clean up eigen storage */
      for(i = 0; i < Nvecs_alloc; i++) free(eigVec[i]);
      free(eigVal); free(eigVec); free(resid);

      ENDTIME("save eigenvectors (if requested)");

#endif
    }

    /* Clean up quark sources, both base and modified */
    for(i = 0; i < param.num_base_source + param.num_modified_source; i++)
      clear_qs(&param.src_qs[i]);

    /* Free links */
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif

/****************************************************************/
/* Compute GB baryon propagators */

#ifdef GB_BARYON

    STARTTIME;
    int iqo0,iqo1,iqo2;
    ks_prop_field *qko0[8];
    ks_prop_field *qko1[8];
    ks_prop_field *qko2[8];
    #ifdef GB_BARYON_MMAP
        int jqo0=0,jqo1=0,jqo2=0;
        mmap_cache *tmp_cache0; /* Pointers to temporarily retain memory */
        mmap_cache *tmp_cache1;
        mmap_cache *tmp_cache2;
        node0_printf("Creating gb baryon cache container\n");
        create_gb_qk_cache(3);
    #endif
        for(i = 0; i < param.num_gb_triplet; i++){

          /* Index for the quarks making up this gb baryon */
          iqo0 = param.qk8triplet[i][0];
          iqo1 = param.qk8triplet[i][1];
          iqo2 = param.qk8triplet[i][2];
          node0_printf("Golterman-Bailey baryon for quark octets %d, %d, and %d\n",
           iqo0,iqo1,iqo2);
          node0_printf("Octet %d :  ",iqo0);
          for(j = 0; j < 8; j++){
            iq0 = param.qk_oct[iqo0][j];
            node0_printf(" %d ",iq0);
            if(iq0 == -1) qko0[j] = NULL;
            else qko0[j] = quark[iq0];
          }
          node0_printf("\nOctet %d :  ",iqo1);
          for(j = 0; j < 8; j++){
            iq1 = param.qk_oct[iqo1][j];
            node0_printf(" %d ",iq1);
            if(iq1 == -1) qko1[j] = NULL;
            else qko1[j] = quark[iq1];

          }
          node0_printf("\nOctet %d :  ",iqo2);
          for(j = 0; j < 8; j++){
            iq2 = param.qk_oct[iqo2][j];
            node0_printf(" %d ",iq2);
            if(iq2 == -1) qko2[j] = NULL;
            else qko2[j] = quark[iq2];
          }
          node0_printf("\n");

#ifdef GB_BARYON_MMAP
          node0_printf("Creating gb baryon cache\n");
	  double gbcachestart = dclock();
          if (i == 0) {
            // all new, create or copy
            create_qk_oct_cache(qko0,0,param.r_offset_gb[i],ape_links);
            if      (iqo1 == iqo0) { copy_qk_oct_cache(1,0); }
            else {
              create_qk_oct_cache(qko1,1,param.r_offset_gb[i],ape_links);
            }
            if      (iqo2 == iqo0) { copy_qk_oct_cache(2,0); }
            else if (iqo2 == iqo1) { copy_qk_oct_cache(2,1); }
            else {
              create_qk_oct_cache(qko2,2,param.r_offset_gb[i],ape_links);
            }
          } else {

            /* If cache is used again, retain memory in temporary pointers
               clear all the mmap pointers
               wipe memory for caches that are no longer necessary */
            if      (jqo0 != iqo0 && jqo0 != iqo1 && jqo0 != iqo2) {
              tmp_cache0 = NULL;
              destroy_qk_oct_cache(0);
              if (jqo1 == jqo0) { tmp_cache1 = NULL; unmap_qk_oct_cache(1); }
              if (jqo2 == jqo0) { tmp_cache2 = NULL; unmap_qk_oct_cache(2); }
            }
            else { /* Retain cache pointer, unmap caches */
              tmp_cache0 = get_qk_cache_pointer(0);
              unmap_qk_oct_cache(0);
              if (jqo1 == jqo0) { tmp_cache1 = tmp_cache0; unmap_qk_oct_cache(1); }
              if (jqo2 == jqo0) { tmp_cache2 = tmp_cache0; unmap_qk_oct_cache(2); }
            }
            if (jqo0 != jqo1) { /* Already unmapped if jqo0 == jqo1 */
              if (jqo1 != iqo0 && jqo1 != iqo1 && jqo1 != iqo2) {
                  tmp_cache1 = NULL;
                  destroy_qk_oct_cache(1);
                  if (jqo2 == jqo1) { tmp_cache2 = NULL; unmap_qk_oct_cache(2); }
              }
              else {
                tmp_cache1 = get_qk_cache_pointer(1);
                unmap_qk_oct_cache(1);
                if (jqo2 == jqo1) { tmp_cache2 = tmp_cache1; unmap_qk_oct_cache(2); }
              }
            }
            if (jqo0 != jqo2 && jqo1 != jqo2) {
              if (jqo2 != iqo0 && jqo2 != iqo1 && jqo2 != iqo2) {
                  tmp_cache2 = NULL;
                  destroy_qk_oct_cache(2);
              }
              else {
                tmp_cache2 = get_qk_cache_pointer(2);
                unmap_qk_oct_cache(2);
              }
            }

            /* Reassign cache pointers for caches that are reused */
            if      (iqo0 == jqo0) { assign_qk_cache_pointer(tmp_cache0,0); }
            else if (iqo0 == jqo1) { assign_qk_cache_pointer(tmp_cache1,0); }
            else if (iqo0 == jqo2) { assign_qk_cache_pointer(tmp_cache2,0); }
            if      (iqo1 == jqo0) { assign_qk_cache_pointer(tmp_cache0,1); }
            else if (iqo1 == jqo1) { assign_qk_cache_pointer(tmp_cache1,1); }
            else if (iqo1 == jqo2) { assign_qk_cache_pointer(tmp_cache2,1); }
            if      (iqo2 == jqo0) { assign_qk_cache_pointer(tmp_cache0,2); }
            else if (iqo2 == jqo1) { assign_qk_cache_pointer(tmp_cache1,2); }
            else if (iqo2 == jqo2) { assign_qk_cache_pointer(tmp_cache2,2); }

            /* Create new caches */
            if (iqo0 != jqo0 && iqo0 != jqo1 && iqo0 != jqo2) {
              create_qk_oct_cache(qko0,0,param.r_offset_gb[i],ape_links);
            }
            if (iqo1 != jqo0 && iqo1 != jqo1 && iqo1 != jqo2) {
              if      (iqo1 == iqo0) { copy_qk_oct_cache(1,0); }
              else {
                create_qk_oct_cache(qko1,1,param.r_offset_gb[i],ape_links);
              }
            }
            if (iqo2 != jqo0 && iqo2 != jqo1 && iqo2 != jqo2) {
              if      (iqo2 == iqo0) { copy_qk_oct_cache(2,0); }
              else if (iqo2 == iqo1) { copy_qk_oct_cache(2,1); }
              else {
                create_qk_oct_cache(qko2,2,param.r_offset_gb[i],ape_links);
              }
            }
          }

          /* Copy octet indices for next iteration */
          jqo0 = iqo0;
          jqo1 = iqo1;
          jqo2 = iqo2;
          endtime=dclock();
          node0_printf("Done creating gb baryon cache\n");
          //node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
          node0_printf("Time to create gb baryon cache %e sec\n",endtime-gbcachestart);
    #endif

          /* Tie together to generate hadron spectrum */
          spectrum_ks_gb_baryon(qko0,qko1,qko2,ape_links,i);

    #ifdef GB_BARYON_MMAP
          if (i == param.num_gb_triplet-1) { /* Done, remove all */
            destroy_qk_oct_cache(0);
            if (iqo0 != iqo1)                 { destroy_qk_oct_cache(1); }
            if (iqo0 != iqo2 && iqo1 != iqo2) { destroy_qk_oct_cache(2); }
          }
    #endif
        }
    #ifdef GB_BARYON_MMAP
        node0_printf("Destroying gb baryon cache container\n");
        destroy_gb_qk_cache();
    #endif
    ENDTIME("tie gb baryon correlators");
    endtime=dclock();

    node0_printf("GB BARYON COMPLETED\n");
    //node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
#endif 

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();
    
    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    starttime = endtime; /* In case we continue looping over readin */
    node0_printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
    printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
    fflush(stdout);
    
    for(i = 0; i < param.num_qk; i++){
      if(quark[i] != NULL)node0_printf("destroy quark[%d]\n",i);
      destroy_ksp_field(quark[i]); quark[i] = NULL;
    }
    
    destroy_ape_links_4D(ape_links);

#if ! defined(HAVE_QUDA) && defined(GAUSS_SMEAR_KS_TWOLINK)
    gauss_smear_delete_2link_cpu();
#endif
    
    /* Destroy fermion links (created in readin() */
    
    fn_links = NULL;
    starttime = endtime;
  } /* readin(prompt) */
  
  free_lattice();

#ifdef HAVE_QUDA
  finalize_quda();
#endif
  
#ifdef HAVE_QPHIX
  finalize_qphix();
#endif

#ifdef HAVE_GRID
  finalize_grid();
#endif

  normal_exit(0);
  return 0;
}
