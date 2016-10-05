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
#endif
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

int main(int argc, char *argv[])
{
  int prompt;
  int i, j, k, iq0, iq1, iq2;
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
  ks_prop_field *prop[MAX_PROP];
  ks_prop_field *quark[MAX_QK];
  int prop_nc[MAX_PROP];
#if EIGMODE == EIGCG || EIGMODE == DEFLATION
  int Nvecs_curr;
  double *resid = NULL;
  imp_ferm_links_t **fn;
#endif
  
  initialize_machine(&argc,&argv);

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
#ifdef HYPISQ_SVD_COUNTER
    hypisq_svd_counter = 0;
#endif
    
#if EIGMODE == DEFLATION
    /**************************************************************/
    /* Compute Dirac eigenpairs           */

    STARTTIME;

    active_parity = EVEN;
    fn = get_fm_links(fn_links);
    Nvecs_curr = Nvecs_tot = param.Nvecs;

    /* compute eigenpairs if requested */
    if(param.ks_eigen_startflag == FRESH){
      int total_R_iters;
      total_R_iters=Kalkreuter(eigVec, eigVal, param.eigenval_tol, param.error_decr,
			       Nvecs_curr, param.MaxIter, param.Restart, param.Kiters, 1);
      node0_printf("total Rayleigh iters = %d\n", total_R_iters);

#if 0 /* If needed for debugging */
      /* (The Kalkreuter routine uses the random number generator to
	 initialize the eigenvector search, so, if you want to compare
	 first results with and without deflation, you need to
	 re-initialize here.) */
      initialize_site_prn_from_seed(iseed);
#endif
    }
    
    /* Calculate and print the residues and norms of the eigenvectors */
    resid = (double *)malloc(Nvecs_curr*sizeof(double));
    check_eigres( resid, eigVec, eigVal, Nvecs_curr, EVEN, fn[0] );

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

    ENDTIME("calculate Dirac eigenpairs");

#endif
    
    /**************************************************************/
    /* Compute chiral condensate and related quantities           */

    STARTTIME;

    /* Make fermion links if not already done */

    for(i = 0; i < param.num_pbp_masses; i++){
#ifdef U1_FIELD
      u1phase_on(param.charge_pbp[i], u1_A);
      invalidate_fermion_links(fn_links);
#endif
      restore_fermion_links_from_site(fn_links, param.qic_pbp[i].prec);
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
#endif
    }

    ENDTIME("calculate pbp, etc");

    /**************************************************************/
    /* Fix the gauge */
    
    if( param.fixflag == COULOMB_GAUGE_FIX)
      {
	if(this_node == 0) 
	  printf("Fixing to Coulomb gauge\n");

	rephase( OFF );
	STARTTIME;
	gaugefix(TUP,(Real)1.8,500,GAUGE_FIX_TOL);
	//gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
	ENDTIME("gauge fix");

	/* (Re)construct APE smeared links after gauge fixing.  
	   No KS phases here! */
	destroy_ape_links_4D(ape_links);
	ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
	apply_apbc( ape_links );

	rephase( ON );
	invalidate_fermion_links(fn_links);

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
#endif


    if(this_node==0)printf("END OF HEADER\n");
    

    /**************************************************************/


    /* Loop over sets of propagators */

    STARTTIME;
    for(k=0; k<param.num_set; k++){
      int num_prop = param.end_prop[k] - param.begin_prop[k] + 1;
      int i0 = param.begin_prop[k];

      for(i=param.begin_prop[k]; i <= param.end_prop[k]; i++){
      
	/**************************************************************/
	/* Read and/or generate quark propagator */
	
	prop_nc[i] = param.src_qs[k].ncolor;
	prop[i] = create_ksp_field(prop_nc[i]);
	
	if(prop[i] == NULL){
	  printf("main(%d): No room for prop\n",this_node);
	  terminate(1);
	}
	
	node0_printf("Mass= %g source %s ",
		     (double)param.ksp[i].mass,
		     param.src_qs[k].descrp);
#ifdef U1_FIELD
	node0_printf("Q %g ",param.charge[k]);
#endif
	node0_printf("residue= %g rel= %g\n",
		     (double)param.qic[i].resid,
		     (double)param.qic[i].relresid);
	
      }
      
      /* We pass the beginning addresses of the set data */
      
      total_iters += solve_ksprop(num_prop,
				  param.startflag_ks + i0,
				  param.startfile_ks + i0,
				  param.saveflag_ks + i0,
				  param.savefile_ks + i0,
				  prop + i0,
				  &param.src_qs[k], 
				  param.qic + i0, 
				  param.ksp + i0,
				  param.charge[k],
				  param.bdry_phase[i0],
				  param.coord_origin,
				  param.check[i0]);
      
      clear_qs(&param.src_qs[k]);
      
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
	    free_ksp_field(prop[oldip0]);  prop[oldip0] = NULL;
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
    
    /* Now destroy all remaining propagator fields */
    
    for(i = 0; i <= param.end_prop[param.num_set-1]; i++){
      if(prop[i] != NULL){
	node0_printf("destroy prop[%d]\n",i);
	destroy_ksp_field(prop[i]);
	prop[i] = NULL;
      }
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

#if EIGMODE == EIGCG || EIGMODE == DEFLATION

    STARTTIME;

    /* save eigenvectors if requested */
    int status = save_ks_eigen(param.ks_eigen_saveflag, param.ks_eigen_savefile,
			       Nvecs_curr, eigVal, eigVec, resid, 1);
    if(status != 0){
      node0_printf("ERROR writing eigenvectors\n");
    }

    /* Clean up eigen storage */
    for(i = 0; i < Nvecs_tot; i++) free(eigVec[i]);
    free(eigVal); free(eigVec); free(resid);

    ENDTIME("save eigenvectors (if requested)");

#endif

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();
    
    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
    printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
    printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
    fflush(stdout);
    
    for(i = 0; i < param.num_qk; i++){
      if(quark[i] != NULL)node0_printf("destroy quark[%d]\n",i);
      destroy_ksp_field(quark[i]); quark[i] = NULL;
    }
    
    destroy_ape_links_3D(ape_links);
    
    /* Destroy fermion links (created in readin() */
    
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
    destroy_fermion_links_hypisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
  } /* readin(prompt) */
  

#ifdef HAVE_QUDA
  qudaFinalize();
#endif
  
#ifdef HAVE_QPHIX
  finalize_qphix();
#endif

  normal_exit(0);
  return 0;
}
