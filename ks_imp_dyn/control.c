/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int meascount,traj_done;
  int prompt;
  int s_iters, avs_iters, avspect_iters, avbcorr_iters;
  double dtime, dclock();
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

//  restore_random_state_scidac_to_site("randsave", F_OFFSET(site_prn));
//  restore_color_vector_scidac_to_site("xxx1save", F_OFFSET(xxx1),1);
//  restore_color_vector_scidac_to_site("xxx2save", F_OFFSET(xxx2),1);

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
    dtime = -dclock();
    for( traj_done=0; traj_done < warms; traj_done++ ){
      update();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avspect_iters = avs_iters = avbcorr_iters = 0;
    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
      
      /* do the trajectories */
      s_iters=update();
      
      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	rephase(OFF);
	g_measure( );
	rephase(ON);

	/* Load fat and long links for fermion measurements */
	load_ferm_links(&fn_links, &ks_act_paths);
#ifdef DM_DU0
	load_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);
#endif

	/* Measure pbp, etc */
#ifdef ONEMASS
	f_meas_imp(F_OFFSET(phi),F_OFFSET(xxx),mass, &fn_links, 
		   &fn_links_dmdu0);
#else
	f_meas_imp( F_OFFSET(phi1), F_OFFSET(xxx1), mass1, 
		    &fn_links, &fn_links_dmdu0);
	f_meas_imp( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
		    &fn_links, &fn_links_dmdu0);
#endif

	/* Measure derivatives wrto chemical potential */
#ifdef D_CHEM_POT
#ifdef ONEMASS
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass,
		  &fn_links, &fn_links_dmdu0);
#else
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass1,
		  &fn_links, &fn_links_dmdu0);
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass2,
		  &fn_links, &fn_links_dmdu0);
#endif
#endif

#ifdef SPECTRUM 
	/* Fix TUP Coulomb gauge - gauge links only*/
	rephase( OFF );
	gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
	rephase( ON );
#ifdef FN
	invalidate_ferm_links(&fn_links);
#ifdef DM_DU0
	invalidate_ferm_links(&fn_links_dmdu0);
#endif
#endif
	/* Load fat and long links for fermion measurements */
	load_ferm_links(&fn_links, &ks_act_paths);
#ifdef DM_DU0
	load_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);
#endif	
	if(strstr(spectrum_request,",spectrum,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum2(mass,F_OFFSET(phi),F_OFFSET(xxx),
				     &fn_links);
#else
	  avspect_iters += spectrum2( mass1, F_OFFSET(phi1),
				      F_OFFSET(xxx1), &fn_links);
	  avspect_iters += spectrum2( mass2, F_OFFSET(phi1),
				      F_OFFSET(xxx1), &fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_point,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_fzw(mass,F_OFFSET(phi),F_OFFSET(xxx),
					&fn_links);
#else
	  avspect_iters += spectrum_fzw( mass1, F_OFFSET(phi1),
					 F_OFFSET(xxx1), &fn_links);
	  avspect_iters += spectrum_fzw( mass2, F_OFFSET(phi1),
					 F_OFFSET(xxx1), &fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",nl_spectrum,") != NULL){
#ifdef ONEMASS
	  avspect_iters += nl_spectrum(mass,F_OFFSET(phi),F_OFFSET(xxx),
				       F_OFFSET(tempmat1),F_OFFSET(staple),
				       &fn_links);
#else
	  avspect_iters += nl_spectrum( mass1, F_OFFSET(phi1), 
		F_OFFSET(xxx1), F_OFFSET(tempmat1),F_OFFSET(staple),
					&fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_mom,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_mom(mass,mass,F_OFFSET(phi),5e-3,
					&fn_links);
#else
	  avspect_iters += spectrum_mom( mass1, mass1, 
					 F_OFFSET(phi1), 1e-1,
					 &fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_multimom,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_multimom(mass,
					     spectrum_multimom_low_mass,
					     spectrum_multimom_mass_step,
					     spectrum_multimom_nmasses,
					     5e-3, &fn_links);
#else
	  avspect_iters += spectrum_multimom(mass1,
					     spectrum_multimom_low_mass,
					     spectrum_multimom_mass_step,
					     spectrum_multimom_nmasses,
					     5e-3, &fn_links);

#endif
	}
	
#ifndef ONEMASS
	if(strstr(spectrum_request,",spectrum_nd,") != NULL){
	  avspect_iters += spectrum_nd( mass1, mass2, 1e-1,
					&fn_links);
	}
#endif
	if(strstr(spectrum_request,",spectrum_nlpi2,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_nlpi2(mass,mass,F_OFFSET(phi),5e-3,
					  &fn_links );
#else
	  avspect_iters += spectrum_nlpi2( mass1, mass1, 
					   F_OFFSET(phi1),1e-1,
					   &fn_links );
	  avspect_iters += spectrum_nlpi2( mass2, mass2, 
					   F_OFFSET(phi1),1e-1,
					   &fn_links );
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_singlets,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_singlets(mass, 5e-3, F_OFFSET(phi),
					     &fn_links);
#else
	  avspect_iters += spectrum_singlets(mass1, 5e-3, F_OFFSET(phi1),
					     &fn_links );
	  avspect_iters += spectrum_singlets(mass2, 5e-3, F_OFFSET(phi1),
					     &fn_links );
#endif
	}

	if(strstr(spectrum_request,",fpi,") != NULL)
	  {
	    avspect_iters += fpi_2( fpi_mass, fpi_nmasses, 2e-3,
				    &fn_links );
	  }
	
#ifdef HYBRIDS
	if(strstr(spectrum_request,",spectrum_hybrids,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_hybrids( mass,F_OFFSET(phi),1e-1,
					     &fn_links);
#else
	  avspect_iters += spectrum_hybrids( mass1, F_OFFSET(phi1), 5e-3,
					     &fn_links);
	  avspect_iters += spectrum_hybrids( mass2, F_OFFSET(phi1), 2e-3,
					     &fn_links);
#endif
	}
#endif
	if(strstr(spectrum_request,",hvy_pot,") != NULL){
	  rephase( OFF );
	  hvy_pot( F_OFFSET(link[XUP]) );
	  rephase( ON );
	}
#endif /* SPECTRUM */
	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
#ifdef SPECTRUM
      node0_printf("average cg iters for spectrum = %e\n",
		   (double)avspect_iters/meascount);
#endif
    }
    
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
    }
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
#ifdef HAVE_QIO
//       save_random_state_scidac_from_site("randsave", "Dummy file XML",
//        "Random number state", QIO_SINGLEFILE, F_OFFSET(site_prn));
//       save_color_vector_scidac_from_site("xxx1save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx1),1);
//       save_color_vector_scidac_from_site("xxx2save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx2),1);
#endif
    }
  }
#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  normal_exit(0);
  return 0;
}
