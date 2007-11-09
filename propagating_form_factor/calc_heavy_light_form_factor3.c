/*********************** calc_heavy_light_form_factor3.c **************/
/* MIMD version 6 */
/*  Routine to calculate all the required correlators,
 *  for both heavy to heavy and heavy to light form factors.
 * 
 *  Both the three and and two point functions are calculated.
 *
 *  There are no subroutine arguments
 *
 *   GLOBAL VARIABLES (used in the site structure)
 *
 *   Spin_wilson_vectors for the various quark propagators 
 *     quark_zonked
 *     quark_sequential
 *     quark_spectator
 *
*/


#include "prop_form_includes.h"
#include <string.h>

#ifdef DEBUGDEF
#include "debug_form.h"
#endif


#ifdef LARGE_MEMORY_SMEAR
#define M_SINK_SMEAR(QUARK,S_FUNC) \
	sink_smear_light_pos_quark_large(F_OFFSET(QUARK), \
				   F_OFFSET(large_fft_one), F_OFFSET(large_fft_two), F_OFFSET(S_FUNC) ) 
#else
#define M_SINK_SMEAR(QUARK,S_FUNC) \
      for(spin = 0 ; spin < 4 ; ++spin)\
	sink_smear_light_pos_quark(F_OFFSET(QUARK.d[spin] ) , \
				   F_OFFSET(mp), F_OFFSET(chi), F_OFFSET(S_FUNC) ) 
#endif

/**********************************************************************/

int file_exists_broadcast(char *file)
{
  FILE *dummy_fp;
  int non_null_fp;

  /* Test whether file exists and is openable and broadcast result */
  if(this_node==0){
    dummy_fp = fopen(file,"rb");
    non_null_fp = (dummy_fp != NULL);
    if(non_null_fp)fclose(dummy_fp);
  }
  broadcast_bytes((char *)&non_null_fp,sizeof(int));
  return non_null_fp;
}

/**********************************************************************/

int find_matching_zonked_light(int k_spectator)
{
  int match;
  int k_zonked_light;
  
  match = 0;
  for(k_zonked_light=0 ; k_zonked_light < no_zonked_light;
      k_zonked_light++)
    {
      /* Identical means, same kappa, same source
	 specifications */
      if(kappa_spectator[k_spectator] == 
	 kappa_zonked_light[k_zonked_light] &&
	 
	 wqs_spectator[k_spectator].type 
	 == wqs_zonked_light[k_zonked_light].type &&
	 
	 wqs_spectator[k_spectator].r0 
	 == wqs_zonked_light[k_zonked_light].r0 &&
	 
	 wqs_spectator[k_spectator].x0 
	 == wqs_zonked_light[k_zonked_light].x0 &&
	 
	 wqs_spectator[k_spectator].y0 
	 == wqs_zonked_light[k_zonked_light].y0 &&
	 
	 wqs_spectator[k_spectator].z0 
	 == wqs_zonked_light[k_zonked_light].z0 &&
	 
	 wqs_spectator[k_spectator].t0 
	 == wqs_zonked_light[k_zonked_light].t0)
	{
	  match = 1;
	  break;
	}
    }
  if(match == 1)
    return k_zonked_light;
  else
    return -1;

} /** end of find_matching_zonked_light **/

/**********************************************************************/

void load_in_spectator(int color, int spin, int k_spectator,
		       field_offset dest)
{
  int MinCG;
  int restart_flag_spectator ; 

  w_prop_file *spectator_fp_in ; /*** Quark propagator IO stuff **/
    
  node0_printf("Loading spectator kappa = %f\n",
	       kappa_spectator[k_spectator]);
  fflush(stdout);
  
  if( startflag_spectator[k_spectator] == FRESH )
    restart_flag_spectator = 0 ;
  else
    restart_flag_spectator = 1 ;   
  
  
  /*** open the spectator light quark file *****/
  kappa = kappa_spectator[k_spectator] ;
  
  spectator_fp_in = r_open_wprop(startflag_spectator[k_spectator], 
				qfile_spectator[k_spectator]);
  /*** Load in the spectator quark propagator ****/
  if(reload_wprop_sc_to_site(startflag_spectator[k_spectator],
		       spectator_fp_in, spin, color, dest, 1)!=0)
    terminate(1);
  
  /**** check the wilson vector loaded in, 
	by using it as a new solution to the inverter *****/
  /* Complete the definition of source structure */
  wqs_spectator[k_spectator].color = color;
  wqs_spectator[k_spectator].spin = spin;
  
  kappa = kappa_spectator[k_spectator] ;
  
  /* If we are starting afresh, we set a minimum number
     of iterations */
  if(startflag_spectator[k_spectator] == FRESH)MinCG = nt; 
  else MinCG = 0;
  
  /* Load inversion control structure */
  qic_spectator.prec = PRECISION;
  qic_spectator.min = MinCG;
  qic_spectator.max = niter_spectator;
  qic_spectator.nrestart = nrestart_spectator;
  qic_spectator.resid = resid_spectator;
  qic_spectator.start_flag = restart_flag_spectator;
  
#ifdef CLOVER
  /* Load Dirac matrix parameters */
  dcp.Kappa = kappa_spectator[k_spectator];
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  
  wilson_invert_site_wqs(F_OFFSET(chi), dest,
		     w_source,&wqs_spectator[k_spectator],
		     bicgilu_cl_site,&qic_spectator,(void *)&dcp);
  
#else
  /* Load Dirac matrix parameters */
  dwp.Kappa = kappa_spectator[k_spectator];
  
  wilson_invert_site_wqs(F_OFFSET(chi), dest,
		     w_source,&wqs_spectator[k_spectator],
		     mrilu_w_site,&qic_spectator,(void *)&dwp);
#endif

  /*** close the spectator light quark file *****/
  r_close_wprop(startflag_spectator[k_spectator],spectator_fp_in);
      
  
}  /***** end of load_in_spectator  ******/

/**********************************************************************/

void restore_local_spectator(int color, int k_spectator)
{
  int k_zonked_light;
  int spin;

  /* If the spectator quark is the same as one of the zonked light
     quarks, copy the preloaded zonked light quark instead */
  
  k_zonked_light = find_matching_zonked_light(k_spectator);
  
  if(k_zonked_light >= 0)
    {
      copy_site_spin_wilson_vector(
              F_OFFSET(quark_zonked_light[k_zonked_light]), 
	      F_OFFSET(quark_spectator)); 
    }
  else
    {
      for(spin=0; spin<4; spin++)
	load_in_spectator(color, spin, k_spectator,
			  F_OFFSET(quark_spectator.d[spin]));
    }
}

/**********************************************************************/

void restore_smeared_spectator(int color, int k_spectator)
{
  int k_zonked_light;
  int spin;

  /* If the spectator quark is the same as one of the zonked light
     quarks, copy the preloaded smeared zonked light quark instead */
  
  k_zonked_light = find_matching_zonked_light(k_spectator);
  
  if(k_zonked_light >= 0)
    {
      copy_site_spin_wilson_vector(
              F_OFFSET(quark_zonked_light[k_zonked_light]), 
	      F_OFFSET(quark_spectator)); 
    }
  else
    {
      for(spin=0; spin<4; spin++)
	load_in_spectator(color, spin, k_spectator,
			  F_OFFSET(quark_spectator.d[spin]));

      /**** sink smear the spectator_shell quark 
	    with the generic smearing function ****/
      
      node0_printf("Smearing light spectator.\n");
      fflush(stdout);
      M_SINK_SMEAR(quark_spectator,heavy_smear_func_mom[0]) ;
    }
}

/**********************************************************************/

void calc_heavy_light_form()
{
  int color , spin ; 
  int k_sequential ; 
  int p_insert ;
  int k_spectator , k_zonked_light , k_zonked_heavy  ;
  int HH3_corr_dim , HH3_corr_stride ;
  int HL3_corr_dim  , HL3_corr_stride ;

  complex *HH3_corr ;      /*** The heavy-heavy 3pt correlators *****/
  complex *HL3_corr ;      /*** The heavy-light 3pt correlators *****/

  complex *HL2_GE_corr ;      /*** The heavy-light 2pt correlators *****/
  complex *LL2_GG_corr ;      /*** The light-light 2pt correlators *****/

  complex *HL2_GG_corr ;      /*** The sequential correlators *****/
  complex *HL2_GL_corr ; /*** The local-sink correlators required for fB *****/

  int HL2_GE_corr_dim ;
  int LL2_GG_corr_dim ; 
  int HL2_GG_corr_dim ; 
  int HL2_GL_corr_dim , HL2_GL_corr_stride ; 

  w_prop_file  *zonked_light_ssink_fp[MAX_KAPPA] ; /*** Quark propagator IO stuff **/
  w_prop_file  *zonked_heavy_fp[MAX_KAPPA] ; /*** Quark propagator IO stuff **/
  w_prop_file  *zonked_heavy_ssink_fp[MAX_KAPPA] ; /*** Quark propagator IO stuff **/

  int t_source = wqs_zonked_heavy[0].t0   ; /** assume all the kappa start from the same point **/
  field_offset null_argument = 0 ; 
  int do_heavy_heavy;
  int change_mom_smear;
  int exists_zonked_light_ssink[MAX_KAPPA];
  int exists_zonked_heavy_ssink[MAX_KAPPA];
  int exists_zonked_heavy[MAX_KAPPA];

  /**********----------------------------------------**********/

  setup_w_meson_store();
  setup_HL3_corr( &HL3_corr, &HL3_corr_dim, &HL3_corr_stride );
  setup_HH3_corr( &HH3_corr, &HH3_corr_dim, &HH3_corr_stride); 

  setup_LL2_corr(&LL2_GG_corr , &LL2_GG_corr_dim );

  setup_HL2_corr(&HL2_GE_corr, &HL2_GE_corr_dim);
  setup_HL2_corr(&HL2_GG_corr, &HL2_GG_corr_dim) ; 
  setup_HL2_corr_with_rotations(&HL2_GL_corr , &HL2_GL_corr_dim , &HL2_GL_corr_stride) ; 

  set_zonked_save_intermediate();


  /*** Generate heavy quark propagators
       Gaussian smear the propagators at the sink
       Calculate the light-light and heavy-light smeared sink
       2 pt function 
       Calculate heavy-heavy local sink 2 pt function 
  *****/

  /* Figure out what needs to be done and open any needed and available
     sink-smeared LIGHT quark propagator output files */

  for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
      ++k_zonked_light)
    {
      exists_zonked_light_ssink[k_zonked_light] = 
	file_exists_broadcast(qfile_zonked_light_ssink[k_zonked_light]);

      kappa = kappa_zonked_light[k_zonked_light];
      /* If it doesn't exist and we are computing two pt functions
         that use it, we will create it, so open the file for writing */
	    
      wqstmp = wqs_zonked_light[k_zonked_light];   /* For clover_info */
      strcat(wqstmp.descrp,"; SMEARED SINK.");
      if(!exists_zonked_light_ssink[k_zonked_light] && 
	 ((saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET) ))
	zonked_light_ssink_fp[k_zonked_light]
	  = w_open_wprop(saveflag_zonked_light_ssink, 
			 qfile_zonked_light_ssink[k_zonked_light],
			 wqs_zonked_light[k_zonked_light].type); 
    }

  /* See if the heavy-heavy correlators need to be computed */
  do_heavy_heavy = !file_exists_broadcast(filename_HH2_GL);

  /* Figure out what needs to be done and open any needed and available
     sink-smeared HEAVY quark propagator output files */

  for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
      ++k_zonked_heavy)
    {
      exists_zonked_heavy[k_zonked_heavy] =
	file_exists_broadcast(qfile_zonked_heavy[k_zonked_heavy]);

      exists_zonked_heavy_ssink[k_zonked_heavy] =
	file_exists_broadcast(qfile_zonked_heavy_ssink[k_zonked_heavy]);

      /* We (re)compute the heavy_heavy 2 pt function if any of the
	 heavy local propagator files is missing */
      if(!exists_zonked_heavy[k_zonked_heavy])do_heavy_heavy = 1;

      kappa = kappa_zonked_heavy[k_zonked_heavy];
      /* If it doesn't exist, we will write, so open the file for writing */
      wqstmp = wqs_zonked_heavy[k_zonked_heavy];   /* For clover_info */
      if(!exists_zonked_heavy[k_zonked_heavy])
	zonked_heavy_fp[k_zonked_heavy] 
	  = w_open_wprop(saveflag_zonked_heavy[k_zonked_heavy], 
			 qfile_zonked_heavy[k_zonked_heavy],
			 wqs_zonked_heavy[k_zonked_heavy].type); 

      /* If it doesn't exist and we need it for the HL2_GG correlator,
         we will write, so open the file for writing */
      wqstmp = wqs_zonked_heavy[k_zonked_heavy];   /* For clover_info */
      strcat(wqstmp.descrp,"; SMEARED SINK.");
      if(!exists_zonked_heavy_ssink[k_zonked_heavy] && 
	 (saveflag_HL2_GG != FORGET))
	zonked_heavy_ssink_fp[k_zonked_heavy]
	  = w_open_wprop(saveflag_zonked_heavy_ssink, 
			 qfile_zonked_heavy_ssink[k_zonked_heavy],
			 wqs_zonked_heavy[k_zonked_heavy].type); 
    }

  /* We will compute the heavy_heavy spectrum if any one of the heavy
     local-sink propagator files is missing or if the correlator file
     doesn't exist.  Call to initialize the calculation */
  if(do_heavy_heavy)
    meson_spectrum(null_argument , t_source, 0 , no_zonked_heavy,  
		   SETUP_CORR,filename_HH2_GL) ; 

  /* Calculate the local and smeared two-point functions HL2_GL HL2_GG
     and LL2_GG */
  for(color=0 ; color < 3 ; ++color)
    {
      node0_printf("\nSTARTING SMEARED COLOR %d\n",color);
      fflush(stdout);
      /*** generate smeared light quark propagators ***/
      /* We don't load the light propagators and don't calculate the
	 smeared sink two-point functions unless a flag is set (see
	 setup_form) */
      if((saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET) ){
	
	for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
	    ++k_zonked_light)
	  {
	    /* If the smeared light quark prop file exists, use it */
	    if(exists_zonked_light_ssink[k_zonked_light])
	      {
		/* Load the smeared sink propagator */
		node0_printf("REUSING previous zonked light tmp file %s\n",
			     qfile_zonked_light_ssink[k_zonked_light]); 
		fflush(stdout);
		for(spin=0; spin < 4 ; ++spin ) 
		  load_in_zonked_light_ssink(color,spin,k_zonked_light,
			F_OFFSET(quark_zonked_light[k_zonked_light].d[spin]));
	      } /* end if exists */
	    else
	      {
		/* Load the local sink propagator and smear at sink */
		for(spin=0; spin < 4 ; ++spin ) 
		  load_in_zonked_light2(color,spin,k_zonked_light,
	              F_OFFSET(quark_zonked_light[k_zonked_light].d[spin]));
		
		/*** smear the light zonked quarks at the sink ****/
		/* NOTE: WE ARE ASSUMING THE SHELL SMEARING FUNCTIONS
		   ARE THE SAME FOR ZONKED AND SPECTATOR QUARKS */
		M_SINK_SMEAR(quark_zonked_light[k_zonked_light],  
			     heavy_smear_func_mom[0] ) ; 
		
		/** write the light ssink quark props to disk ***/
		for(spin=0; spin < 4 ; ++spin ) 
		  save_wprop_sc_from_site(saveflag_zonked_light_ssink, 
                        zonked_light_ssink_fp[k_zonked_light],
			spin, color, 
			F_OFFSET(quark_zonked_light[k_zonked_light].d[spin]),
			1) ;
	      } /* end if file doesn't exist */
	  } /* k_zonked_light */
      } /* end if (saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET) */

      /*** generate the heavy zonked local sink propagators if they
	don't already exist. We'll need them for the 3 pt functions even
	if we aren't calculating the smeared two pt functions. ***/
      
      for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
	  ++k_zonked_heavy)
	{
	  if(!exists_zonked_heavy[k_zonked_heavy])
	    {
	      for(spin = 0 ; spin < 4 ; ++spin)
		generate_heavy_zonked(color,spin,
		   kappa_zonked_heavy[k_zonked_heavy],
		   inverter_type_zonked_heavy[k_zonked_heavy],
		   F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin])) ; 

	      /*** store the heavy quark propagator to disk ****/
	      for(spin=0; spin < 4 ; ++spin ) 
		save_wprop_sc_from_site(saveflag_zonked_heavy[k_zonked_heavy], 
                        zonked_heavy_fp[k_zonked_heavy],
			spin, color, 
			F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]), 
			1) ;
	    } /* end if doesn't exist */
	  else if(do_heavy_heavy)
	    {
	      /* Load heavy zonked local from existing tmp now if we
                 need it for the heavy_heavy prop */
	      for(spin=0; spin < 4 ; ++spin ) 
		load_in_zonked_heavy_local(color,spin,k_zonked_heavy,
		     F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]));
	    } /* end else if exists but do_heavy_heavy */
	} /* k_zonked_heavy */

      /** Calculate the heavy degenerate spectrum (HH2_GL) **/
      if(do_heavy_heavy)
	{
	  for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
	      ++k_zonked_heavy)
	    meson_spectrum(F_OFFSET(quark_zonked_heavy[k_zonked_heavy]),
			   t_source ,k_zonked_heavy ,
			   no_zonked_heavy, CALCULATE_SPECTRUM,
			   filename_HH2_GL) ; 
	}

      /*** generate the smeared zonked local sink propagators if they
	don't already exist and if we need them for two point functions ***/

      if(saveflag_HL2_GG != FORGET){
	for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
	    ++k_zonked_heavy)
	  {
	    /* If the sink smeared heavy quark prop file exists, use it */
	    if(exists_zonked_heavy_ssink[k_zonked_heavy])
	      {
		node0_printf("REUSING previous zonked heavy tmp file %s\n",
			     qfile_zonked_heavy_ssink[k_zonked_heavy]); 
		fflush(stdout);
		for(spin=0; spin < 4 ; ++spin ) 
		  load_in_zonked_heavy_smear(color,spin,k_zonked_heavy,
			F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]));
	      } /* end if exists */
	    else
	      {
		/* Load local heavy quark prop unless we just computed
                   or loaded it */
		if(exists_zonked_heavy[k_zonked_heavy] && 
		   !do_heavy_heavy)
		  {
		    node0_printf("REUSING previous zonked heavy file %s\n",
			       qfile_zonked_heavy[k_zonked_heavy]); 
		    fflush(stdout);
		    for(spin=0; spin < 4 ; ++spin ) 
		      load_in_zonked_heavy_local(color,spin,k_zonked_heavy,
			F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]));
		  }

		/* smear heavy zonked quark at sink */

		M_SINK_SMEAR(quark_zonked_heavy[k_zonked_heavy], 
			     heavy_smear_func_mom[0]); 
		
		/*** store the smeared--smeared quark propagator to disk ***/
		for(spin=0; spin < 4 ; ++spin ) 
		  save_wprop_sc_from_site(saveflag_zonked_heavy_ssink, 
                     zonked_heavy_ssink_fp[k_zonked_heavy],
		     spin, color, 
		     F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]), 1) ;
	      } /* end if doesn't exist */
	    
	  } /* k_zonked_heavy */

      } /* if(saveflag_HL2_GG != FORGET) */
	
      /** Calculate the source- and sink-smeared 2 pt functions **/
      
      if((saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET)) {
	for(k_spectator = 0 ; k_spectator < no_spectator ; ++k_spectator)
	  {
	    
	    /** generate or copy sink smeared spectator quark propagator **/
	    
	    /* If the spectator quark is the same as one of the zonked light
	       quarks, copy the preloaded smeared zonked light quark instead 
	       NOTE: WE ARE ASSUMING THE SMEARING FUNCTIONS ARE THE SAME
	       FOR ZONKED AND SPECTATOR QUARKS */
	    
	    restore_smeared_spectator(color,k_spectator);
	    
	    /****-----  
	      light-light two-pt (LL2_GG) 
	      (symmetric source and sink smearing)
	      -----*****/
	    
	    /* Skip this calculation if flag is set (see setup_form) */
	    if(saveflag_LL2_GG != FORGET){
	      
	      node0_printf("Computing LL2_GG correlator\n");
	      fflush(stdout);
	      for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
		  ++k_zonked_light)
		{
		  
		  /*** calculate the light-light two point functions 
		    (LL2_GG) (shell smearing functions) *****/
		  contract_LL2(LL2_GG_corr,  
		       F_OFFSET(quark_zonked_light[k_zonked_light]) , 
		       F_OFFSET(quark_spectator ) ,
		       k_zonked_light, k_spectator)  ;
		}  /*** end of the loop over the zonked light quark reads  ***/
	      
	    } /* end of if(saveflag_LL2_GG != FORGET) */
	    
	    /****-----  
	      heavy-light two pt (HL2_GG) (shell smeared source and sink)
	      -----*****/
	    
	    /* Skip this calculation if flag is set (see setup_form) */
	    if(saveflag_HL2_GG != FORGET){
	      node0_printf("Computing HL2_GG correlator\n");
	      fflush(stdout);
	      
	      for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
		  ++k_zonked_heavy)
		{
		  
		  /*** contract the heavy-light two point function 
		    (HL2_GG) (BAG SMEARING AT THE SINK) ******/
		  contract_HL2(HL2_GG_corr ,
		      F_OFFSET(quark_zonked_heavy[k_zonked_heavy]) , 
		      F_OFFSET(quark_spectator),
		      k_zonked_heavy, k_spectator)  ; 
		  
		}  /*** end of the loop over the zonked heavy quarks  ***/
	      
	    } /* End of if(saveflag_HL2_GG != FORGET) */

	  } /* k_spectator */
	    
      } /* End of if((saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET) ) */

    } /* color */
  
  /* Write heavy-heavy spectrum */
  if(do_heavy_heavy)
    meson_spectrum(null_argument ,t_source ,0 ,no_zonked_heavy, 
		   WRITE_RESULTS,filename_HH2_GL) ; 
  
  /*** write the sink-smeared 2 pt correlators to disk ****/
  
  finish_LL2_GG_corr(LL2_GG_corr, LL2_GG_corr_dim ) ; 
  finish_HL2_GG_corr(HL2_GG_corr, HL2_GG_corr_dim ) ; 
  
  free(LL2_GG_corr); 
  free(HL2_GG_corr) ; 
  
  
  /* Close the temporary propagator output files if we opened them before */

  for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
      ++k_zonked_light)
    {
      if(!exists_zonked_light_ssink[k_zonked_light] && 
	 ((saveflag_HL2_GG != FORGET) || (saveflag_LL2_GG != FORGET) ))
	w_close_wprop(saveflag_zonked_light_ssink, 
		     zonked_light_ssink_fp[k_zonked_light]); 
    }

  for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
      ++k_zonked_heavy)
    {
      if(!exists_zonked_heavy[k_zonked_heavy])
	w_close_wprop(saveflag_zonked_heavy[k_zonked_heavy],
		     zonked_heavy_fp[k_zonked_heavy]); 
      if(!exists_zonked_heavy_ssink[k_zonked_heavy] && 
	 (saveflag_HL2_GG != FORGET))
	w_close_wprop(saveflag_zonked_heavy_ssink,
		     zonked_heavy_ssink_fp[k_zonked_heavy]); 
    }


  /*** Calculate the heavy to heavy and heavy to light 3 pt functions
       Calculate the heavy-light 2 pt functions requiring local-sink quark
       propagators ***/

  for(color=0 ; color < 3 ; ++color)
  {
    node0_printf("\nSTARTING LOCAL COLOR %d\n",color);
    fflush(stdout);

    if((saveflag_HL2_GL != FORGET) ||
       (saveflag_HL2_GE != FORGET) ||
       (saveflag_HH3 != FORGET) ||
       (saveflag_HL3 != FORGET))
      {
	/*** reload all the light zonked quark propagators for this color, 
	  smeared at the source only **/
	for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
	    ++k_zonked_light)
	  for(spin=0; spin < 4 ; ++spin ) 
	    load_in_zonked_light2(color,spin,k_zonked_light,
		    F_OFFSET(quark_zonked_light[k_zonked_light].d[spin]));
	
	/*** reload all the heavy zonked quark propagators for this color, 
	  smeared at the source only **/
	for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
	    ++k_zonked_heavy)
	  for(spin=0; spin < 4 ; ++spin ) 
	    load_in_zonked_heavy_local(color,spin,k_zonked_heavy,
		     F_OFFSET(quark_zonked_heavy[k_zonked_heavy].d[spin]));
      }
	
    for(k_spectator = 0 ; k_spectator < no_spectator ; ++k_spectator)
    {
      
      /** restore spectator quark propagator to "quark_spectator" **/
      restore_local_spectator(color, k_spectator);

      /****-----  
	   heavy-light two-pt (HL2_GL) (shell smeared sources, local sink)
	   -----*****/

      /* Skip this calculation if flag is set (see setup_form) */
      if(saveflag_HL2_GL != FORGET){
	for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
	  ++k_zonked_heavy)
	{
	  
	  node0_printf("Computing HL2_GL correlator\n");
	  fflush(stdout);

	  /*** contract the heavy-light two point function 
	       (HL2_GL) (smeared-source local-sink) */
	  contract_HL2_with_rotations(HL2_GL_corr, 
	       F_OFFSET(quark_zonked_heavy[k_zonked_heavy]) , 
	       F_OFFSET(quark_spectator),
	       F_OFFSET(quark_rot ),
	       k_zonked_heavy, k_spectator)  ;
	  
	  
	}  /*** end of the loop over the zonked heavy inversions  ***/

      } /* End of if(saveflag_HL2_GL != FORGET) */


      for(p_insert = 0 ; p_insert < no_p_values ; ++p_insert)
	{
	  node0_printf("\nStarting momentum insertion %d\n",p_insert);
	  fflush(stdout);

	  /**** sink smear the spectator quark 
		with the sequential smearing function ****/
	  /* If the smearing function did not change, retain the 
	     previous result */
	  if(p_insert == 0 || 
	     strcmp(seq_smear_file[p_insert],
		    seq_smear_file[p_insert-1]) != 0)
	    {
	      node0_printf("Sequential smearing for momentum insertion %d\n",p_insert);
	      fflush(stdout);

	      restore_local_spectator(color, k_spectator);
	      M_SINK_SMEAR(quark_spectator,seq_smear_func[ p_insert ] ) ;
	      change_mom_smear = 1;
	    }
	  else
	    {
	      node0_printf("REUSING sequentially smeared spectator from previous momentum insertion\n");
	      fflush(stdout);
	      change_mom_smear = 0;
	    }


	  /****-----  
	       heavy-light two pt (HL2_GE) (shell source, relative sink)
	       -----*****/
	  
	  /* Skip this calculation if flag is set (see setup_form) */
	  /* Note: We may, in the future, want to distiguish
	     smearing functions for different B meson momenta.
	     However, at present we have not made any provision
	     for such a distinction in the relative-smeared
	     two-point function.  So we do the calculation only
	     for the first momentum in the list! CD */

	  if((saveflag_HL2_GE != FORGET) && p_insert > 0 && change_mom_smear)
	    node0_printf("WARNING: HL2_GE correlator is computed only for the first smearing function in the list\n");
	  fflush(stdout);
	  
	  if((saveflag_HL2_GE != FORGET) && p_insert == 0){
	    
	    node0_printf("Computing HL2_GE correlator\n");
	    fflush(stdout);
	    for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
		++k_zonked_heavy)
	      {
		
		/*** contract the heavy-light two point functiion 
		  (HL2_GE) (RELATIVE SMEARING AT THE SINK) */
		contract_HL2(
		      HL2_GE_corr, 
		      F_OFFSET(quark_zonked_heavy[k_zonked_heavy]) , 
		      F_OFFSET(quark_spectator),
		      k_zonked_heavy, k_spectator)  ;
		
	      }  /*** end of the loop over the zonked heavy inversions  ***/
	  }


	  if( (saveflag_HH3 != FORGET) || (saveflag_HL3 != FORGET))
	    
	    for(k_sequential = 0 ; k_sequential < no_sequential   ; ++k_sequential  )
	      {

		copy_site_spin_wilson_vector(F_OFFSET(quark_spectator) ,
					     F_OFFSET(quark_sequential) );
		
		if( this_node == 0 ) 
		  printf("Computing form factors for k_sequential = %g k_spectator = %g p = %d,%d,%d \n",
			 kappa_sequential[k_sequential], kappa_spectator[k_spectator] , 
			 p_momstore[p_insert ][0] , p_momstore[p_insert ][1] , p_momstore[p_insert ][2] ) ; 
		fflush(stdout);
		
		
		/*** generate the sequential source propagator ****/
		for(spin = 0 ; spin < 4 ; ++spin )
		  {
		    
		    /*** do the sequential source inversion ****/
		    kappa = kappa_sequential[k_sequential] ;
		    sequential_source(
			    F_OFFSET(quark_sequential.d[spin] ),
			    F_OFFSET(psi), 
			    p_momstore[p_insert ][0] ,
			    p_momstore[p_insert ][1] ,
			    p_momstore[p_insert ][2] ,
			    tf, color,spin,
			    kappa_sequential[k_sequential] ,  
			    inverter_type_sequential[k_sequential],
			    niter_zonked_heavy,  
			    nrestart_zonked_heavy, 
			    resid_zonked_heavy, 
			    p_insert ) ;  
		  } /*** end of the loop over spin ***/
		
		
		
		/****-----  
		  heavy to light form factor 
		  -----*****/
		
		if( (saveflag_HL3 != FORGET) ) {
		  for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; 
		      ++k_zonked_light)
		    {
		      
		      /*** tie the propagators together to calculate 
			the three point function *****/
		      contract_HL3(
			     HL3_corr, 
			     F_OFFSET(quark_zonked_light[k_zonked_light]),
			     F_OFFSET(quark_sequential) ,
			     F_OFFSET(quark_rot ),
			     p_insert, k_sequential, 
			     k_zonked_light, k_spectator )  ;
		      
		      
		    }  /*** end of the loop over the zonked light
		  	 quark reads ***/
		  
		} /* End of if(saveflag_HL3 != FORGET) */
		/****-----  
		  heavy to heavy form factor 
		  -----*****/
		
		if( (saveflag_HH3 != FORGET) ) {
		  for(k_zonked_heavy=0 ; k_zonked_heavy < no_zonked_heavy ; 
		      ++k_zonked_heavy)
		    {
		      
		      /*** tie the propagators together *****/
		      contract_HH3(
			     HH3_corr, 
			     F_OFFSET(quark_zonked_heavy[k_zonked_heavy]),
			     F_OFFSET(quark_sequential  ),
			     F_OFFSET(quark_rot ),
			     p_insert, k_sequential, 
			     k_zonked_heavy, k_spectator )  ;
		      
		    }  /*** end of the loop over the zonked heavy quarks ***/
		} /* End of if(saveflag_HH3 != FORGET) */
		
	      } /*** end of the loop over  p_insert  ****/
	}   /*** end of the loop over  k_sequential   ****/
      
    } /*** end of the loop over   k_spectator  ****/
    
  }  /*** end the loop over the source colour *****/

  /*** write the correlators to disk ****/

  finish_HL3_corr(HL3_corr, HL3_corr_dim , HL3_corr_stride) ; 
  finish_HH3_corr(HH3_corr, HH3_corr_dim , HH3_corr_stride) ;

  finish_HL2_GE_corr(HL2_GE_corr, HL2_GE_corr_dim ) ; 
  finish_HL2_GL_corr(HL2_GL_corr,
      HL2_GL_corr_dim, HL2_GL_corr_stride );
  
  /****  free up the memory used in this code  *****/
  free(HH3_corr) ; 
  free(HL3_corr) ; 
  free(HL2_GE_corr); 
  free(HL2_GL_corr);


}  /***** end of calc_heavy_light_form ****/


/** clean up the local macro definition *****/
#undef M_SINK_SMEAR

/*** -- end of the file -- end of the file -- end of the file -- end of the file -- **/
