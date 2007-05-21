/************* calc_hqet_light_form_factor.c **************************/
/* MIMD version 4  **/

/*   Routine to calculate the heavy to light form factors,
 *   using HQET.
 *
 *  Subroutine arguments
 *
 *   GLOBAL VARIABLES (used in the site structure)
 *
 *    variable name              variable type
 *--------------------------------------------------
 *    quark_spectate        ||   spin_wilson_vector
 *    quark_zonked          ||   spin_wilson_vector
 *    quark_zonked_rot      ||   spin_wilson_vector
 *    quark_sequential      ||   spin_wilson_vector
 *    mp                    ||   wilson_vector
 *    chi                   ||   wilson_vector
 *    heavy_prop            ||   su3_matrix
 *
*/

/* Modifications
   C. McNeile 1997  Original version 
   C. DeTar 5/24/97 Added SW rotatation of zonked clover propagator
 */


#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif

/***********************/

void calc_hqet_light_form(void)
{
  int color , spin ; 
  int v_pt = 0 ; 
  int k_spectator = 0, k_zonked_light = 0  ;
  field_offset rv;

  complex *hl_corr ;      /*** The heavy-light correlators *****/
  complex *two_corr ;      /***  light-light two point function ****/
  complex *seq_corr ;      /***  hqet-light two point function ****/

  int   hl_corr_dim = 0 ; 
  int   two_corr_dim = 0 ; 
  int   seq_corr_dim = 0 ; 

  w_prop_file *spectate_fp_in , *zonked_fp_in ; 

  double t_fin ; 
  double t_qload ; 

  int restart_flag_zonked ; 
  int restart_flag_spectate ; 
  int MinCG;

  /*************************************************************************/

  setup_hlcorr(&hl_corr ,  &hl_corr_dim ) ;
  setup_twocorr(&two_corr,  &two_corr_dim ) ;
  setup_seq_corr(&seq_corr, &seq_corr_dim ) ;


  /*** calculate the correlators *****/

    for(k_spectator = 0 ; k_spectator < no_spectator ; ++k_spectator)
      for(k_zonked_light=0 ; k_zonked_light < no_zonked_light ; ++k_zonked_light)
      {

	if( startflag_spectator[k_spectator] == FRESH )
	  restart_flag_spectate = 0 ;
	else
	  restart_flag_spectate = 1 ;   /**** MORE WORK*******/
	
	
	if( startflag_zonked[k_zonked_light] == FRESH )
	  restart_flag_zonked = 0 ;
	else
	  restart_flag_zonked  = 1 ;   /**** MORE WORK*******/


	/*** open the spectator light quark file *****/
	kappa = kappa_spectator[k_spectator]  ;
	
	spectate_fp_in = r_open_wprop(startflag_spectator[k_spectator], 
				     qfile_spectator[k_spectator]);

	/** open the light quark zonked propagator ***/
	kappa =   kappa_zonked_light[k_spectator] ;

	zonked_fp_in = r_open_wprop(startflag_zonked[ k_zonked_light ],
				   qfile_zonked[ k_zonked_light ]);

	for(color=0 ; color < 3 ; ++color)
	{
	  if( this_node == 0 )
	  {
	    printf("**************************************************\n"); 
	    printf("Starting kappa_spectator = %f  kappa_zonked = %f colour = %d\n",
		   kappa_spectator[k_spectator] , kappa_zonked_light[k_spectator]  , color); 
	    printf("**************************************************\n"); 
	  }

	  t_qload = -dclock(); 

	  for(spin = 0 ; spin < 4 ; ++spin)
	  {

	    IF_MASTER printf("------> starting to load light spectator k = %f <-----\n",kappa_spectator[k_spectator]);
	    /*** Load in the spectator quark propagator ****/
	    kappa = kappa_spectator[k_spectator]  ;
	    if(reload_wprop_sc_to_site( startflag_spectator[k_spectator],
			      spectate_fp_in, spin, color, 
			      F_OFFSET(quark_spectate.d[spin]),1)!=0)
	      terminate(1);

	    /**** check the wilson vector loaded in , by using it as a
                  new solution to MR *****/

	    /* Complete the definition of source structure */
	    wqs_spectator[k_spectator].color = color;
	    wqs_spectator[k_spectator].spin = spin;

	    /* For clover_info if we ever use it */
	    wqstmp = wqs_spectator[k_spectator];
	    
	    /* If we are starting afresh, we set a minimum number
	       of iterations */
	    if(startflag_spectator[k_spectator] == FRESH)MinCG = nt/2; 
	    else MinCG = 0;

	    /* Load inversion control structure */
	    qic_spectator.prec = PRECISION;
	    qic_spectator.min = MinCG;
	    qic_spectator.max = niter_spectator;
	    qic_spectator.nrestart = nrestart_spectator;
	    qic_spectator.resid = resid_spectator;
	    qic_spectator.start_flag = restart_flag_spectate;

#ifdef BICG_CLOVER
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa_spectator[k_spectator];
	    dcp.Clov_c = clov_c;
	    dcp.U0 = u0;

	    rv  = F_OFFSET(tmpb);
	    wilson_invert_site_wqs(F_OFFSET(chi), F_OFFSET(quark_spectate.d[spin]),
			       w_source,&wqs_spectator[k_spectator],
			       bicgilu_cl_site,&qic_spectator,(void *)&dcp);

#else
	    /* Load Dirac matrix parameters */
	    dwp.Kappa = kappa_spectator[k_spectator];
	
	    wilson_invert_site_wqs(F_OFFSET(chi), F_OFFSET(quark_spectate.d[spin]),
			       w_source,&wqs_spectator[k_spectator],
			       mrilu_w_site,&qic_spectator,(void *)&dwp);
#endif

	    IF_MASTER printf("------> starting to load light zonked quark k = %f <-----\n",kappa_zonked_light[k_spectator]);

	    /*** load the light zonked quark propagagor from disk ***/
	    kappa =   kappa_zonked_light[ k_zonked_light ] ;

	    if(reload_wprop_sc_to_site( startflag_zonked[k_zonked_light],
			      zonked_fp_in, spin, color, 
			      F_OFFSET(quark_zonked.d[spin]), 1)!=0)
	      terminate(1);

	    /**** check the wilson vector loaded in , by using it as a
                  new solution to MR *****/

	    /* Complete the definition of source structure */
	    wqs_zonked_light[k_zonked_light].color = color;
	    wqs_zonked_light[k_zonked_light].spin = spin;

	    /* For clover_info if we ever use it */
	    wqstmp = wqs_zonked_light[k_zonked_light];

	    /* If we are starting afresh, we set a minimum number
	       of iterations */
	    if(startflag_zonked[k_zonked_light] == FRESH)MinCG = nt/2; 
	    else MinCG = 0;

	    /* Load inversion control structure */
	    qic_zonked_light.prec = PRECISION;
	    qic_zonked_light.min = MinCG;
	    qic_zonked_light.max = niter_zonked;
	    qic_zonked_light.nrestart = nrestart_zonked;
	    qic_zonked_light.resid = resid_zonked;
	    qic_zonked_light.start_flag = startflag_zonked[k_zonked_light];

#ifdef BICG_CLOVER
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa_zonked_light[k_zonked_light];
	    dcp.Clov_c = clov_c;
	    dcp.U0 = u0;
	    
	    wilson_invert_site_wqs(F_OFFSET(chi), F_OFFSET(quark_zonked.d[spin]),
			       w_source,&wqs_zonked_light[k_zonked_light],
			       bicgilu_cl_site,&qic_zonked_light,(void *)&dcp);

#else
	    /* Load Dirac matrix parameters */
	    dwp.Kappa = kappa_zonked_light[k_zonked_light];
	    
	    wilson_invert_site_wqs(F_OFFSET(chi), F_OFFSET(quark_zonked.d[spin]),
			       w_source,&wqs_zonked_light[k_zonked_light],
			       mrilu_w_site,&qic_zonked_light,(void *)&dwp);
			  
#endif
	  }  /*** end the loop over the spin ****/


	  t_qload += dclock(); 
	  IF_VERBOSE_ON(1)
	    printf("Time to load quark props for color = %d = %f sec\n",color, t_qload) ;

#ifdef FLIP_Q_SOURCE_REP
	    /***** flip the source representation *******/
	    flip_source_re( F_OFFSET(quark_zonked)) ;    
	    flip_source_re( F_OFFSET(quark_spectate));  
#endif




	  clover_rotate( F_OFFSET(quark_zonked), F_OFFSET(quark_zonked_rot) );

	  contract_light_twopt(two_corr, F_OFFSET(quark_zonked ),
			       F_OFFSET(quark_spectate),
			       k_zonked_light, k_spectator )  ; 


	  /****-----  calculate the three point function  -----*****/
	  for(v_pt = 0 ; v_pt < novel ; ++v_pt)
	  {

	    /*** because of teh smearing, quark_spectate gets overwritten for each velocity **/
	    copy_lattice_spin_wilson_vector(F_OFFSET(quark_sequential) , F_OFFSET(quark_spectate) ) ;


	    /*** smear the light quark and calculate the HQET-light correlators ****/
	    twopt_sequential_corr(seq_corr, F_OFFSET(heavy_prop),v_pt, k_spectator, color) ; 

	    for(spin = 0 ; spin < 4 ; ++spin )
	    {

	      /*** do the sequential source HQET inversion ****/
	      smeared_sequential_source(F_OFFSET(quark_sequential.d[spin] ),
					F_OFFSET(heavy_prop),tf, v_pt) ;   

	    } /*** end the loop over spin *****/

	    /*** tie the propagators together *****/
	    contract_hqet_to_light(hl_corr, F_OFFSET(quark_zonked ),
				   F_OFFSET(quark_zonked_rot ),
				   F_OFFSET(quark_sequential),
				   v_pt, k_zonked_light, k_spectator )  ; 


	  }  /*** end of the loop over the velocity  ***/



	}  /*** end the loop over the source colour *****/

	/*** close the spectator light quark file *****/
	r_close_wprop(startflag_spectator[k_spectator],spectate_fp_in);

	/** close light quark zonked propagator ***/
	r_close_wprop(startflag_zonked[ k_zonked_light ], zonked_fp_in);
	
      } /*** end of the loop over k_zonked_light, k_spectator  ****/



  /** write the data to disk ****/
  t_fin = dclock() ; 

  finish_hlcorr(hl_corr, hl_corr_dim ) ;    
  free(hl_corr) ;  


  finish_twocorr(two_corr, two_corr_dim ) ; 
  free(two_corr) ; 

  finish_seqcorr(seq_corr, seq_corr_dim ) ; 
  free(seq_corr) ; 

  IF_VERBOSE_ON(1)
    printf("calc_hqet_light_form::Time to save all the correlators = %g sec\n",dclock() - t_fin) ;


}


