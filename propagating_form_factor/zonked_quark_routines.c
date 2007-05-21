/******************** zonked_quark_routines.c *******************/
/* MIMD version 6 */
/*  A collection of routines to generate and
   smear the zonked quarks.


**/


#include "prop_form_includes.h"

#ifdef DEBUGDEF
#include "debug_form.h"
#endif



#ifdef COMPILEME
/**
 **  Load in the light zonked quark propagator for a 
 **  particular kappa value.
 **
 **/


void load_in_zonked_light(int color, int k_zonked_light)
{
  int spin ; 
  int restart_flag_zonked_light ; 
  int MinCG;
  w_prop_file *zonked_fp_in ; /*** Quark propagator IO stuff **/

  /************************************************************/

      
  /*** load in the zonked quark propagator ********/
  kappa = kappa_zonked_light[k_zonked_light] ;
  
  /** open the light quark zonked propagator ***/
  
  zonked_fp_in = r_open_wprop(startflag_zonked_light[ k_zonked_light ],
			     qfile_zonked_light[ k_zonked_light ]);

  node0_printf("Loading from %s\n",qfile_zonked_light[ k_zonked_light ]);

  for(spin = 0 ; spin < 4 ; ++spin )
  {
    /*** load the light zonked quark propagagor from disk ***/
    if(reload_wprop_sc_to_site( startflag_zonked_light[k_zonked_light],
			 zonked_fp_in, spin, color, 
			 F_OFFSET(quark_zonked.d[spin]), 1)!=0)
      terminate(1);
    
    /**** check the wilson vector loaded in , by using it
      as a new solution to the inverter *****/
    /* Complete the definition of source structure */
    wqs_zonked_light[k_zonked_light].color = color;
    wqs_zonked_light[k_zonked_light].spin = spin;

    /* For clover_info if we ever use it */
    wqstmp = wqs_zonked_light[k_zonked_light];
    
    /* If we are starting afresh, we set a minimum number
       of iterations */
    if(startflag_zonked_light[k_zonked_light] == FRESH)MinCG = nt; 
    else MinCG = 0;

    /* Load inversion control structure */
    qic_zonked_light.prec = PRECISION;
    qic_zonked_light.min = MinCG;
    qic_zonked_light.max = niter_zonked_light;
    qic_zonked_light.nrestart = nrestart_zonked_light;
    qic_zonked_light.resid = resid_zonked_light;
    qic_zonked_light.start_flag = startflag_zonked_light[k_zonked_light];

#ifdef CLOVER
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

  } /*** end the loop over the source spin ****/

  /** close light quark zonked propagator ***/
  r_close_wprop(startflag_zonked_light[ k_zonked_light ], zonked_fp_in);
	


}  /***** end of load_in_zonked_light  ******/
#endif  /* COMPILEME */

void load_in_zonked_light2(int color, int spin, int k_zonked_light,
			   field_offset dest)
{
  int MinCG;
  w_prop_file *zonked_fp_in ; /*** Quark propagator IO stuff **/

  /************************************************************/

      
  /*** load in the zonked quark propagator ********/
  kappa = kappa_zonked_light[k_zonked_light] ;
  
  /** open the light quark zonked propagator ***/
  
  zonked_fp_in = r_open_wprop(startflag_zonked_light[ k_zonked_light ],
			     qfile_zonked_light[ k_zonked_light ]);

  /*** load the light zonked quark propagagor from disk ***/
  node0_printf("Loading from %s\n",qfile_zonked_light[ k_zonked_light ]);
  if(reload_wprop_sc_to_site( startflag_zonked_light[k_zonked_light],
			zonked_fp_in, spin, color, 
			dest, 1)!=0)
    terminate(1);
  
  /**** check the wilson vector loaded in , by using it
	as a new solution to the inverter *****/
  /* Complete the definition of source structure */
  wqs_zonked_light[k_zonked_light].color = color;
  wqs_zonked_light[k_zonked_light].spin = spin;
  
  /* For clover_info if we ever use it */
  wqstmp = wqs_zonked_light[k_zonked_light];
  
  /* If we are starting afresh, we set a minimum number
     of iterations */
  if(startflag_zonked_light[k_zonked_light] == FRESH)MinCG = nt; 
  else MinCG = 0;
  
  /* Load inversion control structure */
  qic_zonked_light.prec = PRECISION;
  qic_zonked_light.min = MinCG;
  qic_zonked_light.max = niter_zonked_light;
  qic_zonked_light.nrestart = nrestart_zonked_light;
  qic_zonked_light.resid = resid_zonked_light;
  qic_zonked_light.start_flag = startflag_zonked_light[k_zonked_light];
  
#ifdef CLOVER
  /* Load Dirac matrix parameters */
  dcp.Kappa = kappa_zonked_light[k_zonked_light];
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  
  wilson_invert_site_wqs(F_OFFSET(chi), dest,
		     w_source,&wqs_zonked_light[k_zonked_light],
		     bicgilu_cl_site,&qic_zonked_light,(void *)&dcp);
#else
  /* Load Dirac matrix parameters */
  dwp.Kappa = kappa_zonked_light[k_zonked_light];
  
  wilson_invert_site_wqs(F_OFFSET(chi), dest,
		     w_source,&wqs_zonked_light[k_zonked_light],
		     mrilu_w_site,&qic_zonked_light,(void *)&dwp);
  
#endif
  
  
  /** close light quark zonked propagator ***/
  r_close_wprop(startflag_zonked_light[ k_zonked_light ], zonked_fp_in);

}  /***** end of load_in_zonked_light2  ******/


void load_in_zonked_heavy_smear(int color, int spin, int k_zonked_heavy,
			   field_offset dest)
{
  w_prop_file *zonked_fp_in ; /*** Quark propagator IO stuff **/

  /************************************************************/

  zonked_fp_in = 
    r_open_wprop(reloadflag_zonked_heavy_ssink, 
		qfile_zonked_heavy_ssink[k_zonked_heavy] ); 

  node0_printf("Loading from %s\n",qfile_zonked_heavy_ssink[ k_zonked_heavy ]);
  if(reload_wprop_sc_to_site(reloadflag_zonked_heavy_ssink, zonked_fp_in,
		       spin, color, dest, 1) != 0 )  
    terminate(1);  

  r_close_wprop(reloadflag_zonked_heavy_ssink,zonked_fp_in); 

}  /***** end of load_in_zonked_heavy_smear  ******/


void load_in_zonked_heavy_local(int color, int spin, int k_zonked_heavy,
			   field_offset dest)
{
  w_prop_file *zonked_fp_in ; /*** Quark propagator IO stuff **/

  /************************************************************/

  zonked_fp_in = 
    r_open_wprop(startflag_zonked_heavy[k_zonked_heavy], 
		qfile_zonked_heavy[k_zonked_heavy] ); 
  node0_printf("Loading from %s\n",qfile_zonked_heavy[ k_zonked_heavy ]);
  if(reload_wprop_sc_to_site(startflag_zonked_heavy[k_zonked_heavy], zonked_fp_in,
		       spin, color, dest, 1) != 0 )  
    terminate(1);  

  r_close_wprop(startflag_zonked_heavy[k_zonked_heavy],zonked_fp_in); 

}  /***** end of load_in_zonked_heavy_local  ******/


void load_in_zonked_light_ssink(int color, int spin, int k_zonked_light,
			   field_offset dest)
{
  w_prop_file *zonked_fp_in ; /*** Quark propagator IO stuff **/

  /************************************************************/

  zonked_fp_in = 
    r_open_wprop(reloadflag_zonked_light_ssink, 
		qfile_zonked_light_ssink[k_zonked_light] ); 
  node0_printf("Loading from %s\n",qfile_zonked_light_ssink[ k_zonked_light ]);
  if(reload_wprop_sc_to_site(reloadflag_zonked_light_ssink, zonked_fp_in,
		       spin, color, dest, 1) != 0 )  
    terminate(1);  

  r_close_wprop(reloadflag_zonked_light_ssink,zonked_fp_in); 

}  /***** end of load_in_zonked_light_ssink  ******/


/**
 **  Generate and store the heavy zonked quark propagators to
 **  temporary disk files.
 **/

void generate_heavy_zonked(int color, int spin, 
			   Real Kappa, int inverter_type, field_offset dest)
{
  int MinCG;

  if( this_node == 0 )  printf("Starting heavy_zonked inversions, kappa = %f\n",
			       Kappa); 

  
  /*** do the heavy quark inversion ****/
  load_wilson_source(F_OFFSET(heavy_smear_func), F_OFFSET(psi), color, spin) ; 
  
  /* To be sure we "hop" far enough in the inversion */
  /* Note: one iteration advances two time slices
     so nt covers the lattice twice */
  MinCG = nt;
  
  /* Load inversion control structure */
  qic_zonked_heavy.prec = PRECISION;
  qic_zonked_heavy.min = MinCG;
  qic_zonked_heavy.max = niter_zonked_heavy;
  qic_zonked_heavy.nrestart = nrestart_zonked_heavy;
  qic_zonked_heavy.resid = resid_zonked_heavy;
  qic_zonked_heavy.start_flag = START_ZERO_GUESS;
  
#ifdef CLOVER
  /* Load Dirac matrix parameters */
  dcp.Kappa = Kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  
  /* Use specified inverter */
  if(inverter_type == HOPILU)
    {
      qic_zonked_heavy.nrestart = 1;   /* No restarts with hopping inverter */
      wilson_invert_site(F_OFFSET(psi), dest,
			 hopilu_cl_site,&qic_zonked_heavy,(void *)&dcp);
    }
  else if(inverter_type == BICGILU)
    wilson_invert_site(F_OFFSET(psi), dest,
		       bicgilu_cl_site,&qic_zonked_heavy,(void *)&dcp);
  else
    {
      printf("generate_heavy_zonked: ERROR inverter type %d unknown\n",
	     inverter_type);
      terminate(1);
    }
  
#else
  /* Load Dirac matrix parameters */
  dwp.Kappa = Kappa;
  
  wilson_invert_site(F_OFFSET(psi), dest,
		     mrilu_w_site,&qic_zonked_heavy,(void *)&dwp);
#endif
  /*** store the quark propagator to disk ****/

} /*** end of generate_and_store_heavy_zonked  ***/





/**  
 ** Set up the IO routines for the intermediate saves
 ** of the zonked quark propagators
 **
 **/


void set_zonked_save_intermediate()
{
  int i ;

  /** create the file names for the temporary files ***/
  for(i = 0 ; i < no_zonked_light ; ++i)
  {
    sprintf(qfile_zonked_light_ssink[i],"%s%s",qfile_zonked_light[i],
	    qfile_suffix_zonked_light) ; 

    wqs_zonked_light_tmp[i] = wqs_zonked_light[i] ; 
  }

  for(i = 0 ; i < no_zonked_heavy ; ++i)
  {
    sprintf(qfile_zonked_heavy_ssink[i],"%s%s",qfile_zonked_heavy[i],
	    qfile_suffix_zonked_heavy) ; 

    wqs_zonked_heavy_tmp[i] = wqs_zonked_heavy[i] ;
  }

}  /** end of the function **/




/**
 **  Set the reload_flag according to the value of 
 **  the saveflag
 **/


int set_reload_flag(int l_save_flag)
{
  int r_flag ; 

  switch (l_save_flag ) 
  {
  case SAVE_ASCII:
    r_flag = RELOAD_ASCII; 
    break ; 
  case SAVE_SERIAL:
    r_flag = RELOAD_SERIAL; 
    break ; 
  case SAVE_PARALLEL :
    r_flag = RELOAD_PARALLEL; 
    break ; 
  case SAVE_CHECKPOINT :
    r_flag = RELOAD_PARALLEL; 
    break ; 
  case SAVE_MULTIDUMP :
    r_flag = RELOAD_MULTIDUMP;
    break ;
  case FORGET :
    r_flag = FRESH;
    break ;
  default :
    IF_MASTER printf("ERROR:set_reload_flag() save_flag = %d is out of range\n",l_save_flag); 
    terminate(1) ;
    break ; 
  }


  return  r_flag ; 
}

