/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//  $Log: setup.c,v $
//  Revision 1.21  2009/04/14 17:09:24  detar
//  AB: Calculate structure of multi_x array in terms of different Naik corrections
//  AB: Preset Naik corrections in eps_naik[] controlled by HISQ_NAIT_ADJUSTABLE
//
//  Revision 1.20  2008/03/28 15:20:22  detar
//  Bazavov changes 3/08.  Debugging and Naik term mass.
//
//  Revision 1.19  2007/12/14 04:51:45  detar
//  Add HISQ code.
//
//  Revision 1.18  2007/11/09 16:06:28  detar
//  Pull FN link calculation out of inverter.
//
//  Revision 1.17  2007/10/09 19:51:02  detar
//  Add ferm_links_t and ks_action_paths structures and pass them as params
//
//  Revision 1.16  2007/05/21 04:23:32  detar
//  Add Precision selection for fermion force in QOP and QDP
//
//  Revision 1.15  2007/03/27 20:48:54  detar
//  Fix n_pseudo check.
//
//  Revision 1.14  2006/12/13 18:41:40  detar
//  Add precision arg to mat_invert_uml and mat_invert_cg
//
//  Revision 1.13  2006/12/09 14:10:57  detar
//  Add mixed precision capability for QDP and QOP inverters
//
//  Revision 1.12  2006/11/16 06:07:50  detar
//  Moved optimization selections to Makefile
//  Put rational function parameters in preamble of input file
//  Recreated fiducial output files and remade errtol files
//  Remade rationals.sample files
//  Removed obsolete files
//
//  Revision 1.11  2006/11/07 02:30:58  detar
//  Fix some omissions to complete the previous update.
//
//  Revision 1.10  2006/11/04 23:35:14  detar
//  Add separate CG control for MD, FA, GR
//  Add nrestart parameter
//  Remove some debugging lines from sample output files and errtol files
//
//  Revision 1.9  2006/10/29 02:35:54  detar
//  Abandon rationals.h and load parameters from file instead.
//
//  Revision 1.8  2006/10/09 03:44:04  detar
//  Move fermion_force_fn to generic_ks/fermion_force_fn_multi.c
//  and path_transport_field to generic/path_transport.c
//  Change ks_multicg selection method to a set_opts call.
//
//  Revision 1.7  2006/10/02 04:13:50  detar
//  Distinguish inverter residuals for molecular dynamics and action
//
//  Revision 1.6  2006/09/19 03:06:58  detar
//  Upgrade for concurrent EOS calculations
//
//  Revision 1.5  2006/08/24 04:29:02  detar
//  Remove or protect unused variable declarations.
//  Remove unwanted Makefile
//
//  Revision 1.4  2006/08/22 19:32:22  detar
//  Upgrade to QDP and train error tolerance file
//
//  Revision 1.3  2006/08/13 15:16:06  detar
//  Train error files.  Trivia.
//
//  Revision 1.2  2006/08/13 04:02:32  detar
//  Switch to function pointer for selecting multicg inverter species
//
//  Revision 1.1  2006/08/09 04:22:19  detar
//  Adding Doug's code to repository
//
//  Revision 1.3  2005/11/10 16:58:44  detar
//  Experimenting with tags and versions
//
//  Revision 1.2  2005/11/10 16:55:09  detar
//  Add CVS version and log tags to file
//
//
*/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include "quark_action.h"
#include "lattice_qdp.h"
#include "../include/su3_mat_op.h"

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

gauge_file *r_parallel_i(char *);
void r_parallel(gauge_file *, field_offset);
void r_parallel_f(gauge_file *);

gauge_file *r_binary_i(char *);
void r_binary(gauge_file *);
void r_binary_f(gauge_file *);
void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int initial_set();
void make_3n_gathers();

int
setup()
{
  int prompt;
#ifdef HAVE_QDP
  int i;
#endif
  
  /* print banner, get volume, seed */
  prompt = initial_set();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);

  /* Mark fermion links as unallocated */
  init_ferm_links(&fn_links);
#ifdef DM_DU0
  init_ferm_links(&fn_links_dmdu0);
#endif
  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  node0_printf("Made 3nn gathers\n"); fflush(stdout);
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  
#ifdef HAVE_QDP
  for(i=0; i<4; ++i) {
    shiftdirs[i] = QDP_neighbor[i];
    shiftdirs[i+4] = neighbor3[i];
  }
  for(i=0; i<8; ++i) {
    shiftfwd[i] = QDP_forward;
    shiftbck[i] = QDP_backward;
  }
#endif

  node0_printf("Finished setup\n"); fflush(stdout);
  return( prompt );
}


/* SETUP ROUTINES */
int 
initial_set()
{
  int prompt,status,i,tmporder;
  Real current_naik_mass;

  /* On node zero, read lattice size, seed, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("MIMD version 7 $Name:  $\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    printf("Rational function hybrid Monte Carlo algorithm\n");
    /* Print list of options selected */
    node0_printf("Options selected...\n");
    show_generic_opts();
    show_generic_ks_opts();
#ifdef INT_ALG
    node0_printf("INT_ALG=%s\n",ks_int_alg_opt_chr());
#endif
#ifdef HISQ_NAIK_ADJUSTABLE
    node0_printf("HISQ_NAIK_ADJUSTABLE (means Naik correction is full epsilon and not just mass)\n");
#endif
#ifdef HISQ_FORCE_FILTER
    node0_printf("HISQ_FORCE_FILTER=%f\n",HISQ_FORCE_FILTER);
#endif
#ifdef HISQ_REUNIT_ALLOW_SVD
    node0_printf("HISQ_REUNIT_ALLOW_SVD\n");
#endif
#ifdef HISQ_REUNIT_SVD_ONLY
    node0_printf("HISQ_REUNIT_SVD_ONLY (used together with HISQ_REUNIT_ALLOW_SVD)\n");
#endif
#ifdef MILC_GLOBAL_DEBUG
    node0_printf("MILC_GLOBAL_DEBUG ***********************\n");
#endif
#ifdef HISQ_REUNITARIZATION_DEBUG
    node0_printf("HISQ_REUNITARIZATION_DEBUG is ON\n");
#endif
#ifdef HISQ_FF_MULTI_WRAPPER
    node0_printf("HISQ_FF_MULTI_WRAPPER is ON\n");
#endif
#ifdef HISQ_FF_DEBUG
    node0_printf("HISQ_FF_DEBUG is ON\n");
#endif
    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   par_buf.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   par_buf.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );
    /* Number of pseudofermions */
    IF_OK status += get_i(stdin, prompt,"n_pseudo", &par_buf.n_pseudo );
    if(par_buf.n_pseudo > MAX_N_PSEUDO){
      printf("Error:  Too many pseudofermion fields.  Recompile. Current max is %d\n"
	     ,MAX_N_PSEUDO);
      terminate(1);
    }
    /* get name of file containing rational function parameters */
    IF_OK status += get_s(stdin, prompt, "load_rhmc_params", 
			  par_buf.rparamfile);
    /* beta, quark masses */
    IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

    IF_OK status += get_i(stdin, prompt,"n_dyn_masses", &par_buf.n_dyn_masses );
    IF_OK status += get_vf(stdin, prompt, "dyn_mass", par_buf.dyn_mass, par_buf.n_dyn_masses);
    IF_OK status += get_vi(stdin, prompt, "dyn_flavors", par_buf.dyn_flavors, par_buf.n_dyn_masses);

    IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(mynode()==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )
    normal_exit(0);
  
  nx        = par_buf.nx;
  ny        = par_buf.ny;
  nz        = par_buf.nz;
  nt        = par_buf.nt;
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = par_buf.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = par_buf.ionode_geometry[i];
#endif
#endif
  iseed     = par_buf.iseed;
  n_pseudo  = par_buf.n_pseudo;
  strcpy(rparamfile,par_buf.rparamfile);
  
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  total_iters=0;

  /* Load rational function parameters */
  rparam = load_rhmc_params(rparamfile, n_pseudo);  
  if(rparam == NULL)terminate(1);

  /* Determine the maximum rational fcn order */
  max_rat_order = 0;
  for(i = 0; i < n_pseudo; i++){
    if(rparam[i].MD.order > max_rat_order)max_rat_order = rparam[i].MD.order;
    if(rparam[i].GR.order > max_rat_order)max_rat_order = rparam[i].GR.order;
    if(rparam[i].FA.order > max_rat_order)max_rat_order = rparam[i].FA.order;
  }
  node0_printf("Maximum rational func order is %d\n",max_rat_order);


  /* Determine the number of different Naik masses */
  current_naik_mass = rparam[0].naik_term_mass;
  tmporder = 0;
  n_naiks = 0;
  n_order_naik_total = 0;
  for( i=0; i<n_pseudo; i++ ) {
    if( rparam[i].naik_term_mass != current_naik_mass ) {
      if( tmporder > 0 ) {
        n_orders_naik[n_naiks] = tmporder;
        masses_naik[n_naiks] = current_naik_mass;
        current_naik_mass = rparam[i].naik_term_mass;
        n_naiks++;
        n_order_naik_total += tmporder;
        tmporder = 0;
      }
    }
    tmporder += rparam[i].MD.order;
    n_pseudo_naik[n_naiks]++;
  }
  if( tmporder > 0 ) {
    n_orders_naik[n_naiks] = tmporder;
    masses_naik[n_naiks] = current_naik_mass;
    n_order_naik_total += tmporder;
    n_naiks++;
  }
  // calculate epsilon corrections for different Naik terms
  if( 0!=masses_naik[0] ) {
    node0_printf("IN HISQ ACTION FIRST SET OF PSEUDO FERMION FIELDS SHOULD HAVE EPSILON CORRECTION TO NAIK TERM ZERO.\n");
    terminate(1);
  }
  eps_naik[0] = 0.0; // first set of X links always has 0 correction
  for( i=1; i<n_naiks; i++ ) {
#ifdef HISQ
#ifdef HISQ_NAIK_ADJUSTABLE
    // value read from rational function file is considered full epsilon correction
    eps_naik[i] = masses_naik[i];
#else
    // value read from rational function file is considered quark mass
    // and epsilon correction is calculated with the second order perturbation theory,
    // HISQ_NAIK_2ND_ORDER is set in the hisq_action.h
    eps_naik[i] = HISQ_NAIK_2ND_ORDER*masses_naik[i]*masses_naik[i];
#endif
#else /* HISQ */
    // IT IS ASSUMED THAT ACTIONS OTHER THAN HISQ DO NOT HAVE
    // ANY EPSILON CORRECTION TERMS
    eps_naik[i] = 0;
#endif /* HISQ */
  }
  node0_printf("Naik term correction structure of multi_x:\n");
  node0_printf("n_naiks %d\n",n_naiks);
  for( i=0; i<n_naiks; i++ ) {
    node0_printf("n_pseudo_naik[%d]=%d\n",i,n_pseudo_naik[i]);
    node0_printf("n_orders_naik[%d]=%d\n",i,n_orders_naik[i]);
    node0_printf("masses_naik[%d]=%f\n",i,masses_naik[i]);
    node0_printf("eps_naik[%d]=%f\n",i,eps_naik[i]);
  }
  node0_printf("n_order_naik_total %d\n",n_order_naik_total);
#ifdef HISQ
  if( n_naiks+1 > MAX_NAIK ) {
    node0_printf("MAX_NAIK=%d < n_naiks+1=%d\n", MAX_NAIK, n_naiks+1 );
    node0_printf("Increase MAX_NAIK\n");
    terminate(1);
  }
#else /* HISQ */
  if( n_naiks>1 ) {
    node0_printf("FOR ACTIONS OTHER THAN HISQ EPSILON CORRECTION IS NOT USED.\n");
    node0_printf("ONLY ONE SET OF X LINKS IS USED.\n");
    node0_printf("SET ALL naik_mass TO 0 IN RATIONAL FUNCTION FILE.\n");
    terminate(1);
  }
#endif /* HISQ */

  beta = par_buf.beta;
  
  n_dyn_masses = par_buf.n_dyn_masses;
  for(i = 0; i < n_dyn_masses; i++){
    dyn_mass[i] = par_buf.dyn_mass[i];
    dyn_flavors[i] = par_buf.dyn_flavors[i];
  }
  u0 = par_buf.u0;

  return(prompt);
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  Real x;
  int i;
#ifdef SPECTRUM
  char request_buf[MAX_SPECTRUM_REQUEST];
#endif
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {
    
    printf("\n\n");
    status=0;
    
    /* warms, trajecs */
    IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
    IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );
    
    /* trajectories between propagator measurements */
    IF_OK status += 
      get_i(stdin, prompt,"traj_between_meas", &par_buf.propinterval );
    
    /* microcanonical time step */
    IF_OK status += 
      get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );
    
    /*microcanonical steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
    
    /* Data for each pseudofermion */

    for(i = 0; i < par_buf.n_pseudo; i++){
      Real tmp[3]; int itmp[3];

      /* Residuals for multicg solves */
      IF_OK status += get_vf(stdin, prompt,"cgresid_md_fa_gr", tmp, 3 );
      /* rsqmin is r**2 in conjugate gradient */
      IF_OK {
	par_buf.rsqmin_md[i] = tmp[0]*tmp[0];
	par_buf.rsqmin_fa[i] = tmp[1]*tmp[1];
	par_buf.rsqmin_gr[i] = tmp[2]*tmp[2];
      }

      /* Max CG iterations for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "max_multicg_md_fa_gr", itmp, 3);
      IF_OK {
	par_buf.niter_md[i] = itmp[0];
	par_buf.niter_fa[i] = itmp[1];
	par_buf.niter_gr[i] = itmp[2];
      }

      /* Precision for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "cgprec_md_fa_gr", itmp, 3);
      IF_OK {
	par_buf.prec_md[i] = itmp[0];
	par_buf.prec_fa[i] = itmp[1];
	par_buf.prec_gr[i] = itmp[2];
      }
    }

    /* Precision for fermion force calculation */
    IF_OK status = get_i(stdin, prompt, "prec_ff", &par_buf.prec_ff);

    /* error for propagator conjugate gradient */
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
    IF_OK par_buf.rsqprop = x*x;
    
    /* maximum no. of conjugate gradient iterations for propagator
       etc. and maximum no. of restarts */
    IF_OK status += get_i(stdin, prompt,"max_cg_prop", &par_buf.niter );
    IF_OK status += get_i(stdin, prompt,"max_cg_prop_restarts", 
			  &par_buf.nrestart );
    
#ifdef NPBP_REPS
    /* number of random sources npbp_reps and precision for inversions */
    IF_OK status += get_i(stdin, prompt,"npbp_reps", &par_buf.npbp_reps_in );
    IF_OK status += get_i(stdin, prompt,"prec_pbp", &par_buf.prec_pbp );
#endif
    
#ifdef SPECTRUM
    /* request list for spectral measurments */
    /* prepend and append a comma for ease in parsing */
    IF_OK status += get_s(stdin, prompt,"spectrum_request", request_buf );
    IF_OK strcpy(par_buf.spectrum_request,",");
    IF_OK strcat(par_buf.spectrum_request,request_buf);
    IF_OK strcat(par_buf.spectrum_request,",");
    
    /* source time slice and increment */
    IF_OK status += get_i(stdin, prompt,"source_start", &par_buf.source_start );
    IF_OK status += get_i(stdin, prompt,"source_inc", &par_buf.source_inc );
    IF_OK status += get_i(stdin, prompt,"n_sources", &par_buf.n_sources );
    
    /* Additional parameters for spectrum_multimom */
    if(strstr(par_buf.spectrum_request,",spectrum_multimom,") != NULL){
      IF_OK status += get_i(stdin, prompt,"spectrum_multimom_nmasses",
			    &par_buf.spectrum_multimom_nmasses );
      IF_OK status += get_f(stdin, prompt,"spectrum_multimom_low_mass",
			    &par_buf.spectrum_multimom_low_mass );
      IF_OK status += get_f(stdin, prompt,"spectrum_multimom_mass_step",
			    &par_buf.spectrum_multimom_mass_step );
    }
    /* Additional parameters for fpi */
    par_buf.fpi_nmasses = 0;
    if(strstr(par_buf.spectrum_request,",fpi,") != NULL){
      IF_OK status += get_i(stdin, prompt,"fpi_nmasses",
			    &par_buf.fpi_nmasses );
      if(par_buf.fpi_nmasses > MAX_FPI_NMASSES){
	printf("Maximum of %d exceeded.\n",MAX_FPI_NMASSES);
	terminate(1);
      }
      for(i = 0; i < par_buf.fpi_nmasses; i++){
	IF_OK status += get_f(stdin, prompt,"fpi_mass",
			      &par_buf.fpi_mass[i]);
      }
    }
    
#endif /*SPECTRUM*/

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
					  par_buf.startfile );
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
					par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				  par_buf.stringLFN );
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )return par_buf.stopflag;
  
  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  propinterval = par_buf.propinterval;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  for(i = 0; i< n_pseudo; i++){
    niter_md[i] = par_buf.niter_md[i];
    niter_fa[i] = par_buf.niter_fa[i];
    niter_gr[i] = par_buf.niter_gr[i];

    rsqmin_md[i] = par_buf.rsqmin_md[i];
    rsqmin_fa[i] = par_buf.rsqmin_fa[i];
    rsqmin_gr[i] = par_buf.rsqmin_gr[i];

    prec_md[i] = par_buf.prec_md[i];
    prec_fa[i] = par_buf.prec_fa[i];
    prec_gr[i] = par_buf.prec_gr[i];
  }
  prec_ff = par_buf.prec_ff;
  npbp_reps_in = par_buf.npbp_reps_in;
  prec_pbp = par_buf.prec_pbp;
  rsqprop = par_buf.rsqprop;
  epsilon = par_buf.epsilon;
  n_pseudo = par_buf.n_pseudo;
#ifdef SPECTRUM
  strcpy(spectrum_request,par_buf.spectrum_request);
  source_start = par_buf.source_start;
  source_inc = par_buf.source_inc;
  n_sources = par_buf.n_sources;
  spectrum_multimom_nmasses = par_buf.spectrum_multimom_nmasses;
  spectrum_multimom_low_mass = par_buf.spectrum_multimom_low_mass;
  spectrum_multimom_mass_step = par_buf.spectrum_multimom_mass_step;
  fpi_nmasses = par_buf.fpi_nmasses;
  for(i = 0; i < fpi_nmasses; i++){
    fpi_mass[i] = par_buf.fpi_mass[i];
  }
#endif /*SPECTRUM*/
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  int isite, idir;
  site *s;
  FORALLSITES(isite,s) {
    for( idir=XUP;idir<=TUP;idir++ ) {
      lattice[isite].on_step_Y[idir] = 0;
      lattice[isite].on_step_W[idir] = 0;
      lattice[isite].on_step_V[idir] = 0;
      /* zero out information from previous time step
         if fresh lattice; keep everything for continuation */
      if( startflag != CONTINUE ){
        lattice[isite].phase_Y_previous[idir] = 0.0;
        lattice[isite].phase_Y[idir] = 0.0;
        lattice[isite].Vdet[idir] = 0.0;
        clear_su3mat( &(lattice[isite].Wlink[idir]) );
        clear_su3mat( &(lattice[isite].Wlink_previous[idir]) );
      }
    }
  }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */

  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  startlat_p = reload_lattice( startflag, startfile );
  /* if a lattice was read in, put in KS phases and AP boundary condition */
#ifdef FN
  invalidate_all_ferm_links(&fn_links);
  invalidate_all_ferm_links(&fn_links_dmdu0);
#endif
  phases_in = OFF;
  rephase( ON );
  
  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();
  /* make tables of coefficients and permutations of paths in quark action */
  init_path_table(&ks_act_paths);
  init_path_table(&ks_act_paths_dmdu0);
  make_path_table(&ks_act_paths, &ks_act_paths_dmdu0);
  
  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
void 
make_3n_gathers()
{
  int i;
#ifdef HAVE_QDP
  int disp[4]={0,0,0,0};
#endif
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X3UP);

#ifdef HAVE_QDP
  for(i=0; i<4; i++) {
    disp[i] = 3;
    neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
  }
#endif
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

void 
third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
	"forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+3)%nx; break;
  case XDOWN: *xp = (x+4*nx-3)%nx; break;
  case YUP: *yp = (y+3)%ny; break;
  case YDOWN: *yp = (y+4*ny-3)%ny; break;
  case ZUP: *zp = (z+3)%nz; break;
  case ZDOWN: *zp = (z+4*nz-3)%nz; break;
  case TUP: *tp = (t+3)%nt; break;
  case TDOWN: *tp = (t+4*nt-3)%nt; break;
  default: printf("third_neighb: bad direction\n"); exit(1);
  }
}
