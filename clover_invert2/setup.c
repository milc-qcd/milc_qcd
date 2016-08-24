/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

* 5/30/07 Created from setup_cl.c */

//  $Log: setup.c,v $
//  Revision 1.13  2013/12/24 17:49:21  detar
//  Add support for COMBO type (linear combinations of props)
//  Make coordinate offset inherited.
//  Support OK action type
//  Support KS0 type (staggered prop originating from extended Dirac)
//
//  Revision 1.12  2012/11/23 23:43:21  detar
//  Support saving propagator as source.  Add OK action (IFLA) from Bugra.  Add stride to covariant smearing.
//
//  Revision 1.11  2012/04/25 03:37:42  detar
//  Initialize boundary_phase
//
//  Revision 1.10  2012/03/06 03:21:23  detar
//  Select GPU inverter precision with input parameter
//
//  Revision 1.9  2011/11/29 22:03:51  detar
//  Support arbitrary colors and momentum and boundary twists.
//
//  Revision 1.8  2009/06/04 16:37:09  detar
//  Make clover term persistent. Accommodate changes to generic_clover/make_clov2.c
//
//  Revision 1.7  2009/05/31 03:26:38  detar
//  Correct fix to "continue" handling
//
//  Revision 1.6  2009/05/31 02:00:56  detar
//  Fix "continue" and NULL startlat_p bug in clover_info.c and setup*.c
//
//  Revision 1.5  2009/04/05 16:43:07  detar
//  Utilities for open meson propagators
//
//  Revision 1.4  2008/04/18 15:12:11  detar
//  Revise metadata format for correlator files.
//
//  Revision 1.3  2008/03/28 15:42:04  detar
//  Switch to indexing by momentum-operator pair.  Add another test case.
//
//  Revision 1.2  2007/11/09 16:15:58  detar
//  Support changes in propagator I/O
//
//  Revision 1.1  2007/10/07 20:02:32  detar
//  Add new application.  Generarlizes clover_invert.
//
//


#include "cl_inv_includes.h"
//#include "lattice_qdp.h"
#include "quark_action.h"
#include <string.h>
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h
static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers();

#include "params.h"

/* Forward declarations */
static int hash_corr_label(char meson_label[MAX_CORR][MAX_MESON_LABEL], 
			   char mom_label[MAX_CORR][MAX_MOM_LABEL],
			   char *meson_label_in, char *mom_label_in, int *n);
static char decode_parity(char *parity_label_in);
static double decode_factor(char *factor_op, double factor);
static void broadcast_heap_params(void);

int setup(void)   {
  int prompt, dir;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;

  /* initialize the node random number generator */
  initialize_prn( &node_prn, param.iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  FORALLUPDIR(dir){
    boundary_phase[dir] = 0.;
  }
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  twist_in = OFF;
  /* Create clover structure */
  gen_clov = create_clov();

  return(prompt);
}


/* SETUP ROUTINES */
static int initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 clover/naive valence fermions\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    show_scidac_opts();
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");
    get_utc_datetime(utc_date_time);
    
    /* Print list of options selected */
    node0_printf("Options selected...\n");
    show_generic_opts();
    show_generic_ks_opts();

#if FERM_ACTION == HISQ
    show_su3_mat_opts();
    show_hisq_links_opts();
#endif
    status = get_prompt(stdin,  &prompt );
    
    IF_OK status += get_i(stdin,prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin,prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin,prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin,prompt,"nt", &param.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   param.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_i(stdin, prompt,"iseed", &param.iseed );
    IF_OK status += get_s(stdin, prompt,"job_id",param.job_id);
    
    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    normal_exit(0);

  nx=param.nx;
  ny=param.ny;
  nz=param.nz;
  nt=param.nt;
  iseed=param.iseed;
  
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif

  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  char savebuf[128];
  int i;
  int ipair;
  int max_cg_iterations, max_cg_restarts;
  Real bdry_phase[4];
#ifdef PRTIME
  double dtime;
#endif

  STARTTIME;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    /*------------------------------------------------------------*/
    /* Gauge configuration section                                */
    /*------------------------------------------------------------*/

    IF_OK status += ask_starting_lattice(stdin,  prompt, &param.startflag,
	param.startfile );
    IF_OK status += get_f(stdin, prompt,"u0", &param.u0 );

    IF_OK if (prompt==1) 
      printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
    IF_OK scanf("%s",savebuf);
    IF_OK printf("%s\n",savebuf);
    IF_OK {
      if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
	param.fixflag = COULOMB_GAUGE_FIX;
      }
      else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
	param.fixflag = NO_GAUGE_FIX;
      }
      else{
	printf("error in input: fixing_command %s is invalid\n",savebuf); status++;
      }
    }
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
			     param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				  param.stringLFN );

    /* APE smearing parameters (if needed) */
    /* Zero suppresses APE smearing */
    IF_OK status += get_f(stdin, prompt, "staple_weight", 
			  &param.staple_weight);
    IF_OK status += get_i(stdin, prompt, "ape_iter",
			  &param.ape_iter);

    /* Coordinate origin for KS phases and antiperiodic boundary condition */
    IF_OK status += get_vi(stdin, prompt, "coordinate_origin", param.coord_origin, 4);
    
    /*------------------------------------------------------------*/
    /* Propagator inversion control                               */
    /*------------------------------------------------------------*/

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &max_cg_iterations );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &max_cg_restarts );
    
    /*------------------------------------------------------------*/
    /* Base sources                                               */
    /*------------------------------------------------------------*/

    IF_OK status += get_i(stdin,prompt,"number_of_base_sources", 
			  &param.num_base_source);
    IF_OK {
      if(param.num_base_source > MAX_SOURCE){
	printf("Exceeded dimension %d\n",MAX_SOURCE);
	status++;
      }
    }

    for(i = 0; i < param.num_base_source; i++){
      
      IF_OK init_qs(&param.base_src_qs[i]);
      IF_OK status += get_wv_quark_source( stdin, prompt, 
					   &param.base_src_qs[i]);
      /* Base sources have no parents or ops */
      IF_OK param.parent_source[i] = BASE_SOURCE_PARENT;
      IF_OK init_qss_op(&param.src_qs_op[i]);
      IF_OK set_qss_op_offset(&param.src_qs_op[i], param.coord_origin);

      /* Get optional file for saving the base source */
      IF_OK {
	int source_type, saveflag_s;
	char descrp[MAXDESCRP];
	char savefile_s[MAXFILENAME];
	status += 
	  ask_output_quark_source_file( stdin, prompt, &saveflag_s,
					&source_type, NULL, descrp,
					savefile_s );
	IF_OK {
	  param.base_src_qs[i].savetype = source_type;
	  param.base_src_qs[i].saveflag = saveflag_s;
	  strcpy(param.base_src_qs[i].save_file, savefile_s);
	  if(saveflag_s != FORGET && source_type != DIRAC_FIELD_FILE
	     && source_type != VECTOR_FIELD_FILE){
	    printf("Unsupported output source type\n");
	    status++;
	  }
	} /* OK */
      } /* OK */
    }

    /*------------------------------------------------------------*/
    /* Modified sources                                           */
    /*------------------------------------------------------------*/

    IF_OK status += get_i(stdin,prompt,"number_of_modified_sources", 
			  &param.num_modified_source);

    IF_OK {
      if(param.num_base_source + param.num_modified_source > MAX_SOURCE){
	printf("Total including base sources exceeds dimension %d\n",
	       MAX_SOURCE);
	status++;
      }
    }

    for(i = 0; i < param.num_modified_source; i++){
      /* We append the modified sources to the list of base sources */
      int is = param.num_base_source + i;
      
      IF_OK status += get_i(stdin,prompt,"source", &param.parent_source[is]);
      
      IF_OK {
	if( param.parent_source[is] >= is){
	  printf("Source index must be less than %d here\n",is);
	  status++;
	}
      }

      IF_OK init_qss_op(&param.src_qs_op[is]);
      set_qss_op_offset(&param.src_qs_op[is], param.coord_origin);

      /* Get source operator attributes */
      IF_OK status += get_wv_field_op( stdin, prompt, &param.src_qs_op[is]);

      /* Copy parent source attributes to the derived source structure */
      IF_OK {
	int p = param.parent_source[is];
	param.base_src_qs[is] = param.base_src_qs[p];
	param.base_src_qs[is].op = copy_qss_op_list(param.base_src_qs[p].op);
	
	/* Add the new operator to the linked list */
	insert_qss_op(&param.base_src_qs[is], &param.src_qs_op[is]);
	
	/* Append the operator info to the description if the operator
	   is nontrivial, but simply copy the label */
	if(param.src_qs_op[is].type != IDENTITY){
	  char *descrp = param.base_src_qs[is].descrp;
	  char *op_descrp = param.src_qs_op[is].descrp;
	  char *label = param.base_src_qs[is].label;
	  char *op_label = param.src_qs_op[is].label;
	  strncat(descrp, "/", MAXDESCRP-strlen(descrp)-1);
	  strncat(descrp, op_descrp, MAXDESCRP-strlen(descrp)-1);
	  strncpy(label,  op_label, MAXSRCLABEL-strlen(label)-1);
	}
      }

      /* Get optional file for saving the modified source */
      IF_OK {
	int source_type, saveflag_s;
	char descrp[MAXDESCRP];
	char savefile_s[MAXFILENAME];
	status += 
	  ask_output_quark_source_file( stdin, prompt, &saveflag_s,
					&source_type, NULL, descrp,
					savefile_s );
	IF_OK {
	    param.base_src_qs[is].savetype = source_type;
	    param.base_src_qs[is].saveflag = saveflag_s;
	    strcpy(param.base_src_qs[is].save_file, savefile_s);
	  if(saveflag_s != FORGET && source_type != DIRAC_FIELD_FILE &&
	     source_type != VECTOR_FIELD_FILE){
	    printf("Unsupported output source type\n");
	    status++;
	  }
	} /* OK */
      } /* OK */
    }
	
    /*------------------------------------------------------------*/
    /* Propagators and their sources                              */
    /*------------------------------------------------------------*/

    /* Number of propagators */
    IF_OK status += get_i(stdin,prompt,"number_of_propagators", 
			  &param.num_prop );
    if( param.num_prop>MAX_PROP ){
      printf("num_prop = %d must be <= %d!\n", param.num_prop, MAX_PROP);
      status++;
    }
    
    /* Get propagator parameters */

    IF_OK for(i = 0; i < param.num_prop; i++){
      
      /* Initialize dependency */
      param.prop_dep_qkno[i] = 0;

      /* Type of propagator */

      IF_OK status += get_s(stdin, prompt,"propagator_type", savebuf );
      IF_OK {
	/* Standard clover */
	if(strcmp(savebuf,"clover") == 0)param.prop_type[i] = CLOVER_TYPE;
	/* Standard staggered (asqtad or HISQ) to be converted to naive */
	else if(strcmp(savebuf,"KS") == 0)param.prop_type[i] = KS_TYPE;
	/* Same as standard staggered, but conversion to naive is based on hypercube offset 0 */
	else if(strcmp(savebuf,"KS0") == 0)param.prop_type[i] = KS0_TYPE;
	/* Staggered propagator originating from an extended Dirac source */
	else if(strcmp(savebuf,"KS4") == 0)param.prop_type[i] = KS4_TYPE;
	/* Improved fermion lattice action (OK action) */
	else if(strcmp(savebuf,"ifla") == 0 )param.prop_type[i] = IFLA_TYPE;
	else {
	  printf("Unknown quark type %s\n",savebuf);
	  status++;
	}
      }

      /* Mass parameters, etc */

      if(param.prop_type[i] == CLOVER_TYPE){

	IF_OK status += get_s(stdin, prompt,"kappa", param.kappa_label[i]);
	IF_OK param.dcp[i].Kappa = atof(param.kappa_label[i]);
	IF_OK status += get_f(stdin, prompt,"clov_c", &param.dcp[i].Clov_c );
	param.dcp[i].U0 = param.u0;

      } else if(param.prop_type[i] == IFLA_TYPE) { 

	printf("Ifla Type Fermion\n");
#ifndef HAVE_QOP
	printf("Compilation with the QOP package is required for this fermion type\n");
	terminate(1);
#endif
	
	IF_OK status += get_s(stdin,prompt,"kapifla",param.kappa_label[i]);
	IF_OK param.nap[i].kapifla = atof(param.kappa_label[i]);
	IF_OK status += get_f(stdin, prompt, "kappa_s", &param.nap[i].kappa_s);
	IF_OK status += get_f(stdin, prompt, "kappa_t", &param.nap[i].kappa_t);
	IF_OK status += get_f(stdin, prompt, "r_s",     &param.nap[i].r_s);
	IF_OK status += get_f(stdin, prompt, "r_t",     &param.nap[i].r_t);
	IF_OK status += get_f(stdin, prompt, "zeta",    &param.nap[i].zeta);
	IF_OK status += get_f(stdin, prompt, "c_E",     &param.nap[i].c_E);
	IF_OK status += get_f(stdin, prompt, "c_B",     &param.nap[i].c_B);
	IF_OK status += get_f(stdin, prompt, "c_1",     &param.nap[i].c_1);
	IF_OK status += get_f(stdin, prompt, "c_2",     &param.nap[i].c_2);
	IF_OK status += get_f(stdin, prompt, "c_3",     &param.nap[i].c_3);
	IF_OK status += get_f(stdin, prompt, "c_4",     &param.nap[i].c_4);
	IF_OK status += get_f(stdin, prompt, "c_5",     &param.nap[i].c_5);
	IF_OK status += get_f(stdin, prompt, "c_EE",    &param.nap[i].c_EE);
	param.nap[i].u0 = param.u0;
	
      } else {  /* KS_TYPE || KS0_TYPE || KS4_TYPE */
	
	IF_OK status += get_s(stdin, prompt,"mass", param.mass_label[i] );
	IF_OK param.ksp[i].mass = atof(param.mass_label[i]);
#if FERM_ACTION == HISQ
	IF_OK status += get_f(stdin, prompt, "naik_term_epsilon", 
			      &param.ksp[i].naik_term_epsilon);
#else
	IF_OK param.ksp[i].naik_term_epsilon = 0.0;
#endif
      }

      /* Should we solve for the propagator? */
      IF_OK status += get_s(stdin, prompt,"check", savebuf);
      IF_OK {
	/* Should we be checking the propagator by running the solver? */
	if(strcmp(savebuf,"no") == 0)param.check[i] = CHECK_NO;
	else if(strcmp(savebuf,"yes") == 0)
	  param.check[i] = CHECK_YES;
	else if(strcmp(savebuf,"sourceonly") == 0)
	  param.check[i] = CHECK_SOURCE_ONLY;
	else{
	  printf("Unrecognized 'check' option. Wanted 'no', 'yes', or 'sourceonly'\n");
	  status++;
	}
      }

      /* Error for propagator conjugate gradient or bicg */
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &param.qic[i].resid );
      IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			    &param.qic[i].relresid );
      IF_OK status += get_i(stdin, prompt,"precision", &param.qic[i].prec );
#if ! defined(HAVE_QOP) && ! defined(USE_CG_GPU)
      if(param.qic[i].prec != PRECISION){
	node0_printf("WARNING: Compiled precision %d overrides request\n",PRECISION);
	node0_printf("QOP or CG_GPU compilation is required for mixed precision\n");
	param.qic[i].prec = PRECISION;
      }
#endif
      param.qic[i].max = max_cg_iterations;
      param.qic[i].nrestart = max_cg_restarts;
      param.qic[i].parity = EVENANDODD;
      param.qic[i].min = 0;
      param.qic[i].start_flag = 0;
      param.qic[i].nsrc = 1;
      
      /* Momentum twist and time boundary condition */

      IF_OK status += get_vf(stdin, prompt, "momentum_twist",
			     bdry_phase, 3);
      
      IF_OK status += get_s(stdin, prompt,"time_bc", savebuf);

      if(param.prop_type[i] == CLOVER_TYPE || param.prop_type[i] == IFLA_TYPE){

	/* NOTE: The Dirac built-in bc is periodic. */
	IF_OK {
	  if(strcmp(savebuf,"antiperiodic") == 0)bdry_phase[3] = 1;
	  else if(strcmp(savebuf,"periodic") == 0)bdry_phase[3] = 0;
	  else{
	    node0_printf("Expecting 'periodic' or 'antiperiodic' but found %s\n",
			 savebuf);
	    status++;
	  }
	}

      } else {  /* KS_TYPE || KS0_TYPE || KS4_TYPE */

	/* NOTE: The staggered built-in bc is antiperiodic.  We use
	   bdry_phase to alter the built-in convention. */
	IF_OK {
	  if(strcmp(savebuf,"antiperiodic") == 0)bdry_phase[3] = 0;
	  else if(strcmp(savebuf,"periodic") == 0)bdry_phase[3] = 1;
	  else{
	    node0_printf("Expecting 'periodic' or 'antiperiodic' but found %s\n",
			 savebuf);
	    status++;
	  }
	}
      }

      IF_OK {
	int dir;
	FORALLUPDIR(dir)param.bdry_phase[i][dir] = bdry_phase[dir];
      }

      /* Get source index for this set */
      IF_OK status += get_i(stdin,prompt,"source", &param.source[i]);

      IF_OK {
	int ns = param.num_base_source + param.num_modified_source;
	if( param.source[i] >= ns){
	  printf("Source index must be less than %d here\n",ns);
	  status++;
	}
      }

      /* Copy the base source to the propagator source structure */
      
      IF_OK {
	int p = param.source[i];
	init_qs(&param.src_qs[i]);
	param.src_qs[i] = param.base_src_qs[p];
	param.src_qs[i].op = copy_qss_op_list(param.base_src_qs[p].op);
      }

      /* Consistency check */
      IF_OK {
	if((param.prop_type[i] == KS_TYPE ||
	    param.prop_type[i] == KS0_TYPE )
	    && is_dirac_source(param.src_qs[i].type)){
	  node0_printf("ERROR: Can't build a KS propagator from a Dirac source.\n");
	  node0_printf("For an extended Dirac source, use propagator_type KS4.\n");
	  status++;
	}
      }

      /* Names of propagator files input and output */

      if(param.prop_type[i] == CLOVER_TYPE){

	IF_OK status += ask_starting_wprop( stdin, prompt, 
					    &param.startflag_w[i],
					    param.startfile_w[i]);
	
	IF_OK status += ask_ending_wprop( stdin, prompt, &param.saveflag_w[i],
					  param.savefile_w[i]);

      } else {  /* KS_TYPE || KS0_TYPE || KS4_TYPE */
	
	/* Get propagator file names */
	
	if(param.prop_type[i] == KS_TYPE || param.prop_type[i] == KS0_TYPE){
	  
	  IF_OK status += ask_starting_ksprop( stdin, prompt, 
					       &param.startflag_ks[i],
					       param.startfile_ks[i]);
	  
	  IF_OK status += ask_ending_ksprop( stdin, prompt, 
					     &param.saveflag_ks[i],
					     param.savefile_ks[i]);
	  
	} else { 	  /* KS4_TYPE */


	  /* We expect the standard start/end file specifications, but for
	     this type there is no standard propagator input file
	     format with four staggered fermion propagators.  Once the
	     propagator has been computed, it is a Dirac propagator,
	     so it is written as such */

	  IF_OK status += ask_starting_wprop( stdin, prompt, 
					      &param.startflag_w[i],
					      param.startfile_w[i]);
	  
	  IF_OK status += ask_ending_wprop( stdin, prompt, 
					    &param.saveflag_w[i],
					    param.savefile_w[i]);

	  if(param.startflag_w[i] != FRESH)
	    printf("WARNING: Input file specification changed to 'fresh_wprop'\n");
	}
      }
    }

    /*------------------------------------------------------------*/
    /* Quarks                                                     */
    /*------------------------------------------------------------*/

    /* Number of quarks */
    IF_OK status += get_i(stdin,prompt,"number_of_quarks", 
			  &param.num_qk );
    if( param.num_qk>MAX_QK ){
      printf("num_qk = %d must be <= %d!\n", param.num_qk, MAX_QK);
      status++;
    }

    IF_OK for(i = 0; i < param.num_qk; i++){
      char *check_tag;

      /* Initialize dependency */
      param.quark_dep_qkno[i] = 0;

      /* Get the propagator that we act on with the sink operator to
	 form the "quark" field used in the correlator.  It might be a
	 raw "propagator" or it might be a previously constructed
	 "quark" field */
      /* First we look for the type */
      IF_OK {
	if(prompt==1)printf("enter 'propagator' or 'quark' or 'combine'\n");
	check_tag = get_next_tag(stdin, "propagator or quark or combine", "readin");
	printf("%s ",check_tag);
	if(strcmp(check_tag,"propagator") == 0)
	  param.parent_type[i] = PROP_TYPE;
	else if(strcmp(check_tag,"quark") == 0)
	  param.parent_type[i] = QUARK_TYPE;
	else if(strcmp(check_tag,"combine") == 0)
	  param.parent_type[i] = COMBO_TYPE;
	else{
	  printf("\nError: expected 'propagator' or 'quark' or 'combine'\n");
	  status++;
	}
      }
      if( param.parent_type[i] == PROP_TYPE ||  param.parent_type[i] == QUARK_TYPE ){
	/* Next we get its index */
	IF_OK {
	  if(prompt==1)printf("enter the index\n");
	  if(scanf("%d",&param.prop_for_qk[i]) != 1){
	    printf("\nFormat error reading index\n");
	    status++;
	  } else{
	    printf("%d\n",param.prop_for_qk[i]);
	    if(param.parent_type[i] == PROP_TYPE){
	      if(param.prop_for_qk[i] >= param.num_prop){
		printf("Propagator index must be less than %d\n",
		       param.num_prop);
		status++;
	      }
	      /* Update dependency */
	      param.prop_dep_qkno[param.prop_for_qk[i]] = i;
	    }
	    else if(param.parent_type[i] == QUARK_TYPE){
	      if(param.prop_for_qk[i] >= i){
		printf("Quark index must be less than %d here\n",i);
		status++;
	      }
	      /* Update depedency */
	      param.quark_dep_qkno[param.prop_for_qk[i]] = i;
	    }
	  }
	} 
      } else { /* COMBO_TYPE */
	/* Next we get its index */
	IF_OK {
	  if(prompt==1)printf("enter the number or quarks to combine\n");
	  if(scanf("%d",&param.num_combo[i]) != 1){
	    printf("\nFormat error reading index\n");
	    status++;
	  } else {
	    int j;
	    printf("%d\n",param.num_combo[i]);
	    if(param.num_combo[i] > MAX_COMBO){
	      printf("\nError: number to combine must not exceed %d\n", MAX_COMBO);
	      status++;
	    }
	    IF_OK status += get_vi(stdin, prompt, "quarks", param.combo_qk_index[i], param.num_combo[i]);
	    IF_OK for(j = 0; j < param.num_combo[i]; j++){
	      if(param.combo_qk_index[i][j] >= i){
		printf("Error: quark index %d must be less than %d here\n", param.combo_qk_index[i][j], i);
		status++;
	      }
	      /* Update depedency */
	      param.quark_dep_qkno[param.combo_qk_index[i][j]] = i;
	    }
	    IF_OK status += get_vf(stdin, prompt, "coeffs", param.combo_coeff[i], param.num_combo[i]);
	  }
	  /* Provisional: For constructing the ancestry for the correlator file. */
	  param.prop_for_qk[i] = param.combo_qk_index[i][0];
	}
      }
      /* Get sink operator attributes */
      IF_OK init_qss_op(&param.snk_qs_op[i]);
      IF_OK status += get_wv_field_op( stdin, prompt, &param.snk_qs_op[i]);
      /* Get sink quark save attributes */
      IF_OK {
	char descrp[MAXDESCRP];
	status += 
	  ask_ending_wprop_or_wsource( stdin, prompt, &param.saveflag_q[i], 
				       &param.savetype_q[i], NULL, descrp,
				       param.savefile_q[i]);
      }
    }
    
    /*------------------------------------------------------------*/
    /* Meson correlators                                          */
    /*------------------------------------------------------------*/
    
    /* Number of quark pairs */
    IF_OK status += get_i(stdin,prompt,"number_of_pairings", &param.num_pair );
    if( param.num_pair>MAX_PAIR ){
      printf("num_pair = %d must be <= %d!\n", param.num_pair, MAX_PAIR);
      status++;
    }
    
    IF_OK for(ipair = 0; ipair < param.num_pair; ipair++){
      char request_buf[MAX_SPECTRUM_REQUEST];
      
      /* Which quarks in the pair? */
      
      IF_OK status += get_vi(stdin, prompt, "pair", param.qkpair[ipair], 2);
      
      IF_OK {
	int j;
	for(j = 0; j < 2; j++){
	  if(param.qkpair[ipair][j] < 0 || param.qkpair[ipair][j] >= param.num_qk){
	    printf("Quark index %d must be in [0,%d]\n",
		   param.qkpair[ipair][j],param.num_qk-1);
	    status++;
	  }
	  /* Don't destroy this quark until the hadrons are constructed. */
	  param.quark_dep_qkno[param.qkpair[ipair][j]] = param.num_qk;
	}
      }
      
      /* What spectrum calculation to do? */
      /* prepend and append a comma for ease in parsing */
      IF_OK request_buf[0] = '\0';
      IF_OK strcpy(request_buf,",");
      IF_OK status += get_s(stdin, prompt,"spectrum_request", request_buf+1 );
      IF_OK strcat(request_buf,",");
      
      /* Parse the spectrum request */
      IF_OK {
	/* Mesons */
	if(strstr(request_buf,",meson,") != NULL)
	  param.do_meson_spect[ipair] = 1;
	else param.do_meson_spect[ipair] = 0;
	
	/* Baryons */
	if(strstr(request_buf,",baryon,") != NULL)
	  param.do_baryon_spect[ipair] = 1;
	else param.do_baryon_spect[ipair] = 0;

	/* Any meson or baryon */
	param.do_closed_hadron[ipair]  = 
	  param.do_meson_spect[ipair] |
	  param.do_baryon_spect[ipair];

	/* Open meson */
	if(strstr(request_buf,",open_meson,") != NULL)
	  param.do_open_meson[ipair] = 1;
	else param.do_open_meson[ipair] = 0;

	/* We need some correlator specification */
	if(!param.do_open_meson[ipair] && !param.do_closed_hadron[ipair]){
	  printf("Unrecognized spectrum request\n");
	  status++;
	}

	/* But we can't do both closed and open hadrons in the same pairing */
	if(param.do_open_meson[ipair] && param.do_closed_hadron[ipair]){
	  printf("Can't combine open_meson with other options ");
	  printf("in the same pairing\n");
	  status++;
	}
      }
      
      /* What file for the resulting correlators? */
      
      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_c[ipair],
				     param.savefile_c[ipair]);
      /* Correlator values are printed with t relative to the t_offset */
      /* FT phases are computed with x,y,z relative to r_offset */
      IF_OK {
	int r[4];
	status += get_vi(stdin,prompt, "r_offset", r, 4);
	param.r_offset[ipair][0] = r[0];
	param.r_offset[ipair][1] = r[1];
	param.r_offset[ipair][2] = r[2];
	param.t_offset[ipair]    = r[3];
      }
      
      /*------------------------------------------------------------*/
      /* Table of correlator combinations actually needed           */
      /*------------------------------------------------------------*/
      
      /* Number of correlators */
      IF_OK status += get_i(stdin,prompt,"number_of_correlators", 
			    &param.num_corr[ipair] );
      IF_OK {
	if( param.num_corr[ipair]>MAX_CORR ){
	  printf("num_corr = %d must be <= %d!\n", param.num_corr[ipair], 
		 MAX_CORR);
	  status++;
	}
      }
      
      /* Sample format for correlator line:
	 correlator  A1_P5 p200 -i * 1 G5 G5X 2 0 0 E E E */
      
      param.num_corr_report[ipair] = 0;
      IF_OK for(i = 0; i < param.num_corr[ipair]; i++){
	int ok,m;
	char meson_label_in[MAX_MESON_LABEL], mom_label_in[MAX_MOM_LABEL],
	  gam_src_lab[MAXGAMMA], gam_snk_lab[MAXGAMMA], phase_lab[4],
	  factor_op[2], parity_x_in[3], parity_y_in[3], parity_z_in[3];
	double factor;
	
	/* Read the meson label */
	IF_OK status += get_sn(stdin, prompt, "correlator", meson_label_in);
	
	/* Read the momentum label */
	IF_OK {
	  ok = scanf("%s", mom_label_in);
	  if(ok != 1){
	    printf("Error reading momentum label\n");
	    status++;
	  }
	}
	IF_OK printf(" %s", mom_label_in);
	
	/* Find or add to hash table */
	/* The tuple (meson_label, mom_label) is the correlator label */
	/* The hash table lists the unique tuples. m is the index for the tuple. */
	m = hash_corr_label(param.meson_label[ipair], param.mom_label[ipair],
			    meson_label_in, mom_label_in,
			    &param.num_corr_report[ipair]);
	param.corr_index[ipair][i] = m;
	if(m < 0)status++;
	
	/* phase, op, factor, and gamma matrix specification */
	IF_OK {
	  ok = scanf("%s %s %lf %s %s\n",phase_lab,factor_op,&factor,
		     gam_src_lab,gam_snk_lab);
	  if(ok != 5){
	    printf("\nError reading phase, factor, and gammas\n");
	    status++;
	  }
	  else {
	    printf(" %3s %1s %6f %3s %3s",phase_lab,factor_op,factor,
		   gam_src_lab,gam_snk_lab);
	  }
	}
	
	/* decode phase for correlator */
	IF_OK {
	  param.corr_phase[ipair][i] = decode_phase(phase_lab);
	  if(param.corr_phase[ipair][i] == -1){
	    printf("\n%s is not a valid phase label\n",phase_lab);
	    status ++;
	  }
	}
	
	/* decode real factor for correlator */
	/* Permitted syntax is /ddd or *ddd where ddd is a real number */
	IF_OK {
	  param.corr_factor[ipair][i] = decode_factor(factor_op,factor);
	  if(param.corr_factor[ipair][i] == 0 ){
	    printf("\n: Decoding factor %s %g results in zero.\n",
		   factor_op,factor);
	    status++;
	  }
	}
	
	/* decode gamma matrix labels for source and sink */
	IF_OK {
	  param.gam_src[ipair][i] = gamma_index(gam_src_lab);
	  if(param.gam_src[ipair][i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",gam_src_lab);
	    status += 1;
	  }
	  param.gam_snk[ipair][i] = gamma_index(gam_snk_lab);
	  if(param.gam_snk[ipair][i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",gam_snk_lab);
	    status ++;
	  }
	}
	
	/* momentum indices */
	IF_OK {
	  if(scanf("%d %d %d",&param.corr_mom[ipair][i][0],
		   &param.corr_mom[ipair][i][1],
		   &param.corr_mom[ipair][i][2]) != 3){
	    printf("\nFormat error on momentum line\n");
	    status++;
	  }
	  
	  IF_OK printf(" %2d %2d %2d",param.corr_mom[ipair][i][0],
		       param.corr_mom[ipair][i][1],
		       param.corr_mom[ipair][i][2]);
	}
	
	/* parity for each momentum component */
	
	IF_OK {
	  ok = scanf("%s %s %s", parity_x_in, parity_y_in, parity_z_in);
	  if(ok != 3){
	    printf("\nError reading parity label\n");
	    status++;
	  }
	  IF_OK printf(" %2s %2s %2s\n", parity_x_in, parity_y_in, parity_z_in);
	}
	
	/* decode parity labels */
	IF_OK {
	  param.corr_parity[ipair][i][0] = decode_parity(parity_x_in);
	  param.corr_parity[ipair][i][1] = decode_parity(parity_y_in);
	  param.corr_parity[ipair][i][2] = decode_parity(parity_z_in);
	  if(param.corr_parity[ipair][i][0] == 0xf)status++;
	  if(param.corr_parity[ipair][i][1] == 0xf)status++;
	  if(param.corr_parity[ipair][i][2] == 0xf)status++;
	}
      } /* correlators for this pair */
    } /* pairs */
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    

  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0;
  if( param.stopflag != 0 )return param.stopflag;

  if(prompt==2)return 0;

  /* Broadcast parameter values kept on the heap */
  broadcast_heap_params();

  /* Construct the eps_naik table of unique Naik epsilon
     coefficients.  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* First term is always zero */
  start_eps_naik(eps_naik, &n_naiks);
  
  /* Contribution from the propagator epsilons */
  
  for(i = 0; i < param.num_prop; i++){
    if(param.prop_type[i] == CLOVER_TYPE)continue;
    param.ksp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, 
		    &n_naiks, param.ksp[i].naik_term_epsilon);
  }

  /* Assign Naik term indices to quarks based on inheritance from
     propagators */

  for(i = 0; i < param.num_qk; i++){
    if(param.parent_type[i] == PROP_TYPE)
      param.naik_index[i] = param.ksp[param.prop_for_qk[i]].naik_term_epsilon_index;
    else
      param.naik_index[i] = param.naik_index[param.prop_for_qk[i]];
  }

  /* Do whatever is needed to get lattice */
  if( param.startflag != CONTINUE ){
    startlat_p = reload_lattice( param.startflag, param.startfile );
    invalidate_this_clov(gen_clov);
  }

  phases_in = OFF;

  /* Construct APE smeared links */
  ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );

  /* Set options for fermion links in case we use them */
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif
  
  /* We start with no fermion links.
     They are created in make_prop if we need them */

  fn_links = NULL;

  ENDTIME("readin");

  return(0);
}

/* Broadcast operator parameter values.  They are on the heap on node 0. */

static void broadcast_heap_params(void){
  int i;

  for(i = 0; i < param.num_base_source + param.num_modified_source; i++){
    broadcast_quark_source_sink_op_recursive(&param.base_src_qs[i].op);
    broadcast_quark_source_sink_op_recursive(&param.src_qs_op[i].op);
  }

  for(i = 0; i < param.num_prop; i++)
    broadcast_quark_source_sink_op_recursive(&param.src_qs[i].op);

  for(i = 0; i < param.num_qk; i++)
    broadcast_quark_source_sink_op_recursive(&param.snk_qs_op[i].op);

}

/* decode parity label */
static char decode_parity(char *parity_label_in){
  if(strcmp("E",parity_label_in) == 0)return EVEN;
  else if(strcmp("O",parity_label_in) == 0)return ODD;
  else if(strcmp("EO",parity_label_in) == 0)return EVENANDODD;

  printf("\nUnrecognized parity label %s\n",parity_label_in);
  return 0xf;
}
    
/* decode real factor label */
/* Permitted syntax is *ddd or /ddd where ddd is a real number */
static double decode_factor(char *factor_op, double factor){

  if(factor == 0)
    return factor;
  else if(factor_op[0] == '*')
    return factor;
  else if(factor_op[0] == '/')
    return 1/factor;
  else
    printf("Incorrect operator for normalization factor\n");
  return 0.;
}
    

/* Manage hash table for correlator labels */

/* Look the labels up in the table.  If found, return its
   index.  Otherwise, add it to the table and assign it the new index */
static int 
hash_corr_label(char meson_label[MAX_CORR][MAX_MESON_LABEL], 
		char mom_label[MAX_CORR][MAX_MOM_LABEL],
		char *meson_label_in, char *mom_label_in, int *n){
  int i;

  for(i = 0; i < *n; i++){
    if(strcmp(meson_label_in, meson_label[i]) == 0 &&
       strcmp(mom_label_in,   mom_label[i]  ) == 0)
      return i;
  }
  if(*n >= MAX_CORR){
    printf("Too many correlators labels %d\n",*n);
    return -1;
  }
  strcpy(meson_label[*n],meson_label_in);
  strcpy(mom_label[*n],mom_label_in);
  (*n)++;
  return *n-1;
}


/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
static void 
make_3n_gathers()
{
  int i;
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X3UP);

}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

static void 
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
#if 0
    /*------------------------------------------------------------*/
    /* Gamma matrix pairings                                      */
    /*------------------------------------------------------------*/
    
    /* Number of gamma matrix pairings */
    IF_OK status += get_i(stdin,prompt,"number_of_mesons", &param.num_meson );
    IF_OK {
      if( param.num_meson>MAX_MESON ){
	printf("num_meson = %d must be <= %d!\n", param.num_meson, MAX_MESON);
	status++;
      }
      param.num_meson_report = 0;
    }

    IF_OK for(i = 0; i < param.num_meson; i++){
      char meson_label_in[MAX_MESON_LABEL];
      char gam_src_lab[MAXGAMMA], gam_snk_lab[MAXGAMMA];
      char phase_lab[4];

      /* meson label */
      IF_OK status += get_sn(stdin, prompt, "meson_source_sink", meson_label_in);
      IF_OK {
	param.meson_index[i] = 
	  hash_meson_label(param.meson_label, meson_label_in);
	if(param.meson_index[i] < 0)status++;
      }

      /* gamma matrix labels and phase label */
      IF_OK scanf("%s %s %s\n",gam_src_lab,gam_snk_lab,phase_lab);
      IF_OK printf(" \t%s\t%s\t%s\n",gam_src_lab,gam_snk_lab,phase_lab);

      /* decode gamma matrix labels for source and sink */
      IF_OK {
	param.gam_src[i] = gamma_index(gam_src_lab);
	if(param.gam_src[i] < 0){
	  printf("%s is not a valid gamma matrix label\n",gam_src_lab);
	  status += 1;
	}
	param.gam_snk[i] = gamma_index(gam_snk_lab);
	if(param.gam_snk[i] < 0){
	  printf("%s is not a valid gamma matrix label\n",gam_snk_lab);
	  status ++;
	}
      }

      /* decode phase for correlator */
      IF_OK {
	param.meson_phase[i] = decode_phase(phase_lab);
	if(param.meson_phase[i].real == 0 &&
	   param.meson_phase[i].imag == 0 ){
	  printf("%s is not a valid phase label\n",phase_lab);
	  status ++;
	}
      }
    }
    
    /*------------------------------------------------------------*/
    /* Momentum table                                             */
    /*------------------------------------------------------------*/

    /* Number of momenta */
    IF_OK status += get_i(stdin,prompt,"number_of_meson_momenta", 
			  &param.num_mom );
    IF_OK {
      if( param.num_mom>MAX_MESON_MOMENTUM ){
	printf("num_mom = %d must be <= %d!\n", 
	       param.num_mom, MAX_MESON_MOMENTUM);
	status++;
      }
      param.num_mom_report = 0;
    }
    
    /* momentum label */

    IF_OK for(i = 0; i < param.num_mom; i++){
      char mom_label_in[MAX_MOM_LABEL];
      IF_OK status += get_sn(stdin, prompt, "momentum", mom_label_in);

      IF_OK {
	param.mom_index[i] = 
	  hash_momentum_label(param.mom_label, mom_label_in);
	if(param.mom_index[i] < 0)status++;
      }
    
      /* momentum indices */
      IF_OK {
	if(scanf("%d %d %d",&param.meson_mom[i][0],
		 &param.meson_mom[i][1],&param.meson_mom[i][2]) != 3){
	  printf("Format error on momentum line\n");
	  status++;
	}
	
	IF_OK printf(" %d %d %d\n",param.meson_mom[i][0],
		     param.meson_mom[i][1],param.meson_mom[i][2]);
      }
    }
#endif
