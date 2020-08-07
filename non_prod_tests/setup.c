/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ... */

//  $Log: setup.c,v $
//  Revision 1.7  2013/12/24 05:32:40  detar
//  Add combo type.  Support embedded inverse.
//
//  Revision 1.6  2012/11/24 05:14:20  detar
//  Add support for U(1) fields and for future HYPISQ action
//
//  Revision 1.5  2012/04/25 03:21:29  detar
//  Initialize boundary phase
//
//  Revision 1.4  2012/03/06 03:23:26  detar
//  Set GPU inverter precision through input parameter
//
//  Revision 1.3  2012/01/21 21:35:08  detar
//  Support general spin_taste interpolating operators for mesons.
//
//  Revision 1.2  2011/12/03 03:44:36  detar
//  Fix support for mu_eos
//
//  Revision 1.1  2011/11/30 22:11:40  detar
//  Add
//
//


#include "ani_non_prod_tests_includes.h"
#include "lattice_qdp.h"
#include <string.h>
#include <ctype.h>
#include "params.h"
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#endif

/* Forward declarations */
#ifdef ANISOTROPY
static int dirstring2index (char savebuf[], int *status);
#endif
static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers(void);

static int hash_corr_label(char meson_label[MAX_CORR][MAX_MESON_LABEL], 
			   char mom_label[MAX_CORR][MAX_MOM_LABEL],
			   char *meson_label_in, char *mom_label_in, int *n);
static char decode_parity(char *parity_label_in);
static double decode_factor(char *factor_op, double factor);
static void broadcast_heap_params(void);


int setup()   {
  int prompt, dir;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, param.iseed, volume+mynode() );
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
#ifdef U1_FIELD
  u1_A = create_u1_A_field();
#endif
  FORALLUPDIR(dir){
    boundary_phase[dir] = 0.;
  }
  /* Initialize fermion links as unallocated */
  //  init_ferm_links(&fn_links, &ks_act_paths);
  //  init_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* set up K-S phase vectors, antiperiodic boundary conditions */
  phaseset();

  return(prompt);
}


static int n_naiks = 1;
static double eps_naik[MAX_NAIK];

/* SETUP ROUTINES */
static int initial_set(void){
  int prompt=0,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 staggered valence fermions\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
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
#elif FERM_ACTION == HYPISQ
    show_su3_mat_opts();
    show_hypisq_links_opts();
#endif

    status = get_prompt(stdin,  &prompt );
    
    IF_OK status += get_i(stdin,prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin,prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin,prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin,prompt,"nt", &param.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_SUBNODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "subnode_geometry", 
			   param.subnode_geometry, 4);
#endif
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

  if(prompt==2)return prompt;

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

  number_of_nodes = numnodes();
  volume=(size_t)(nx*ny*nz*nt);

  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  char savebuf[128];
  int i,k,nprop;
  int ipair, itriplet;
#ifdef PRTIME
  double dtime;
#endif

  STARTTIME;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    /* Coordinate origin for KS phases and antiperiodic boundary condition */
    IF_OK status += get_vi(stdin, prompt, "coordinate_origin", param.coord_origin, 4);
    IF_OK status += get_s(stdin, prompt, "time_bc", savebuf);
    IF_OK {
      /* NOTE: The staggered default time bc is antiperiodic. */
      if(strcmp(savebuf,"antiperiodic") == 0)param.time_bc = 0;
      else if(strcmp(savebuf,"periodic") == 0)param.time_bc = 1;
      else{
	node0_printf("Expecting 'periodic' or 'antiperiodic' but found %s\n", savebuf);
	status++;
      }
    }

#ifdef ANISOTROPY
    /* Direction of anisotropy */
    IF_OK status += get_s(stdin, prompt,"ani_dir",savebuf);
    IF_OK param.ani_dir = dirstring2index( savebuf, &status);
    IF_OK ani_dir = param.ani_dir;
    /* Bare fermion anisotropy */
    IF_OK status += get_f(stdin, prompt, "ani_xiq", &param.ani_xiq);
    IF_OK ani_xiq = param.ani_xiq;
#ifdef ONEDIM_ANISO_TEST
    /* bare quark isotropic link factor for debugging */
    IF_OK status += get_f(stdin, prompt, "iso_xiq", &param.iso_xiq);
    IF_OK iso_xiq = param.iso_xiq;
#endif
#endif

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
      
      IF_OK init_qs(&param.src_qs[i]);
      IF_OK status += get_v_quark_source( stdin, prompt, 
					  &param.src_qs[i]);
      /* Base sources have no parents or ops */
      IF_OK param.parent_source[i] = BASE_SOURCE_PARENT;
      IF_OK init_qss_op(&param.src_qs_op[i]);
      /* Enforce a uniform boundary condition */
      IF_OK set_qss_op_offset(&param.src_qs_op[i], param.coord_origin);
      /* NOTE: The KS built-in bc is antiperiodic. */
      IF_OK param.src_qs_op[i].bp[3] = param.time_bc;    
      
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
	  param.src_qs[i].savetype = source_type;
	  param.src_qs[i].saveflag = saveflag_s;
	  strcpy(param.src_qs[i].save_file, savefile_s);
	  if(saveflag_s != FORGET && source_type != VECTOR_FIELD_FILE){
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
     
      /* Get source operator attributes */
      IF_OK status += get_v_field_op( stdin, prompt, &param.src_qs_op[is]);
      /* Enforce a uniform boundary condition */
      set_qss_op_offset(&param.src_qs_op[is], param.coord_origin);
      /* NOTE: The KS built-in bc is antiperiodic. */
      param.src_qs_op[is].bp[3] = param.time_bc;    
     
      /* Copy parent source attributes to the derived source structure */
      IF_OK {
        int p = param.parent_source[is];
        param.src_qs[is] = param.src_qs[p];
        param.src_qs[is].op = copy_qss_op_list(param.src_qs[p].op);
     
        /* Add the new operator to the linked list */
        insert_qss_op(&param.src_qs[is], &param.src_qs_op[is]);
     
        /* Append the operator info to the description if the operator
           is nontrivial, but simply copy the label */
        if(param.src_qs_op[is].type != IDENTITY){
          char *descrp = param.src_qs[is].descrp;
          char *op_descrp = param.src_qs_op[is].descrp;
          char *label = param.src_qs[is].label;
          char *op_label = param.src_qs_op[is].label;
          strncat(descrp, "/", MAXDESCRP-strlen(descrp)-1);
          strncat(descrp, op_descrp, MAXDESCRP-strlen(descrp)-1);
          /* it is not possible to conform with previous test cases for correlator keys 
           * and retain a coherent strategy for labeling extended sources; 
           * using strcpy without the reduced length instead 
           * otherwise use  strcat as commented out here 
           * (with large enough MAXSRCLABEL defined i n../include/generic_quark_types.h)
           * strncat(label,  op_label, MAXSRCLABEL-strlen(label)-1);
           *********************************************** */
          strncpy(label,  op_label, MAXSRCLABEL);
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
          param.src_qs[is].savetype = source_type;
          param.src_qs[is].saveflag = saveflag_s;
          strcpy(param.src_qs[is].save_file, savefile_s);
          if(saveflag_s != FORGET && source_type != VECTOR_FIELD_FILE){
            printf("Unsupported output source type\n");
            status++;
          }
        } /* OK */
      } /* OK */
    }

    /*------------------------------------------------------------*/
    /* Propagators and their sources                              */
    /*------------------------------------------------------------*/
    
    /* Number of sets grouped for multimass inversion */
    
    IF_OK status += get_i(stdin,prompt,"number_of_sets", &param.num_set);
    if( param.num_set>MAX_SET ){
      printf("num_set = %d must be <= %d!\n", param.num_set, MAX_SET);
      status++;
    }

    nprop = 0;
    IF_OK for(k = 0; k < param.num_set; k++){
      int max_cg_iterations, max_cg_restarts;
      int check = CHECK_NO;

#ifdef MULTISOURCE
      IF_OK status += get_s(stdin, prompt, "set_type", savebuf);
      IF_OK {
	if(strcmp(savebuf,"multimass") == 0)
	  param.set_type[k] = MULTIMASS_SET;
	else if(strcmp(savebuf,"multisource") == 0)
	  param.set_type[k] = MULTISOURCE_SET;
	else {
	  printf("Unrecognized set type %s\n",savebuf);
	  printf("Choices are 'multimass', 'multisource'\n");
	  status++;
	}
      }
#else
	  param.set_type[k] = MULTIMASS_SET;
#endif

      /* maximum no. of conjugate gradient iterations */
      IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			    &max_cg_iterations );
      
      /* maximum no. of conjugate gradient restarts */
      IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			    &max_cg_restarts );

      /* Should we be checking (computing) the propagator by running
	 the solver? */

      IF_OK status += get_s(stdin, prompt,"check", savebuf);
      IF_OK {
	/* Should we be checking the propagator by running the solver? */
	if(strcmp(savebuf,"no") == 0)check = CHECK_NO;
	else if(strcmp(savebuf,"yes") == 0)
	  check = CHECK_YES;
	else if(strcmp(savebuf,"sourceonly") == 0)
	  check = CHECK_SOURCE_ONLY;
	else {
	  node0_printf("Unrecognized check command %s\n",savebuf);
	  node0_printf("Choices are 'yes', 'no', 'sourceonly'\n");
	  status++;
	}
      }

      /* Boundary conditions (including twist) apply to all members of the set */
      /* They specify the phase that is added to the default, which is
	 periodic in x, y, z */
      /* So 0 0 0 gets the default, the one we commonly use. */
      /* The values are entered as fractions of pi */
      /* So 0 0 0.5 inserts a phase exp(i 0.5 pi) */
      
      Real bdry_phase[4];
      IF_OK status += get_vf(stdin, prompt, "momentum_twist",
			     bdry_phase, 3);
      bdry_phase[3] = param.time_bc;  /* Enforce a uniform boundary condition */

      IF_OK {
	IF_OK status += get_i(stdin, prompt,"precision", &param.qic[0].prec );
#if ! defined(HAVE_QOP) && ! defined(USE_CG_GPU) && !defined(HAVE_QPHIX)
	IF_OK if(param.qic[0].prec != MILC_PRECISION){
	  node0_printf("WARNING: Compiled precision %d overrides request\n",MILC_PRECISION);
	  node0_printf("QOP or CG_GPU or QPHIX compilation is required for mixed precision\n");
	  param.qic[0].prec = MILC_PRECISION;   /* Same for all members of a set*/
	}
#endif
      }

#ifdef U1_FIELD
      /* Charge, if U(1) */
      IF_OK status += get_s(stdin, prompt, "charge", param.charge_label[k]);
      IF_OK param.charge[k] = atof(param.charge_label[k]);
#endif

      int tmp_src;
      Real tmp_naik;
      
      IF_OK {

	if(param.set_type[k] == MULTIMASS_SET){
	  
	  /* Get source index common to this set */
	  IF_OK status += get_i(stdin,prompt,"source", &tmp_src);
	} else {
	  
	  /* Get mass label common to this set */
	  IF_OK status += get_s(stdin,prompt,"mass", savebuf);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	  IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
				&tmp_naik);
#else
	  tmp_naik = 0.0;
#endif
	}
	
      }

      /* Number of propagators in this set */
      IF_OK status += get_i(stdin,prompt,"number_of_propagators", 
			    &param.num_prop[k]);
      if( param.num_prop[k]>MAX_PROP ){
	printf("num_prop = %d must be <= %d!\n", param.num_prop[k], MAX_PROP);
	status++;
      }

      /* Indexing range for set */
      param.begin_prop[k] = nprop;
      param.end_prop[k] = nprop + param.num_prop[k] - 1;
      if(param.end_prop[k] > MAX_PROP){
	printf("Total number_of_propagators must be <= %d!\n", MAX_PROP);
	status++;
      }

      IF_OK for(i = 0; i < param.num_prop[k]; i++){
    
	/* Propagator parameters */

	IF_OK {
	  
	  if(param.set_type[k]  == MULTIMASS_SET){
	    
	    /* Get mass label common to this set */
	    IF_OK status += get_s(stdin,prompt,"mass", param.mass_label[nprop]);
	    
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	    IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
				  &param.ksp[nprop].naik_term_epsilon);
#else
	    param.ksp[nprop].naik_term_epsilon = 0.0;
#endif
	    param.source[nprop] = tmp_src;
	    
	  } else {
	    
	    /* Get source index common to this set */
	    IF_OK status += get_i(stdin,prompt,"source", &(param.source[nprop]));
	    strcpy(param.mass_label[nprop], savebuf);
	    param.ksp[nprop].naik_term_epsilon = tmp_naik;
	  }
	}

	IF_OK param.ksp[nprop].mass = atof(param.mass_label[nprop]);

	IF_OK {
	  int dir;
	  FORALLUPDIR(dir)param.bdry_phase[nprop][dir] = bdry_phase[dir];
	}
	IF_OK param.check[nprop] = check;  /* Same for all members of a set */
	IF_OK param.set[nprop] = k;  /* The set to which this prop belongs */

	/*------------------------------------------------------------*/
	/* Propagator inversion control                               */
	/*------------------------------------------------------------*/
	
	/* maximum no. of conjugate gradient iterations */
	param.qic[nprop].max = max_cg_iterations;
      
	/* maximum no. of conjugate gradient restarts */
	param.qic[nprop].nrestart = max_cg_restarts;
      
	/* Should we be deflating? */
	param.qic[nprop].deflate = 0;
	IF_OK {
	  if(param.eigen_param.Nvecs > 0){  /* Need eigenvectors to deflate */
	    IF_OK status += get_s(stdin, prompt,"deflate", savebuf);
	    IF_OK {
	      if(strcmp(savebuf,"yes") == 0)param.qic[nprop].deflate = 1;
	    }
	  }
	}

	/* error for clover propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			      &param.qic[nprop].resid );
	IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			      &param.qic[nprop].relresid );
#if defined(HALF_MIXED) && defined(HAVE_QOP)
	/* Parameter used by QOPQDP inverter for mixed-precision solves */
	IF_OK status += get_f(stdin, prompt, "mixed_rsq", &param.qic[nprop].mixed_rsq );
#endif
	/* Precision for all members of the set must be the same */
	param.qic[nprop].prec = param.qic[0].prec;
	
	/* Parity is always EVENANDODD for spectroscopy */
	param.qic[nprop].parity = EVENANDODD;
    
	IF_OK status += ask_starting_ksprop( stdin, prompt, 
					     &param.startflag_ks[nprop],
					     param.startfile_ks[nprop]);
	
	IF_OK status += ask_ending_ksprop( stdin, prompt, 
					   &param.saveflag_ks[nprop],
					   param.savefile_ks[nprop]);
	
	nprop++;
      }
    } // END IF_OK for(k = 0; k < param.num_set; k++)


#ifdef FUNNYLINKS
    IF_OK status += get_vf(stdin, prompt, "umu", param.umu, 4);
#endif

    /* warms, trajecs */
    IF_OK status += get_i(stdin, prompt,"warms", &param.warms );
    IF_OK status += get_i(stdin, prompt,"trajecs", &param.trajecs );

    if( status > 0) return status;
    /* trajectories between propagator measurements */
    IF_OK status +=
        get_i(stdin, prompt,"traj_between_meas", &param.propinterval );
    
#ifndef ANISOTROPY
    /* beta */
    IF_OK status += get_f(stdin, prompt,"beta", &param.beta );
#else
    /* beta[0] - space, beta[1] - time */
    IF_OK status += get_vf(stdin, prompt,"beta", param.beta, 2 );
#endif
    /* no dynamical masses for pure gauge */
    n_dyn_masses = 0;
    dyn_flavors[0] = 0;

#ifdef HMC_ALGORITHM
    /* steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &param.steps );
    /* microcanonical time step */
    IF_OK status +=
        get_f(stdin, prompt,"microcanonical_time_step", &param.epsilon );
#endif  

#ifdef ORA_ALGORITHM
    /*overrelaxed steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &param.steps );
    IF_OK status += get_i(stdin, prompt,"qhb_steps", &param.stepsQ );
#endif   

    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
  
  
  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0 = 1.;
  if( param.stopflag != 0 )return param.stopflag;

  if(prompt==2)return 0;

  /* Broadcast parameter values kept on the heap */
  broadcast_heap_params();

  /* Construct the eps_naik table of unique Naik epsilon coefficients.
     Also build the hash table for mapping a mass term to its Naik
     epsilon index.  We need to list any required Naik epsilons. */

  /* First term is always zero */
  start_eps_naik(eps_naik, &n_naiks);
  
  /* Contribution from the chiral condensate epsilons */
  for(i = 0; i < param.num_pbp_masses; i++){
    param.ksp_pbp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, 
	  &n_naiks, param.ksp_pbp[i].naik_term_epsilon);
  }

  /* Contribution from the propagator epsilons */
  nprop = param.end_prop[param.num_set-1] + 1;
  for(i = 0; i < nprop; i++)
    param.ksp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, 
		    &n_naiks, param.ksp[i].naik_term_epsilon);

  /* Requests from any embedded inverse and hopping operators in the
     modified ops */
  for(i = 0; i < param.num_modified_source; i++){
    Real eps = 0.;
    int is = param.num_base_source + i;
    /* If the operator uses Dslash, get the requested Naik epsilon */
    if(get_qss_eps_naik(&eps, &param.src_qs_op[is])){
      insert_qss_eps_naik_index(fill_eps_naik(eps_naik, &n_naiks, eps), &param.src_qs_op[is]);
    }
  }

  /* Requests from any embedded inverse and hopping operators in
     the sink ops */
  for(i = 0; i < param.num_qk; i++){
    Real eps = 0.;
    /* If the operator uses Dslash, get the requested Naik epsilon */
    if(get_qss_eps_naik(&eps, &param.snk_qs_op[i])){
      insert_qss_eps_naik_index(fill_eps_naik(eps_naik, &n_naiks, eps), &param.snk_qs_op[i]);
    }
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
  if( param.startflag == CONTINUE ){
    rephase( OFF );
  }
  if( param.startflag != CONTINUE ){
    startlat_p = reload_lattice( FRESH, param.startfile );
  }
  /* if a lattice was read in, put in KS phases and AP boundary condition */
  phases_in = OFF;
  rephase( ON );

  /* Set options for fermion links */
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif

  /* Don't need to save HISQ auxiliary links */
  fermion_links_want_aux(0);
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )

#ifdef DM_DEPS
  fermion_links_want_deps(1);
#endif

  fn_links = create_fermion_links_from_site(MILC_PRECISION, n_naiks, eps_naik);

#else

#ifdef DM_DU0
  fermion_links_want_du0(1);
#endif

  fn_links = create_fermion_links_from_site(MILC_PRECISION, 0, NULL);

#endif

  /* Construct APE smeared links without KS phases, but with
     conventional antiperiodic bc.  This is the same initial
     setup as the gauge field itself.  Later the phases are
     adjusted according to boundary phases and momentum twists. */
/*
  rephase( OFF );
  ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
  if(param.time_bc == 0)apply_apbc( ape_links, param.coord_origin[3] );
  rephase( ON );
*/

  warms = param.warms;
  trajecs = param.trajecs;
  steps = param.steps;
  stepsQ = param.stepsQ;
  propinterval = param.propinterval;
  epsilon = param.epsilon;

#ifndef ANISOTROPY
  beta = par_buf.beta;
#else
  beta[0] = param.beta[0];
  beta[1] = param.beta[1];
#endif

  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();

#ifdef ANISOTROPY
  /* figure out which loops are temporal and which are spatial */
  path_determine_ani();
#endif

  ENDTIME("readin");

  return 0;
}

/* Broadcast operator parameter values.  They are on the heap on node 0. */

static void broadcast_heap_params(void){
  int i, k;

  for(i = 0; i < param.num_base_source + param.num_modified_source; i++){
    broadcast_quark_source_sink_op_recursive(&param.src_qs[i].op);
    broadcast_quark_source_sink_op_recursive(&param.src_qs_op[i].op);
  }

//  for(k = 0; k < param.num_set; k++)
//    broadcast_quark_source_sink_op_recursive(&param.src_qs[k].op);

  for(i = 0; i < param.num_qk; i++)
    broadcast_quark_source_sink_op_recursive(&param.snk_qs_op[i].op);
}

/* decode parity label */
char decode_parity(char *parity_label_in){
  if(strcmp("E",parity_label_in) == 0)return EVEN;
  else if(strcmp("O",parity_label_in) == 0)return ODD;
  else if(strcmp("EO",parity_label_in) == 0)return EVENANDODD;

  printf("\nUnrecognized parity label %s\n",parity_label_in);
  return 0xf;
}
    
/* decode real factor label */
/* Permitted syntax is *ddd or /ddd where ddd is a real number */
static double 
decode_factor(char *factor_op, double factor){

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

#ifdef ANISOTROPY
static int dirstring2index (char savebuf[], int *status) {
  short mydir;

  if ( savebuf[0] >= 'A' && savebuf[0] <= 'Z' )
    mydir = tolower(savebuf[0]);
  else
    mydir = savebuf[0];
  switch(mydir) {
    case XUP:
    case '0':
    case 'x': mydir = XUP; break;
    case YUP:
    case '1':
    case 'y': mydir = YUP; break;
    case ZUP:
    case '2':
    case 'z': mydir = ZUP; break;
    case TUP:
    case '3':
    case 't': mydir = TUP; break;
    default:
      node0_printf("Expecting direction \
as x,y,z,t, X,Y,Z,T, or 0,1,2,3;  instead %c\n", mydir);
     (*status)++;
  }
  return mydir;
}
#endif


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
make_3n_gathers(void)
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
