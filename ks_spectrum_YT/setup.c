/******** setup.c *********/
/* MIMD version 7 */

#define _POSIX_C_SOURCE 200112L /* for gethostname */
#define IF_OK if(status==0)

#include "ks_spectrum_includes.h"
#include "lattice_qdp.h"
#include <string.h>
#include "params.h"
#include <unistd.h>
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#endif

/* Forward declarations */

static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers(void);
static void second_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
			    int *xp, int *yp, int *zp, int *tp);
static void make_2n_gathers(void);

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
  ape_links = NULL;
  FORALLUPDIR(dir){
    boundary_phase[dir] = 0.;
  }
  /* Initialize fermion links as unallocated */
  //  init_ferm_links(&fn_links, &ks_act_paths);
  //  init_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);

  /* The following make gathers must appear in this order in order
     to be compatible with the include/dirs.h macros */
  /* They must also be called before any make_gathers call */
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* make 2nd neighbor gathers for two-link gauss smearing */
#ifdef GAUSS_SMEAR_KS_TWOLINK
  make_2n_gathers();
#endif
  /* set up K-S phase vectors, antiperiodic boundary conditions */
  phaseset();

  return(prompt);
}


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
    IF_OK {
      int iseed_in;
      status += get_i(stdin, prompt,"iseed", &iseed_in);
      param.iseed = iseed_in;
    }
    IF_OK status += get_s(stdin, prompt,"job_id",param.job_id);

    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  fflush(stdout);
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    terminate(1);

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
  volume=(size_t)nx*ny*nz*nt;

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
#ifdef GB_BARYON
  int ioctet, iqk;
#endif
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
    IF_OK status += scanf("%s",savebuf)==1 ? 0 : 1;
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
#ifdef U1_FIELD
    /* what kind of starting U(1) lattice to use, read filename */
    IF_OK status+=ask_starting_u1_lattice(stdin,prompt,
					  &param.start_u1flag, param.start_u1file );
    IF_OK status+=ask_ending_u1_lattice(stdin,prompt,
					&param.save_u1flag, param.save_u1file );
#endif
    /* Provision is made to build covariant sources from smeared
       links */
    /* APE smearing parameters (if needed) */
    /* Zero suppresses APE smearing */
#ifdef APE_LINKS_FILE
    IF_OK status += ask_starting_apelinks(stdin, prompt, &param.start_ape_flag, param.start_ape_file);
    if(param.start_ape_flag == FRESH){
      IF_OK status += get_f(stdin, prompt, "staple_weight",
			    &param.staple_weight);
      IF_OK status += get_i(stdin, prompt, "ape_iter",
			    &param.ape_iter);
    }
    IF_OK status += ask_ending_apelinks(stdin, prompt, &param.save_ape_flag, param.save_ape_file);
#else
    IF_OK status += get_f(stdin, prompt, "staple_weight",
			  &param.staple_weight);
    IF_OK status += get_i(stdin, prompt, "ape_iter",
			  &param.ape_iter);
#endif
    
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
    IF_OK param.ani_dir = dirchar2index( savebuf[0], &status);
    IF_OK status += ( param.ani_dir*param.ani_dir > TUP*TUP );
    /* Bare fermion anisotropy */
    IF_OK status += get_f(stdin, prompt, "ani_xiq", &param.ani_xiq);
#  ifdef ONEDIM_ANISO_TEST
    /* bare quark isotropic link factor for debugging */
    IF_OK status += get_f(stdin, prompt, "iso_xiq", &param.iso_xiq);
#  endif
#endif


    /* number of eigenpairs */
    IF_OK status += get_i(stdin, prompt,"max_number_of_eigenpairs", &param.eigen_param.Nvecs);

    IF_OK if(param.eigen_param.Nvecs > 0){

      /* eigenvector input */
      IF_OK status += ask_starting_ks_eigen(stdin, prompt, &param.ks_eigen_startflag,
					    param.ks_eigen_startfile);

      /* eigenvector output */
      IF_OK status += ask_ending_ks_eigen(stdin, prompt, &param.ks_eigen_saveflag,
					  param.ks_eigen_savefile);

      /* If we are reading in eigenpairs, we don't regenerate them */

#if EIGMODE != EIGCG
      if(param.ks_eigen_startflag == FRESH){

	/*------------------------------------------------------------*/
	/* Dirac eigenpair parameters                                 */
	/*------------------------------------------------------------*/

	status += read_ks_eigen_param(&param.eigen_param, status, prompt);

      }
#else

      /* for eigcg */

      /* maximum number of eigenvectors */
      param.eigcgp.Nvecs_max =  param.eigen_param.Nvecs;

      /* If we are reading in eigenpairs, we don't regenerate them */

      if(param.ks_eigen_startflag == FRESH){


	/* restart for Lanczos */
	IF_OK status += get_i(stdin, prompt,"restart_lanczos", &param.eigcgp.m);

	/* number of eigenvectors per inversion */
	IF_OK status += get_i(stdin, prompt,"Number_of_eigenvals_per_inversion", &param.eigcgp.Nvecs);

	IF_OK {
	  if(param.eigcgp.m <= 2*param.eigcgp.Nvecs){
	    printf("restart_lanczos should be larger than 2*Number_of_eigenvals!\n");
	    status++;
	  }
	}
      } else {
	param.eigcgp.m = 0;
	param.eigcgp.Nvecs = 0;
      }

      param.eigcgp.Nvecs_curr = 0;
      param.eigcgp.H = NULL;
#endif

    }

    /*------------------------------------------------------------*/
    /* Chiral condensate and related quantities                   */
    /*------------------------------------------------------------*/

    IF_OK status += get_i(stdin, prompt, "number_of_pbp_masses",
			  &param.num_pbp_masses);
    if(param.num_pbp_masses > MAX_MASS_PBP){
      printf("Number of masses exceeds dimension %d\n",MAX_MASS_PBP);
      status++;
    }
    IF_OK if(param.num_pbp_masses > 0){
      IF_OK status += get_i(stdin, prompt, "max_cg_iterations",
			    &param.qic_pbp[0].max);
      IF_OK status += get_i(stdin, prompt, "max_cg_restarts",
			    &param.qic_pbp[0].nrestart);
      IF_OK status += get_i(stdin, prompt, "npbp_reps", &param.npbp_reps );
      IF_OK status += get_i(stdin, prompt, "prec_pbp", &param.qic_pbp[0].prec);
      IF_OK for(i = 0; i < param.num_pbp_masses; i++){
	IF_OK status += get_f(stdin, prompt, "mass", &param.ksp_pbp[i].mass);
#if ( FERM_ACTION == HISQ )
	IF_OK status += get_f(stdin, prompt, "naik_term_epsilon", 
			      &param.ksp_pbp[i].naik_term_epsilon);
#else
	param.ksp_pbp[i].naik_term_epsilon = 0.0;
#endif
#ifdef U1_FIELD
	IF_OK status += get_f(stdin, prompt, "charge", &param.charge_pbp[i]);
#endif
	param.qic_pbp[i].min = 0;
	param.qic_pbp[i].nsrc = 1;
	param.qic_pbp[i].max = param.qic_pbp[0].max;
	param.qic_pbp[i].nrestart = param.qic_pbp[0].nrestart;
	param.qic_pbp[i].prec = param.qic_pbp[0].prec;
	/* Should we be deflating? */
	param.qic_pbp[i].deflate = 0;
	IF_OK {
	  if(param.eigen_param.Nvecs > 0){  /* Need eigenvectors to deflate */
	    IF_OK status += get_s(stdin, prompt,"deflate", savebuf);
	    IF_OK {
	      if(strcmp(savebuf,"yes") == 0)param.qic_pbp[i].deflate = 1;
	    }
	  }
	}
	IF_OK status += get_f(stdin, prompt, "error_for_propagator", &param.qic_pbp[i].resid);
	IF_OK status += get_f(stdin, prompt, "rel_error_for_propagator", &param.qic_pbp[i].relresid );
#ifdef HALF_MIXED
	IF_OK status += get_f(stdin, prompt, "mixed_rsq", &param.qic_pbp[i].mixed_rsq );
#endif
	IF_OK param.qic_pbp[i].inv_type = UMLTYPE;
      }
    }

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
      int max_cg_iterations, max_inner_cg_iterations, max_cg_restarts;
      enum check_type check = CHECK_NO;
      char mgparamfile[MAXFILENAME] = "";

      IF_OK status += get_s(stdin, prompt, "set_type", savebuf);
      IF_OK {
	if(strcmp(savebuf,"multimass") == 0)
	  param.set_type[k] = MULTIMASS_SET;
	else if(strcmp(savebuf,"multisource") == 0)
	  param.set_type[k] = MULTISOURCE_SET;
	else if(strcmp(savebuf,"single") == 0)
	  param.set_type[k] = SINGLES_SET;
	else if(strcmp(savebuf,"multicolorsource") == 0)
	  param.set_type[k] = MULTICOLORSOURCE_SET;
	else {
	  printf("Unrecognized set type %s\n",savebuf);
	  printf("Choices are 'single', 'multimass', 'multisource', 'multicolorsource'\n");
	  status++;
	}
      }

      IF_OK status += get_s(stdin, prompt, "inv_type", savebuf);
      IF_OK {
	if(strcmp(savebuf,"MG") == 0)
	  param.inv_type[k] = MGTYPE;
	else if(strcmp(savebuf,"CG") == 0)
	  param.inv_type[k] = CGTYPE;
	else if(strcmp(savebuf,"CGZ") == 0)
	  param.inv_type[k] = CGZTYPE;
	else if(strcmp(savebuf,"UML") == 0)
	  param.inv_type[k] = UMLTYPE;
	else {
	  printf("Unrecognized inverter type %s\n",savebuf);
	  printf("Choices are 'CG', 'CGZ', 'MG', 'UML'\n");
	  status++;
	}
      }
      
      IF_OK {
        if (param.inv_type[k] == MGTYPE) {
          IF_OK status += get_s(stdin, prompt, "MGparams", mgparamfile);
        }

	/* maximum no. of conjugate gradient iterations */
        IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			      &max_cg_iterations );

	/* maximum no. of conjugate gradient restarts */
        IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			      &max_cg_restarts );

#if (defined(HALF_MIXED) || defined(MAX_MIXED)) && ! defined(HAVE_QUDA) && defined(HAVE_GRID)
	/* (QUDA sets its own value).  We need this value for GRID mixed precision */
        IF_OK status += get_i(stdin,prompt,"max_inner_cg_iterations", 
			      &max_inner_cg_iterations );
#else
	max_inner_cg_iterations = 0;
#endif
      }
      
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
	  printf("Unrecognized check command %s\n",savebuf);
	  printf("Choices are 'yes', 'no', 'sourceonly'\n");
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
#if ! defined(HAVE_QOP) && ! defined(USE_CG_GPU) && !defined(HAVE_QPHIX) && ! defined(HAVE_GRID)
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

	if(param.set_type[k] == MULTISOURCE_SET ||
	   param.set_type[k] == MULTICOLORSOURCE_SET){
	  /* Get mass label common to this set */
	  IF_OK status += get_s(stdin,prompt,"mass", savebuf);
#if ( FERM_ACTION == HISQ )
	  IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
				&tmp_naik);
#else
	  tmp_naik = 0.0;
#endif
	} else {
	  /* MULTIMASS_SET or SINGLES_SET */
	  /* Get source index common to this set */
	  IF_OK status += get_i(stdin,prompt,"source", &tmp_src);
	}
      }

      /* Number of propagators in this set */
      IF_OK status += get_i(stdin,prompt,"number_of_propagators",
			    &param.num_prop[k]);
      if( param.num_prop[k]>MAX_PROP ){
	printf("num_prop = %d must be <= %d!\n", param.num_prop[k], MAX_PROP);
	status++;
      }

      if( param.inv_type[k] == MGTYPE && param.set_type[k] == MULTIMASS_SET
	  && param.num_prop[k] > 1){
	node0_printf("WARNING: Multigrid support for multimass is currently emulated via separate inversions\n");
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
	  
	  if(param.set_type[k]  == MULTISOURCE_SET ||
	     param.set_type[k]  == MULTICOLORSOURCE_SET){

	    /* Get source index common to this set */
	    IF_OK status += get_i(stdin,prompt,"source", &param.source[nprop]);
	    strcpy(param.mass_label[nprop], savebuf);
	    param.ksp[nprop].naik_term_epsilon = tmp_naik;
	    
	  } else {

	    /* MULTIMASS_SET or SINGLES_SET */
	    /* Get mass label common to this set */
	    IF_OK status += get_s(stdin,prompt,"mass", param.mass_label[nprop]);
	    
#if ( FERM_ACTION == HISQ )
	    IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
				  &param.ksp[nprop].naik_term_epsilon);
#else
	    param.ksp[nprop].naik_term_epsilon = 0.0;
#endif
	    param.source[nprop] = tmp_src;
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
	
        /* inversion type */
        param.qic[nprop].inv_type = param.inv_type[k];

	/* maximum no. of conjugate gradient iterations */
	param.qic[nprop].max = max_cg_iterations;

	/* maximum no. of inner conjugate gradient iterations */
	/* (QUDA has its own internal definition) */
	param.qic[nprop].max_inner = max_inner_cg_iterations;

	/* maximum no. of conjugate gradient restarts */
	param.qic[nprop].nrestart = max_cg_restarts;

	/* multigrid parameter file name */
	strncpy(param.qic[nprop].mgparamfile, mgparamfile, MAXFILENAME);
      
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
	/* Parameter used by QOPQDP inverter for mixed-precision solves ?? */
	IF_OK status += get_f(stdin, prompt, "mixed_rsq", &param.qic[nprop].mixed_rsq );
#endif

#ifdef MULTIGRID
  /* parameter within MG solve to specify how to refresh the coarse op */

	IF_OK {
	  if (param.inv_type[k] == MGTYPE) {
	    IF_OK status += get_s(stdin, prompt, "rebuild_type", savebuf);
	    IF_OK {
	      if(strcmp(savebuf,"FULL") == 0)
		param.qic[nprop].mg_rebuild_type = FULLREBUILD;
	      else if(strcmp(savebuf,"THIN") == 0)
		param.qic[nprop].mg_rebuild_type = THINREBUILD;
	      else if(strcmp(savebuf,"CG") == 0)
		param.qic[nprop].mg_rebuild_type = CGREBUILD;
	      else {
		printf("Unrecognized rebuild type %s\n",savebuf);
		printf("Choices are 'FULL', 'THIN', 'CG'\n");
		status++;
	      }
	    }
	  }
	}
#else
  param.qic[nprop].mg_rebuild_type = CGREBUILD;
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
      const char *check_tag;
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
	  printf("\nError: expected 'propagator' or 'quark'\n");
	  status++;
	}
      }

      IF_OK init_qss_op(&param.snk_qs_op[i]);
      /* Enforce a uniform boundary condition */
      set_qss_op_offset(&param.snk_qs_op[i], param.coord_origin);
      /* NOTE: The KS built-in bc is antiperiodic. */
      param.snk_qs_op[i].bp[3] = param.time_bc;

      if( param.parent_type[i] == PROP_TYPE ||  param.parent_type[i] == QUARK_TYPE ){
	/* Next we get its index */
	IF_OK {
	  if(prompt==1)printf("enter the index\n");
	  if(scanf("%d",&param.prop_for_qk[i]) != 1){
	    printf("\nFormat error reading index\n");
	    status++;
	  }
	  else{
	    printf("%d\n",param.prop_for_qk[i]);
	    if(param.parent_type[i] == PROP_TYPE &&
	       param.prop_for_qk[i] >= nprop){
	      printf("Propagator index must be less than %d\n", nprop);
	      status++;
	    }
	    else if(param.parent_type[i] == QUARK_TYPE &&
		    param.prop_for_qk[i] >= i){
	      printf("Quark index must be less than %d here\n",i);
	      status++;
	    }
	  }
	}

	/* Get sink operator attributes */
	IF_OK status += get_v_field_op( stdin, prompt, &param.snk_qs_op[i]);

      } else { /* COMBO_TYPE */
	strcpy(param.snk_qs_op[i].descrp, "combination");
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
	    }
	    IF_OK status += get_vf(stdin, prompt, "coeffs", param.combo_coeff[i], param.num_combo[i]);
	  }
	  /* Provisional: For constructing the ancestry for the correlator file. */
	  param.prop_for_qk[i] = param.combo_qk_index[i][0];
	}
      }

      /* Will we save this propagator? */
      IF_OK status += ask_ending_ksprop( stdin, prompt, &param.saveflag_q[i],
					 param.savefile_q[i]);

    }
    

    /*------------------------------------------------------------*/
    /* Meson correlators                                          */
    /*------------------------------------------------------------*/

    /* Number of quark pairs */
    IF_OK status += get_i(stdin,prompt,"number_of_mesons", &param.num_pair );
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
	/* Mesons. The only choice for now. */
	if(strstr(request_buf,",meson,") != NULL)
	  param.do_meson_spect[ipair] = 1;
	else param.do_meson_spect[ipair] = 0;

	/* We need some correlator specification */
	if(!param.do_meson_spect[ipair]){
	  printf("Unrecognized spectrum request\n");
	  status++;
	}
      }

      /* What file for the resulting correlators? */

      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_m[ipair],
				     param.savefile_m[ipair]);
      /* Correlator values are printed with t relative to the t_offset_m */
      /* FT phases are computed with x,y,z relative to r_offset_m */
      IF_OK {
	int r[4];
	status += get_vi(stdin,prompt, "r_offset", r, 4);
	param.r_offset_m[ipair][0] = r[0];
	param.r_offset_m[ipair][1] = r[1];
	param.r_offset_m[ipair][2] = r[2];
	param.r_offset_m[ipair][3] = r[3];
      }

      /*------------------------------------------------------------*/
      /* Table of correlator combinations actually needed           */
      /*------------------------------------------------------------*/

      /* Number of correlators */
      IF_OK status += get_i(stdin,prompt,"number_of_correlators",
			    &param.num_corr_m[ipair] );
      IF_OK {
	if( param.num_corr_m[ipair]>MAX_CORR ){
	  printf("num_corr = %d must be <= %d!\n", param.num_corr_m[ipair],
		 MAX_CORR);
	  status++;
	}
      }

      /* Sample format for correlator line:
	 correlator P5-P5_V1-S_T13_m0.5744 p211 -1 / 4608.0 GX-G1 -2 -1 -1 EO EO EO */

      param.num_corr_report[ipair] = 0;
      IF_OK for(i = 0; i < param.num_corr_m[ipair]; i++){
	int ok,m;
	char meson_label_in[MAX_MESON_LABEL], mom_label_in[MAX_MOM_LABEL],
	  spin_taste_string[16], phase_lab[4], 
	  factor_op[2], parity_x_in[3], parity_y_in[3], parity_z_in[3];
	double factor;

	IF_OK status += get_sn(stdin, prompt, "correlator", meson_label_in);

	/* Read the momentum label next */
	IF_OK {
	  ok = scanf("%s", mom_label_in);
	  if(ok != 1){
	    printf("Error reading momentum label\n");
	    status++;
	  }
	}
	IF_OK printf(" %s", mom_label_in);

	/* Add to hash table */
	m = hash_corr_label(param.meson_label[ipair], param.mom_label[ipair],
			    meson_label_in, mom_label_in,
			    &param.num_corr_report[ipair]);
	param.corr_index[ipair][i] = m;
	if(m < 0)status++;

	/* phase, op, factor, and spin-taste label */
	IF_OK {
	  ok = scanf("%s %s %lf %s",phase_lab,factor_op,&factor,
		     spin_taste_string);
	  if(ok != 4){
	    printf("\nError reading phase, factor, and spin-taste\n");
	    status++;
	  }
	  else {
	    printf(" %3s %1s %g %6s",phase_lab,factor_op,factor,
		   spin_taste_string);
	  }
	}

	/* decode phase for correlator */
	IF_OK {
	  param.meson_phase[ipair][i] = decode_phase(phase_lab);
	  if(param.meson_phase[ipair][i] == -1){
	    printf("\n%s is not a valid phase label\n",phase_lab);
	    status ++;
	  }
	}

	/* decode real factor for correlator */
	/* Permitted syntax is /ddd or *ddd where ddd is a real number */
	IF_OK {
	  param.meson_factor[ipair][i] = decode_factor(factor_op,factor);
	  if(param.meson_factor[ipair][i] == 0 ){
	    printf("\n: Decoding factor %s %g results in zero.\n",
		   factor_op,factor);
	    status++;
	  }
	}

	/* decode spin-taste sink label */
	IF_OK {
	  param.spin_taste_snk[ipair][i] = spin_taste_index(spin_taste_string);
	  if(param.spin_taste_snk[ipair][i] < 0 ){
	    printf("\n: Unrecognized spin-taste label %s.\n",
		   spin_taste_string);
	    status++;
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
    

    /*------------------------------------------------------------*/
    /* Baryon correlators                                          */
    /*------------------------------------------------------------*/

    /* Number of quark triplets */
    IF_OK status += get_i(stdin,prompt,"number_of_baryons", &param.num_triplet );
    if( param.num_triplet>MAX_TRIPLET ){
      printf("number_of_baryons = %d must be <= %d!\n", param.num_triplet, MAX_TRIPLET);
      status++;
    }

    IF_OK for(itriplet = 0; itriplet < param.num_triplet; itriplet++){
      char request_buf[MAX_SPECTRUM_REQUEST];

      /* Which quarks in the triplet? */

      IF_OK status += get_vi(stdin, prompt, "triplet", param.qktriplet[itriplet], 3);

      IF_OK {
	int j;
	for(j = 0; j < 2; j++){
	  if(param.qktriplet[itriplet][j] < 0 || param.qktriplet[itriplet][j] >= param.num_qk){
	    printf("Quark index %d must be in [0,%d]\n",
		   param.qktriplet[itriplet][j],param.num_qk-1);
	    status++;
	  }
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
	param.do_baryon_spect[itriplet] = 0;
	/* Nucleons (The only choice for now) */
	if(strstr(request_buf,",baryon,") != NULL)
	  param.do_baryon_spect[itriplet] = 1;
	if(!param.do_baryon_spect[itriplet]){
	  printf("Unrecognized spectrum request\n");
	  status++;
	}
      }

      /* What file for the resulting correlators? */

      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_b[itriplet],
				     param.savefile_b[itriplet]);
      /* Correlator values are printed with t relative to the t_offset_b */
      /* FT phases are computed with x,y,z relative to r_offset_b */
      IF_OK {
	int r[4];
	status += get_vi(stdin,prompt, "r_offset", r, 4);
	param.r_offset_b[itriplet][0] = r[0];
	param.r_offset_b[itriplet][1] = r[1];
	param.r_offset_b[itriplet][2] = r[2];
	param.r_offset_b[itriplet][3] = r[3];
      }

      /*------------------------------------------------------------*/
      /* Table of correlator combinations actually needed           */
      /*------------------------------------------------------------*/

      /* Number of correlators */

      IF_OK status += get_i(stdin,prompt,"number_of_correlators",
			    &param.num_corr_b[itriplet] );
      IF_OK {
	if( param.num_corr_b[itriplet]>MAX_CORR ){
	  printf("num_corr = %d must be <= %d!\n", param.num_corr_b[itriplet],
		 MAX_CORR);
	  status++;
	}
      }

      /* Sample format for baryon correlator line:
	 correlator  NUCLEON -i * 1 nucleon */
      IF_OK for(i = 0; i < param.num_corr_b[itriplet]; i++){
	int ok;
	char baryon_label_in[MAX_BARYON_LABEL], phase_lab[4], baryon_type_label[32],
	  factor_op[2];
	double factor;

	IF_OK status += get_sn(stdin, prompt, "correlator", baryon_label_in);

	IF_OK strcpy(param.baryon_label[itriplet][i], baryon_label_in);

	/* phase, op, factor, and baryon type label */
	IF_OK {
	  ok = scanf("%s %s %lf %s",phase_lab,factor_op,&factor,
		     baryon_type_label);
	  if(ok != 4){
	    printf("\nError reading phase, factor, and baryon label\n");
	    status++;
	  }
	  else {
	    printf(" %3s %1s %6f %8s\n",phase_lab,factor_op,factor,
		   baryon_type_label);
	  }
	}

	/* decode phase for correlator */
	IF_OK {
	  param.baryon_phase[itriplet][i] = decode_phase(phase_lab);
	  if(param.baryon_phase[itriplet][i] == -1){
	    printf("\n%s is not a valid phase label\n",phase_lab);
	    status ++;
	  }
	}

	/* decode real factor for correlator */
	/* Permitted syntax is /ddd or *ddd where ddd is a real number */
	IF_OK {
	  param.baryon_factor[itriplet][i] = decode_factor(factor_op,factor);
	  if(param.baryon_factor[itriplet][i] == 0 ){
	    printf("\n: Decoding factor %s %g results in zero.\n",
		   factor_op,factor);
	    status++;
	  }
	}

	/* decode baryon type label */
	IF_OK {
	  param.baryon_type_snk[itriplet][i] = baryon_type_index(baryon_type_label);
	  if(param.baryon_type_snk[itriplet][i] < 0 ){
	    printf("\n: Unrecognized baryon-type label %s.\n",
		   baryon_type_label);
	    status++;
	  }
	}
      } /* correlators for this baryon */
    } /* baryons */

		/*------------------------------------------------------------*/
		/*                        GB BARYONS                          */
		/*------------------------------------------------------------*/
		#ifdef GB_BARYON

		    IF_OK status += get_i(stdin,prompt,"number_of_quark_octets", &param.num_oct);
		    if( param.num_oct > MAX_OCTET ){
		      printf("number_of_octets = %d must be <= %d!\n",
		        param.num_oct , MAX_OCTET );
		      status++;
		    }

		    IF_OK if(prompt==1) {
		      printf("Enter indices of quarks in octets\n");
		      printf("Unneeded quarks can substitute index -1\n");
		      printf("-> 0 x y xy z zx yz xyz\n");
		    }
		    IF_OK for(ioctet = 0; ioctet < param.num_oct; ioctet++){

		      IF_OK status += get_vi(stdin, prompt, "octet", param.qk_oct[ioctet], 8);

		      for (iqk = 0; iqk < 8; iqk++) {
		        // index of -1 for unneeded quarks
		        if(param.qk_oct[ioctet][iqk] < -1 || param.qk_oct[ioctet][iqk] >= param.num_qk){
		          printf("Quark index %d must be in [-1,%d]\n",
		          param.qk_oct[itriplet][iqk],param.num_qk-1);
		          status++;
		        }
		      }

		    }

		    /* Number of quark triplets (triplets of octets)*/
		    IF_OK status += get_i(stdin,prompt,"number_of_gb_baryons", &param.num_gb_triplet );
		    if( param.num_gb_triplet > MAX_TRIPLET ){
		      printf("number_of_gb_baryons = %d must be <= %d!\n",
		        param.num_gb_triplet, MAX_TRIPLET);
		      status++;
		    }

		    IF_OK for(itriplet = 0; itriplet < param.num_gb_triplet; itriplet++){
		      char request_buf[MAX_SPECTRUM_REQUEST];

		      /* Which quark OCTETS in the triplet? */

		      IF_OK status += get_vi(stdin, prompt, "triplet", param.qk8triplet[itriplet], 3);

		      IF_OK {
			int j;
			for(j = 0; j < 2; j++){
			  if(param.qk8triplet[itriplet][j] < 0 || param.qk8triplet[itriplet][j] >= param.num_oct){
			    printf("Quark octet index %d must be in [0,%d]\n",
				   param.qk8triplet[itriplet][j],param.num_oct-1);
			    status++;
			  }
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
					param.do_gbbaryon_spect[itriplet] = 0;
					if(strstr(request_buf,",gb_baryon,") != NULL)
					  param.do_gbbaryon_spect[itriplet] = 1;
					if(!param.do_gbbaryon_spect[itriplet]){
					  printf("Unrecognized spectrum request\n");
					  status++;
					}
	      }

	      /* What is the triplet quark content? */
	      char qk_content_in[4];
	      IF_OK status += get_sn(stdin,prompt,"quark_content",qk_content_in);
	      IF_OK {
	        param.qk8num_d[itriplet] = 0;
	        param.qk8num_s[itriplet] = 0;
	        int j;
	        for (j = 0; j < 3; j++) {
	          switch (qk_content_in[j])
	          {
	            case 'u': continue; break;
	            case 'd': param.qk8num_d[itriplet]++;
	                      continue; break;
	            case 's': param.qk8num_s[itriplet]++;
	                      continue; break;
	            default:  printf("error in input: unknown quark content %c\n",qk_content_in[j]);
	                      status++; break;
	          }
	        }
	      }
	      IF_OK printf("\n");

	      /* What file for the resulting correlators? */
	      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_gb[itriplet],
					     param.savefile_gb[itriplet]);
	      /* Correlator values are printed with t relative to the t_offset_b */
	      /* FT phases are computed with x,y,z relative to r_offset_b */
	      IF_OK {
					int r[4];
					status += get_vi(stdin,prompt, "r_offset", r, 4);
					param.r_offset_gb[itriplet][0] = r[0];
					param.r_offset_gb[itriplet][1] = r[1];
					param.r_offset_gb[itriplet][2] = r[2];
					param.r_offset_gb[itriplet][3] = r[3];
			  }

	      IF_OK status += get_vi(stdin, prompt, "momentum", param.snkmom_gb[itriplet], 3);

	      /* Get momentum to include in sink*/
	      IF_OK {
	        printf(" %2d %2d %2d",
	          param.snkmom_gb[itriplet][0],
	      	  param.snkmom_gb[itriplet][1],
	      	  param.snkmom_gb[itriplet][2]);
	      }

	      /* Parity for each momentum component */
	      int ok;
	      char parity_x_in[3], parity_y_in[3], parity_z_in[3];

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
	        param.snkpar_gb[itriplet][0] = decode_parity(parity_x_in);
	        param.snkpar_gb[itriplet][1] = decode_parity(parity_y_in);
	        param.snkpar_gb[itriplet][2] = decode_parity(parity_z_in);
	        if(param.snkpar_gb[itriplet][0] == 0xf)status++;
	        if(param.snkpar_gb[itriplet][1] == 0xf)status++;
	        if(param.snkpar_gb[itriplet][2] == 0xf)status++;
	      }

	      char gb_spin_taste_in[8];
	      IF_OK status += get_sn(stdin,prompt,"spin_taste",gb_spin_taste_in);

	      /* decode spin-taste sink label */
	      IF_OK {
	        if(strcmp("2point",gb_spin_taste_in) == 0){
	         /* use scalar-scalar spin taste */
	         //param.gb_spintaste[itriplet] = spin_taste_index("G1-G1"); // scalar-scalar
	         param.gb_spintaste[itriplet] = GB_2POINT_BACKPROP; // will use scalar-scalar
	        } else {
	         param.gb_spintaste[itriplet]
	           = spin_taste_index(gb_spin_taste_in);
	        }
	        if(param.gb_spintaste[itriplet] < 0 ){
	          printf("\n: Unrecognized spin-taste label %s.\n",
	                 gb_spin_taste_in);
	          status++;
		        }
		      }

		      /*------------------------------------------------------------*/
		      /* Table of correlator combinations actually needed           */
		      /*------------------------------------------------------------*/

		      /* Number of correlators */

		      IF_OK status += get_i(stdin,prompt,"number_of_correlators",
					    &param.num_corr_gb[itriplet] );
		      IF_OK {
							if( param.num_corr_gb[itriplet]>MAX_CORR ){
							  printf("num_corr = %d must be <= %d!\n",
						            param.num_corr_gb[itriplet], MAX_CORR);
							  status++;
							}
		      }

			/* Sample format for baryon correlator line:
			 correlator  NUCLEON -i * 1 16+ M0 2 16+ M0 42 cube
			 correlator  NUCLEON -i * 1 16+ M0 2 M0 42 cube  <- old format x1
			 correlator  NUCLEON -i * 1 16+ M0 2 M0 42 <- old format x2
				 correlator  NUCLEON -i * 1 nucleon  <- old format x3
			*/
		  IF_OK for(i = 0; i < param.num_corr_gb[itriplet]; i++){
			int ok;
			char baryon_label_in[MAX_BARYON_LABEL],
		             phase_lab[4], factor_op[2], //gts_label[4],
		             gts_label_src[4], gts_label_snk[4],
		             symiso_label_src[6], symiso_label_snk[7], sink_tie[7],
					 cube_pos[7];
			double factor;
		        int gb_class_src, gb_class_snk;

			IF_OK status += get_sn(stdin, prompt, "correlator", baryon_label_in);

			IF_OK strcpy(param.gbbaryon_label[itriplet][i], baryon_label_in);

			/* phase, op, factor, GTSirrep, symiso, class, sink tieup type, construction*/
			IF_OK {
			  ok = scanf("%s %s %lf %s %s %i %s %s %i %s %s",
		                     phase_lab,factor_op,&factor,
		                     gts_label_src,symiso_label_src,&gb_class_src,
		                     gts_label_snk,symiso_label_snk,&gb_class_snk,
							 sink_tie, cube_pos);

			  if(ok != 11){
			    printf("\nError reading Golterman-Bailey baryon parameters\n");
			    status++;
			  }
			  else {
			    printf(" %3s %1s %6f %4s %6s %2i %4s %7s %2i %6s %6s\n",
		                   phase_lab,factor_op,factor,
		                   gts_label_src,symiso_label_src,gb_class_src,
		                   gts_label_snk,symiso_label_snk,gb_class_snk,
						   sink_tie, cube_pos);
			  }
			}

			/* decode phase for correlator */
			IF_OK {
			  param.gbbaryon_phase[itriplet][i] = decode_phase(phase_lab);
			  if(param.gbbaryon_phase[itriplet][i] == -1){
			    printf("\n%s is not a valid phase label\n",phase_lab);
			    status ++;
			  }
			}

			/* Decode real factor for correlator */
			/* Permitted syntax is /ddd or *ddd where ddd is a real number */
			IF_OK {
			  param.gbbaryon_factor[itriplet][i] = decode_factor(factor_op,factor);
			  if(param.gbbaryon_factor[itriplet][i] == 0 ){
			    printf("\n: Decoding factor %s %g results in zero.\n",
				   factor_op,factor);
			    status++;
			  }
			}

		  /* Check input symmetry and isospin */
			IF_OK {
		          if(strcmp("S",     symiso_label_src) == 0 ||
		             strcmp("A",     symiso_label_src) == 0 ||
		             strcmp("M1/2",  symiso_label_src) == 0   ) {
			       } /* match is okay */
			  else {
		            printf("error in input: unknown operator symmetry %s\n",symiso_label_src);
			    status++;
		          }
		        }
			IF_OK {
		          if(strcmp("S",     symiso_label_snk) == 0 ||
		             strcmp("A",     symiso_label_snk) == 0 ||
		             strcmp("M1/2",  symiso_label_snk) == 0 ||
		             strcmp("S*",    symiso_label_snk) == 0 ||
		             strcmp("A*",    symiso_label_snk) == 0 ||
		             strcmp("M1/2*", symiso_label_snk) == 0   ) {
			       } /* match is okay */
			  else  {
		            printf("error in input: unknown operator symmetry %s\n",symiso_label_snk);
			    status++;
		          }
		        }

			IF_OK {
		          if(strcmp("S*",    symiso_label_snk) == 0 ||
		             strcmp("A*",    symiso_label_snk) == 0 ||
		             strcmp("M1/2*", symiso_label_snk) == 0   ){
		            param.gb_parity_flip_snk[itriplet][i] = 1;
		          } else {
		            param.gb_parity_flip_snk[itriplet][i] = 0;
		          }
		        }

		  /* Check input GTS irrep */
			IF_OK {
		          if (strcmp("8",  gts_label_src) == 0 ||
		              strcmp("8'", gts_label_src) == 0 ||
		              strcmp("16+",gts_label_src) == 0 ||
		              strcmp("16-",gts_label_src) == 0) {
			       } /* match is okay */
			  else {
			    printf("error in input: unknown GTS irrep %s\n",gts_label_src);
			    status++;
		          }
		        }

		  /* Check input GTS irrep */
	IF_OK {
          if (strcmp("8",  gts_label_snk) == 0 ||
              strcmp("8'", gts_label_snk) == 0 ||
              strcmp("16+",gts_label_snk) == 0 ||
              strcmp("16-",gts_label_snk) == 0) {
	      }  /* match is okay */
	  else {
	      printf("error in input: unknown GTS irrep %s\n",gts_label_snk);
	      status++; }
	}

      /* Check input source class */
      IF_OK {
        if (gb_class_src != 1  &&
            gb_class_src != 2  &&
            gb_class_src != 3  &&
            gb_class_src != 41 &&
            gb_class_src != 42 &&
            gb_class_src != 5  &&
            gb_class_src != 61 &&
            gb_class_src != 62 &&
            gb_class_src != 7){
          printf("error in input: unknown source class %i\n",gb_class_src);
          status++;
        			}
      }

      /* Check input sink class */
      IF_OK {
        if (gb_class_snk != 1  &&
            gb_class_snk != 2  &&
            gb_class_snk != 3  &&
            gb_class_snk != 41 &&
            gb_class_snk != 42 &&
            gb_class_snk != 5  &&
            gb_class_snk != 61 &&
            gb_class_snk != 62 &&
            gb_class_snk != 7){
          	printf("error in input: unknown sink class %i\n",gb_class_snk);
          	status++;
	    }
       }

      IF_OK {
        param.gbbaryon_src[itriplet][i] =
          decode_gb_op(symiso_label_src,gts_label_src,param.qk8num_s[itriplet],gb_class_src);
        if (param.gbbaryon_src[itriplet][i] == GB_UNDEFINED){
          printf("error: undefined baryon source\n");
          status++;
			        }
			      }

      IF_OK {
        param.gbbaryon_snk[itriplet][i] =
          decode_gb_op(symiso_label_snk,gts_label_snk,param.qk8num_s[itriplet],gb_class_snk);
        if (param.gbbaryon_snk[itriplet][i] == GB_UNDEFINED){
          printf("error: undefined baryon sink\n");
          status++;
				      }
				    }
      /* Decode sink tieup options for correlator */
      /* Permitted choices are 'point' or 'wall' */
      IF_OK {
        if     (strcmp("point", sink_tie) == 0 ){ param.gb_wall[itriplet][i] = 0x0; }
        else if(strcmp("wall", sink_tie) == 0 ){ param.gb_wall[itriplet][i] = 0x1; }
        else {
          printf("\n: Invalid operator construction %s.\n", sink_tie);
          status++;
				      }
				    }

      /* Decode cube construction for correlator */
      /* Permitted choices are 'corner' or 'cube' */
      IF_OK {
        if     (strcmp("corner",cube_pos) == 0 ){ param.gb_corner[itriplet][i] = 0x0; }
        else if(strcmp("cube",  cube_pos) == 0 ){ param.gb_corner[itriplet][i] = 0x1; }
        else {
          printf("\n: Invalid operator construction %s.\n",cube_pos);
          status++;
				      }
				    }

		      } /* Correlators for this baryon */
		    } /* gb baryons */
		#endif /* GB_BARYON */

    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
  
  
  fflush(stdout);
  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0;
  if( param.stopflag != 0 )return param.stopflag;

  if(prompt==2)return 0;

  /* Broadcast parameter values kept on the heap */
  broadcast_heap_params();
  fflush(stdout);

#ifdef ANISOTROPY
    /* Direction of anisotropy */
    ani_dir = param.ani_dir;
    /* Bare fermion anisotropy */
    ani_xiq = param.ani_xiq;
#ifdef ONEDIM_ANISO_TEST
    /* bare quark isotropic link factor for debugging */
    iso_xiq = param.iso_xiq;
#endif
#endif
#ifdef CORDIR
  cor_dir = param.cor_dir;
#endif


  /* Construct the eps_naik table of unique Naik epsilon coefficients.
     Also build the hash table for mapping a mass term to its Naik
     epsilon index.  We need to list any required Naik epsilons. */

  /* First term is always zero */
  start_eps_naik(param.eps_naik, &param.n_naiks);

  /* Contribution from the chiral condensate epsilons */
  for(i = 0; i < param.num_pbp_masses; i++){
    param.ksp_pbp[i].naik_term_epsilon_index =
      fill_eps_naik(param.eps_naik,
	  &param.n_naiks, param.ksp_pbp[i].naik_term_epsilon);
  }

  /* Contribution from the propagator epsilons */
  if(param.num_set > 0){
    nprop = param.end_prop[param.num_set-1] + 1;
    for(i = 0; i < nprop; i++)
      param.ksp[i].naik_term_epsilon_index = 
	fill_eps_naik(param.eps_naik, 
		      &param.n_naiks, param.ksp[i].naik_term_epsilon);
  }

  /* Requests from any embedded inverse and hopping operators in the
     modified ops */
  for(i = 0; i < param.num_modified_source; i++){
    Real eps = 0.;
    int is = param.num_base_source + i;
    /* If the operator uses Dslash, get the requested Naik epsilon */
    if(get_qss_eps_naik(&eps, &param.src_qs_op[is])){
      insert_qss_eps_naik_index(fill_eps_naik(param.eps_naik, &param.n_naiks, eps), &param.src_qs_op[is]);
    }
  }

  /* Requests from any embedded inverse and hopping operators in
     the sink ops */
  for(i = 0; i < param.num_qk; i++){
    Real eps = 0.;
    /* If the operator uses Dslash, get the requested Naik epsilon */
    if(get_qss_eps_naik(&eps, &param.snk_qs_op[i])){
      insert_qss_eps_naik_index(fill_eps_naik(param.eps_naik, &param.n_naiks, eps), &param.snk_qs_op[i]);
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


  ENDTIME("read parameters");

  return 0;
}

/* Broadcast operator parameter values.  They are on the heap on node 0. */

static void broadcast_heap_params(void){
  int i;

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

/*------------------------------------------------------------*/
/* Set up comlink structures for the 2nd nearest gather pattern; 
   make_lattice(), make_nn_gathers() and make_3n_gathers() must be called before.
*/
/* Hwancheol Jeong 4/2024 */

static void
make_2n_gathers(void)
{
  int i;
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(second_neighbor, &i, WANT_INVERSE,
                ALLOW_EVEN_ODD, SAME_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X2UP, X2DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X2UP);
}

/*------------------------------------------------------------*/
/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 2nd nearest neighbor in that direction */
/* Hwancheol Jeong 4/2024 */

static void
second_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
                int *xp, int *yp, int *zp, int *tp)
/* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
   "forwards/backwards"  */
/* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+2)%nx; break;
  case XDOWN: *xp = (x+2*nx-2)%nx; break;
  case YUP: *yp = (y+2)%ny; break;
  case YDOWN: *yp = (y+2*ny-2)%ny; break;
  case ZUP: *zp = (z+2)%nz; break;
  case ZDOWN: *zp = (z+2*nz-2)%nz; break;
  case TUP: *tp = (t+2)%nt; break;
  case TDOWN: *tp = (t+2*nt-2)%nt; break;
  default: printf("second_neighb: bad direction\n"); exit(1);
  }
}

