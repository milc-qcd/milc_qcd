/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "ks_measure_includes.h"
#include <string.h>
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h
#include "params.h"
#ifdef U1_FIELD
#include "../include/generic_u1.h"
#include "../include/io_u1lat.h"
#endif

/* Forward declarations */

static int initial_set();
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers();


int setup()   {
  int prompt;

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
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  return(prompt);
}

static int n_naiks = 1;
static double eps_naik[MAX_NAIK];
static double charge[MAX_CHARGE];

/* SETUP ROUTINES */
static int 
initial_set(){
  int prompt=0,status;

  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 staggered fermion measurements\n");
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
  for(int i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(int i = 0; i < 4; i++)
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
  int i, k, npbp_masses = 0;
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
    IF_OK status += get_f(stdin, prompt, "staple_weight", 
			  &param.staple_weight);
    IF_OK status += get_i(stdin, prompt, "ape_iter",
			  &param.ape_iter);

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

#if EIGMODE == EIGCG
      /* for eigcg */

      /* maximum number of eigenvectors */
      param.eigcgp.Nvecs_max =  param.eigen_param.Nvecs;

      /* If we are reading in eigenpairs, we don't regenerate them */

      if(param.ks_eigen_startflag == FRESH){
	
	/* restart for Lanczos */
	IF_OK status += get_i(stdin, prompt,"restart_lanczos", &param.eigcgp.m);
	
	/* number of eigenvectors per inversion */
	IF_OK status += get_i(stdin, prompt,"Number_of_eigenvals", &param.eigcgp.Nvecs);
	
	if(param.eigcgp.m <= 2*param.eigcgp.Nvecs){
	  printf("restart_lanczos should be larger than 2*Number_of_eigenvals!\n");
	  status++;
	}
      } else {
	param.eigcgp.m = 0;
	param.eigcgp.Nvecs = 0;
      }
      
      param.eigcgp.Nvecs_curr = 0;
      param.eigcgp.H = NULL;
#else

      /*------------------------------------------------------------*/
      /* Dirac eigenpair calculation                                */
      /*------------------------------------------------------------*/
      
      if(param.ks_eigen_startflag == FRESH){
	
	status += read_ks_eigen_param(&param.eigen_param, status, prompt);

      }
    }

#endif

    /*------------------------------------------------------------*/
    /* Chiral condensate and related quantities                   */
    /*------------------------------------------------------------*/
    
    IF_OK status += get_i(stdin,prompt,"number_of_sets", &param.num_set);
    if( param.num_set>MAX_SET ){
      printf("num_set = %d must be <= %d!\n", param.num_set, MAX_SET);
      status++;
    }

    npbp_masses = 0;
    IF_OK for(k = 0; k < param.num_set; k++){
      int max_cg_iterations, max_cg_restarts, prec_pbp;
      Real error_for_propagator, rel_error_for_propagator;
#ifdef CURRENT_DISC
      int max_cg_iterations_sloppy, max_cg_restarts_sloppy, prec_pbp_sloppy;
      Real error_for_propagator_sloppy, rel_error_for_propagator_sloppy;
#endif

      /* Number of stochastic sources */
      IF_OK status += get_i(stdin, prompt, "npbp_reps", &param.npbp_reps[k] );

#ifdef CURRENT_DISC
      /* For some applications.  Random source count between writes */
      IF_OK status += get_i(stdin, prompt, "source_spacing", &param.thinning[k] );
      IF_OK {
	if(param.thinning[k] < 2){
	  printf("Source spacing must be at least 2\n");
	  status++;
	}
      }
      /* For truncated solver Take difference of sloppy and precise?*/
      char savebuf[128];
      IF_OK status += get_s(stdin, prompt, "take_truncate_diff", savebuf);
      IF_OK {
	if(strcmp(savebuf,"no") == 0)param.truncate_diff[k] = 0;
	else if(strcmp(savebuf,"yes") == 0)param.truncate_diff[k] = 1;
	else {
	  printf("Unrecognized response %s\n",savebuf);
	  printf("Choices are 'yes' and 'no'\n");
	  status++;
	}
      }
#endif

      /* The following parameters are common to the set and will be copied
	 to each member */

      /* maximum no. of conjugate gradient iterations */
      IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			    &max_cg_iterations );
      
      /* maximum no. of conjugate gradient restarts */
      IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			    &max_cg_restarts );

      IF_OK status += get_i(stdin, prompt, "prec_pbp", 
			    &prec_pbp );

#ifdef CURRENT_DISC
      /* If we are taking the difference between a sloppy and a precise solve,
	 get the sloppy solve parameters */
      if(param.truncate_diff[k]){
	Real error_for_propagator_sloppy, rel_error_for_propagator_sloppy;

	/* The following parameters are common to the set and will be copied
	   to each member */
	
	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin,prompt,"max_cg_iterations_sloppy", 
			      &max_cg_iterations_sloppy );
	
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin,prompt,"max_cg_restarts_sloppy", 
			      &max_cg_restarts_sloppy );
	
	IF_OK status += get_i(stdin, prompt, "prec_pbp_sloppy", 
			      &prec_pbp_sloppy );
	
      }
#endif

      /* Number of pbp masses in this set */
      IF_OK status += get_i(stdin, prompt, "number_of_pbp_masses",
			    &param.num_pbp_masses[k]);

      /* Indexing range for set */
      param.begin_pbp_masses[k] = npbp_masses;
      param.end_pbp_masses[k] = npbp_masses + param.num_pbp_masses[k] - 1;
      if(param.end_pbp_masses[k] > MAX_MASS_PBP){
	printf("Total number of masses must be <= %d!\n", MAX_MASS_PBP);
	status++;
      }

      IF_OK for(int i = 0; i < param.num_pbp_masses[k]; i++){
    
	/* PBP mass parameters */
	
	IF_OK status += get_s(stdin, prompt,"mass", param.mass_label[npbp_masses] );
	IF_OK param.ksp_pbp[npbp_masses].mass = atof(param.mass_label[npbp_masses]);
#if ( FERM_ACTION == HISQ )
	IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
			      &param.ksp_pbp[npbp_masses].naik_term_epsilon );
#else
	IF_OK param.ksp_pbp[npbp_masses].naik_term_epsilon = 0.0;
#endif
#ifdef U1_FIELD
	IF_OK status += get_s(stdin, prompt,"charge", param.charge_label[npbp_masses] );
	IF_OK param.ksp_pbp[npbp_masses].charge = atof(param.charge_label[npbp_masses]);
#else
	IF_OK strcpy(param.charge_label[npbp_masses],"0.");
	IF_OK param.ksp_pbp[npbp_masses].charge = 0.;
#endif
	/* error for staggered propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			      &error_for_propagator );
	IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			      &rel_error_for_propagator );
	
#ifdef CURRENT_DISC
	if(param.truncate_diff[k]){
	  /* error for staggered propagator conjugate gradient */
	  IF_OK status += get_f(stdin, prompt,"error_for_propagator_sloppy", 
				&error_for_propagator_sloppy );
	  IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator_sloppy", 
				&rel_error_for_propagator_sloppy );
	}

#if 0 	/* Deprecate saving the entire site-wise density  */
	IF_OK status += get_s(stdin, prompt, "save_file", param.pbp_filenames[npbp_masses] );
#endif
#endif

	/* The set to which this pbp_mass belongs */
	IF_OK param.set[npbp_masses] = k;

        /* inversion type */
        param.qic_pbp[npbp_masses].inv_type = UMLTYPE;

	/* maximum no. of conjugate gradient iterations */
	param.qic_pbp[npbp_masses].max = max_cg_iterations;
      
	/* maximum no. of conjugate gradient restarts */
	param.qic_pbp[npbp_masses].nrestart = max_cg_restarts;
      
	/* precision */
	param.qic_pbp[npbp_masses].prec = prec_pbp;

	/* errors */
	param.qic_pbp[npbp_masses].resid = error_for_propagator;
	param.qic_pbp[npbp_masses].relresid = rel_error_for_propagator;

#ifdef CURRENT_DISC
	param.qic_pbp[npbp_masses].parity = EVEN;
#else
	param.qic_pbp[npbp_masses].parity = EVENANDODD;
#endif
	param.qic_pbp[npbp_masses].min = 0;
	param.qic_pbp[npbp_masses].start_flag = 0;
	param.qic_pbp[npbp_masses].nsrc = 1;

	/* Should we be deflating? */
	param.qic_pbp[npbp_masses].deflate = 0;
	IF_OK {
	  /* Always deflate if we have eigenvectors */
	  if(param.eigen_param.Nvecs > 0)param.qic_pbp[npbp_masses].deflate = 1;
	}

#ifdef CURRENT_DISC
      /* If we are taking the difference between a sloppy and a precise solve,
	 get the sloppy solve parameters */
	if(param.truncate_diff[k]){
	  
	  /* maximum no. of conjugate gradient iterations */
	  param.qic_pbp_sloppy[npbp_masses].max = max_cg_iterations_sloppy;
	  
	  /* maximum no. of conjugate gradient restarts */
	  param.qic_pbp_sloppy[npbp_masses].nrestart = max_cg_restarts_sloppy;
	  
	  /* precision */
	  param.qic_pbp_sloppy[npbp_masses].prec = prec_pbp_sloppy;
	  
	  /* errors */
	  param.qic_pbp_sloppy[npbp_masses].resid = error_for_propagator_sloppy;
	  param.qic_pbp_sloppy[npbp_masses].relresid = rel_error_for_propagator_sloppy;
	  
	  param.qic_pbp_sloppy[npbp_masses].parity = EVENANDODD;
	  param.qic_pbp_sloppy[npbp_masses].min = 0;
	  param.qic_pbp_sloppy[npbp_masses].start_flag = 0;
	  param.qic_pbp_sloppy[npbp_masses].nsrc = 1;
	  
	}
	
#endif
	
	npbp_masses++;
      } /* i */
    } /* k */
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    

  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0;
  if( param.stopflag != 0 )return param.stopflag;

  if(prompt==2)return 0;

  npbp_masses = param.end_pbp_masses[param.num_set-1] + 1;

  /* Construct the eps_naik table of unique Naik epsilon
     coefficients.  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* First term is always zero */
  start_eps_naik(eps_naik, &n_naiks);
  
  /* Contribution from the chiral condensate epsilons */
  for(int i = 0; i < npbp_masses; i++){
    param.ksp_pbp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, 
		    &n_naiks, param.ksp_pbp[i].naik_term_epsilon);
  }

  /* Construct a table of quark charges and build a hash table
     for mapping a charge to its charge index */

  /* First term is always zero */
  start_charge(charge, &n_charges);
  
  /* Contribution from the chiral condensate charges */
  for(int i = 0; i < npbp_masses; i++){
    param.ksp_pbp[i].charge_index = 
      fill_charge(charge, &n_charges, param.ksp_pbp[i].charge);
  }

  /* Keep the current lattice or load a new one */
  if( param.startflag == CONTINUE ){
    rephase( OFF );
  } else {
    startlat_p = reload_lattice( param.startflag, param.startfile );
    phases_in = OFF;
  }

#if 0
  su3_matrix *G = create_random_m_field();
  gauge_transform_links(G);
  d_plaquette(&g_ssplaq,&g_stplaq);
  d_linktrsum(&linktrsum);
  nersc_checksum = nersc_cksum();
  node0_printf("CHECK PLAQ: %.16e %.16e\n",g_ssplaq,g_stplaq);
  node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
	       linktrsum.real/3.,nersc_checksum);
#endif

  /* By defailt KS and BC phases are in the gauge links */
  rephase( ON );

#ifdef U1_FIELD
  /* Read the U(1) gauge field, if wanted */
  start_u1lat_p = reload_u1_lattice( param.start_u1flag, param.start_u1file);
#endif

  /* Set options for fermion links */
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif

#ifdef DM_DU0
  fermion_links_want_du0(1);
#endif
  
  /* Don't need to save HISQ auxiliary links */
  fermion_links_want_aux(0);

#if FERM_ACTION == HISQ && defined(DM_DEPS)
  fermion_links_want_deps(1);
#endif
  
  /* Create an array of fermion links structures for all unique Naik epsilon values and charges */
  fn_links_charge = (fermion_links_t **)malloc(sizeof(fermion_links_t *) * n_charges);
  for(int i = 0; i < n_charges; i++){
#ifdef U1_FIELD
    if(charge[i] != 0)u1phase_on(charge[i], u1_A);
#endif
    /* Create a set of fermion links for all eps_naik */
    fn_links_charge[i] = create_fermion_links_from_site(MILC_PRECISION, n_naiks, eps_naik);
#ifdef U1_FIELD
    if(charge[i] != 0)u1phase_off();
#endif
  }
  /* For compatibility. The first charge is always zero */
  fn_links = fn_links_charge[0];

  /* Construct APE smeared links, but without KS phases */
  rephase( OFF );
  ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
  /* We put in antiperiodic bc to the APE links to match what we did to the gauge field */
  apply_apbc( ape_links, 0 );
  rephase( ON );
  refresh_ape_links = 1;
  /* Put the KS phases into APE links to match what we did to the gauge field */
  ape_links_ks_phases = OFF;
  rephase_field_offset( ape_links, ON, &ape_links_ks_phases, param.coord_origin );

#if EIGMODE == EIGCG
  int Nvecs_max = param.eigcgp.Nvecs_max;
  if(param.ks_eigen_startflag == FRESH)
    //    Nvecs_tot = ((Nvecs_max - 1)/param.eigcgp.Nvecs)*param.eigcgp.Nvecs
    //      + param.eigcgp.m;
    Nvecs_tot = Nvecs_max + param.eigcgp.m - 1;
  else
    Nvecs_tot = Nvecs_max;

  eigVal = (double *)malloc(Nvecs_tot*sizeof(double));
  eigVec = (su3_vector **)malloc(Nvecs_tot*sizeof(su3_vector *));
  node0_printf("Allocating space for %d eigenvectors\n", Nvecs_tot);
  for(int i = 0; i < Nvecs_tot; i++)
    eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /* Do whatever is needed to get eigenpairs -- assumed charge 0 */
  imp_ferm_links_t **fn = get_fm_links(fn_links);
  status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
			   &Nvecs_tot, eigVal, eigVec, fn[0], 1);
  if(status != 0) terminate(1);
  //  if(param.fixflag != NO_GAUGE_FIX){
  //    node0_printf("WARNING: Gauge fixing does not readjust the eigenvectors\n");
  //  }

  if(param.ks_eigen_startflag != FRESH){
    param.eigcgp.Nvecs = 0;
    param.eigcgp.Nvecs_curr = Nvecs_tot;
    param.eigcgp.H = (double_complex *)malloc(Nvecs_max*Nvecs_max
					      *sizeof(double_complex));
    for(int i = 0; i < Nvecs_max; i++){
      for(k = 0; k < i; k++)
	param.eigcgp.H[k + Nvecs_max*i] = dcmplx((double)0.0, (double)0.0);
      param.eigcgp.H[(Nvecs_max+1)*i] = dcmplx(eigVal[i], (double)0.0);
    }
  }
#endif
  
#if EIGMODE != EIGCG
  if(param.eigen_param.Nvecs > 0){
    /* malloc for eigenpairs */
    eigVal = (double *)malloc(param.eigen_param.Nvecs*sizeof(double));
    eigVec = (su3_vector **)malloc(param.eigen_param.Nvecs*sizeof(su3_vector *));
    for(i=0; i < param.eigen_param.Nvecs; i++){
      eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      if(eigVec[i] == NULL){
	printf("No room for eigenvector\n");
	terminate(1);
      }
    }
    
    /* Do whatever is needed to get eigenpairs -- assumed charge 0 */
    imp_ferm_links_t **fn = get_fm_links(fn_links);
    status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
			     &param.eigen_param.Nvecs, eigVal, eigVec, fn[0], 1);
    if(status != 0)terminate(1);
#if 0
    for(int j = 0; j < param.eigen_param.Nvecs; j++){
      gauge_transform_v_field(eigVec[j], G);
    }
    destroy_m_field(G);
#endif
  }
#endif

  ENDTIME("readin");
  fflush(stdout);

  return 0;
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
