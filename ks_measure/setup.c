/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ... */

//  $Log: setup.c,v $
//  Revision 1.2  2012/11/24 05:21:42  detar
//  Add support for future HYPISQ action
//
//  Revision 1.1  2011/12/02 04:38:15  detar
//  Add
//
//


#include "ks_measure_includes.h"
#include <string.h>
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h
#include "params.h"

/* Forward declarations */

static int initial_set();
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers();


int setup()   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;

  /* initialize the node random number generator */
  initialize_prn( &node_prn, param.iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* Initialize fermion links as unallocated */
//  init_ferm_links(&fn_links, &ks_act_paths);
//  init_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);
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

/* SETUP ROUTINES */
static int 
initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
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
  int i, k, npbp_masses;
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

    /* Provision is made to build covariant sources from smeared
       links */
    /* APE smearing parameters (if needed) */
    /* Zero suppresses APE smearing */
    IF_OK status += get_f(stdin, prompt, "staple_weight", 
			  &param.staple_weight);
    IF_OK status += get_i(stdin, prompt, "ape_iter",
			  &param.ape_iter);

#if EIGMODE == EIGCG
    /* for eigcg */
    /* restart for Lanczos */
    IF_OK status += get_i(stdin, prompt,"restart_lanczos", &param.eigcgp.m);

    /* number of eigenvectors per inversion */
    IF_OK status += get_i(stdin, prompt,"Number_of_eigenvals", &param.eigcgp.Nvecs);

    if(param.eigcgp.m <= 2*param.eigcgp.Nvecs){
      printf("restart_lanczos should be larger than 2*Number_of_eigenvals!\n");
      status++;
    }

    /* maximum number of eigenvectors */
    IF_OK status += get_i(stdin, prompt,"Max_Number_of_eigenvals",
			  &param.eigcgp.Nvecs_max);

    /* eigenvector input */
    IF_OK status += ask_starting_ks_eigen(stdin, prompt, &param.ks_eigen_startflag,
					  param.ks_eigen_startfile);

    /* eigenvector output */
    IF_OK status += ask_ending_ks_eigen(stdin, prompt, &param.ks_eigen_saveflag,
					param.ks_eigen_savefile);

    param.eigcgp.Nvecs_curr = 0;
    param.eigcgp.H = NULL;
#endif

#if EIGMODE == DEFLATION
    /*------------------------------------------------------------*/
    /* Dirac eigenpair calculation                                */
    /*------------------------------------------------------------*/

    /* number of eigenvectors */
    IF_OK status += get_i(stdin, prompt,"Number_of_eigenvals", &param.Nvecs);

    /* max  Rayleigh iterations */
    IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &param.MaxIter);

    /* Restart  Rayleigh every so many iterations */
    IF_OK status += get_i(stdin, prompt,"Restart_Rayleigh", &param.Restart);

    /* Kalkreuter iterations */
    IF_OK status += get_i(stdin, prompt,"Kalkreuter_iters", &param.Kiters);

     /* Tolerance for the eigenvalue computation */
    IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", &param.eigenval_tol);

     /* error decrease per Rayleigh minimization */
    IF_OK status += get_f(stdin, prompt,"error_decrease", &param.error_decr);

    /* eigenvector input */
    IF_OK status += ask_starting_ks_eigen(stdin, prompt, &param.ks_eigen_startflag,
					  param.ks_eigen_startfile);

    /* eigenvector output */
    IF_OK status += ask_ending_ks_eigen(stdin, prompt, &param.ks_eigen_saveflag,
					param.ks_eigen_savefile);
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
      IF_OK status += get_i(stdin, prompt, "nwrite", &param.nwrite[k] );
      IF_OK status += get_i(stdin, prompt, "source_spacing", &param.thinning[k] );
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
      if(param.num_pbp_masses[k] > MAX_MASS_PBP){
	printf("Number of masses exceeds dimension %d\n",MAX_MASS_PBP);
	status++;
      }

      /* Indexing range for set */
      param.begin_pbp_masses[k] = npbp_masses;
      param.end_pbp_masses[k] = npbp_masses + param.num_pbp_masses[k] - 1;
      if(param.end_pbp_masses[k] > MAX_PBP_MASSES){
	printf("Total number of masses must be <= %d!\n", MAX_PBP_MASSES);
	status++;
      }

      IF_OK for(i = 0; i < param.num_pbp_masses[k]; i++){
    
	/* PBP mass parameters */
	
	IF_OK status += get_s(stdin, prompt,"mass", param.mass_label[npbp_masses] );
	IF_OK param.ksp_pbp[npbp_masses].mass = atof(param.mass_label[npbp_masses]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", 
			      &param.ksp_pbp[npbp_masses].naik_term_epsilon );
#else
	IF_OK param.ksp_pbp[npbp_masses].naik_term_epsilon = 0.0;
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

	IF_OK status += get_s(stdin, prompt, "save_file", param.pbp_filenames[npbp_masses] );
#endif

	/* The set to which this pbp_mass belongs */
	IF_OK param.set[npbp_masses] = k;

	/* maximum no. of conjugate gradient iterations */
	param.qic_pbp[npbp_masses].max = max_cg_iterations;
      
	/* maximum no. of conjugate gradient restarts */
	param.qic_pbp[npbp_masses].nrestart = max_cg_restarts;
      
	/* precision */
	param.qic_pbp[npbp_masses].prec = prec_pbp;

	/* errors */
	param.qic_pbp[npbp_masses].resid = error_for_propagator;
	param.qic_pbp[npbp_masses].relresid = rel_error_for_propagator;

	param.qic_pbp[npbp_masses].parity = EVENANDODD;
	param.qic_pbp[npbp_masses].min = 0;
	param.qic_pbp[npbp_masses].start_flag = 0;
	param.qic_pbp[npbp_masses].nsrc = 1;

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
      }
    }
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    

  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0;
  if( param.stopflag != 0 )return param.stopflag;

  if(prompt==2)return 0;

  /* Construct the eps_naik table of unique Naik epsilon
     coefficients.  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* First term is always zero */
  start_eps_naik(eps_naik, &n_naiks);
  
  /* Contribution from the chiral condensate epsilons */
  for(k = 0; k < param.num_set; k++)
    for(i = param.begin_pbp_masses[k]; i < param.num_pbp_masses[k]; i++){
      param.ksp_pbp[i].naik_term_epsilon_index = 
	fill_eps_naik(eps_naik, 
		      &n_naiks, param.ksp_pbp[i].naik_term_epsilon);
    }
  
  /* Do whatever is needed to get lattice */
  if( param.startflag == CONTINUE ){
    rephase( OFF );
  }
  if( param.startflag != CONTINUE ){
    startlat_p = reload_lattice( param.startflag, param.startfile );
  }

  /* if a lattice was read in, put in KS phases and AP boundary condition */
  phases_in = OFF;
  rephase( ON );

  /* Set options for fermion links */
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif

#ifdef DM_DU0
  fermion_links_want_du0(1);
#endif
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ ) &  defined(DM_DEPS)
  fermion_links_want_deps(1);
#endif
  
  fn_links = create_fermion_links_from_site(PRECISION, n_naiks, eps_naik);

  /* Construct APE smeared links, but without KS phases */
  rephase( OFF );
  ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
  rephase( ON );

/* We put in antiperiodic bc to match what we did to the gauge field */
  apply_apbc( ape_links );

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
  for(i = 0; i < Nvecs_tot; i++)
    eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /* Do whatever is needed to get eigenpairs */
  status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
			   &Nvecs_tot, eigVal, eigVec, 1);
  if(status != 0) terminate(1);

  if(param.ks_eigen_startflag != FRESH){
    param.eigcgp.Nvecs = 0;
    param.eigcgp.Nvecs_curr = Nvecs_tot;
    param.eigcgp.H = (double_complex *)malloc(Nvecs_max*Nvecs_max
					      *sizeof(double_complex));
    for(i = 0; i < Nvecs_max; i++){
      for(k = 0; k < i; k++)
	param.eigcgp.H[k + Nvecs_max*i] = dcmplx((double)0.0, (double)0.0);
      param.eigcgp.H[(Nvecs_max+1)*i] = dcmplx(eigVal[i], (double)0.0);
    }
  }
#endif

#if EIGMODE == DEFLATION
  /* malloc for eigenpairs */
  eigVal = (double *)malloc(param.Nvecs*sizeof(double));
  eigVec = (su3_vector **)malloc(param.Nvecs*sizeof(su3_vector *));
  for(i=0; i < param.Nvecs; i++)
    eigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /* Do whatever is needed to get eigenpairs */
  node0_printf("Reading %d eigenvectors\n", param.Nvecs); fflush(stdout);
  status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
			   &param.Nvecs, eigVal, eigVec, 1);
  if(status != 0)terminate(1);
#endif

  ENDTIME("readin");

  return(0);
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
