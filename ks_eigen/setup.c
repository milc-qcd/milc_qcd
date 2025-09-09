/************************** setup.c ***************************/
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
//
*/

#define IF_OK if(status==0)

#include "ks_eig_includes.h"	/* definitions files and prototypes */
//#include "lattice_qdp.h"
#include <string.h>
#include "params.h"
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h
#ifdef U1_FIELD
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"
#endif

//#ifdef HAVE_QOP
//#include "../include/generic_qop.h"
//#endif

/* Forward declarations */

static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers(void);

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
  /* set up neighbor pointers and comlink structures
     code for this routine is in com_machine.c  */
  make_nn_gathers();
#ifdef FN
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
#endif
  /* set up K-S phase vectors, antiperiodic boundary conditions */
  phaseset();

  return( prompt );
}


/* SETUP ROUTINES */
static int initial_set(){
  int prompt=0,status;
  /* On node zero, read lattice size, seed, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Eigenvalues and eigenvectors\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes(ranks)\n",machine_type(),numnodes());
    
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

    status = get_prompt(stdin, &prompt);
    
    IF_OK status += get_i(stdin, prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &param.nt );
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
  /* read in parameters for the calculation	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  char savebuf[128];
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

#ifdef U1_FIELD
    /* what kind of starting U(1) lattice to use, read filename */
    IF_OK status+=ask_starting_u1_lattice(stdin,prompt,
					  &param.start_u1flag, param.start_u1file );
    IF_OK status+=ask_ending_u1_lattice(stdin,prompt,
					&param.save_u1flag, param.save_u1file );
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
    
#ifdef U1_FIELD
      /* Charge, if U(1) */
      IF_OK status += get_s(stdin, prompt, "charge", param.charge_label);
      IF_OK param.charge = atof(param.charge_label);
#endif

    /*------------------------------------------------------------*/
    /* Dirac eigenpair calculation                                */
    /*------------------------------------------------------------*/

    /* number of eigenpairs */
    IF_OK status += get_i(stdin, prompt,"max_number_of_eigenpairs", &param.eigen_param.Nvecs);

    /* eigenvector input */
    IF_OK status += ask_starting_ks_eigen(stdin, prompt, &param.ks_eigen_startflag,
					  param.ks_eigen_startfile);

    IF_OK {
      if(param.ks_eigen_startflag == FRESH) param.eigen_param.Nvecs_in = 0;
      /* Start by assuming we have the max number as input -- changed later, if not. */
      else param.eigen_param.Nvecs_in = param.eigen_param.Nvecs;
    }

    /* eigenvector output */
    IF_OK status += ask_ending_ks_eigen(stdin, prompt, &param.ks_eigen_saveflag,
					param.ks_eigen_savefile);
    /*------------------------------------------------------------*/
    /* Dirac eigenpair parameters                                 */
    /*------------------------------------------------------------*/
    
    status += read_ks_eigen_param(&param.eigen_param, status, prompt);
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  u0 = param.u0;
  if( param.stopflag != 0 )return param.stopflag;
  
  if(prompt==2)return 0;
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  n_naiks = 1;
  eps_naik[0] = 0.0;
#endif
  
  /* Do whatever is needed to get lattice */
  if( param.startflag == CONTINUE ){
    rephase( OFF );
  }
  if( param.startflag != CONTINUE )
    startlat_p = reload_lattice( param.startflag, param.startfile );
  /* if a lattice was read in, put in KS phases and AP boundary condition */
  phases_in = OFF;
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
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fn_links = create_fermion_links_from_site(MILC_PRECISION, n_naiks, eps_naik);
#else
  fn_links = create_fermion_links_from_site(MILC_PRECISION, 0, NULL);
#endif
  
  eigVal = (double *)malloc(param.eigen_param.Nvecs*sizeof(double));
  eigVec = (su3_vector **)malloc(param.eigen_param.Nvecs*sizeof(su3_vector*));
  for(int i=0;i<param.eigen_param.Nvecs;i++)
    eigVec[i]= (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));
  
  /* Do whatever is needed to get eigenpairs */
  imp_ferm_links_t *fn = get_fm_links(fn_links, 0);
  status = reload_ks_eigen(param.ks_eigen_startflag, param.ks_eigen_startfile, 
			   &param.eigen_param.Nvecs_in, eigVal, eigVec, fn, 1);
#if 1
  /* Calculate and print the residues and norms of the eigenvectors */
  double *resid = (double *)malloc(param.eigen_param.Nvecs_in*sizeof(double));
  node0_printf("Even site residuals\n");
  check_eigres( resid, eigVec, eigVal, param.eigen_param.Nvecs_in, EVEN, fn );
  free(resid);
  destroy_fn_links(fn);
#endif
  
  if(param.fixflag != NO_GAUGE_FIX){
    node0_printf("WARNING: Gauge fixing does not readjust the eigenvectors");
  }
  if(status != 0) normal_exit(0);

  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
 */

static void 
make_3n_gathers(){
   int i;
 
   for(i=XUP;i<=TUP;i++) {
      make_gather(third_neighbor,&i,WANT_INVERSE,
		  ALLOW_EVEN_ODD,SWITCH_PARITY);
   }
   
    /* Sort into the order we want for nearest neighbor gathers,
       so you can use X3UP, X3DOWN, etc. as argument in calling them. */

   sort_eight_gathers(X3UP);
}
 

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

static void 
third_neighbor(int x,int y,int z,int t,int *dirpt,int FB,int *xp,int *yp,int *zp,int *tp)
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
