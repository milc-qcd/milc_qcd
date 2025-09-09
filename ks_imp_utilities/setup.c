/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik improved
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//
*/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#include <lattice_qdp.h>

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

/* Each node has a params structure for passing simulation parameters */
#include "params.h"

/* Forward declarations */

int initial_set();
void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
void make_3n_gathers();

int
setup()
{
  int prompt;

  /* print banner, get volume, nflavors1,nflavors2, seed */
  prompt = initial_set();
  if(prompt == 2)return prompt;
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);
  /* Mark t_longlink and t_fatlink unallocted */
  // init_ferm_links(&fn_links, &ks_act_paths);
  /* set up neighbor pointers and comlink structures
     code for this routine is in com_machine.c  */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  node0_printf("Made 3nn gathers\n"); fflush(stdout);
  /* set up K-S phase vectors, boundary conditions */
  phaseset();

  node0_printf("Finished setup\n"); fflush(stdout);
  return  prompt;
}

static int n_naiks = 1;
static double eps_naik[MAX_NAIK];

/* SETUP ROUTINES */
int
initial_set()
{
  int prompt=0,status;
  /* On node zero, read lattice size, seed and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
#if defined(CHECK_INVERT)
    printf("Inversion checking\n");
#elif defined(FERMION_FORCE)
    printf("Fermion-force checking\n");
#elif defined(LINK_FATTENING)
    printf("Creating FN link files\n");
#elif defined(CHECK_LINK_FATTENING)
    printf("Checking FN link fattening files\n");
#elif defined(REUNIT)
    printf("Reunitarization checking\n");
#else
#error "Must specify what is being checked"
#endif

    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

    status=get_prompt(stdin, &prompt);
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
    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    return param.stopflag;

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
#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = 0;
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = 0;
#endif

  return prompt;
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int i, status, current_index;
  char savebuf[128];
#ifdef CHECK_INVERT
  char invert_string[16];
#endif

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {

    printf("\n\n");
    status=0;

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(param.startflag),
					  param.startfile );

    IF_OK status += get_f(stdin, prompt,"u0", &param.u0 );

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
					param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				  param.stringLFN );

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
    
    /* find out what to do with longlinks at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.savelongflag),
				       param.savelongfile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.savelongflag,
				  param.stringLFNlong );
    /* find out what to do with fatlinks at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.savefatflag),
				       param.savefatfile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.savefatflag,
				  param.stringLFNfat );
    IF_OK status += get_i(stdin, prompt,"withKSphases", &param.withKSphases );

    /* Eigenpairs not supported */
    param.eigen_param.Nvecs = 0;

#if defined(CHECK_INVERT) || defined(FERMION_FORCE)    
    /* Inversion parameters */
    IF_OK status += get_i(stdin, prompt,"number_of_masses", &param.nmass );

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &param.qic[0].max );

    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &param.qic[0].nrestart );
#endif    
#ifdef CHECK_INVERT
    /* find out what kind of color vector source to use */
    IF_OK status += ask_color_vector( prompt, &(param.srcflag[0]),
				      param.srcfile[0] );
#endif
#if defined(CHECK_INVERT) || defined(FERMION_FORCE)    
    IF_OK for(i = 0; i < param.nmass; i++){
    
#ifdef FERMION_FORCE
      IF_OK status += ask_color_vector( prompt, &(param.srcflag[i]),
					param.srcfile[i] );
      IF_OK if(param.srcflag[i] != param.srcflag[0]){
	node0_printf("Must reload all or save all alike.");
	status++;
      }
#endif
      IF_OK status += get_f(stdin, prompt,"mass", &param.ksp[i].mass );
#if FERM_ACTION == HISQ
      IF_OK status += get_f(stdin, prompt, "naik_term_epsilon", 
			    &param.ksp[i].naik_term_epsilon);
      IF_OK {
	if(i == 0){
	  if(param.ksp[i].naik_term_epsilon != 0.0){
	    node0_printf("First Naik term epsilon must be zero.");
	    status++;
	  }
	} else {
	  if(param.ksp[i].naik_term_epsilon > param.ksp[i-1].naik_term_epsilon){
	    node0_printf("Naik term epsilons must be in descending order.");  // Why?
	    status++;
	  }
	}
      }

#else
      param.ksp[i].naik_term_epsilon = 0.0;
#endif

      param.qic[i].min = 0;
      param.qic[i].start_flag = 0;
      param.qic[i].nsrc = 1;
      param.qic[i].max = param.qic[0].max;
      param.qic[i].nrestart = param.qic[0].nrestart;
      param.qic[i].prec = MILC_PRECISION;
      param.qic[i].parity = EVENANDODD;
      /* error for propagator conjugate gradient */
      IF_OK status += get_f(stdin, prompt, "error_for_propagator", 
			    &param.qic[i].resid);
      IF_OK status += get_f(stdin, prompt, "rel_error_for_propagator", 
			    &param.qic[i].relresid );
#ifdef CHECK_INVERT
      /* find out what kind of color vector result to use */
      IF_OK status += ask_color_vector( prompt, &(param.ansflag[i]),
					param.ansfile[i] );
      IF_OK if(param.ansflag[i] != param.ansflag[0]){
	node0_printf("Must reload all or save all alike.");
	status++;
      }
#endif
    }
#endif // CHECK_INVERT or FERMION_FORCE
#if defined(FERMION_FORCE) || defined(CHECK_FATTENING)
    /* Optional answer for fat links or fermion force */
    IF_OK status += ask_color_matrix( prompt, &(param.ansflag[0]),
				      param.ansfile[0] );
#ifdef CHECK_FATTENING
    /* Optional answer for fat links or long links */
    IF_OK status += ask_color_matrix( prompt, &(param.ansflag[1]),
				      param.ansfile[1] );
#endif
#endif


#ifdef CHECK_INVERT
    /* find out which inversion to check */
    IF_OK status += get_s(stdin,  prompt, "invert", invert_string);
    if(status == 0){
      if(strcmp(invert_string,"M")==0)
	param.inverttype = INVERT_M;
      else if(strcmp(invert_string,"MdaggerM")==0)
	param.inverttype = INVERT_MdaggerM;
      else{
	printf("Unrecognized invert string %s\n",invert_string);
	status++;
      }
    }
#endif

    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 ){
    return param.stopflag;
  }

  if(prompt==2)return 0;

  nmass = param.nmass;
  niter = param.qic[0].max;
  nrestart = param.qic[0].nrestart;
  u0 = param.u0;
  startflag = param.startflag;
  saveflag = param.saveflag;
  savelongflag = param.savelongflag;
  savefatflag = param.savefatflag;
  ansflag = param.ansflag[0];  /* Must all be the same */
  srcflag = param.srcflag[0];
  strcpy(startfile,param.startfile);
  strcpy(savefile,param.savefile);
  strcpy(stringLFN, param.stringLFN);
  strcpy(savelongfile,param.savelongfile);
  strcpy(stringLFNlong, param.stringLFNlong);
  strcpy(savefatfile,param.savefatfile);
  strcpy(stringLFNfat, param.stringLFNfat);
  inverttype = param.inverttype;

  /* Construct the eps_naik table of unique Naik epsilon
     coefficients.  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* We require the first naik_term_epsilon to be zero 
     and the remaining naik_term_epsilons to be in descending order */
  start_eps_naik(eps_naik, &n_naiks);
  current_index = 0;
  n_orders_naik[current_index] = 0;
  n_order_naik_total = nmass;
  
  for(i = 0; i < param.nmass; i++){
    int next_index = 
      fill_eps_naik(eps_naik, &n_naiks, param.ksp[i].naik_term_epsilon);
    param.ksp[i].naik_term_epsilon_index = next_index;
    if(next_index == current_index)
      n_orders_naik[current_index]++;
    else {
      current_index++;
      n_orders_naik[current_index] = 1;
    }
  }

  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile );

  /* if a lattice was read in, put in KS phases and AP boundary condition */
  phases_in = OFF;
  rephase( ON );

#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif

#if FERM_ACTION == HISQ
  fn_links = create_fermion_links_from_site(MILC_PRECISION, n_naiks, eps_naik);
#else
  fn_links = create_fermion_links_from_site(MILC_PRECISION, 0, NULL);
#endif

  return 0;
}

/* Set up comlink structures for 3rd nearest gather pattern;
   make_lattice() and  make_nn_gathers() must be called first,
   preferably just before calling make_3n_gathers().
*/
void
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
