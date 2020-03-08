/************************ setup.c ****************************/
/* MIMD version 7 */

#define IF_OK if(status==0)

#include "gauge_replicate_includes.h"	/* definitions files and prototypes */

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

int
setup()
{
  int initial_set();
  int prompt;
  
  /* print banner, get initial parameters */
  prompt = initial_set();
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();

  node0_printf("Made lattice\n"); fflush(stdout);
  /* set up neighbor pointers and comlink structures
     code for this routine is in com_machine.c  */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);
  
  node0_printf("Finished setup\n"); fflush(stdout);
  return prompt;
}


/* SETUP ROUTINES */
int 
initial_set()
{
  int prompt=0,status;
  int i;

  /* On node zero, read lattice size, seed, nflavors1, nflavors2,
     nflavors, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 gauge utility program\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");

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
    param.iseed = 1234;  /* Not used */
    
    /* beta, quark masses */
    IF_OK status += get_f(stdin, prompt,"beta", &param.beta );

    IF_OK status += get_i(stdin, prompt,"n_dyn_masses", &param.n_dyn_masses );
    IF_OK status += get_vf(stdin, prompt, "dyn_mass", param.dyn_mass, param.n_dyn_masses);
    IF_OK status += get_vi(stdin, prompt, "dyn_flavors", param.dyn_flavors, param.n_dyn_masses);
    /* Get tadpole factor */
    IF_OK status += get_f(stdin, prompt, "u0", &param.u0);

    /* Get translation vector */

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
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif
  iseed=param.iseed;
  
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  beta = param.beta;
  
  n_dyn_masses = param.n_dyn_masses;
  for(i = 0; i < n_dyn_masses; i++){
    dyn_mass[i] = param.dyn_mass[i];
    dyn_flavors[i] = param.dyn_flavors[i];
  }
  u0 = param.u0;
  return prompt;
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {
    
    printf("\n\n");
    status=0;

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(param.startflag),
					 param.startfile );

    /* Get translation vector */
    IF_OK status += get_vi(stdin, prompt, "replicate", param.reps, 4);

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
				       param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				 param.stringLFN );
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));
  
  if( param.stopflag != 0 )return param.stopflag;
  if(prompt==2)return 0;
  
  /* Get lattice no KS phases in this program */
  if( param.startflag != CONTINUE )
    startlat_p = reload_lattice( param.startflag, param.startfile );
  return 0;
}
