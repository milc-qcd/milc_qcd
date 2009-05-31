/******** setup.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)

#include "schroed_pg_includes.h"
int initial_set();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int prompt;

        /* print banner, get volume, nflavors, seed */
    prompt=initial_set();
        /* initialize the node random number generator */
    initialize_prn(&node_prn,iseed,volume+mynode());
        /* Initialize the layout functions, which decide where sites live */
    setup_layout();
        /* allocate space for lattice, set up coordinate fields */
    make_lattice();
        /* set up neighbor pointers and comlink structures */
    make_nn_gathers();

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
        /* print banner */
        printf("Schroedinger functional for pure gauge SU3\n");
#ifdef RMD_ALGORITHM
        printf("Microcanonical simulation with refreshing\n");
#endif
#ifdef HMC_ALGORITHM
        printf("Microcanonical simulation with refreshing\n");
#endif
        printf("MIMD version 6\n");
        printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
        printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef ORA_ALGORITHM
        printf("Overrelaxed/quasi-heat bath algorithm\n");
#endif
        status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
	IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );

	if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(mynode()==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    nx=par_buf.nx;
    ny=par_buf.ny;
    nz=par_buf.nz;
    nt=par_buf.nt;
    iseed=par_buf.iseed;
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=nx*ny*nz*nt;
    return(prompt);
}

/* read in parameters and coupling constants    */
int readin(int prompt) {
/* read in parameters for su3 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input */

int status;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

        printf("\n\n");
	status=0;
    
	/* warms, trajecs */
	IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
	IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );
    
	/* trajectories between propagator measurements */
	IF_OK status += 
	    get_i(stdin, prompt,"traj_between_meas", &par_buf.propinterval );
    
	/* get couplings and broadcast to nodes	*/
	/* beta */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

	/* boundary condition flag */
	IF_OK status += get_i(stdin, prompt,"bc_flag", &par_buf.bc_flag );

#if ( defined HMC_ALGORITHM || defined RMD_ALGORITHM )
        /* microcanonical time step */
	IF_OK status +=
            get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );
#endif  
        /*microcanonical steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
    
#ifdef ORA_ALGORITHM
        /*qhb steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"qhb_steps", &par_buf.stepsQ );
#endif   

        /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

        /* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );
	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				      par_buf.stringLFN );
 
        /* send parameter structure */
	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    warms = par_buf.warms;
    trajecs = par_buf.trajecs;
    steps = par_buf.steps;
    stepsQ = par_buf.stepsQ;
    propinterval = par_buf.propinterval;
    startflag = par_buf.startflag;
    saveflag = par_buf.saveflag;
    epsilon = par_buf.epsilon;
    beta = par_buf.beta;
    bc_flag = par_buf.bc_flag;
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);

    c_t11 = 0;

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );

    return(0);
} /*readin()*/
