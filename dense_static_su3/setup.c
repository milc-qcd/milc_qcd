/******** setup.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)

#include "su3_dense_includes.h"

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
   int initial_set();
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
#ifdef LIGHT_PBP
    phaseset();
#endif

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
        /* print banner */
        printf("Dense \"almost quenched\" SU(3)\n");
        printf("Overrelaxed/Metropolis algorithm\n");
	time_stamp("start");
        status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(prompt,"nt", &par_buf.nt );
	IF_OK status += get_i(prompt,"iseed", &par_buf.iseed );

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
int readin( int prompt){
/* read in parameters for su3 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input  */

int status,status2;
Real x;
char savebuf[128];

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

        printf("\n\n");
	status=0;
    
	/* warms, trajecs */
	IF_OK status += get_i(prompt,"warms", &par_buf.warms );
	IF_OK status += get_i(prompt,"trajecs", &par_buf.trajecs );
    
	/* trajectories between propagator measurements */
	IF_OK status += 
	    get_i(prompt,"traj_between_meas", &par_buf.propinterval );
    
	/* get couplings and broadcast to nodes	*/
	/* beta, mass */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

        /* get and broadcast C */
	IF_OK status += get_f(stdin, prompt,"C", &par_buf.C );

        /* overrelaxed steps per trajectory */
	IF_OK status += get_i(prompt,"overrelax_steps", &par_buf.steps_over );
    
        /*metropolis steps per time link */
	IF_OK status += get_i(prompt,"update_steps", &par_buf.steps_update );
    
	/* find out what kind of starting lattice to use */
	IF_OK ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );
 
#ifdef ORA_ALGORITHM
       IF_OK if (prompt==1) printf(
           "enter '[n]o_gauge_fix', or 'coulomb_gauge_fix'\n");
       status2=scanf("%s",savebuf);
       IF_OK if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
          par_buf.fixflag = COULOMB_GAUGE_FIX;
          if(this_node==0)printf("fixing to coulomb gauge\n");
       }
       else if((strcmp("no_gauge_fix",savebuf) == 0) || 
	       (strcmp("n",savebuf) == 0)) {
          par_buf.fixflag = NO_GAUGE_FIX;
          if(this_node==0)printf("NOT fixing the gauge\n");
       }
       else{
           printf("error in input: fixing_command is invalid\n"); status++;
       }
#endif
 
	/* find out what to do with lattice at end */
	IF_OK ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );
	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				      par_buf.stringLFN );

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    warms = par_buf.warms;
    trajecs = par_buf.trajecs;
    steps_over = par_buf.steps_over;
    steps_update = par_buf.steps_update;
    propinterval = par_buf.propinterval;
    startflag = par_buf.startflag;
    fixflag = par_buf.fixflag;
    saveflag = par_buf.saveflag;
    epsilon = par_buf.epsilon;
    beta = par_buf.beta;
    C = par_buf.C;
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );
    phases_in = OFF;
    return(0);
}





