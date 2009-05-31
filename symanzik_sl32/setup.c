/******** setup.c *********/
/* MIMD version 7 */
/* Version for Symanzik improved type pure gauge actions */
#define IF_OK if(status==0)

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#include "symanzik_sl32_includes.h"
#include <string.h>

void make_sublattices();

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
    make_lattice();     /* Makes standard lattice */
    make_sublattices(); /* Redefines "parity" and sets up "neighsub" */
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
	printf("Pure gauge SU3\n");
#ifdef HMC_ALGORITHM
        printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef ORA_ALGORITHM
        printf("Overrelaxed/quasi-heat bath algorithm\n");
	/**	if(NREPS > 1){
	    printf("Dummy! Use overrelaxed/quasi-heat bath only for NREPS = 1\n");
	    printf("But this executable has NREPS = %d\n", NREPS);
	    terminate(-1);
	    }**/
	/* Test is done by assertion in dsdu_qhb_subl */
#endif
	printf("MIMD version 7\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

	status = get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin,  prompt, "nx", &par_buf.nx );
	IF_OK status += get_i(stdin,  prompt, "ny", &par_buf.ny );
	IF_OK status += get_i(stdin,  prompt, "nz", &par_buf.nz );
	IF_OK status += get_i(stdin,  prompt, "nt", &par_buf.nt );
	IF_OK status += get_i(stdin,  prompt, "iseed", &par_buf.iseed );
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

void make_sublattices(){
  register int i,j;		/* scratch */
  int x,y,z,t;		/* coordinates */
  int x2,y2,z2,t2,subl,cb;
  
  /* Make list of neighbouring sublattices */
  for(subl=0;subl<N_SUBL32;subl++){
    x2 = subl%2;
    j = subl/2;
    y2 = j%2;
    j = j/2;
    z2 = j%2;
    j = j/2;
    t2 = j%2;
    cb = j/2;
    if( x2 == 0 ){
      neighsubl[subl][XUP] = 1 + 2*y2 + 4*z2 + 8*t2 + 16*cb;
      neighsubl[subl][XDOWN] = 1 + 2*y2 + 4*z2 + 8*t2 + 16*((cb+1)%2);
    }
    else{
      neighsubl[subl][XUP] = 2*y2 + 4*z2 + 8*t2 + 16*((cb+1)%2);
      neighsubl[subl][XDOWN] = 2*y2 + 4*z2 + 8*t2 + 16*cb;
    }
    if( y2 == 0 ){
      neighsubl[subl][YUP] = x2 + 2 + 4*z2 + 8*t2 + 16*cb;
      neighsubl[subl][YDOWN] = x2 + 2 + 4*z2 + 8*t2 + 16*((cb+1)%2);
    }
    else{
      neighsubl[subl][YUP] = x2 + 4*z2 + 8*t2 + 16*((cb+1)%2);
      neighsubl[subl][YDOWN] = x2 + 4*z2 + 8*t2 + 16*cb;
    }
    if( z2 == 0 ){
      neighsubl[subl][ZUP] = x2 + 2*y2 + 4 + 8*t2 + 16*cb;
      neighsubl[subl][ZDOWN] = x2 + 2*y2 + 4 + 8*t2 + 16*((cb+1)%2);
    }
    else{
      neighsubl[subl][ZUP] = x2 + 2*y2 + 8*t2 + 16*((cb+1)%2);
      neighsubl[subl][ZDOWN] = x2 + 2*y2 + 8*t2 + 16*cb;
    }
    if( t2 == 0 ){
      neighsubl[subl][TUP] = x2 + 2*y2 + 4*z2 + 8 + 16*cb;
      neighsubl[subl][TDOWN] = x2 + 2*y2 + 4*z2 + 8 + 16*((cb+1)%2);
    }
    else{
      neighsubl[subl][TUP] = x2 + 2*y2 + 4*z2 + 16*((cb+1)%2);
      neighsubl[subl][TDOWN] = x2 + 2*y2 + 4*z2 + 16*cb;
    }
  }
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    if(node_number(x,y,z,t)==mynode()){
      i=node_index(x,y,z,t);
      x2 = x/2; y2 = y/2; z2 = z/2; t2 = t/2;
      subl = (x%2) + 2*(y%2) + 4*(z%2) + 8*(t%2);
      if( (x2+y2+z2+t2)%2 == 0)lattice[i].parity=subl;
      else		     lattice[i].parity=subl+16;
    }
  }
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
/* read in parameters for su3 monte carlo	*/
/* argument "prompt" is 1 if prompts are to be given for input	*/

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

	/* no dynamical masses for pure gauge */
	n_dyn_masses = 0;
	dyn_flavors[0] = 0;

	/* u0 */
	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

#ifdef HMC_ALGORITHM
	/* steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
	/* microcanonical time step */
	IF_OK status +=
	    get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );
#endif  

#ifdef ORA_ALGORITHM
	/*overrelaxed steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
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
    u0 = par_buf.u0;
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );

    return(0);
}
