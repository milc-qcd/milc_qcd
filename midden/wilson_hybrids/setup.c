/******** setup.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)

#include "wi_hyb_includes.h"

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int initial_set();
int prompt;

	/* print banner, get volume, seed */
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
    /* On node zero, read lattice size, seed, and send to others */
    if(mynode()==0){
	/* print banner */
	printf("SU3 spectrum with Wilson fermions\n");
	printf("Microcanonical simulation with refreshing\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");

	status=get_prompt(&prompt);

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
    total_iters=0;
    return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
/* read in parameters for su3 monte carlo	*/
/* argument "prompt" is 1 if prompts are to be given for input	*/

int status;
Real x;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;
    
	/* get couplings and broadcast to nodes	*/
	/* beta, kappa */
	IF_OK status += get_f(prompt,"beta", &par_buf.beta );
	IF_OK status += get_f(prompt,"kappa", &par_buf.kappa );

	/* source time slice and increment */
	IF_OK status += get_i(prompt,"source_start",&par_buf.source_start);
	IF_OK status += get_i(prompt,"source_inc",&par_buf.source_inc);
	IF_OK status += get_i(prompt,"n_sources",&par_buf.n_sources);
    
	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(prompt,"max_cg_iterations", &par_buf.niter );
    
	/* error for propagator conjugate gradient */
	IF_OK status += get_f(prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;
    
        /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice( prompt, &(par_buf.startflag),
	    par_buf.startfile );

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    startflag = par_buf.startflag;
    niter = par_buf.niter;
    rsqprop = par_buf.rsqprop;
    beta = par_buf.beta;
    kappa = par_buf.kappa;
    source_start = par_buf.source_start;
    source_inc = par_buf.source_inc;
    n_sources = par_buf.n_sources;
    strcpy(startfile,par_buf.startfile);

    /* Do whatever is needed to get lattice */
    startlat_p = reload_lattice( startflag, startfile );

    return(0);
}
