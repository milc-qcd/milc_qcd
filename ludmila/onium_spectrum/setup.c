/******** setup.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)


#include "ks_dyn_includes.h"
#include <string.h>

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int initial_set();
void make_gen_pt(),setup_layout();
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
	/* set up K-S phase vectors, boundary conditions */
    phaseset();
    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
	/* print banner */
	printf("Heavy-light spectroscpy with Kogut-Susskind light fermions\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");
	status = get_prompt( &prompt );
	IF_OK status += get_i(prompt,"nflavors", &par_buf.nflavors );
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
    nflavors=par_buf.nflavors;
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=nx*ny*nz*nt;
    total_iters=0;
    return(prompt);
}

/* read in parameters and coupling constants	*/

int readin(int prompt)  {

/* argument "prompt" is 1 if prompts are to be given for input	*/

int status,i;
float x;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;

	/* Number of kappas */
	IF_OK status += get_i(prompt,"number_of_kappas", &par_buf.num_kap );
	if( par_buf.num_kap>MAX_KAP ){
	  printf("num_kap = %d must be <= %d!\n", par_buf.num_kap, MAX_KAP);
	  status++;
	}
    
	/* To be safe initialize the following to zero */
	for(i=0;i<MAX_KAP;i++){
	  kap[i] = 0.0;
	}
    
	for(i=0;i<par_buf.num_kap;i++){
	  IF_OK status += get_f(prompt,"kappa", &par_buf.kap[i] );
	}
	IF_OK status += get_f(prompt,"rot_param", &par_buf.d1 );
        
	/* get couplings and broadcast to nodes	*/
	/* beta, mass */
	IF_OK status += get_f(prompt,"beta", &par_buf.beta );
	IF_OK status += get_f(prompt,"mass", &par_buf.mass );
    
	    
	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(prompt,"max_cg_iterations", &par_buf.niter );
    
	/* error per site for conjugate gradient */
	IF_OK status += get_f(prompt,"error_per_site", &x );
	IF_OK par_buf.rsqmin = x*x;   /* rsqmin is r**2 in conjugate gradient */
	    /* New conjugate gradient normalizes rsqmin by norm of source */
    
	/* error for propagator conjugate gradient */
	IF_OK status += get_f(prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;
    
	/* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice( prompt, &(par_buf.startflag),
	    par_buf.startfile );
    	
	/* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice( prompt, &(par_buf.saveflag),
	    par_buf.savefile );
    	
	/* Get ensemble values for NERSC archive */
	IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
		  par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
	  status += get_s( prompt,"ensemble_id", par_buf.ensemble_id );
	IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
		  par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
	  status += get_i( prompt,"sequence_number", 
			   &par_buf.sequence_number );
	IF_OK status += ask_starting_prop( prompt,&par_buf.startflag_w[0],
					   par_buf.startfile_w[0]);
	/* what to do with computed wilson propagator */
	IF_OK status += ask_ending_prop( prompt,&par_buf.saveflag_w[0],
				 par_buf.savefile_w[0]);
	IF_OK status += get_s(prompt,"smear_func_file", par_buf.smearfile[0]);
	IF_OK status += get_s(prompt,"start_ks_prop_file", par_buf.start_ks_prop_file);
	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

   
    startflag = par_buf.startflag;
    saveflag = par_buf.saveflag;
    niter = par_buf.niter;
    rsqmin = par_buf.rsqmin;
    rsqprop = par_buf.rsqprop;
    beta = par_buf.beta;
    mass = par_buf.mass;
    d1 = par_buf.d1;
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(startfile_w[0],par_buf.startfile_w[0]);
    strcpy(savefile_w[0],par_buf.savefile_w[0]);
    strcpy(ensemble_id,par_buf.ensemble_id);
    sequence_number = par_buf.sequence_number;

    /* Do whatever is needed to get lattice */
    if( startflag == CONTINUE ){
        rephase( OFF );
    }
    startlat_p = reload_lattice( startflag, startfile );
    /* if a lattice was read in, put in KS phases and AP boundary condition */
    phases_in = OFF;
    rephase( ON );

    return(0);
}
