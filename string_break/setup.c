/******** setup.c *********/
/* MIMD version 6 */
/* original code by UMH */

#define IF_OK if(status==0)

#include "string_break_includes.h"
#include <string.h>
int initial_set();
void phaseset();
void rephase( int flag );


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
        printf("Pure gauge SU3 Wilson loops and static light mesons\n");
        printf("MIMD version 6\n");
        printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

        status=get_prompt(stdin, &prompt);
        IF_OK status += get_i(stdin, prompt,"nx",&par_buf.nx);
        IF_OK status += get_i(stdin, prompt,"ny",&par_buf.ny);
        IF_OK status += get_i(stdin, prompt,"nz",&par_buf.nz);
        IF_OK status += get_i(stdin, prompt,"nt",&par_buf.nt);
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


/********* phaseset() - set up phase vectors **********/
/** PERIODIC BOUNDARY CONDITIONS IN ALL DIRS **/
void phaseset() {
register site *sit; /* pointer to current site */
register int i;
    /*  phase of link(i,mu) = sum(nu<mu) { -1^i[nu] }           */
    /*  We use peridic boundary conditions in all directions!   */

    FORALLSITES(i,sit){
	sit->phase[TUP] = 1.0;
	if( (sit->t)%2 == 1) sit->phase[XUP] = -1.0;
	  else sit->phase[XUP] = 1.0;
	if( (sit->x)%2 == 1) sit->phase[YUP] = -sit->phase[XUP];
	  else sit->phase[YUP] = sit->phase[XUP];
	if( (sit->y)%2 == 1) sit->phase[ZUP] = -sit->phase[YUP];
	  else sit->phase[ZUP] = sit->phase[YUP];
    }
}

/************************** rephase() ******************************/
/* put Kogut-Sussind and boundary condition phase factors into or
   out of lattice */
void rephase( int flag ){
register int i,j,k,dir;
register site *s;
    /* Check to make sure we are going in expected direction */
    if( !( (flag==ON && phases_in==OFF) || (flag==OFF && phases_in==ON) ) ){
	node0_printf("DUMMY: you fouled up the phases\n");
        terminate(1);
    }
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		s->link[dir].e[j][k].real *= s->phase[dir];
		s->link[dir].e[j][k].imag *= s->phase[dir];
	    }
	}
    }
    phases_in = flag;
}

/* read in parameters and coupling constants    */
int readin(int prompt){
/* read in parameters for su3 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input  */

int status,i;
Real x;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

        printf("\n\n");
	status = 0;
    
	/* get couplings and broadcast to nodes */
	/* smearing iterations and factor */
	IF_OK status += get_i(stdin, prompt,"no_smear_level", 
			      &par_buf.no_smear_level);
	if(par_buf.no_smear_level > MAX_LEVEL){
	    printf("no_smear_level = %d must be <= %d!\n",
		par_buf.no_smear_level, MAX_LEVEL);
	    status++;
	}

	/* To be safe initialize the following to zero */
	for(i=0;i<MAX_LEVEL;i++){
	    smear_num[i] = 0;
	}

	for(i=0;i<par_buf.no_smear_level;i++){
	    IF_OK status += get_i(stdin, prompt,"smear_num", &par_buf.smear_num[i]);
	}

	IF_OK status += get_f(stdin, prompt,"smear_fac", &par_buf.smear_fac);

	/* off_axis_flag : do off-axis Wilson loops? */
	IF_OK status += get_i(stdin, prompt,"off_axis_flag", &par_buf.off_axis_flag);
    
	/* mass */
	IF_OK status += get_f(stdin, prompt,"mass", &par_buf.mass );

	/* r0, for "fuzzy" light quark source/sink */
	IF_OK status += get_i(stdin, prompt,"r0", &par_buf.r0);

	/* number of Gaussian random sources */
	IF_OK status += get_i(stdin, prompt,"num_src", &par_buf.num_src);
	if(par_buf.num_src > MAX_SRC){
	    printf("num_src = %d must be <= %d!\n", par_buf.num_src, MAX_SRC);
	    status++;
	}

	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
    
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );

	/* error per site for conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_per_site", &x );
	IF_OK par_buf.rsqmin = x*x;   /* rsqmin is r**2 in conjugate gradient */
	    /* New conjugate gradient normalizes rsqmin by norm of source */

	/* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

	/* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );
	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				      par_buf.stringLFN );

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);


	no_smear_level = par_buf.no_smear_level;
	for(i=0;i<par_buf.no_smear_level;i++){
	    smear_num[i] = par_buf.smear_num[i];
	}
	smear_fac = par_buf.smear_fac;
	startflag = par_buf.startflag;
	off_axis_flag = par_buf.off_axis_flag;
	mass = par_buf.mass;
	r0 = par_buf.r0;
	num_src = par_buf.num_src;
	niter = par_buf.niter;
	nrestart = par_buf.nrestart;
	rsqmin = par_buf.rsqmin;
	saveflag = par_buf.saveflag;
	strcpy(startfile,par_buf.startfile);
	strcpy(savefile,par_buf.savefile);
	strcpy(stringLFN, par_buf.stringLFN);

    /* Do whatever is needed to get lattice */
	if( startflag != CONTINUE )
	  startlat_p = reload_lattice( startflag, startfile );

    return(0);
}
