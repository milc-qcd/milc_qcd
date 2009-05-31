/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "schroed_ks_includes.h"
#include <string.h>

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
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up K-S phase vectors, boundary conditions */
  phaseset_sf();

  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
	/* print banner */
	printf("SU3 Schroedinger functional with Kogut-Susskind fermions\n");
	printf("Microcanonical simulation with refreshing\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
	printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
	printf("PHI algorithm\n");
#else
	printf("R algorithm\n");
#endif
#ifdef SEXT_WEIN
	printf("with Sexton-Weingarten updating\n");
#endif
	status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin, prompt,"nflavors", &par_buf.nflavors );
#ifdef PHI_ALGORITHM
	IF_OK if(par_buf.nflavors != 4){
	    printf("Dummy! Use phi algorithm only for four flavors\n");
	    status++;
	}
#endif
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
    nflavors=par_buf.nflavors;

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
#ifdef FERM_PHASES
 int i;
#endif

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;

	/* warms, trajecs */
	IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
	IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );

	/* get couplings and broadcast to nodes	*/
	/* beta, mass */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );
	IF_OK status += get_f(stdin, prompt,"mass", &par_buf.mass );

#ifdef REWEIGH
	IF_OK status += get_f(stdin, prompt,"gamma", &par_buf.gamma_rv );
#endif

	/* boundary condition flag */
	IF_OK status += get_i(stdin, prompt,"bc_flag", &par_buf.bc_flag );

#ifdef FERM_PHASES
	/* fermion phase factors */
	IF_OK status += get_f(stdin, prompt,"ferm_phases[0]", &par_buf.ferm_phas[0] );
	IF_OK status += get_f(stdin, prompt,"ferm_phases[1]", &par_buf.ferm_phas[1] );
	IF_OK status += get_f(stdin, prompt,"ferm_phases[2]", &par_buf.ferm_phas[2] );
#endif

	/* microcanonical time step */
	IF_OK status +=
	    get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );

#ifdef SEXT_WEIN
	/* n for Sexton-Weingarten update */
	IF_OK status += get_i(stdin, prompt,"sexton_weingarten_n", &par_buf.n_sxw );
#endif

	/*microcanonical steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );

	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
    
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );

	/* error per site for conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_per_site", &x );
	IF_OK par_buf.rsqmin = x*x;	/* rsqmin is r**2 in conjugate gradient */

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
    startflag = par_buf.startflag;
    saveflag = par_buf.saveflag;
    niter = par_buf.niter;
    nrestart = par_buf.nrestart;
    rsqmin = par_buf.rsqmin;
    epsilon = par_buf.epsilon;
#ifdef SEXT_WEIN
    n_sxw = par_buf.n_sxw;
#endif
    beta = par_buf.beta;
    mass = par_buf.mass;
#ifdef REWEIGH
    gamma_rv = par_buf.gamma_rv;
#endif
    bc_flag = par_buf.bc_flag;
#ifdef FERM_PHASES
    for(i=0;i<3;i++){
            ferm_phases[i] = par_buf.ferm_phas[i];
    }
#endif
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);

    c_t11 = -(Real)(nflavors)*0.02841;

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );
    /* put in KS phases and fermion phases, if desired */
    phases_in = OFF;
    rephase_sf( ON );

    return(0);
}

/********* phaseset_sf() - set up phase vectors **********/
void phaseset_sf() {
register site *sit; /* pointer to current site */
register int i;
    /*	phase of link(i,mu) = sum(nu<mu) { -1^i[nu] }		*/
    /* Without antiperiodic bc here */

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

/************************** rephase_sf() ******************************/
/* put Kogut-Sussind phase factors into or out of lattice,
   as well as the "fermonion phases", if desired */
void rephase_sf( int flag ){
  register int i,j,k,dir;
  register site *s;
#ifdef FERM_PHASES
  Real ddr, ddi, ethr, ethi;
#endif

    /* Check to make sure we are going in expected direction */
    if( !( (flag==ON && phases_in==OFF) || (flag==OFF && phases_in==ON) ) ){
        node0_printf("DUMMY: you fouled up the phases\n");
        terminate(0);
    }

#ifndef FERM_PHASES
    FORALLSITES(i,s){
	/* The phases for the time direction are all 1! */
	for(dir=XUP;dir<TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		s->link[dir].e[j][k].real *= s->phase[dir];
		s->link[dir].e[j][k].imag *= s->phase[dir];
	    }
	}
	/* The boundary links for t=nt are stored in boundary at t=nt-1 */
	if(s->t == (nt-1)) for(dir=XUP;dir<TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		s->boundary[dir].e[j][k].real *= -s->phase[dir];
		s->boundary[dir].e[j][k].imag *= -s->phase[dir];
	    }
	}
	/* The derivatives for t=nt are stored in boundary at t=nt-2 */
	if(s->t == (nt-2)) for(dir=XUP;dir<TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		s->boundary[dir].e[j][k].real *= s->phase[dir];
		s->boundary[dir].e[j][k].imag *= s->phase[dir];
	    }
	}
	/* The derivatives for t=0 are stored in boundary at t=0 */
	if(s->t == 0) for(dir=XUP;dir<TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		s->boundary[dir].e[j][k].real *= s->phase[dir];
		s->boundary[dir].e[j][k].imag *= s->phase[dir];
	    }
	}
    }
#else
    /* The phases for the time direction are all 1! */
    for(dir=XUP;dir<=ZUP;dir++){
	switch(dir){
		case(XUP):	ddr = ferm_phases[dir] / (Real)nx;	break;
		case(YUP):	ddr = ferm_phases[dir] / (Real)ny;	break;
		case(ZUP):	ddr = ferm_phases[dir] / (Real)nz;	break;
	}
	ethr = cos((double)ddr);
        if( flag == ON )
	    ethi = sin((double)ddr);
	else
	    ethi = -sin((double)ddr);

	FORALLSITES(i,s){
	    if( s->t != (nt-1) ){
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    ddr = s->link[dir].e[j][k].real * s->phase[dir];
		    ddi = s->link[dir].e[j][k].imag * s->phase[dir];
		    s->link[dir].e[j][k].real = ethr*ddr - ethi*ddi;
		    s->link[dir].e[j][k].imag = ethr*ddi + ethi*ddr;
		    ddr = s->boundary[dir].e[j][k].real * s->phase[dir];
		    ddi = s->boundary[dir].e[j][k].imag * s->phase[dir];
		    s->boundary[dir].e[j][k].real = ethr*ddr - ethi*ddi;
		    s->boundary[dir].e[j][k].imag = ethr*ddi + ethi*ddr;
		}
	    }
	    else{
		/* boundary links for t=nt are stored in boundary at t=nt-1 */
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    ddr = s->link[dir].e[j][k].real * s->phase[dir];
		    ddi = s->link[dir].e[j][k].imag * s->phase[dir];
		    s->link[dir].e[j][k].real = ethr*ddr - ethi*ddi;
		    s->link[dir].e[j][k].imag = ethr*ddi + ethi*ddr;
		    ddr = - s->boundary[dir].e[j][k].real * s->phase[dir];
		    ddi = - s->boundary[dir].e[j][k].imag * s->phase[dir];
		    s->boundary[dir].e[j][k].real = ethr*ddr - ethi*ddi;
		    s->boundary[dir].e[j][k].imag = ethr*ddi + ethi*ddr;
		}
	    }
	}
    }
#endif

    phases_in = flag;
}

