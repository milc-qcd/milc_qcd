/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "gluon_prop_includes.h"

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int initial_set();
void make_3n_gathers();
void periodic_bc();
 int prompt;

        /* print banner, get volume, nflavors, seed */
    prompt=initial_set();
        /* Initialize the layout functions, which decide where sites live */
    setup_layout();
        /* allocate space for lattice, set up coordinate fields */
    make_lattice();
    node0_printf("Made lattice\n"); fflush(stdout);
        /* set up neighbor pointers and comlink structures */
    make_nn_gathers();
    node0_printf("Made nn gathers\n"); fflush(stdout);
#ifdef QUARK_PROP
     /* Mark t_longlink and t_fatlink as unallocated */
    //    init_ferm_links(&fn_links, &ks_act_paths);
        /* set up 3rd nearest neighbor pointers and comlink structures
           code for this routine is below  */
    make_3n_gathers();
    node0_printf("Made 3nn gathers\n"); fflush(stdout);
        /* set up K-S phase vectors, boundary conditions */
    phaseset();
	/* make boundary condition in time periodic */
    periodic_bc();

#endif

    node0_printf("Finished setup\n"); fflush(stdout);
    return(prompt);
}

#if defined(QUARK_PROP) && FERM_ACTION == HISQ
static int n_naiks = 1;
static double eps_naik[MAX_NAIK];
#endif

/* SETUP ROUTINES */
int initial_set(){
int prompt,status;

    /* On node zero, read lattice size and send to others */
    if(mynode()==0){
        /* print banner */
#ifdef GLUON_PROP
        printf("SU3 gluon propagator\n");
#endif
#ifdef QUARK_PROP
        printf("SU3 quark propagator\n");
#ifdef NAIK
        printf("With improved fermions with Naik term\n");
#endif
#endif
        printf("MIMD version 7 $Name:  $\n");
        printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");
        status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );

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
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=nx*ny*nz*nt;
    return(prompt);
}

/* read in parameters and coupling constants    */
int readin(int prompt) {
/* read in parameters for su3 gluon/quark propagator       */
/* argument "prompt" is 1 if prompts are to be given for input */

  int status,status2;
#ifdef QUARK_PROP
  int i;
  Real x;
#endif
  char savebuf[128];

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;

	/* get couplings and broadcast to nodes	*/
	/* beta, mass */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

#ifdef QUARK_PROP
	IF_OK status += get_i(stdin, prompt,"number_of_masses", &par_buf.num_mass );
	if (par_buf.num_mass > MAX_NUM_MASS){
	    printf("num_mass = %d must be <= %d!\n",
		par_buf.num_mass, MAX_NUM_MASS);
	    status++;
	}

	/* To be safe initialize the following to zero */
	for(i=0;i<MAX_NUM_MASS;i++){
	    mass[i] = 0.0;
	}
	for(i=0;i<par_buf.num_mass;i++){
	    IF_OK status += get_f(stdin, prompt,"mass", &par_buf.mass[i] );
	}

	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );

	/* error for propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;
#else
#ifdef IMP_GFIX
	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );
#endif
#endif

	/* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

 
	IF_OK if (prompt==1) printf(
	    "enter 'no_gauge_fix', 'landau_gauge_fix', or 'coulomb_gauge_fix'\n");
	IF_OK status2=scanf("%s",savebuf);
	IF_OK {
	    if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
		par_buf.fixflag = COULOMB_GAUGE_FIX;
		printf("fixing to coulomb gauge\n");
	    }
	    else if(strcmp("landau_gauge_fix",savebuf) == 0 ) {
		par_buf.fixflag = LANDAU_GAUGE_FIX;
		printf("fixing to landau gauge\n");
	    }
	    else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
		par_buf.fixflag = NO_GAUGE_FIX;
		printf("NOT fixing the gauge\n");
	    }
	    else{
		printf("error in input: fixing_command %s is invalid\n",savebuf); 
		status++;
	    }
#ifndef GFIX
	    if(par_buf.fixflag != NO_GAUGE_FIX) {
		printf("Gauge fixing not allowed in this version\n");
		status++;
	    }
#endif
	}
 
	IF_OK if (prompt==1) printf(
	    "enter 'no_gauge_fix', 'landau_gauge_fix', or 'coulomb_gauge_fix'\n");
	IF_OK status2=scanf("%s",savebuf);
	IF_OK {
	    if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
		par_buf.fixflag_ft = COULOMB_GAUGE_FIX;
		printf("FFT fixing to coulomb gauge\n");
	    }
	    else if(strcmp("landau_gauge_fix",savebuf) == 0 ) {
		par_buf.fixflag_ft = LANDAU_GAUGE_FIX;
		printf("FFT fixing to landau gauge\n");
	    }
	    else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
		par_buf.fixflag_ft = NO_GAUGE_FIX;
		printf("NOT FFT fixing the gauge\n");
	    }
	    else{
		printf("error in input: fixing_command %s is invalid\n",savebuf); 
		status++;
	    }
#ifndef GFIX
	    if(par_buf.fixflag_ft != NO_GAUGE_FIX) {
		printf("FFT gauge fixing not allowed in this version\n");
		status++;
	    }
#endif
	}
 
	/* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );
	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
	    par_buf.stringLFN );
 
#ifdef QUARK_PROP
	/* find out starting propagator */
	for(i=0;i<par_buf.num_mass;i++){
	  IF_OK status += ask_starting_ksprop( stdin, prompt,&par_buf.ksstartflag[i],
				par_buf.ksstartfile[i]);
	}

	/* find out whether to do inversion */
	for(i=0;i<par_buf.num_mass;i++){
	    IF_OK status2 = scanf("%s",savebuf);
	    IF_OK{
		if(strcmp("run_inversion",savebuf) == 0 ){
		    par_buf.run_CG_flag[i] = 1;
		    printf("Inversion to be done\n");
		}
		else if(strcmp("dont_run_inversion",savebuf) == 0 ){
		    if(par_buf.ksstartflag[i] != FRESH ){
			par_buf.run_CG_flag[i] = 0;
			printf("No inversion to be done\n");
		    }
		    else{
			/* With fresh start need inversion! */
			par_buf.run_CG_flag[i] = 1;
			printf("Inversion to be done\n");
		    }
		}
		else{
		    printf("error in input: invert %s is invalid\n",savebuf); 
		    status++;
		}
	    }
	}

	/* what to do with computed propagator */
	for(i=0;i<par_buf.num_mass;i++){
	  IF_OK status += ask_ending_ksprop( stdin, prompt,&par_buf.kssaveflag[i],
				par_buf.kssavefile[i]);
	}
#endif

	/* send parameter structure */
	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 ) return par_buf.stopflag;

    startflag = par_buf.startflag;
    fixflag = par_buf.fixflag;
    fixflag_ft = par_buf.fixflag_ft;
    saveflag = par_buf.saveflag;
    beta = par_buf.beta;
#ifdef QUARK_PROP
    num_mass = par_buf.num_mass;
    for(i=0;i<par_buf.num_mass;i++){
	mass[i] = par_buf.mass[i];
	run_CG_flag[i] = par_buf.run_CG_flag[i];
	ksstartflag[i] = par_buf.ksstartflag[i];
	kssaveflag[i] = par_buf.kssaveflag[i];
	strcpy(ksstartfile[i],par_buf.ksstartfile[i]);
	strcpy(kssavefile[i],par_buf.kssavefile[i]);
    }
    u0 = par_buf.u0;
    niter = par_buf.niter;
    nrestart = par_buf.nrestart;
    rsqprop = par_buf.rsqprop;
#else
#ifdef IMP_GFIX
    u0 = par_buf.u0;
#endif
#endif
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );
#ifdef QUARK_PROP
    phases_in = OFF;

    /* make table of coefficients and permutations of paths in quark action */
    //init_path_table(fn_links.ap);
    //make_path_table(fn_links.ap, NULL);
    rephase( ON );
  /* Set uptions for fermion links */
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif
  
#if FERM_ACTION == HISQ
  /* WARNING: Not fully supported.  We need to read the Naik epsilons first */
  fn_links = create_fermion_links_from_site(PRECISION, n_naiks, eps_naik);
#else
  fn_links = create_fermion_links_from_site(PRECISION, 0, NULL);
#endif

    rephase( OFF );
    
#endif

    /* For archive writing, for now */
    strcpy(ensemble_id,"none");
    sequence_number = 0;

    return(0);
} /*readin()*/


#ifdef QUARK_PROP

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
void make_3n_gathers(){
   int i;
   void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
 
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

void third_neighbor(int x,int y,int z,int t,int *dirpt,int FB,int *xp,int *yp,int *zp,int *tp)
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

/* Interchange periodic and anti-periodic bc in time direction */
void periodic_bc(){
register site *s;
register int i;

    FORALLSITES(i,s) if(s->t == nt-1){
	s->phase[TUP] = -s->phase[TUP];
    }
}
#endif	/* #ifdef QUARK_PROP */
