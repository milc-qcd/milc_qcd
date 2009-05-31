/*************************** setup.c *******************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//
*/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include <string.h>
#include <lattice_qdp.h>

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

gauge_file *r_parallel_i(char *);
void r_parallel(gauge_file *, field_offset);
void r_parallel_f(gauge_file *);

gauge_file *r_binary_i(char *);
void r_binary(gauge_file *);
void r_binary_f(gauge_file *);

void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
int initial_set();
void make_3n_gathers();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
#ifdef HAVE_QDP
  int i;
#endif
    int prompt;

	/* print banner, get volume, seed */
    prompt=initial_set();
   	/* initialize the node random number generator */
    initialize_prn( &node_prn, iseed, volume+mynode() );
	/* Initialize the layout functions, which decide where sites live */
    setup_layout();
	/* allocate space for lattice, set up coordinate fields */
    make_lattice();
    node0_printf("Made lattice\n"); fflush(stdout);
    init_ferm_links(&fn_links);
#ifdef DM_DU0
    init_ferm_links(&fn_links_dmdu0);
#endif

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

#ifdef HAVE_QDP
  for(i=0; i<4; ++i) {
    shiftdirs[i] = QDP_neighbor[i];
    shiftdirs[i+4] = neighbor3[i];
  }
  for(i=0; i<8; ++i) {
    shiftfwd[i] = QDP_forward;
    shiftbck[i] = QDP_backward;
  }
#endif

node0_printf("Finished setup\n"); fflush(stdout);
    return( prompt );
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
    /* On node zero, read lattice size, seed, and send to others */
    if(mynode()==0){
	/* print banner */
	printf("SU3 with improved KS action\n");
	printf("MIMD version 7\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	printf("Multimass inverter\n");

	status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin, prompt,"nflavors1", &par_buf.nflavors1 );
	IF_OK status += get_i(stdin, prompt,"nflavors2", &par_buf.nflavors2 );
	IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
#ifdef FIX_NODE_GEOM
	IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			       par_buf.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
	IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			       par_buf.ionode_geometry, 4);
#endif
#endif
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
#ifdef FIX_NODE_GEOM
    for(i = 0; i < 4; i++)
      node_geometry[i] = par_buf.node_geometry[i];
#ifdef FIX_IONODE_GEOM
    for(i = 0; i < 4; i++)
      ionode_geometry[i] = par_buf.ionode_geometry[i];
#endif
#endif
    iseed=par_buf.iseed;
    nflavors1=par_buf.nflavors1;
    nflavors2=par_buf.nflavors2;
    dyn_flavors[0] = nflavors1;
    dyn_flavors[1] = nflavors2;
    
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
     char savebuf[128];
     int i;
     Real x;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;
    
    
	/* get couplings and broadcast to nodes	*/
	/* beta, mass */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );
	IF_OK status += get_f(stdin, prompt,"mass1", &par_buf.mass1);
	IF_OK status += get_f(stdin, prompt,"mass2", &par_buf.mass2 );
	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
    
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );
    
	/* error for propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;

	/* Parameters for fpi */
        /* source time slice and increment */
	IF_OK status += get_i(stdin, prompt,"source_start", &par_buf.source_start );
	IF_OK status += get_i(stdin, prompt,"source_inc", &par_buf.source_inc );
	IF_OK status += get_i(stdin, prompt,"n_sources", &par_buf.n_sources );
	
	par_buf.fpi_nmasses = 0;
	IF_OK status += get_i(stdin, prompt,"nmasses",
			      &par_buf.fpi_nmasses );
	if(par_buf.fpi_nmasses > MAX_FPI_NMASSES){
	  printf("Maximum of %d exceeded.\n",MAX_FPI_NMASSES);
	  terminate(1);
	}
	for(i = 0; i < par_buf.fpi_nmasses; i++){
	  IF_OK status += get_f(stdin, prompt,"mass",
				&par_buf.fpi_mass[i]);
	}

	/* Parameters for multimass_inverter. */

        /* point source locations */
	IF_OK status += get_i(stdin, prompt,"n_sources_mminv", 
			      &par_buf.mminv.n_sources );
	for(i = 0; i < par_buf.mminv.n_sources; i++ ){
	  IF_OK status += get_vi(stdin, prompt,"r0", par_buf.mminv.r0[i],4);
	  /* We want an even source */
	  if((par_buf.mminv.r0[i][0] + par_buf.mminv.r0[i][1] + 
	      par_buf.mminv.r0[i][2] + par_buf.mminv.r0[i][3]) % 2 != 0){
	    printf("ERROR: Source coordinate should be even\n");
	    status = 1;
	  }
	}
	IF_OK status += get_i(stdin, prompt,"nmasses_mminv", &par_buf.mminv.nmasses );
	if(par_buf.mminv.nmasses > MAX_MMINV_NMASSES){
	  printf("Maximum of %d exceeded.\n",MAX_FPI_NMASSES);
	  terminate(1);
	}
	for(i = 0; i < par_buf.mminv.nmasses; i++){
	  IF_OK status += get_f(stdin, prompt,"mass",
				&par_buf.mminv.masses[i]);
	}

	/* error for propagator multimass conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator_mminv", &x );
	IF_OK par_buf.mminv.rsqprop = x*x;


        /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

	/* decide about gauge fixing */
    	IF_OK if (prompt!=0)
      		printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
    	IF_OK scanf("%s",savebuf);
    	IF_OK printf("%s\n",savebuf);
    	IF_OK {
      		if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
        		par_buf.fixflag = COULOMB_GAUGE_FIX;
      		}
      	else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
        		par_buf.fixflag = NO_GAUGE_FIX;
      		}
      	else{
        		printf("error in input: fixing_command is invalid\n"); status++;
      		}
    	}


        /* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );

	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				      par_buf.stringLFN );

        /* find out whether or not to save propagator at end */
        IF_OK status += ask_ending_ksprop( stdin, prompt, &(par_buf.kssaveflag),
                                           par_buf.kssavefile );


	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    niter = par_buf.niter;
    nrestart = par_buf.nrestart;
    rsqprop = par_buf.rsqprop;
    beta = par_buf.beta;
    mass1 = par_buf.mass1;
    mass2 = par_buf.mass2;
    n_dyn_masses = 2;
    u0 = par_buf.u0;
    source_start = par_buf.source_start;
    source_inc = par_buf.source_inc;
    n_sources = par_buf.n_sources;
    fpi_nmasses = par_buf.fpi_nmasses;
    for(i = 0; i < fpi_nmasses; i++){
      fpi_mass[i] = par_buf.fpi_mass[i];
    }
    memcpy(&mminv,&par_buf.mminv,sizeof(params_mminv)); 

    startflag = par_buf.startflag;
    saveflag = par_buf.saveflag;
    strcpy(startfile,par_buf.startfile);
    fixflag = par_buf.fixflag;
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);
    kssaveflag = par_buf.kssaveflag;
    strcpy(kssavefile,par_buf.kssavefile);

    /* Do whatever is needed to get lattice */
    if( startflag == CONTINUE ){
        rephase( OFF );
    }
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );
    /* if a lattice was read in, put in KS phases and AP boundary condition */
#ifdef FN
    invalidate_all_ferm_links(&fn_links);
#ifdef DM_DU0
    invalidate_all_ferm_links(&fn_links_dmdu0);
#endif
#endif
    phases_in = OFF;
    rephase( ON );

    /* make table of coefficients and permutations of loops in gauge action */
    make_loop_table();
    /* make table of coefficients and permutations of paths in quark action */
    init_path_table(&ks_act_paths);
    make_path_table(&ks_act_paths, NULL);

    return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
 */
void make_3n_gathers(){
   int i;
#ifdef HAVE_QDP
   int disp[4]={0,0,0,0};
#endif
 
   for(i=XUP;i<=TUP;i++) {
      make_gather(third_neighbor,&i,WANT_INVERSE,
		  ALLOW_EVEN_ODD,SWITCH_PARITY);
   }
   
    /* Sort into the order we want for nearest neighbor gathers,
       so you can use X3UP, X3DOWN, etc. as argument in calling them. */

   sort_eight_neighborlists(X3UP);
#ifdef HAVE_QDP

   for(i=0; i<4; i++) {
     disp[i] = 3;
     neighbor3[i] = QDP_create_shift(disp);
     disp[i] = 0;
   }
#endif
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
