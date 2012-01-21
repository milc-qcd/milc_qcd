/************************** setup.c ***************************/
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

#define IF_OK if(status==0)

#include "ks_eig_includes.h"	/* definitions files and prototypes */
//#include "lattice_qdp.h"
#include <string.h>
#ifdef HAVE_QOP
#include "../include/generic_qop.h"
#endif

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

gauge_file *r_parallel_i(char *);
void r_parallel(gauge_file *, field_offset);
void r_parallel_f(gauge_file *);

gauge_file *r_binary_i(char *);
void r_binary(gauge_file *);
void r_binary_f(gauge_file *);

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

//#ifdef HAVE_QDP
//void
//initial_li(QLA_Int *li, int coords[])
//{
//  int i,t;
// 
//  t = coords[0];
//  for(i=1; i<4; i++) {
//    t = t*QDP_coord_size(i) + coords[i];
//  }
//  *li = t;
//}
// 
//void
//make_rand_seed(void)
//{
//  QDP_Int *li;
//
//  rand_state = QDP_create_S();
//  li = QDP_create_I();
//
//  QDP_I_eq_func(li, initial_li, QDP_all);
//  QDP_S_eq_seed_i_I(rand_state, iseed, li, QDP_all);
//
//  QDP_destroy_I(li);
//}
//#endif

int  setup()   {
    int initial_set();
    void make_gen_pt();
    void make_3n_gathers(),
        setup_layout();
    int prompt;
//#ifdef HAVE_QDP
//    int i;
//#endif

	/* print banner, get volume, seed */
    prompt=initial_set();
   	/* initialize the node random number generator */
    initialize_prn( &node_prn, iseed, volume+mynode() );
	/* Initialize the layout functions, which decide where sites live */
    setup_layout();
	/* allocate space for lattice, set up coordinate fields */
    make_lattice();
    node0_printf("Made lattice\n"); fflush(stdout);
    //init_ferm_links(&fn_links, &ks_act_paths);
	/* set up neighbor pointers and comlink structures
	   code for this routine is in com_machine.c  */
    make_nn_gathers();
node0_printf("Made nn gathers\n"); fflush(stdout);
#ifdef FN
	/* set up 3rd nearest neighbor pointers and comlink structures
	   code for this routine is below  */
    make_3n_gathers();
node0_printf("Made 3nn gathers\n"); fflush(stdout);
#endif
	/* set up K-S phase vectors, boundary conditions */
    phaseset();

#if HAVE_QOP
  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    node0_printf("setup: Error initializing QOP\n");
    terminate(1);
  }
#endif

//#ifdef HAVE_QDP
//    make_rand_seed();
//node0_printf("Made random seed\n"); fflush(stdout);
//
//  for(i=0; i<4; ++i) {
//    shiftdirs[i] = QDP_neighbor[i];
//    shiftdirs[i+4] = neighbor3[i];
//  }
//  for(i=0; i<8; ++i) {
//    shiftfwd[i] = QDP_forward;
//    shiftbck[i] = QDP_backward;
//  }
//#endif

node0_printf("Finished setup\n"); fflush(stdout);
    return( prompt );
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, and send to others */
    if(mynode()==0){
	/* print banner */
	printf("SU3 with improved KS action\n");
	printf("Eigenvalues and eigenvectors\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

	gethostname(hostname, 128);
	printf("Host(0) = %s\n",hostname);
	printf("Username = %s\n", getenv("USER"));
	time_stamp("start");
	get_utc_datetime(utc_date_time);
	
	/* Print list of options selected */
	node0_printf("Options selected...\n");
	show_generic_opts();
	show_generic_ks_opts();
	
#if FERM_ACTION == HISQ
	show_su3_mat_opts();
	show_hisq_links_opts();
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
	/* beta, mass */
	IF_OK status += get_f(stdin, prompt,"mass", &par_buf.mass );
	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

	/* maximum no. of conjugate gradient iterations */
	IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
	/* maximum no. of conjugate gradient restarts */
	IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );
    
	/* error per site for conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_per_site", &x );
	IF_OK par_buf.rsqmin = x*x;   /* rsqmin is r**2 in conjugate gradient */
	    /* New conjugate gradient normalizes rsqmin by norm of source */
    
	/* error for propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;
	IF_OK status += get_i(stdin, prompt,"Number_of_eigenvals", &par_buf.Nvecs );
	IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &par_buf.MaxIter );
	IF_OK status += get_i(stdin, prompt,"Restart_Rayleigh", &par_buf.Restart );
	IF_OK status += get_i(stdin, prompt,"Kalkreuter_iters", &par_buf.Kiters );
	IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", 
			      &par_buf.eigenval_tol );
	IF_OK status += get_f(stdin, prompt,"error_decrease", &par_buf.error_decr);

        /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )return par_buf.stopflag;

    niter = par_buf.niter;
    nrestart = par_buf.nrestart;
    rsqmin = par_buf.rsqmin;
    rsqprop = par_buf.rsqprop;
    mass = par_buf.mass;
    u0 = par_buf.u0;
    Nvecs = par_buf.Nvecs ;
    MaxIter = par_buf.MaxIter ;
    Restart = par_buf.Restart ;
    Kiters = par_buf.Kiters ;
    eigenval_tol = par_buf.eigenval_tol ;
    error_decr = par_buf.error_decr ;
    startflag = par_buf.startflag;
    strcpy(startfile,par_buf.startfile);

#if FERM_ACTION == HISQ
    n_naiks = 1;
    eps_naik[0] = 0.0;
#endif

    /* Do whatever is needed to get lattice */
    if( startflag == CONTINUE ){
        rephase( OFF );
    }
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );
    /* if a lattice was read in, put in KS phases and AP boundary condition */
    phases_in = OFF;
    rephase( ON );

  /* Set uptions for fermion links */
    
#ifdef DBLSTORE_FN
    /* We want to double-store the links for optimization */
    fermion_links_want_back(1);
#endif
    
#if FERM_ACTION == HISQ
    fn_links = create_fermion_links_from_site(PRECISION, n_naiks, eps_naik);
#else
    fn_links = create_fermion_links_from_site(PRECISION, 0, NULL);
#endif
    
//    node0_printf("Calling for path table\n");fflush(stdout);
//    /* make table of coefficients and permutations of paths in quark action */
//    init_path_table(fn_links.ap);
//    make_path_table(fn_links.ap, NULL);
//    node0_printf("Done with path table\n");fflush(stdout);

    return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
 */

void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);

void make_3n_gathers(){
   int i;
//#ifdef HAVE_QDP
//   int disp[4]={0,0,0,0};
//#endif
 
   for(i=XUP;i<=TUP;i++) {
      make_gather(third_neighbor,&i,WANT_INVERSE,
		  ALLOW_EVEN_ODD,SWITCH_PARITY);
   }
   
    /* Sort into the order we want for nearest neighbor gathers,
       so you can use X3UP, X3DOWN, etc. as argument in calling them. */

   sort_eight_gathers(X3UP);

//#ifdef HAVE_QDP
//  for(i=0; i<4; i++) {
//    disp[i] = 3;
//    neighbor3[i] = QDP_create_shift(disp);
//    disp[i] = 0;
//  }
//#endif
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
