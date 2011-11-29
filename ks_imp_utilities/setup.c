/************************ setup.c ****************************/
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
#include <lattice_qdp.h>

EXTERN gauge_header start_lat_hdr;
gauge_file *gf;

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;
void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
void make_3n_gathers();

int
setup()
{
  int initial_set();
  int prompt;

  /* print banner, get volume, nflavors1,nflavors2, seed */
  prompt = initial_set();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);
  /* Mark t_longlink and t_fatlink unallocted */
  // init_ferm_links(&fn_links, &ks_act_paths);
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

  node0_printf("Finished setup\n"); fflush(stdout);
  return( prompt );
}

/* SETUP ROUTINES */
int
initial_set()
{
  int prompt,status;
  /* On node zero, read lattice size, seed and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Inversion checking\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

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

  if( par_buf.stopflag != 0 ) return par_buf.stopflag;

  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
  iseed=par_buf.iseed;

  dyn_flavors[0] = 1;
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  return(prompt);
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  Real x;
#ifdef CHECK_INVERT
  char invert_string[16];
#endif

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {

    printf("\n\n");
    status=0;

    /* get couplings and broadcast to nodes	*/
    /* mass, u0 */
    IF_OK status += get_f(stdin, prompt,"mass", &par_buf.mass );
    IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );

    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );
    
    /* error for propagator conjugate gradient */
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
    IF_OK par_buf.rsqprop = x*x;

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
					  par_buf.startfile );

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
					par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				  par_buf.stringLFN );

    /* find out what to do with longlinks at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.savelongflag),
					par_buf.savelongfile );

    /* find out what to do with fatlinks at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.savefatflag),
					par_buf.savefatfile );

    /* find out what kind of color vector source to use */
    IF_OK status += ask_color_vector( prompt, &(par_buf.srcflag),
				      par_buf.srcfile );
#ifdef CHECK_INVERT
    /* find out what kind of color vector result to use */
    IF_OK status += ask_color_vector( prompt, &(par_buf.ansflag),
				      par_buf.ansfile );
#else
    /* find out what kind of color matrix result to use */
    IF_OK status += ask_color_matrix( prompt, &(par_buf.ansflag),
				      par_buf.ansfile );
#endif

#ifdef CHECK_INVERT
    /* find out which inversion to check */
    IF_OK status += get_s(stdin,  prompt, "invert", invert_string);
    if(status == 0){
      if(strcmp(invert_string,"M")==0)
	par_buf.inverttype = INVERT_M;
      else if(strcmp(invert_string,"MdaggerM")==0)
	par_buf.inverttype = INVERT_MdaggerM;
      else{
	printf("Unrecognized invert string %s\n",invert_string);
	status++;
      }
    }
#endif

    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));

  if( par_buf.stopflag != 0 ) return par_buf.stopflag;

  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqprop = par_buf.rsqprop;
  mass = par_buf.mass;
  n_dyn_masses = 1;
  u0 = par_buf.u0;
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  savelongflag = par_buf.savelongflag;
  savefatflag = par_buf.savefatflag;
  srcflag = par_buf.srcflag;
  ansflag = par_buf.ansflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);
  strcpy(savelongfile,par_buf.savelongfile);
  strcpy(savefatfile,par_buf.savefatfile);
  strcpy(srcfile,par_buf.srcfile);
  strcpy(ansfile,par_buf.ansfile);
  inverttype = par_buf.inverttype;

  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile );

  /* if a lattice was read in, put in KS phases and AP boundary condition */
  phases_in = OFF;
  rephase( ON );

#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif

  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();

  fn_links = create_fermion_links_from_site(PRECISION, 0, NULL);

  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern;
   make_lattice() and  make_nn_gathers() must be called first,
   preferably just before calling make_3n_gathers().
*/
void
make_3n_gathers()
{
  int i;

  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }

  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */

  sort_eight_gathers(X3UP);

}


/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

void
third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
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
