/* ************************************************************	*/
/*								*/
/* 			     SETUP.C		   		*/
/*								*/
/* This program setup the lattice in various stages through  	*/
/* the routine 'setup'. It reads various lattice parameters     */
/* through  'initial_set' and allocate  spaces for  lattice     */
/* and other storage matrices via 'make_lattice'.               */
/*								*/
/* 1. setup()							*/
/* 2. initial_set()						*/
/* 3. readin()							*/
/*								*/
/* Last updated on 10.30.07					*/
/*								*/
/* ************************************************************	*/

#define IF_OK if(status==0)

#include "include_u1g.h"
#include "params.h"
#include <string.h>

params par_buf;

/* Forward declarations */
static int initial_set(void);

/* setup: setting up lattice */
int setup(void)   
{

  int prompt;

  /* print banner, get lattice parameters */
  prompt=initial_set();

  /* initialize the layout -- lattice across the nodes */
  setup_layout();
  this_node=mynode();

  /* initialize the node random number generator */
  initialize_prn(&node_prn,iseed,volume+mynode());

  /* allocate space for lattice, set up coordinates */
  make_lattice();
  node0_printf("Made lattice!\n"); fflush(stdout);
  u1_A = create_u1_A_field();

  /* set up nearest neighbor gathers */
  make_nn_gathers();
  node0_printf("Made nn_gathers!\n"); fflush(stdout);

  node0_printf("Finished setup!\n"); fflush(stdout);
  return(prompt);

}  /* end of setup() */

/* basic setup */
static int initial_set(void)
{

  int prompt=0,status;

  if(mynode()==0)
    {
    /* print banner */
    printf("U(1) Gauge field generation ... \n");
    printf("stored as A(mu,x), convert to U(mu,x) to couple to fermions.\n");
    printf("Machine = %s, with %d nodes!\n",machine_type(),numnodes());
    node0_printf("Options selected...\n");
    show_u1_gauge_field_opts(GAUGE, SCHEME);
    show_generic_opts();
    time_stamp("Start");
    status=get_prompt(stdin,&prompt);

    /* get lattice size */
    IF_OK status+=get_i(stdin,prompt,"nx",&par_buf.nx);
    IF_OK status+=get_i(stdin,prompt,"ny",&par_buf.ny);
    IF_OK status+=get_i(stdin,prompt,"nz",&par_buf.nz);
    IF_OK status+=get_i(stdin,prompt,"nt",&par_buf.nt);

    /* initializing random number */
    IF_OK status+=get_i(stdin,prompt,"iseed",&par_buf.iseed);

    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* if(mynode()==0)-ends */

  /* broadcasts parameters from node0 to all nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if(par_buf.stopflag!=0) normal_exit(0);

  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
  iseed=par_buf.iseed;

  number_of_nodes=numnodes();
  volume=(size_t)nx*ny*nz*nt;

  return(prompt);

}  /* end of initial_set() */

/* read in parameters and couplings */
int readin(int prompt)
{

  int status;

  if(this_node==0)
    {
    status=0;
    printf("\n\n");

    /* get couplings */
    IF_OK status+=get_f(stdin,prompt,"electron_charge",&par_buf.echarge);

    /* what kind of starting U(1) lattice to use, read filename */
    IF_OK status+=ask_starting_u1_lattice(stdin,prompt,
			&(par_buf.start_u1flag),par_buf.start_u1file);

    /* what to do with U(1) lattice at end, read filename */
    IF_OK status+=ask_ending_u1_lattice(stdin,prompt,
			&(par_buf.save_u1flag),par_buf.save_u1file);

    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* if(this_node==0)-ends */

  /* broadcasts parameters from node0 to all nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if(par_buf.stopflag!=0) return(par_buf.stopflag);

  echarge=par_buf.echarge;
  
  start_u1flag=par_buf.start_u1flag;
  
  strcpy(start_u1file,par_buf.start_u1file);
  
  save_u1flag=par_buf.save_u1flag;
  
  strcpy(save_u1file,par_buf.save_u1file);
  
  start_u1lat_p=reload_u1_lattice(start_u1flag,start_u1file);
  
  return(0);
} /* end of readin() */

/* ************************************************************	*/
