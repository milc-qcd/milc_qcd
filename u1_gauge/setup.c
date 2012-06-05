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
/* 4. set_mom_table()						*/
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
static void set_mom_table(void);
static int initial_set(void);

/* setup: setting up lattice */
int setup(void)   
{

  int prompt;

  /* print banner, get lattice parameters */
  prompt=initial_set();

  /* initialize the node random number generator */
  initialize_prn(&node_prn,iseed,volume+mynode());

  /* initialize the layout -- lattice across the nodes */
  setup_layout();

  /* allocate space for lattice, set up coordinates */
  make_lattice();
  node0_printf("Made lattice!\n"); fflush(stdout);
  u1_A = create_u1_A_field();

  /* set up nearest neighbor gathers */
  make_nn_gathers();
  node0_printf("Made nn_gathers!\n"); fflush(stdout);

  /* set up momentum table */
  set_mom_table();
  node0_printf("Momentum table prepared!\n"); fflush(stdout);

  node0_printf("Finished setup!\n"); fflush(stdout);
  return(prompt);

}  /* end of setup() */

/* basic setup */
static int initial_set(void)
{

  int prompt,status;

  if(mynode()==0)
    {
    /* print banner */
    printf("U(1) [Coulomb gauge-fixed] gauge field generation ... \n");
    printf("stored as A(mu,x), convert to U(mu,x) to couple to fermions.\n");
    printf("Machine = %s, with %d nodes!\n",machine_type(),numnodes());
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

  this_node=mynode();
  number_of_nodes=numnodes();
  volume=nx*ny*nz*nt;

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

/*------------------------------------------------------------------*/
/* Added to support MPP gauge field generation that gives the same
   random field as the original single-process field generation
   code.  This is the single-process node_index function. */

static size_t scalar_node_index(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];
  int sum = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
    sum += coords[d];
  }

  if(sum % 2 == 0)return rank/2;
  else return (rank + sites_on_node)/2;

}

/* ************************************************************	*/

/* setup momentum table */
static void set_mom_table(void)
{

  int x,y,z,t;
  int i,dir;
  size_t j,k;
  int dimn[4] = {nx, ny, nz, nt};
  site *s;

  /* initialization */
  junk_id=-9999;
  
  latin=(int *)malloc(sites_on_node*sizeof(int));

  FORALLSITES(i,s){
    latin[i] = junk_id;
  }

  /* This algorithm is needed to achieve backward compatibility with
     the single-process code */

  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++){
	  int st[4] = { x, y, z, t };
	  int stin[4];
	  
	  /* The parity-inverted coordinate */
	  FORALLUPDIR(dir){
	    stin[dir] = (dimn[dir] - st[dir]) % dimn[dir];
	  }
	    
	  /* If the current node has the inverse j, set latin */
	  j = scalar_node_index(stin, 4, dimn);
	  k = scalar_node_index(st, 4, dimn);
	  if(j <= k){
	    if(node_number(stin[XUP], stin[YUP], stin[ZUP], stin[TUP]) == this_node){
	      latin[j]=k;
	    }
	  }
	}

  // DEBUG
  //  FORALLSITES(i,s){
  //    printf("%d -> %d , %d %d %d %d\n",i, latin[i], s->x, s->y ,s->z, s->t);
  //  }

} /* end of set_mom_table() */

Real sqr(Real val)
{

  return(val*val);

} /* end of sqr() */

/* ************************************************************	*/

