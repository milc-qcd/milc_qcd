/************************ setup.c ****************************/
/* Handles lattice and parameter setup                       */ 

#define IF_OK if(status==0)

/* definitions, files, and prototypes */
#include "wilson_flow_includes.h"
#include <string.h>	

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int initial_set();

int
setup()
{
  int prompt;
  
  /* print banner, get initial parameters */
  prompt = initial_set();

  /* Initialize the layout functions, which decide where sites live */
  setup_layout();

  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);

  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);

  node0_printf("Finished setup\n"); fflush(stdout);
  return prompt;
}


/* SETUP ROUTINES */
int 
initial_set()
{
  int prompt, status;

  if(mynode()==0){
    /* print banner */
    printf("Wilson/Symanzik Flow application\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());

    /* Read prompt type and lattice dimensions */
    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );

    if(status>0) 
      par_buf.stopflag=1; 
    else 
      par_buf.stopflag=0;
  } /* end if(mynode()==0) */
  
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )
    normal_exit(0);
  
  /* Update global variables with parameters */
  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
  
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  
  return prompt;
}

/* Read in configuration specific parameters */
int
readin(int prompt)
{
  int status;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {
    
    printf("\n\n");
    status=0;

    /* Identify the starting configuration */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag), 
                                         par_buf.startfile);

    /* Get flow parameters */
    IF_OK {
      /* Get flow description */
      if (prompt==1)
        printf("enter 'wilson' or 'symanzik'\n");
      scanf("%s", par_buf.flow_description);
      printf("%s\n", par_buf.flow_description);

      /* Check flow description and set staple flag */
      if( strcmp("wilson", par_buf.flow_description) == 0 ) {
        par_buf.stapleflag = WILSON;
        printf("set staple to wilson\n");
      }
      else if( strcmp("symanzik", par_buf.flow_description) == 0 ) {
        par_buf.stapleflag = SYMANZIK;
        printf("set staple to symanzik\n");
      }
      else {
        printf("Error: flow_description %s is invalid\n", 
               par_buf.flow_description);
        status++;
      }
    } /*end: flow_description IF_OK */

    IF_OK status += get_i(stdin, prompt, "exp_order", &par_buf.exp_order);
    IF_OK status += get_f(stdin, prompt, "stepsize", &par_buf.stepsize);
    IF_OK status += get_f(stdin, prompt, "stoptime", &par_buf.stoptime);

    /* Determine what to do with the final configuration */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag), 
                                       par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin, prompt, par_buf.saveflag, 
                                 par_buf.stringLFN );

    if( status > 0) 
      par_buf.stopflag=1; 
    else
      par_buf.stopflag=0;

  } /* end if(this_node==0) */
  
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));

  if( par_buf.stopflag != 0 )
    return par_buf.stopflag;

  /* Update global variables with parameter buffer */
  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);

  strcpy(flow_description, par_buf.flow_description);
  stapleflag = par_buf.stapleflag;
  exp_order = par_buf.exp_order;
  stepsize = par_buf.stepsize;
  stoptime = par_buf.stoptime;

  saveflag = par_buf.saveflag;
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

  /* Load configuration (no phases in this program) */
  if( par_buf.startflag != CONTINUE )
    startlat_p = reload_lattice( par_buf.startflag, par_buf.startfile );

  return 0;
}
