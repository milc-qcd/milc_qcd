/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routine for checking the derivative
// of the reunitarization calculation
*/

#define IF_OK if(status==0)

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#include <lattice_qdp.h>

/* Each node has a params structure for passing simulation parameters */
#include "params.h"

/* Forward declarations */

int initial_set();

int
setup()
{
  int prompt;

  /* print banner, get volume, nflavors1,nflavors2, seed */
  prompt = initial_set();
  if(prompt == 2)return prompt;
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);

  node0_printf("Finished setup\n"); fflush(stdout);
  return  prompt;
}

/* SETUP ROUTINES */
int
initial_set()
{
  int prompt=0,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
  /* On node zero, read lattice size, seed and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Reunitarization checking\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());

    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &param.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   param.ionode_geometry, 4);
#endif
#endif
    IF_OK {
      int iseed_in;
      status += get_i(stdin, prompt,"iseed", &iseed_in);
      param.iseed = iseed_in;
    }
    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    return param.stopflag;

  if(prompt==2)return prompt;

  nx=param.nx;
  ny=param.ny;
  nz=param.nz;
  nt=param.nt;
  iseed=param.iseed;

#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif

  number_of_nodes = numnodes();
  volume=(size_t)nx*ny*nz*nt;
#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = 0;
#endif
      
  return prompt;
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int i, status, current_index;
  char savebuf[128];

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {

    printf("\n\n");
    status=0;

    /* Eigenpairs not supported */
    param.eigen_param.Nvecs = 0;

    /* Is there a file with the fiducial result? */
    IF_OK status += ask_color_matrix( prompt, &(param.ansflag[0]),
				      param.ansfile[0] );
    
    /* Is there a file with the fiducial result? */
    IF_OK status += ask_color_matrix( prompt, &(param.ansflag[1]),
				      param.ansfile[1] );
    
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 ){
    return param.stopflag;
  }

  if(prompt==2)return 0;

  return 0;
}

