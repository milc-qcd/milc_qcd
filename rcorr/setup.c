/******** setup.c *********/

/* Read parameters and initialize the lattice */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#define IF_OK if(status==0)

/* Modifications ... */

#include "rcorr_includes.h"

/* Forward declarations */

static int initial_set(void);

int setup()   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;

  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();

  /* Initialize fermion links as unallocated */
  //  init_ferm_links(&fn_links, &ks_act_paths);
  //  init_ferm_links(&fn_links_dmdu0, &ks_act_paths_dmdu0);
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  return(prompt);
}

/* SETUP ROUTINES */
static int initial_set(void){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif

  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("Radial correlator code\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");
    get_utc_datetime(utc_date_time);

    status = get_prompt(stdin,  &prompt );
    IF_OK status += get_i(stdin,prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin,prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin,prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin,prompt,"nt", &param.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   param.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_s(stdin, prompt,"job_id",param.job_id);
    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    normal_exit(0);

  if(prompt==2)return prompt;

  nx=param.nx;
  ny=param.ny;
  nz=param.nz;
  nt=param.nt;

#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif

  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;

  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status, jflav;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    IF_OK status += get_i(stdin, prompt, "number_of_random_sources",
			  &param.nrand);

   if(param.nrand < 2){
      fprintf(stderr, "ERROR: need more than 1 random source to compute correlations\n");
      status++;
    }

    IF_OK status += get_i(stdin, prompt, "number_of_flavors",
			  &param.nflav);

    for(jflav = 0; jflav < param.nflav; jflav++){
      IF_OK status += get_f(stdin, prompt, "charge", &param.charges[jflav]);
      IF_OK status += get_s(stdin, prompt, "file", param.fname[jflav]);
    }

    IF_OK status += get_s(stdin, prompt, "save_corr", param.corrfile);
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
  
  
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )return param.stopflag;

  return 0;

} /* setup.c */

