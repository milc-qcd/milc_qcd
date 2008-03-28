/************************* setup_baryon.c *********************************/
/* MIMD version 7 */

#define IF_OK if(status==0)
#include "ks_hl_spectrum_includes.h"
#include <string.h>
#include <time.h>

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int initial_set();
void make_gen_pt(),setup_layout();
int prompt;

	/* print banner, get volume, nflavors, seed */
    prompt=initial_set();
   	/* Initialize the layout functions, which decide where sites live */
    setup_layout();
	/* allocate space for lattice, set up coordinate fields */
    make_lattice();
	/* set up neighbor pointers and comlink structures */
    make_nn_gathers();

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, and send to others */
    if(mynode()==0){
	/* print banner */
	printf("Baryon spectroscpy with Kogut-Susskind light fermions and Wilson heavy quarks\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");
	status = get_prompt(stdin, &prompt );
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
    total_iters=0;
    return(prompt);
}



/* read in parameters and coupling constants for singly heavy baryon*/

int readin(int prompt)  {
  
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status,i;
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    status=0;
    
    /* Unless calculate rotation of fields, we don't need lattice configurations. -H */
    /* find out what kind of starting lattice to use */
    /*	IF_OK status += ask_starting_lattice( prompt, &(par_buf.startflag),
	par_buf.startfile ); */
    
    /* find out what to do with lattice at end */
    /*IF_OK status += ask_ending_lattice( prompt, &(par_buf.saveflag),
      par_buf.savefile );*/
    
    /* Get ensemble values for NERSC archive */
    /*IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
      par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
      status += get_s( prompt,"ensemble_id", par_buf.ensemble_id );
      IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
      par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
      status += get_i( prompt,"sequence_number", 
      &par_buf.sequence_number );*/
    /* Number of kappas */
    
    IF_OK status += get_i(stdin, prompt, "number_of_kappas", 
			  &par_buf.num_kap );
    if( par_buf.num_kap>MAX_KAP ){
      printf("num_kap = %d must be <= %d!\n", par_buf.num_kap, MAX_KAP);
      status++;
    }
    
    /* To be safe initialize the following to zero */
    for(i=0;i<MAX_KAP;i++){
      kap[i] = 0.0;
    }
    
    for(i=0;i<par_buf.num_kap;i++){
      IF_OK status += ask_starting_wprop(stdin, prompt, 
					 &par_buf.startflag_w[i], 
					 par_buf.startfile_w[i]);
      IF_OK status += get_f(stdin, prompt, "kappa", &par_buf.kap[i] );       
    }
    
    /*get information for multiple light quarks*/
    IF_OK status += get_i(stdin, prompt, "number_of_light", 
			  &par_buf.num_light );
    if( par_buf.num_light>MAX_LIGHT ){
      printf("num_light = %d must be <= %d!\n", par_buf.num_light, MAX_LIGHT);
      status++;
    }
    
    for(i=0;i<par_buf.num_light;i++){
      IF_OK status += ask_starting_ksprop(stdin, prompt, 
					  &par_buf.start_ks_light_flag[i], 
					  par_buf.start_ks_light_file[i]);
      IF_OK status += get_f(stdin, prompt,"m_light", &par_buf.m_light[i] );
    }
    
    /*get information for multiple strange quarks*/
    IF_OK status += get_i(stdin, prompt, "number_of_strange", 
			  &par_buf.num_strange );
    if( par_buf.num_strange>MAX_STRANGE ){
      printf("num_strange = %d must be <= %d!\n", par_buf.num_strange, 
	     MAX_STRANGE);
      status++;
    }
    
    for(i=0;i<par_buf.num_strange;i++){
      IF_OK status += ask_starting_ksprop(stdin, prompt, 
					  &par_buf.start_ks_strange_flag[i], 
					  par_buf.start_ks_strange_file[i]);
      IF_OK status += get_f(stdin, prompt, "m_strange", 
			    &par_buf.m_strange[i] );
    }
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    
  } /* end if(this_node==0) */
  
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 ) return par_buf.stopflag;
  
  
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  
  
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  num_kap = par_buf.num_kap;
  for(i=0;i<par_buf.num_kap;i++){
    kap[i] = par_buf.kap[i]; 
    strcpy(startfile_w[i],par_buf.startfile_w[i]);
    strcpy(savefile_w[i],par_buf.savefile_w[i]);
    startflag_w[i] = par_buf.startflag_w[i];
  }
  
  num_light = par_buf.num_light;
  for(i=0;i<par_buf.num_light;i++){
    m_light[i] = par_buf.m_light[i];
    strcpy(start_ks_light_file[i],par_buf.start_ks_light_file[i]);
    start_ks_light_flag[i] = par_buf.start_ks_light_flag[i];
  }
  
  num_strange = par_buf.num_strange;
  for(i=0;i<par_buf.num_strange;i++){
    m_strange[i] = par_buf.m_strange[i];
    strcpy(start_ks_strange_file[i],par_buf.start_ks_strange_file[i]);
    start_ks_strange_flag[i] = par_buf.start_ks_strange_flag[i];
  }
  
  return(0);
}

