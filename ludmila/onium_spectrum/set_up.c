/******** setup.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)


#include <onium_generic.h>
#include <string.h>

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
	printf("Heavy-Heavy spectroscopy\n");
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");
	status = get_prompt( &prompt );
	IF_OK status += get_i(prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(prompt,"nt", &par_buf.nt );

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

/* read in parameters and coupling constants	*/

int readin(int prompt)  {

/* argument "prompt" is 1 if prompts are to be given for input	*/

int status,i;
float x;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	status=0;
        
	/* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice( prompt, &(par_buf.startflag),
	    par_buf.startfile );
    	
	/* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice( prompt, &(par_buf.saveflag),
	    par_buf.savefile );
    	
	/* Get ensemble values for NERSC archive */
	IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
		  par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
	  status += get_s( prompt,"ensemble_id", par_buf.ensemble_id );
	IF_OK if (par_buf.saveflag == SAVE_SERIAL_ARCHIVE ||
		  par_buf.saveflag == SAVE_PARALLEL_ARCHIVE)
	  status += get_i( prompt,"sequence_number", 
			   &par_buf.sequence_number );

	IF_OK status += get_s(prompt,"antiquark_prop_file", par_buf.a_startfile_w);
	IF_OK status += get_f(prompt,"d1", &par_buf.a_d1 );
	IF_OK status += get_i(prompt,"format", &(par_buf.a_format));
	/* Number of kappas */

	IF_OK status += get_i(prompt,"number_of_kappas", &par_buf.num_kap);
	if( par_buf.num_kap>MAX_KAP ){
	  printf("num_kap = %d must be <= %d!\n", par_buf.num_kap, MAX_KAP);
	  status++;
	}
    
	/* To be safe initialize the following to zero */
	for(i=0;i<MAX_KAP;i++){
	  kap[i] = 0.0;
	}
    
	for(i=0;i<par_buf.num_kap;i++){
	  IF_OK status += get_s(prompt,"wilson_prop_file", par_buf.startfile_w[i]);
	  IF_OK status += get_f(prompt,"kappa", &par_buf.kap[i] );       
	  IF_OK status += get_f(prompt,"d1", &par_buf.d1[i] );
	  IF_OK status += get_i(prompt,"format", &(par_buf.format[i]));
	}
	IF_OK status += get_i(prompt,"number_of_smearings", &par_buf.num_smear );
	for(i=0;i<par_buf.num_smear;i++)
	  IF_OK status += get_s(prompt,"smear_func_file", par_buf.smearfile[i]);

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

   
    startflag = par_buf.startflag;
    saveflag = par_buf.saveflag;
   
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    a_d1 = par_buf.a_d1;
    num_kap = par_buf.num_kap;
    for(i=0;i<par_buf.num_kap;i++){
      kap[i] = par_buf.kap[i]; 
      d1[i] = par_buf.d1[i];
      strcpy(startfile_w[i],par_buf.startfile_w[i]);
      strcpy(savefile_w[i],par_buf.savefile_w[i]);
      format[i] = par_buf.format[i];
    }
    for(i=0;i<par_buf.num_smear;i++)
      strcpy(smearfile[i],par_buf.smearfile[i]);

    strcpy(ensemble_id,par_buf.ensemble_id);

    strcpy(a_startfile_w, par_buf.a_startfile_w); 
    sequence_number = par_buf.sequence_number;
    num_smear = par_buf.num_smear;
    a_format = par_buf.a_format;
    startlat_p = reload_lattice( startflag, startfile );

    return(0);
}
