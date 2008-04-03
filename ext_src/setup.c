/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

* 5/30/07 Created from setup_cl.c */

//  $Log: setup.c,v $
//  Revision 1.2  2008/04/03 11:43:28  detar
//  Fix precision bug: If there is a blank after 1 or 2, the libraries get the wrong name.
//
//  Revision 1.1  2008/03/28 15:24:10  detar
//  Add
//
//  Revision 1.2  2007/11/09 16:15:58  detar
//  Support changes in propagator I/O
//
//  Revision 1.1  2007/10/07 20:02:32  detar
//  Add new application.  Generarlizes clover_invert.
//
//


#include "ext_src_includes.h"
#include "lattice_qdp.h"
#include <string.h>
int initial_set();

#include "params.h"

int setup()   {
  int prompt;
#ifdef HAVE_QDP
  int i;
#endif

  /* print banner, get volume */
  prompt=initial_set();
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();

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
  
  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 clover valence fermions\n");
    printf("MIMD version 7 $Name:  $\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    time_stamp("start");
    
    status = get_prompt(stdin,  &prompt );
    
    IF_OK status += get_i(stdin,prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin,prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin,prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin,prompt,"nt", &param.nt );
    IF_OK status += get_s(stdin, prompt,"job_id",param.job_id);
    
    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    normal_exit(0);

  nx=param.nx;
  ny=param.ny;
  nz=param.nz;
  nt=param.nt;
  
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  char savebuf[128];
  int i;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    /*------------------------------------------------------------*/
    /* Quarks and their sources                                   */
    /*------------------------------------------------------------*/

    /* Number of quarks */
    IF_OK status += get_i(stdin,prompt,"number_of_quarks", &param.num_qk );
    if( param.num_qk>MAX_QK ){
      printf("num_qk = %d must be <= %d!\n", param.num_qk, MAX_QK);
      status++;
    }

    IF_OK for(i = 0; i < param.num_qk; i++){
    
      /* Quark parameters */
      IF_OK status += get_s(stdin, prompt,"quark_type", savebuf );
      IF_OK {
	if(strcmp(savebuf,"clover") == 0)param.qk_type[i] = CLOVER_TYPE;
	else if(strcmp(savebuf,"KS") == 0)param.qk_type[i] = KS_TYPE;
	else {
	  printf("Unknown quark type %s\n",savebuf);
	  status++;
	}
      }
      if(param.qk_type[i] == CLOVER_TYPE){
	int source_type, t0, saveflag_s;
	char descrp[MAXDESCRP];
	char savefile_s[MAXFILENAME];

	IF_OK status += ask_starting_wprop( stdin, prompt, 
					    &param.startflag_w[i],
					    param.startfile_w[i]);
	
	IF_OK status += 
	  ask_output_w_quark_source_file( stdin, prompt, &saveflag_s,
					  &source_type, &t0, descrp,
					  savefile_s );
	
	IF_OK {
	  if(source_type == DIRAC_FIELD_FILE){
	    init_wqs(&param.dst_wqs[i]);
	    param.dst_wqs[i].type = source_type;
	    param.dst_wqs[i].t0   = t0;
	    param.dst_wqs[i].flag = saveflag_s;
	    strcpy(param.dst_wqs[i].descrp, descrp);
	    strcpy(param.dst_wqs[i].source_file, savefile_s);
	    param.dst_type[i] = CLOVER_TYPE;

	    init_wqs(&param.snk_wqs[i]);
	    status += get_w_quark_sink(stdin, prompt, &param.snk_wqs[i]);
	    param.snk_wqs[i].t0   = t0;
	  } else {
	    printf("Unsupported output source type\n");
	    status++;
	  }
	}
	
      } else {  /* KS_TYPE */
	int source_type, t0, saveflag_s;
	char descrp[MAXDESCRP];
	char savefile_s[MAXFILENAME];

	IF_OK status += ask_starting_ksprop( stdin, prompt, 
					     &param.startflag_ks[i],
					     param.startfile_ks[i]);
	
	IF_OK status += 
	  ask_output_ks_quark_source_file( stdin, prompt, &saveflag_s,
					   &source_type, &t0, descrp,
					   savefile_s );
	
	/* We could generate either a staggered extended source or a
	   naive (Dirac) extended source */
	IF_OK {
	  if(source_type == DIRAC_FIELD_FILE){
	    init_wqs(&param.dst_wqs[i]);
	    param.dst_wqs[i].type = source_type;
	    param.dst_wqs[i].t0   = t0;
	    param.dst_wqs[i].flag = saveflag_s;
	    strcpy(param.dst_wqs[i].descrp, descrp);
	    strcpy(param.dst_wqs[i].source_file, savefile_s);
	    param.dst_type[i] = CLOVER_TYPE;

	    init_wqs(&param.snk_wqs[i]);
	    status += get_w_quark_sink(stdin, prompt, &param.snk_wqs[i]);
	    param.snk_wqs[i].t0   = t0;
	  }
	  else if(source_type == VECTOR_FIELD_FILE){
	    init_ksqs(&param.dst_ksqs[i]);
	    param.dst_ksqs[i].type = source_type;
	    param.dst_ksqs[i].t0   = t0;
	    param.dst_ksqs[i].flag = saveflag_s;
	    strcpy(param.dst_ksqs[i].descrp, descrp);
	    strcpy(param.dst_ksqs[i].source_file, savefile_s);
	    param.dst_type[i] = KS_TYPE;

	    init_ksqs(&param.snk_ksqs[i]);
	    status += get_ks_quark_sink(stdin, prompt, &param.snk_ksqs[i]);
	    param.snk_ksqs[i].t0   = t0;
	  } else {
	    printf("Unsupported output source type\n");
	    status++;
	  }
	}
      } /* else KS type */

      /* Get the sink gamma matrix */

      if(param.dst_type[i] == CLOVER_TYPE){
	IF_OK get_s(stdin, prompt, "sink_gamma", savebuf);
	IF_OK {
	  param.snk_gam[i] = gamma_index(savebuf);
	  if(param.snk_gam[i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",savebuf);
	    status ++;
	  }
	}
      } else {
	param.snk_gam[i] = -999; /* Illegal if we ever try to use it */
      }
    }

    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    
  broadcast_bytes((char *)&param,sizeof(param));
  if( param.stopflag != 0 )
    normal_exit(0);

  return 0;
}

