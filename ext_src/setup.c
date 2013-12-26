/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

* 5/30/07 Created from setup_cl.c */

//  $Log: setup.c,v $
//  Revision 1.6  2013/12/26 16:02:19  detar
//  Fix error handling abort.
//
//  Revision 1.5  2011/11/29 22:10:38  detar
//  New KS4 type for extended dirac propagators.  New source structure.
//
//  Revision 1.4  2009/06/11 16:24:22  detar
//  Allow writing multiple sources from the same file.  Changes parameter inputs.
//
//  Revision 1.3  2009/05/31 02:00:57  detar
//  Fix "continue" and NULL startlat_p bug in clover_info.c and setup*.c
//
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
#include <string.h>

static int initial_set(void);
static void broadcast_heap_params(void);

#include "params.h"

int setup(void)   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;
  /* initialize the node random number generator */
  initialize_prn( &node_prn, param.iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();

  return(prompt);
}


/* SETUP ROUTINES */
static int initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 clover, staggered and naive valence fermions\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");
    
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
  int i,j;

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

    /* Quark parameters */
    IF_OK for(i = 0; i < param.num_qk; i++){

      /* Input propagator type */
      IF_OK status += get_s(stdin, prompt,"quark_type", savebuf );
      IF_OK {
	if(strcmp(savebuf,"clover") == 0)param.qk_type[i] = CLOVER_TYPE;
	else if(strcmp(savebuf,"KS") == 0)param.qk_type[i] = KS_TYPE;
	else {
	  printf("Unknown quark type %s\n",savebuf);
	  status++;
	}
      }

      /* Output propagator type */
      /* The KS4_TYPE is designed for extended naive propagators that
         follow the spin conventions of the clover_invert2 code. */
      IF_OK status += get_s(stdin, prompt,"output_type", savebuf );
      IF_OK {
	if(strcmp(savebuf,"clover") == 0)param.dst_type[i] = CLOVER_TYPE;
	else if(strcmp(savebuf,"KS") == 0)param.dst_type[i] = KS_TYPE;
	else if(strcmp(savebuf,"KS4") == 0)param.dst_type[i] = KS4_TYPE;
	else {
	  printf("Unknown output type %s\n",savebuf);
	  status++;
	}
      }

      /* Input Dirac propagator */
      if(param.qk_type[i] == CLOVER_TYPE){
	/* Name of starting Dirac propagator file */
	IF_OK status += ask_starting_wprop( stdin, prompt, 
					    &param.startflag_w[i],
					    param.startfile_w[i]);
	IF_OK status += get_i(stdin,prompt,"ncolor", &param.ncolor[i] );
	
	/* Dirac sink operator */
	IF_OK {
	  if(param.dst_type[i] == CLOVER_TYPE ||
	     param.dst_type[i] == KS4_TYPE){
	    IF_OK init_qss_op(&param.snk_qs_op[i]);
	    IF_OK status += get_wv_field_op( stdin, prompt, &param.snk_qs_op[i]);
	    
	  } else {
	    printf("Unsupported output source type\n");
	    status++;
	  }
	}
	
      } else {  /* Input KS propagator */
	
        /* Name of starting KS propagator file */
	IF_OK status += ask_starting_ksprop( stdin, prompt, 
					     &param.startflag_ks[i],
					     param.startfile_ks[i]);

	IF_OK status += get_i(stdin,prompt,"ncolor", &param.ncolor[i] );
	
	/* We could generate either a staggered extended source or a
	   naive (Dirac) extended source */
	IF_OK {
	  if(param.dst_type[i] == CLOVER_TYPE ||
	     param.dst_type[i] == KS4_TYPE){
	    IF_OK init_qss_op(&param.snk_qs_op[i]);
	    IF_OK status += get_wv_field_op( stdin, prompt, &param.snk_qs_op[i]);
	  }
	  else if(param.dst_type[i] == KS_TYPE){
	    IF_OK init_qss_op(&param.snk_qs_op[i]);
	    IF_OK status += get_v_field_op( stdin, prompt, &param.snk_qs_op[i]);
	  } else {
	    printf("Unsupported output source type\n");
	    status++;
	  }
	}
      } /* else KS type */

      /* Get the sink gamma matrix */

      IF_OK get_s(stdin, prompt, "sink_gamma", savebuf);

      IF_OK {
	if(param.dst_type[i] == CLOVER_TYPE ||
	   param.dst_type[i] == KS4_TYPE){
	  param.snk_gam[i] = gamma_index(savebuf);
	  if(param.snk_gam[i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",savebuf);
	    status ++;
	  }
	} else {
	  /* For staggered quarks we use spin_taste operators with no
	     displacement in time */
	  param.snk_gam[i] = spin_taste_index(savebuf);
	  if(param.snk_gam[i] < 0){
	    printf("\n%s is not a valid spin_taste label.\n",savebuf);
	    status ++;
	  }
	}
      }

      /* FT and KS phases are computed with x,y,z,t relative to r_offset */
      IF_OK {
	int r[4];
	status += get_vi(stdin,prompt, "r_offset", r, 4);
	param.r_offset[i][0] = r[0];
	param.r_offset[i][1] = r[1];
	param.r_offset[i][2] = r[2];
	param.r_offset[i][3] = r[3];
      }

      /* Set the operator coordinate offset (for some operators) */
      
      set_qss_op_offset(&param.snk_qs_op[i], &param.r_offset[i][0]);

      /* Get the time slice and output source file */

      IF_OK status += get_i(stdin,prompt,"number_of_time_slices", 
			    &param.num_t0[i] );

      IF_OK for(j = 0; j < param.num_t0[i]; j++){
	if(param.qk_type[i] == CLOVER_TYPE){
	  int save_type, t0, saveflag_s;
	  char descrp[MAXDESCRP];
	  char savefile_s[MAXFILENAME];
	  
	  IF_OK status += 
	    ask_output_quark_source_file( stdin, prompt, &saveflag_s,
					  &save_type, &t0, descrp,
					  savefile_s );
	
	  IF_OK {
	    if(t0 >= nt){
	      printf("Source time slice must be less than nt = %d\n", nt);
	      status++;
	    }
	  }

	  IF_OK {
	    if(save_type == DIRAC_FIELD_FILE){
	      if(param.dst_type[i] != CLOVER_TYPE &&
		 param.dst_type[i] != KS4_TYPE){
		printf("This file type requires a Dirac field.\n");
		status++;
	      }
	      /* Set values according to the input propagator */
	      init_qs(&param.dst_qs[i][j]);
	      param.dst_qs[i][j].savetype = save_type;
	      param.dst_qs[i][j].subset = FULL;
	      param.dst_qs[i][j].t0   = t0;
	      param.dst_qs[i][j].saveflag = saveflag_s;
	      strcpy(param.dst_qs[i][j].descrp, descrp);
	      strcpy(param.dst_qs[i][j].save_file, savefile_s);

	    } else {
	      printf("Unsupported output source type\n");
	      status++;
	    }
	  }
	} else {  /* KS_TYPE */
	  int save_type, t0, saveflag_s;
	  char descrp[MAXDESCRP];
	  char savefile_s[MAXFILENAME];
	  
	  IF_OK status += 
	    ask_output_quark_source_file( stdin, prompt, &saveflag_s,
					  &save_type, &t0, descrp,
					  savefile_s );

	  IF_OK {
	    if(t0 >= nt){
	      printf("Source time slice must be less than nt = %d\n", nt);
	      status++;
	    }
	  }
	  
	  /* We could generate either a staggered extended source or a
	     naive (Dirac) extended source */
	  IF_OK {
	    init_qs(&param.dst_qs[i][j]);
	    /* Set values according to the iput propagator */
	    param.dst_qs[i][j].savetype = save_type;
	    param.dst_qs[i][j].subset = FULL;
	    param.dst_qs[i][j].t0   = t0;
	    param.dst_qs[i][j].saveflag = saveflag_s;
	    strcpy(param.dst_qs[i][j].descrp, descrp);
	    strcpy(param.dst_qs[i][j].save_file, savefile_s);

	    if(save_type == DIRAC_FIELD_FILE){
	      if(param.dst_type[i] != CLOVER_TYPE &&
		 param.dst_type[i] != KS4_TYPE){
		printf("This file type requires a Dirac field.\n");
		status++;
	      }
	    }
	    else if(save_type == VECTOR_FIELD_FILE){
	      if(param.dst_type[i] != KS_TYPE){
		printf("This file type requires a color vector field.\n");
		status++;
	      }
	    } else {
	      printf("Unsupported output source type\n");
	      status++;
	    }
	  } /* OK */
	} /* if KS_TYPE */
      } /* j time slices */
    } /* i quarks */

    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    
  broadcast_bytes((char *)&param,sizeof(param));
  if( param.stopflag != 0 )
    normal_exit(0);

  /* Broadcast parameter values kept on the heap */
  broadcast_heap_params();

  return 0;
}

/* Broadcast operator parameter values.  They are on the heap on node 0. */

static void broadcast_heap_params(void){
  int i, j;

  for(i = 0; i < param.num_qk; i++){
    for(j = 0; j < param.num_t0[i]; j++)
      broadcast_quark_source_sink_op_recursive(&param.dst_qs[i][j].op);
    broadcast_quark_source_sink_op_recursive(&param.snk_qs_op[i].op);
  }
}

