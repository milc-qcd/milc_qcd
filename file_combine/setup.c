/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

 * 3/25/12 CD Created  */

#include "file_combine_includes.h"
#include <string.h>

static int initial_set(void);

#include "params.h"

int setup(void)   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  if(prompt == 2)return prompt;
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();

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
    printf("SU3 linear combinations of sources and propagators\n");
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
  
  int status;
  int ifile;
  char *savebuf;
  char descrp[128];
  char myname[] = "readin";

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    /*------------------------------------------------------------*/
    /* File type                                                  */
    /*------------------------------------------------------------*/

    if(prompt == 1){
      printf("enter 'vector_field_files', 'dirac_field_files', ");
      printf("'vector_propagator_files', 'dirac_propagator_files'\n");
    }
    savebuf = get_next_tag(stdin, "file number and type", myname);
    printf("%s ",savebuf);
    if(savebuf == NULL)status++;

    IF_OK {
      if(strcmp("vector_field_files", savebuf)==0)
	param.file_type = VECTOR_FIELD_FILE;
      else if(strcmp("dirac_field_files", savebuf)==0)
	param.file_type = DIRAC_FIELD_FILE;
      else if(strcmp("vector_propagator_files", savebuf)==0)
	param.file_type = VECTOR_PROPAGATOR_FILE;
      else if(strcmp("dirac_propagator_files", savebuf)==0)
	param.file_type = DIRAC_PROPAGATOR_FILE;
      else{
	printf("Unrecognized file type\n");
	status++;
      }
    }

    /*------------------------------------------------------------*/
    /* Number of files                                            */
    /*------------------------------------------------------------*/

    IF_OK {
      if(scanf("%d", &param.nfile) != 1)status++;
      if(status > 0)printf("\n%s: bad number of file\n", myname);
      else printf("%d\n", param.nfile);
      if(param.nfile > MAX_FILES){
	printf("The number of files must be less than or equal to %d\n",
	       MAX_FILES);
	status++;
      }
    }

    /*------------------------------------------------------------*/
    /* Number of colors                                           */
    /*------------------------------------------------------------*/

    IF_OK get_i(stdin, prompt, "ncolor", &param.ncolor);
    param.nspin = 4;  /* Fixed for now */
 
    /*------------------------------------------------------------*/
    /* Time slice restriction                                     */
    /*------------------------------------------------------------*/

    IF_OK get_i(stdin, prompt, "t0", &param.t0);
    /* A negative value means to write the full field */
    if(param.t0 < 0)param.t0 = ALL_T_SLICES;

   /*------------------------------------------------------------*/
    /* Names of files                                            */
    /*------------------------------------------------------------*/

    IF_OK for(ifile = 0; ifile < param.nfile; ifile++){
      if(param.file_type == VECTOR_FIELD_FILE || param.file_type == DIRAC_FIELD_FILE){
	IF_OK status += ask_starting_source(stdin, prompt, &param.startflag[ifile], 
					    param.startfile[ifile]);
      }
      else if(param.file_type == VECTOR_PROPAGATOR_FILE){
	IF_OK status += ask_starting_ksprop( stdin, prompt, &param.startflag[ifile], 
			     param.startfile[ifile]);
      }
      else if(param.file_type == DIRAC_PROPAGATOR_FILE){
	IF_OK status += ask_starting_wprop( stdin, prompt, &param.startflag[ifile], 
			    param.startfile[ifile]);
      }
    }

    /*------------------------------------------------------------*/
    /* Coefficients                                               */
    /*------------------------------------------------------------*/

    IF_OK get_vf(stdin, prompt, "coeffs", param.coeff, param.nfile);

    /*------------------------------------------------------------*/
    /* Output file                                                */
    /*------------------------------------------------------------*/

    if(param.file_type == VECTOR_FIELD_FILE || param.file_type == DIRAC_FIELD_FILE){
      IF_OK status += 
	ask_output_quark_source_file(stdin, prompt, &param.saveflag,
				     &param.savetype, NULL,
				     descrp, param.savefile);
    }
    else if(param.file_type == VECTOR_PROPAGATOR_FILE){
      IF_OK status += ask_ending_ksprop( stdin, prompt, &param.saveflag,
					 param.savefile);
    }
    else if(param.file_type == DIRAC_PROPAGATOR_FILE){
      IF_OK status += ask_ending_wprop( stdin, prompt, &param.saveflag,
					param.savefile);
    }

    /*------------------------------------------------------------*/
    /* End of input fields */
    /*------------------------------------------------------------*/

    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    
  broadcast_bytes((char *)&param,sizeof(param));
  if( param.stopflag != 0 )
    normal_exit(0);

  return 0;
}


