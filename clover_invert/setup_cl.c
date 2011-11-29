/******** setup_cl.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

   8/01/98 Provision for serial write/read of scratch files C.D.
   8/30/96 Added reload_parallel for gauge fields C.D.
   8/15/96 Added prompts and param member for variable scratch file name C.D. 
   8/15/96 Added unitarity checking C.D.
   8/10/96 Revised propagator IO prompts and param file names C.D. */

//  $Log: setup_cl.c,v $
//  Revision 1.17  2011/11/29 04:06:26  detar
//  Move regression tests to "test" directory.  Support discontinued in favor of clover_invert2
//
//  Revision 1.16  2009/06/04 16:37:08  detar
//  Make clover term persistent. Accommodate changes to generic_clover/make_clov2.c
//
//  Revision 1.15  2009/05/31 02:00:56  detar
//  Fix "continue" and NULL startlat_p bug in clover_info.c and setup*.c
//
//  Revision 1.14  2009/04/05 16:49:22  detar
//  Add fixed node geometry and wave function file source
//
//  Revision 1.13  2008/04/18 15:14:46  detar
//  Fix printed comment
//
//  Revision 1.12  2008/03/28 15:37:53  detar
//  Fix heavy-light code and update sample input
//
//  Revision 1.11  2007/10/09 21:10:05  detar
//  Support new wprop options
//
//  Revision 1.10  2007/06/01 22:57:25  detar
//  Upgrade prompts for lattice and prop file names
//
//  Revision 1.9  2007/05/23 02:51:09  detar
//  Support sink wave functions in io_source_w
//  Move quark_source to generic_wilson.h
//
//  Revision 1.8  2007/05/21 04:38:11  detar
//  Reorganize spectrum computation, add QOP support, systematize inverter selection, add relative residue, add new source options
//
//  Revision 1.7  2006/11/07 05:28:10  detar
//  Train error files
//
//  Revision 1.6  2006/11/07 02:30:53  detar
//  Fix some omissions to complete the previous update.
//
//  Revision 1.5  2006/11/04 23:50:03  detar
//  Add file pointer to io_helpers utilities.
//
//  Revision 1.4  2006/02/20 23:47:19  detar
//  Update to version 7 and add hopping parameter inversion option
//


#include "cl_inv_includes.h"
#include <string.h>
int initial_set();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup_cl()   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, par_buf.iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* Create clover structure */
  gen_clov = create_clov();

  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 Wilson valence fermions\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    time_stamp("start");
    
    status = get_prompt(stdin,  &prompt );
    
    IF_OK status += get_i(stdin,prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin,prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin,prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin,prompt,"nt", &par_buf.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   par_buf.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   par_buf.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );
    
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
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = par_buf.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = par_buf.ionode_geometry[i];
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
  
  int status,status2;
  char savebuf[128];
  char save_w[128];
  int i, source_type;
  Real source_r0 = 0;
  int source_loc[4] = { 0,0,0,0 };
  int source_iters = 0;
  char source_file[MAXFILENAME] = "";
  char request_buf[MAX_SPECTRUM_REQUEST];

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;
    
    /* Number of kappas */
    IF_OK status += get_i(stdin,prompt,"number_of_kappas", &par_buf.num_kap );
    if( par_buf.num_kap>MAX_KAP ){
      printf("num_kap = %d must be <= %d!\n", par_buf.num_kap, MAX_KAP);
      status++;
    }
    
    /* To be safe initialize the following to zero */
    for(i=0;i<MAX_KAP;i++){
      kap[i] = 0.0;
      resid[i] = 0.0;
      relresid[i] = 0.0;
    }
    
    for(i=0;i<par_buf.num_kap;i++){
      IF_OK status += get_f(stdin, prompt,"kappa", &par_buf.kap[i] );
    }
    
    /* Clover coefficient, u0 */
    IF_OK status += get_f(stdin, prompt,"clov_c", &par_buf.clov_c );
    IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );
    
    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", &par_buf.niter );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", &par_buf.nrestart );
    
    /* error for propagator conjugate gradient */
    for(i=0;i<par_buf.num_kap;i++){
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid[i] );
      IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			    &par_buf.relresid[i] );
    }
    
    /* request list for spectral measurments */
    /* prepend and append a comma for ease in parsing */
    IF_OK status += get_s(stdin, prompt,"spectrum_request", request_buf );
    IF_OK strcpy(par_buf.spectrum_request,",");
    IF_OK strcat(par_buf.spectrum_request,request_buf);
    IF_OK strcat(par_buf.spectrum_request,",");
    
    /* Get source type */
    init_qs(&par_buf.wqs);
    IF_OK status += ask_w_quark_source(stdin,prompt,&source_type,
				       par_buf.wqs.descrp);
    IF_OK par_buf.wqs.type  = source_type;

    IF_OK {
      if ( source_type == GAUSSIAN ){
	IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
	/* width: psi=exp(-(r/r0)^2) */
	IF_OK status += get_f(stdin, prompt,"r0", &source_r0 );
      }
      else if ( source_type == POINT ){
	IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      }
      else if ( source_type == COVARIANT_GAUSSIAN ){
	IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
	IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
	IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);
      }
      else if ( source_type == COMPLEX_FIELD_FILE ){
	IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
	IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      }
      else if ( source_type == COMPLEX_FIELD_FM_FILE ){
	IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
	IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      }
      else if ( source_type == DIRAC_FIELD_FILE ){
	IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
	IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      }
      else if ( source_type == DIRAC_FIELD_FM_FILE ){
	IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
	IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      }
      else if ( source_type == WAVEFUNCTION_FILE ){
	IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
	IF_OK status += get_s(stdin, prompt, "load_source", source_file);
	IF_OK status += get_f(stdin, prompt, "a", &par_buf.wqs.a);
	IF_OK status += get_vi(stdin, prompt, "momentum", par_buf.wqs.mom, 3);
      }
      else {
	printf("Source type not supported in this application\n");
	status++;
      }
    }

    par_buf.wqs.r0    = source_r0;
    par_buf.wqs.x0    = source_loc[0];
    par_buf.wqs.y0    = source_loc[1];
    par_buf.wqs.z0    = source_loc[2];
    par_buf.wqs.t0    = source_loc[3];
    par_buf.wqs.iters = source_iters;
    strcpy(par_buf.wqs.source_file,source_file);
    
    /* Additional parameters for spectrum_multimom */
    if(strstr(par_buf.spectrum_request,",sink_smear,") != NULL){
      IF_OK status += get_f(stdin, prompt,"sink_r0",
			    &par_buf.sink_r0 );
    }

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
	par_buf.startfile );

    IF_OK if (prompt==1) 
      printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
    IF_OK scanf("%s",savebuf);
    IF_OK printf("%s\n",savebuf);
    IF_OK {
      if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
	par_buf.fixflag = COULOMB_GAUGE_FIX;
      }
      else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
	par_buf.fixflag = NO_GAUGE_FIX;
      }
      else{
	printf("error in input: fixing_command is invalid\n"); status++;
      }
    }
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
			     par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				  par_buf.stringLFN );
    
    /* find out starting propagator */
    IF_OK for(i=0;i<par_buf.num_kap;i++)
      status += ask_starting_wprop( stdin, prompt,&par_buf.startflag_w[i],
			par_buf.startfile_w[i]);
    
    /* what to do with computed propagator */
    IF_OK for(i=0;i<par_buf.num_kap;i++)
      status += ask_ending_wprop( stdin, prompt,&par_buf.saveflag_w[i],
		      par_buf.savefile_w[i]);
    
    IF_OK if(prompt==1) 
      printf("propagator scratch file:\n enter 'serial_scratch_wprop', 'parallel_scratch_wprop' or 'multidump_scratch_wprop'\n");
    IF_OK status2=scanf("%s",save_w);
    IF_OK printf("%s ",save_w);
    IF_OK 
      {
	if(strcmp("serial_scratch_wprop",save_w) == 0 )
	  par_buf.scratchflag = SAVE_SERIAL_SCIDAC;
	else if(strcmp("parallel_scratch_wprop",save_w) == 0 )
	  par_buf.scratchflag = SAVE_PARALLEL_SCIDAC;
	else if(strcmp("multifile_scratch_wprop",save_w) == 0 )
	  par_buf.scratchflag = SAVE_MULTIFILE_SCIDAC;
	else
	  {
	    printf("error in input: %s is not a scratch file command\n",save_w);
	    status++;
	  }
	IF_OK
	  {
	    /*read name of file and load it */
	    if(prompt==1)printf("enter name of scratch file stem for props\n");
	    status2=scanf("%s",par_buf.scratchstem_w);
	    if(status2 !=1) {
	      printf("error in input: scratch file stem name\n"); status++;
	    }
	    printf("%s\n",par_buf.scratchstem_w);
	  }
      }
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if( par_buf.stopflag != 0 )
    normal_exit(0);

  startflag = par_buf.startflag;
  fixflag = par_buf.fixflag;
  saveflag = par_buf.saveflag;
  for(i=0;i<par_buf.num_kap;i++){
    startflag_w[i] = par_buf.startflag_w[i];
    saveflag_w[i] = par_buf.saveflag_w[i];
  }
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  num_kap = par_buf.num_kap;
  clov_c = par_buf.clov_c;
  u0 = par_buf.u0;
  strcpy(spectrum_request,par_buf.spectrum_request);
  sink_r0 = par_buf.sink_r0;
  for(i=0;i<par_buf.num_kap;i++){
    kap[i] = par_buf.kap[i];
    resid[i] = par_buf.resid[i];
    relresid[i] = par_buf.relresid[i];
  }
  wqs = par_buf.wqs;
  wqs.type = par_buf.wqs.type;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  for(i=0;i<par_buf.num_kap;i++){
    strcpy(startfile_w[i],par_buf.startfile_w[i]);
    strcpy(savefile_w[i],par_buf.savefile_w[i]);
  }
  strcpy(scratchstem_w,par_buf.scratchstem_w);
  scratchflag = par_buf.scratchflag;
  
  /* Do whatever is needed to get lattice */
  if( startflag != CONTINUE ){
    startlat_p = reload_lattice( startflag, startfile );
    invalidate_this_clov(gen_clov);
  }
  return(0);
}

