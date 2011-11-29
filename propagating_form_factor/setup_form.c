/******** setup_form.c *********/
/* MIMD version 6 */
/*  set tabstop=2   for easy reading of this file */
/* $Header: /lqcdproj/detar/cvsroot/milc_qcd/propagating_form_factor/setup_form.c,v 1.8 2011/11/29 18:49:04 detar Exp $  ***/
/* MIMD code version 4 */

#include "prop_form_includes.h"
#include <string.h>

#define IF_OK if(status==0)

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;
int prompt;

/*
 *  Driver routines for the setup functions.
 *  The function is called in "main()"
 *
 */
int  setup_h()   {
int initial_set();

	/* print banner, get volume */
    prompt=initial_set();
	/* Initialize the layout functions, which decide where sites live */
    setup_layout();
	/* allocate space for lattice, set up coordinate fields */
    make_lattice();
	/* set up nearest neighbor gathers */
    make_nn_gathers();
	/* start interrupt handler for field_pointer() requests */

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set()
{
  int prompt,status;
  /* On node zero, read lattice size and send to others */
  if(mynode()==0)
  {
    
    /* print banner */
    printf("SU3 Wilson valence fermions;  PROPAGATING form factor code \n");
    printf("Heavy--> Heavy and Heavy-->light form factors\n"); 
    printf("MIMD version 6\n");
#ifdef CLOVER  
    printf("Using CLOVER fermions with the BI-CG and hopping algorithms\n") ; 
#else
    printf("Using WILSON fermions with the minimal residual algorithm\n") ; 
#endif
    
    
    
    printf("Machine = %s, with %d CPUs\n",machine_type(),numnodes());
    time_stamp("start");
    
    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
    
    if(par_buf.nt%2 !=0) 
      {
	printf("nt must be even!! \n"); 
	terminate(1);
      }
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

  return(prompt);
}

/* 
   Node 0 reads in all parameters controlling the
   code, except the dimensions of the lattice (which are
   read in using the function int initial_set()).

   The parameters are then send to all the other nodes.


   argument "prompt" is 1 if prompts are to be given for input	

*/

/* This test is for node 0 information only */
int file_exists(char *file)
{
  FILE *dummy_fp;
  int non_null_fp;

  /* Test whether file exists and is openable result */
  dummy_fp = fopen(file,"rb");
  non_null_fp = (dummy_fp != NULL);
  if(non_null_fp)fclose(dummy_fp);
  return non_null_fp;
}

int ask_ending_correlator( int prompt, int *flag, char *id, 
			  char *filename ){
    char savebuf[256];
    int status;

    if (prompt==1) printf(
        "'forget' %s correlator at end,  'save_ascii', or 'save_binary'\n",
			  id);
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_ending_correlator: ERROR IN INPUT: expecting save command for %s correlator\n",id);
        return(1);
    }
    printf("%s ",savebuf);
    if(strcmp("save_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_binary",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("forget",savebuf) == 0 ) {
        *flag=FORGET;
	printf("\n");
    }
    else {
      printf("ask_ending_correlator: ERROR IN INPUT: %s is not a save correlator command\n",savebuf);
      return(1);
    }


    if(prompt==0){
      status=scanf("%s ",savebuf);
      if(status!=1)
	{
	  printf("ask_ending_correlator: ERROR IN INPUT: Expecting correlator ID %s\n",id);
	  return(1);
	}
      printf("%s ",savebuf);
      if(strcmp(savebuf,id) != 0)
	{
	  printf("ask_ending_correlator: ERROR IN INPUT: Expecting %s but found %s\n",id,savebuf);
	  return(1);
	}
    }

    if( *flag != FORGET ){
        if(prompt==1)printf("enter filename for %s correlator\n",id);
        status=scanf("%s",filename);
        if(status !=1){
    	    printf("ask_ending_correlator: ERROR IN INPUT: expecting filename\n"); return(1);
        }
	printf("%s\n",filename);
    }
    return(0);
}

int ask_quark_inverter( int prompt, int *inverter_type)
{
  char savebuf[256];
  int status;

  if (prompt==1)
    printf("enter 'hop', or 'bicg' for heavy quark inverter type\n");
  status = scanf("%s",savebuf);
  if(status !=1) {
    printf("ask_quark_inverter: ERROR IN INPUT: inverter type command\n");
    return 1;
  }
  if(strcmp("hopilu",savebuf) == 0 ){
    *inverter_type = HOPILU;
  }
  else if(strcmp("bicgilu",savebuf) == 0 ) {
    *inverter_type = BICGILU;
  }
  else{
    printf("ask_quark_inverter: ERROR IN INPUT: inverter command is invalid\n"); 
    return 1;
  }

  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_inverter */

int readin(int prompt) 
{
  register int i,ii;  /******loop variable*****/
  int seen ; 
  int ikappa , imom ; 
  int status;
  int ismear ;
  char descrp[30];

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0)
    {
      
      printf("\n\n");
      status=0;
      
      IF_OK status += get_i(stdin, prompt,"verbose_flag",&par_buf.verbose_flag);
#ifdef CLOVER  
      /* Clover coefficient, u0 */
      IF_OK status += get_f(stdin, prompt,"clov_c",&par_buf.clov_c);
      IF_OK status += get_f(stdin, prompt,"u0",&par_buf.u0);
#endif

      /* find out what kind of starting lattice to use */
      IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
				 par_buf.startfile );

      /* get number of values of SPECTATOR kappa to be run */
      IF_OK status += get_i(stdin, prompt,"nkap_spectator",&par_buf.no_spectator);
      if(par_buf.no_spectator  >MAX_KAPPA) 
	{
	  printf("no_spectator = %d cannot be larger than MAX_KAPPA =%d! \n",
		 par_buf.no_spectator,MAX_KAPPA );
	  status++;
	}
      
      /* maximum no. of spectator conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter_spectator );
      
      /* maximum no. of spectator conjugate gradient restarts */
      IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart_spectator );
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid_spectator );
      
      /*** Load in the spectator kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_spectator  ;i++) 
	  { 
	    IF_OK status += get_f(stdin, prompt,"kappa_spectator",
				  &par_buf.kappa_spectator[i]);
	  }
      }
      
      /* Get source type */
      IF_OK status += ask_w_quark_source( stdin, prompt, &wallflag,descrp);
      
      /* Load in the source widths and names of the files for the
       *   spectator propagators */
      
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK if (prompt==1) 
	printf("enter width(s) r0 as in: source=exp(-(r/r0)^2)\n");

      for(i=0;i<par_buf.no_spectator;i++){
	IF_OK status += get_f(stdin, prompt,"r0", &par_buf.wqs_spectator[i].r0 );
	/* (Same source type for each spectator) */
	IF_OK par_buf.wqs_spectator[i].type = wallflag;
	IF_OK strcpy(par_buf.wqs_spectator[i].descrp,descrp);
	/* (Hardwired source location for each spectator) */
	IF_OK {
	  par_buf.wqs_spectator[i].x0 = source_loc[0];
	  par_buf.wqs_spectator[i].y0 = source_loc[1];
	  par_buf.wqs_spectator[i].z0 = source_loc[2];
	  par_buf.wqs_spectator[i].t0 = source_loc[3];
	}
      }
      
      /***  load the names of the  spectator quarks  ***/
      IF_OK {
	for(i=0;i< par_buf.no_spectator ;i++)
	  {
	    IF_OK status += ask_starting_wprop(stdin, prompt,
					      &par_buf.startflag_spectator[i],
					      par_buf.qfile_spectator[i]);
	  } /****** end of loop over spectator kappa values ******/
      }
      
      /* get number of values of LIGHT zonked kappa to be run */
      IF_OK status += get_i(stdin, prompt,"nkap_light_zonked",
			    &par_buf.no_zonked_light);
      if(par_buf.no_zonked_light  >MAX_ZONKED_LIGHT) 
	{
	  printf("no_light_zonked_light = %d cannot be larger than MAX_KAPPA =%d!\n",
		 par_buf.no_zonked_light,MAX_ZONKED_LIGHT );
	  status++;
	}
      
      /* maximum no. of zonked conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter_zonked_light );
      
      /* maximum no. of zonked conjugate gradient restarts */
      IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart_zonked_light );
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid_zonked_light );
      
      /*** Load in the zonked light kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_light  ;i++) 
	{ 
	  IF_OK status+= get_f(stdin, prompt,"kappa_zonked_light",
			       &par_buf.kappa_zonked_light[i]);
	  /* (Same wallflag for each zonked_light) */
	  IF_OK par_buf.wqs_zonked_light[i].type = wallflag;
	  IF_OK strcpy(par_buf.wqs_zonked_light[i].descrp,descrp);
	  /* (Hardwired source location for each zonked_light) */
	  IF_OK {
	    par_buf.wqs_zonked_light[i].x0 = source_loc[0];
	    par_buf.wqs_zonked_light[i].y0 = source_loc[1];
	    par_buf.wqs_zonked_light[i].z0 = source_loc[2];
	    par_buf.wqs_zonked_light[i].t0 = source_loc[3];
	  }
	}
      }

      /* Get source type */
      IF_OK status += ask_w_quark_source( stdin, prompt, &wallflag, descrp);

      /* Load in the source widths and names of the files for the
       *   zonked propagators */
      
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK if (prompt==1) 
	printf("enter width(s) r0 as in: source=exp(-(r/r0)^2)\n");

      for(i=0;i<par_buf.no_zonked_light;i++){
	IF_OK status += get_f(stdin, prompt,"r0", &par_buf.wqs_zonked_light[i].r0 );
      }
      
      /***  load the names of the  LIGHT zonked quarks  ***/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_light ;i++)
	  {
	    IF_OK status += ask_starting_wprop(stdin, prompt,
				       &par_buf.startflag_zonked_light[i],
				       par_buf.qfile_zonked_light[i]);

	    if( par_buf.startflag_zonked_light[i] == FRESH )
	      sprintf(par_buf.qfile_zonked_light[i],
		      "fresh_light_zonked_k%g",par_buf.kappa_zonked_light[i]); 

	  }
      }

      /*** load the save flag for zonked light sink-smeared quarks ***/
      IF_OK status += ask_ending_wprop(stdin, prompt,
				       &par_buf.saveflag_zonked_light_ssink,
				       par_buf.qfile_suffix_zonked_light);
	
      IF_OK if(par_buf.saveflag_zonked_light_ssink == FORGET)
	sprintf(par_buf.qfile_suffix_zonked_light,"_fresh");
      
      /* get number of values of HEAVY zonked kappa to be run */
      IF_OK status += get_i(stdin, prompt,"nkap_heavy_zonked",
			    &par_buf.no_zonked_heavy);
      if(par_buf.no_zonked_heavy  >MAX_ZONKED_HEAVY) 
	{
	  printf("no_heavy_zonked = %d cannot be larger than MAX_KAPPA =%d!\n",
		 par_buf.no_zonked_heavy,MAX_ZONKED_HEAVY );
	  status++;
	}
      
      /* maximum no. of heqvy zonked conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter_zonked_heavy );
      
      /* maximum no. of heqvy zonked conjugate gradient restarts */
      IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart_zonked_heavy );
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid_zonked_heavy );
      
      /*** Load in the heavy zonked kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_heavy  ;i++) 
	  { 
	    IF_OK status += get_f(stdin, prompt,"kappa_zonked_heavy",
				  &par_buf.kappa_zonked_heavy[i]);
	  }
      }
      
      /*** Get inversion algorithm for heavy zonked quarks ****/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_heavy  ;i++) 
	  { 
	    IF_OK status += ask_quark_inverter(prompt,&par_buf.inverter_type_zonked_heavy[i]);
	  }
      }
      
      /* Get source type for heavy_zonked quark */
      IF_OK status += ask_w_quark_source( stdin, prompt, &wallflag, descrp);

      /* Load in the source widths for the heavy zonked propagators */
      
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK if (prompt==1) 
	printf("enter width(s) r0 as in: source=exp(-(r/r0)^2)\n");

      for(i=0;i<par_buf.no_zonked_heavy;i++){
	IF_OK status += get_f(stdin, prompt,"r0", &par_buf.wqs_zonked_heavy[i].r0 );
	/* (Same source type for each spectator) */
	IF_OK par_buf.wqs_zonked_heavy[i].type = wallflag;
	IF_OK strcpy(par_buf.wqs_zonked_heavy[i].descrp,descrp);
	/* (Hardwired source location for each zonked_heavy) */
	IF_OK {
	  par_buf.wqs_zonked_heavy[i].x0 = source_loc[0];
	  par_buf.wqs_zonked_heavy[i].y0 = source_loc[1];
	  par_buf.wqs_zonked_heavy[i].z0 = source_loc[2];
	  par_buf.wqs_zonked_heavy[i].t0 = source_loc[3];
	}
      }

      /***  load the names of the HEAVY zonked quarks  ***/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_heavy ;i++)
	  {
	    IF_OK status += ask_ending_wprop(stdin, prompt,
				    &par_buf.saveflag_zonked_heavy[i],
				    par_buf.qfile_zonked_heavy[i]);
	  }
      }

      /*** load the save flag for zonked heavy sink-smeared quarks ***/
      IF_OK status += ask_ending_wprop(stdin, prompt,
			       &par_buf.saveflag_zonked_heavy_ssink,
			       par_buf.qfile_suffix_zonked_heavy);
	
      IF_OK if(par_buf.saveflag_zonked_heavy_ssink == FORGET)
	sprintf(par_buf.qfile_suffix_zonked_heavy,"_fresh");
      
      /* get number of values of SEQUENTIAL kappas to be run */
      IF_OK status += get_i(stdin, prompt,"nkap_sequential",&par_buf.no_sequential);
      if(par_buf.no_sequential >MAX_KAPPA) 
	{
	  printf("no_sequential = %d cannot be larger than MAX_KAPPA =%d! \n",
		 par_buf.no_sequential,MAX_KAPPA );
	  status++;
	}

      /*** Load in the sequential kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_sequential  ;i++) 
	  { 
	    IF_OK status += get_f(stdin, prompt,"kappa_seq",
				  &par_buf.kappa_sequential[i]);
	  }	
      }

      /*** the sequential kappa values should be a subset of the zonked
	kappa values ****/
      
      
      IF_OK {
	for(i=0;i< par_buf.no_sequential  ;i++) 
	  { 
	    seen = 0 ; 
	    for(ii=0;ii< par_buf.no_zonked_heavy  ;ii++) 
	      { 
		if( par_buf.kappa_sequential[i] == 
		   par_buf.kappa_zonked_heavy[ii] ) seen  = 1 ;
	      }
	    if( seen == 0 )
	      {
		printf("ERROR sequential kappa value = %g is not in the kappa_zonked_heavy SET\n", 
		       par_buf.kappa_sequential[i]  ); 
		status++;
	      }
	    
	  }
      }
    
      /*** Get inversion algorithm for sequential quarks ****/
      IF_OK {
	for(i=0;i< par_buf.no_sequential  ;i++) 
	  { 
	    IF_OK ask_quark_inverter(prompt,&par_buf.inverter_type_sequential[i]);
	  }
      }
      
      /* get the position in time of the final meson ****/
      IF_OK status += get_i(stdin, prompt,"final_time",&par_buf.tf);
      if( par_buf.tf < 0 || par_buf.tf >= nt )
      {
	printf("ERROR: tf = %d is outside the length of the lattice\n",par_buf.tf ); 
        ++status ;
      }
      
      /* get i/o info */

      /*** load in the B MESON insertion momenta   ****/
      
      IF_OK status += 
	load_momentum(prompt,"p",&par_buf.no_p_values,par_buf.p_momstore , MAXPMOM) ;

      /*** load in the Vertex momenta   ****/
      
      IF_OK status += 
	load_momentum(prompt,"q",&par_buf.no_q_values,par_buf. q_momstore , MAXMOM) ;
      /*** load in the 2-pt MESON momenta   ****/

      IF_OK status += 
	load_momentum(prompt,"k",&par_buf.no_k_values,par_buf. k_momstore , MAXMOM) ;

      /** 
	load in the names of the file to read the sequential smearing functions from 
	***/
      IF_OK {
	for(ismear = 0 ; ismear < par_buf.no_p_values ; ++ismear)
	  {
	    IF_OK status += get_s(stdin, prompt,"seq_smear_func",
				   par_buf.seq_smear_file[ismear]);
	  }
      }
      
      
      /** 
	get output file name and type for the  heavy --> heavy form factors
	***/
      
      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HH3,
			      "HH3",par_buf.filename_HH3);

      IF_OK{
	if(file_exists(par_buf.filename_HH3)){
	  printf("SKIPPING heavy heavy form factor.  File %s exists.\n",
		     par_buf.filename_HH3);
	  par_buf.saveflag_HH3 = FORGET;
	}
      }

	/** 
	  get output file name and type for the heavy --> light form factors
	  ***/
      
      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HL3,
			      "HL3",par_buf.filename_HL3);

      IF_OK{
	if(file_exists(par_buf.filename_HL3)){
	  printf("SKIPPING heavy light form factor.  File %s exists.\n",
		     par_buf.filename_HL3);
	  par_buf.saveflag_HL3 = FORGET;
	}
      }


	/** 
	  get names and types of the files for various 2 pt functions
	  ***/

      /* Heavy-heavy two-pt */
      
      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HH2_GL,
			      "HH2_GL",par_buf.filename_HH2_GL);

      IF_OK {
	if(file_exists(par_buf.filename_HH2_GL)){
	  printf("SKIPPING HH2_GL correlator.  File %s exists.\n",
		     par_buf.filename_HH2_GL);
	  par_buf.saveflag_HH2_GL = FORGET;
	}
	if(par_buf.saveflag_HH2_GL == SAVE_SERIAL)
	  {
	    printf("WARNING: No binary format available for the HH2_GL correlator\n");
	    par_buf.saveflag_HH2_GL = SAVE_ASCII;
	  }
      }


      /* Light-light two-pt */

      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_LL2_GG,
			      "LL2_GG",par_buf.filename_LL2_GG);

      IF_OK {
	if(file_exists(par_buf.filename_LL2_GG)){
	  printf("SKIPPING LL2_GG correlator.  File %s exists.\n",
		     par_buf.filename_LL2_GG);
	  par_buf.saveflag_LL2_GG = FORGET;
	}
      }


      /* Heavy-light two-pt (smeared source, smeared sink) */

      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HL2_GG,
			      "HL2_GG",par_buf.filename_HL2_GG);

      IF_OK {
	if(file_exists(par_buf.filename_HL2_GG)){
	  printf("SKIPPING HL2_GG correlator.  File %s exists.\n",
		     par_buf.filename_HL2_GG);
	  par_buf.saveflag_HL2_GG = FORGET;
	}
      }

      /* Heavy-light two-pt (smeared source, relative wf sink) */

      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HL2_GE,
			      "HL2_GE",par_buf.filename_HL2_GE);

      IF_OK {
	if(file_exists(par_buf.filename_HL2_GE)){
	  printf("SKIPPING HL2_GE correlator.  File %s exists.\n",
		     par_buf.filename_HL2_GE);
	  par_buf.saveflag_HL2_GE = FORGET;
	}
      }

      /* Heavy-light two-pt (smeared source, local sink) */

      IF_OK status += 
	ask_ending_correlator(prompt, &par_buf.saveflag_HL2_GL,
			      "HL2_GL",par_buf.filename_HL2_GL);

      IF_OK {
	if(file_exists(par_buf.filename_HL2_GL)){
	  printf("SKIPPING HL2_GL correlator.  File %s exists.\n",
		     par_buf.filename_HL2_GL);
	 par_buf.saveflag_HL2_GL = FORGET; 
	}
      }


      if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */
  
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  

  if( par_buf.stopflag != 0 )
    normal_exit(0);



  /**** unpack the parameters *****/

  verbose_flag  = par_buf.verbose_flag;
  clov_c = par_buf.clov_c;
  u0 = par_buf.u0;
  startflag = par_buf.startflag  ;
  strcpy(startfile,par_buf.startfile) ; 
  
  no_spectator = par_buf.no_spectator    ;
  niter_spectator = par_buf.niter_spectator ;
  nrestart_spectator = par_buf.nrestart_spectator ;
  resid_spectator = par_buf.resid_spectator ;

  no_zonked_light  = par_buf.no_zonked_light   ;
  niter_zonked_light = par_buf.niter_zonked_light ;
  nrestart_zonked_light = par_buf.nrestart_zonked_light ;
  resid_zonked_light = par_buf.resid_zonked_light ;

  no_zonked_heavy = par_buf.no_zonked_heavy  ;
  niter_zonked_heavy = par_buf.niter_zonked_heavy ;
  nrestart_zonked_heavy = par_buf.nrestart_zonked_heavy ;
  resid_zonked_heavy = par_buf.resid_zonked_heavy ;

  no_sequential = par_buf.no_sequential    ;
  
  for(ikappa =0 ; ikappa < MAX_KAPPA ; ++ikappa )
    {
      kappa_spectator[ ikappa ]  =   par_buf.kappa_spectator[ ikappa ] ; 
      kappa_zonked_light[ ikappa ] = par_buf.kappa_zonked_light[ ikappa ] ;
      kappa_zonked_heavy[ ikappa ] = par_buf.kappa_zonked_heavy[ ikappa ] ;
      kappa_sequential[ ikappa ]  =   par_buf.kappa_sequential[ ikappa ] ; 

      inverter_type_zonked_heavy[ ikappa] = 
	par_buf.inverter_type_zonked_heavy[ ikappa ] ;
      inverter_type_sequential[ ikappa] = 
	par_buf.inverter_type_sequential[ ikappa ] ;

      wqs_spectator[ ikappa ] = par_buf.wqs_spectator [ ikappa ];
      init_qs(&wqs_spectator[ ikappa ]);
      wqs_spectator[ ikappa ].type = par_buf.wqs_spectator [ ikappa ].type;
      wqs_zonked_light[ ikappa ] = par_buf.wqs_zonked_light[ ikappa ];
      init_qs(&wqs_zonked_light[ ikappa ]);
      wqs_zonked_light[ ikappa ].type = par_buf.wqs_zonked_light[ ikappa ].type;
      wqs_zonked_heavy[ ikappa ] = par_buf.wqs_zonked_heavy[ ikappa ];
      init_qs(&wqs_zonked_heavy[ ikappa ]);
      wqs_zonked_heavy[ ikappa ].type = par_buf.wqs_zonked_heavy[ ikappa ].type;

      startflag_spectator[ ikappa ] = par_buf.startflag_spectator[ ikappa ]   ;
      startflag_zonked_light[ ikappa ] =  par_buf.startflag_zonked_light[ ikappa ];
      saveflag_zonked_heavy[ ikappa ] =  par_buf.saveflag_zonked_heavy[ ikappa ];
      
      strcpy(qfile_spectator[ikappa] , par_buf.qfile_spectator[ikappa]);
      strcpy(qfile_zonked_light[ikappa] , par_buf.qfile_zonked_light[ikappa]);
      strcpy(qfile_zonked_heavy[ikappa] , par_buf.qfile_zonked_heavy[ikappa]);
    }

  strcpy(qfile_suffix_zonked_light, par_buf.qfile_suffix_zonked_light);
  strcpy(qfile_suffix_zonked_heavy, par_buf.qfile_suffix_zonked_heavy);
  
  tf = par_buf.tf  ;
  saveflag_zonked_light_ssink   = par_buf.saveflag_zonked_light_ssink ;  
  saveflag_zonked_heavy_ssink   = par_buf.saveflag_zonked_heavy_ssink ;  

  /* Set start and reload flags for propagator files */
  reloadflag_zonked_light_ssink = set_reload_flag(saveflag_zonked_light_ssink);
  reloadflag_zonked_heavy_ssink = set_reload_flag(saveflag_zonked_heavy_ssink);
  for(ikappa = 0; ikappa < par_buf.no_zonked_heavy; ikappa++)
    startflag_zonked_heavy[ikappa] = 
      set_reload_flag(saveflag_zonked_heavy[ikappa]);

  /**** names of the files for the three point function ********************/
  strcpy(filename_HH3,par_buf.filename_HH3);
  strcpy(filename_HL3,par_buf.filename_HL3);

  saveflag_HH3             = par_buf.saveflag_HH3;
  saveflag_HL3             = par_buf.saveflag_HL3;
  
  /**** names of the files for the two point function ********************/
  strcpy(filename_HH2_GL     , par_buf.filename_HH2_GL) ; 
  strcpy(filename_LL2_GG     , par_buf.filename_LL2_GG) ; 
  strcpy(filename_HL2_GG , par_buf.filename_HL2_GG) ;
  strcpy(filename_HL2_GE , par_buf.filename_HL2_GE); 
  strcpy(filename_HL2_GL , par_buf.filename_HL2_GL); 

  saveflag_HH2_GL = par_buf.saveflag_HH2_GL;
  saveflag_LL2_GG = par_buf.saveflag_LL2_GG;
  saveflag_HL2_GG = par_buf.saveflag_HL2_GG;
  saveflag_HL2_GE = par_buf.saveflag_HL2_GE;
  saveflag_HL2_GL = par_buf.saveflag_HL2_GL;
  
  no_p_values  = par_buf.no_p_values ;
  no_q_values = par_buf.no_q_values ;
  no_k_values = par_buf.no_k_values ;
  
  for(i=0 ; i < no_p_values ; ++i )
    {
      strcpy(seq_smear_file[i]  , par_buf.seq_smear_file[i]);
    }

  
  for(imom = 0 ; imom < MAXMOM ; ++imom)
    {
      for(ii = 0 ; ii < 3 ; ++ii)
      {
	q_momstore[imom][ii] = par_buf.q_momstore[imom][ii] ;  
        
	k_momstore[imom][ii] = par_buf.k_momstore[imom][ii] ;  
      }      
    }
  
  for(imom = 0 ; imom < MAXPMOM ; ++imom)
    {
      for(ii = 0 ; ii < 3 ; ++ii)
	p_momstore[imom][ii] = par_buf.p_momstore[imom][ii] ;  
    }
  

  /*** create/load the gauge configuration ****/
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile )  ;
  /* We don't do any gauge fixing here */
  fixflag = NO_GAUGE_FIX;
  
  return(0);
}

