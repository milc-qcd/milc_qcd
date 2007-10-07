/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

* 5/30/07 Created from setup_cl.c */

//  $Log: setup.c,v $
//  Revision 1.1  2007/10/07 20:02:32  detar
//  Add new application.  Generarlizes clover_invert.
//
//


#include "cl_inv_includes.h"
#include <string.h>
int initial_set();

#include "params.h"

int setup()   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  
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

/* Forward declarations */
static int hash_meson_label(char meson_label[MAX_MESON][MAX_MESON_LABEL], 
			    char *meson_label_in);
static int 
hash_momentum_label(char mom_label[MAX_MESON_MOMENTUM][MAX_MOM_LABEL], 
		    char *mom_label_in);

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
    /* Gauge configuration section                                */
    /*------------------------------------------------------------*/

    IF_OK status += ask_starting_lattice(stdin,  prompt, &param.startflag,
	param.startfile );

    IF_OK if (prompt!=0) 
      printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
    IF_OK scanf("%s",savebuf);
    IF_OK printf("%s\n",savebuf);
    IF_OK {
      if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
	param.fixflag = COULOMB_GAUGE_FIX;
      }
      else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
	param.fixflag = NO_GAUGE_FIX;
      }
      else{
	printf("error in input: fixing_command is invalid\n"); status++;
      }
    }
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
			     param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				  param.stringLFN );
    
#if 0
    /*------------------------------------------------------------*/
    /* Scratch file stem for saved propagators                    */
    /*------------------------------------------------------------*/

    IF_OK status += ask_ending_wprop( stdin, prompt,&param.scratchflag,
				      param.scratchstem_w);
#endif
    
    /*------------------------------------------------------------*/
    /* Propagator inversion control                               */
    /*------------------------------------------------------------*/

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &param.qic.max );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &param.qic.nrestart );
    
    /* error for propagator conjugate gradient */
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			  &param.qic.resid );
    IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			  &param.qic.relresid );
    /* Precision fixed to prevailing precision for now */
    param.qic.prec = PRECISION;
    param.qic.parity = EVENANDODD;
    

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
      
      IF_OK status += get_f(stdin, prompt,"kappa", &param.dcp[i].Kappa );
      IF_OK status += get_f(stdin, prompt,"clov_c", &param.dcp[i].Clov_c );
      IF_OK status += get_f(stdin, prompt,"u0", &param.dcp[i].U0 );
      IF_OK status += get_s(stdin, prompt,"check", savebuf);
      IF_OK {
	/* Should we be checking the propagator by running the solver? */
	if(strcmp(savebuf,"no") == 0)param.check[i] = 0;
	else param.check[i] = 1;
      }

      IF_OK status += ask_starting_wprop( stdin, prompt, &param.startflag_w[i],
					  param.startfile_w[i]);
      
      IF_OK status += ask_ending_wprop( stdin, prompt, &param.saveflag_w[i],
					param.savefile_w[i]);
      
      /* Get source type */
      IF_OK init_wqs(&param.src_wqs[i]);
      IF_OK status += get_quark_source( stdin, prompt, &param.src_wqs[i]);
      
      /* Get rotation parameter */
      IF_OK status += get_f( stdin, prompt, "d1", &param.d1[i]);
    }
      
    /*------------------------------------------------------------*/
    /* Meson correlators                                          */
    /*------------------------------------------------------------*/
    
    /* Number of quark pairs */
    IF_OK status += get_i(stdin,prompt,"number_of_pairings", &param.num_pair );
    if( param.num_pair>MAX_PAIR ){
      printf("num_pair = %d must be <= %d!\n", param.num_pair, MAX_PAIR);
      status++;
    }

    IF_OK for(i = 0; i < param.num_pair; i++){
      char request_buf[MAX_SPECTRUM_REQUEST];

      /* Which quarks in the pair? */
      IF_OK get_vi(stdin, prompt, "pair", param.qkpair[i], 2);
      IF_OK {
	int j;
	for(j = 0; j < 2; j++){
	  if(param.qkpair[i][j] < 0 || param.qkpair[i][j] >= param.num_qk){
	    printf("Quark index %d must be in [0,%d]\n",
		   param.qkpair[i][j],param.num_qk);
	    status++;
	  }
	}
      }

      /* What spectrum calculation to do? */
      /* prepend and append a comma for ease in parsing */
      IF_OK request_buf[0] = '\0';
      IF_OK strcpy(request_buf,",");
      IF_OK status += get_s(stdin, prompt,"spectrum_request", request_buf+1 );
      IF_OK strcat(request_buf,",");

      /* Parse the spectrum request */
      IF_OK {
	/* Point sink meson */
	if(strstr(request_buf,",point,") != NULL)
	  param.do_point_meson_spect[i] = 1;
	else param.do_point_meson_spect[i] = 0;

	/* Smeared sink meson */
	if(strstr(request_buf,",smeared,") != NULL)
	  param.do_smear_meson_spect[i] = 1;
	else param.do_smear_meson_spect[i] = 0;

	/* Rotated point sink meson */
	if(strstr(request_buf,",3drotated,") != NULL)
	  param.do_rot_meson_spect[i] = 1;
	else param.do_rot_meson_spect[i] = 0;
	
	/* Baryons */
	if(strstr(request_buf,",baryon,") != NULL)
	  param.do_baryon_spect[i] = 1;
	else param.do_baryon_spect[i] = 0;
      }

      /* What sink smearing wave function? (Only if smearing) */

      IF_OK {
	if(param.do_smear_meson_spect[i]){
	  IF_OK init_wqs(&param.snk_wqs[i]);
	  IF_OK status += get_quark_sink( stdin, prompt, &param.snk_wqs[i]);
	  IF_OK status += get_s(stdin, prompt, "wave_func_label",
				param.snk_label[i]);
	}
      }

      /* What file for the resulting correlators? */

      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_c[i],
				     param.savefile_c[i]);
    }

    /*------------------------------------------------------------*/
    /* Gamma matrix pairings                                      */
    /*------------------------------------------------------------*/
    
    /* Number of quark pairs */
    IF_OK status += get_i(stdin,prompt,"number_of_mesons", &param.num_meson );
    IF_OK {
      if( param.num_meson>MAX_MESON ){
	printf("num_meson = %d must be <= %d!\n", param.num_meson, MAX_MESON);
	status++;
      }
      param.num_meson_report = 0;
    }

    IF_OK for(i = 0; i < param.num_meson; i++){
      char meson_label_in[MAX_MESON_LABEL];
      char gam_src_lab[MAXGAMMA], gam_snk_lab[MAXGAMMA];
      char phase_lab[4];

      /* meson label */
      IF_OK status += get_sn(stdin, prompt, "meson_source_sink", meson_label_in);
      IF_OK {
	param.meson_index[i] = 
	  hash_meson_label(param.meson_label, meson_label_in);
	if(param.meson_index[i] < 0)status++;
      }

      /* gamma matrix labels and phase label */
      IF_OK scanf("%s %s %s\n",gam_src_lab,gam_snk_lab,phase_lab);
      IF_OK printf(" \t%s\t%s\t%s\n",gam_src_lab,gam_snk_lab,phase_lab);

      /* decode gamma matrix labels for source and sink */
      IF_OK {
	param.gam_src[i] = gamma_index(gam_src_lab);
	if(param.gam_src[i] < 0){
	  printf("%s is not a valid gamma matrix label\n",gam_src_lab);
	  status += 1;
	}
	param.gam_snk[i] = gamma_index(gam_snk_lab);
	if(param.gam_snk[i] < 0){
	  printf("%s is not a valid gamma matrix label\n",gam_snk_lab);
	  status ++;
	}
      }

      /* decode phase for correlator */
      IF_OK {
	param.meson_phase[i] = decode_phase(phase_lab);
	if(param.meson_phase[i].real == 0 &&
	   param.meson_phase[i].imag == 0 ){
	  printf("%s is not a valid phase label\n",phase_lab);
	  status ++;
	}
      }
    }
    
    /*------------------------------------------------------------*/
    /* Momentum table                                             */
    /*------------------------------------------------------------*/

    /* Number of momenta */
    IF_OK status += get_i(stdin,prompt,"number_of_meson_momenta", 
			  &param.num_mom );
    IF_OK {
      if( param.num_mom>MAX_MESON_MOMENTUM ){
	printf("num_mom = %d must be <= %d!\n", 
	       param.num_mom, MAX_MESON_MOMENTUM);
	status++;
      }
      param.num_mom_report = 0;
    }
    
    /* momentum label */

    IF_OK for(i = 0; i < param.num_mom; i++){
      char mom_label_in[MAX_MOM_LABEL];
      IF_OK status += get_sn(stdin, prompt, "momentum", mom_label_in);

      IF_OK {
	param.mom_index[i] = 
	  hash_momentum_label(param.mom_label, mom_label_in);
	if(param.mom_index[i] < 0)status++;
      }
    
      /* momentum indices */
      IF_OK {
	if(scanf("%d %d %d",&param.meson_mom[i][0],
		 &param.meson_mom[i][1],&param.meson_mom[i][2]) != 3){
	  printf("Format error on momentum line\n");
	  status++;
	}
	
	IF_OK printf(" %d %d %d\n",param.meson_mom[i][0],
		     param.meson_mom[i][1],param.meson_mom[i][2]);
      }
    }
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    

  broadcast_bytes((char *)&param,sizeof(param));
  if( param.stopflag != 0 )
    normal_exit(0);

  /* Do whatever is needed to get lattice */
  startlat_p = reload_lattice( param.startflag, param.startfile );

  return(0);
}

/* Manage hash table for meson labels */

/* A label can be repeated, in which case the meson propagators will
   be averaged */

/* Look for meson_label_in in the meson_label table.  If found, return its
   index.  Otherwise, add it to the table and assign it the new index */
static int hash_meson_label(char meson_label[MAX_MESON][MAX_MESON_LABEL], 
			    char *meson_label_in){
  int i;

  for(i = 0; i <= param.num_meson_report; i++){
    if(strcmp(meson_label_in, meson_label[i]) == 0)
      return i;
  }
  if(param.num_meson_report >= MAX_MESON){
    printf("Too many mesons labels %d\n",param.num_meson_report);
    return -1;
  }
  strcpy(meson_label[param.num_meson_report],meson_label_in);
  param.num_meson_report++;
  return param.num_meson_report-1;
}

/* Manage hash table for meson momentum labels */

/* A label can be repeated, in which case the propagators will
   be averaged */

/* Look for mom_label_in in the mom_label table.  If found, return its
   index.  Otherwise, add it to the table and assign it the new index */
static int 
hash_momentum_label(char mom_label[MAX_MESON_MOMENTUM][MAX_MOM_LABEL], 
		    char *mom_label_in){
  int i;

  for(i = 0; i <= param.num_mom_report; i++){
    if(strcmp(mom_label_in, mom_label[i]) == 0)
      return i;
  }
  if(param.num_mom_report >= MAX_MESON_MOMENTUM){
    printf("Too many momentum labels %d\n",param.num_mom_report);
    return -1;
  }
  strcpy(mom_label[param.num_mom_report],mom_label_in);
  param.num_mom_report++;
  return param.num_mom_report-1;
}

