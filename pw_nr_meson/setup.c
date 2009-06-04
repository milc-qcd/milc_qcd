/*************************** setup.c *************************************/
/* MIMD version 7 */
#define IF_OK if(status==0)


#include "pw_nr_meson_includes.h"

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
 /* Create clover structure */
 gen_clov = create_clov();

 return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, and send to others */
    if(mynode()==0){
	/* print banner */
	printf("Heavy-Heavy p-wave spectroscopy\n");
	printf("MIMD version 7\n");
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

/* read in parameters and coupling constants	*/

int readin(int prompt)  {
  
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status,i,j;
  char savebuf[128];
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    status=0;
    
    /* find out what kind of starting lattice to use */
    IF_OK status += get_i(stdin, prompt, "sequence", &par_buf.sequence);
    IF_OK status += ask_starting_lattice( stdin, prompt, &(par_buf.startflag),
					  par_buf.startfile );
    
    IF_OK if (prompt!=0) 
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
    IF_OK status += ask_ending_lattice( stdin, prompt, &(par_buf.saveflag),  	    par_buf.savefile );
    
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				 par_buf.stringLFN );

    /* Inversion control */
    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &par_buf.qic.max );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &par_buf.qic.nrestart );
    
    /* error for propagator conjugate gradient */
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			  &par_buf.qic.resid );
    IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			  &par_buf.qic.relresid );
    /* Precision fixed to prevailing precision for now */
    par_buf.qic.prec = PRECISION;
    par_buf.qic.parity = EVENANDODD;
    

    /* Quark parameters */

    IF_OK status += get_f(stdin, prompt,"kappa", &par_buf.dcp.Kappa );
    IF_OK status += get_f(stdin, prompt,"clov_c", &par_buf.dcp.Clov_c );
    IF_OK status += get_f(stdin, prompt,"u0", &par_buf.dcp.U0 );

    IF_OK status += ask_starting_wprop( stdin, prompt, &par_buf.a_startflag_w,
					par_buf.a_startfile_w);
    
    IF_OK status += ask_ending_wprop( stdin, prompt, &par_buf.a_saveflag_w,
				      par_buf.a_savefile_w);
    
    /* Get antiquark source type */
    IF_OK status += get_w_quark_source( stdin, prompt, &par_buf.a_wqs);
    
    /* Number of smearings */
    
    IF_OK status += get_i(stdin, prompt,"number_of_smearings", &par_buf.num_smear);
    if( par_buf.num_smear>NSM ){
      printf("num_smear = %d must be <= %d!\n", par_buf.num_smear, NSM);
      status++;
    }
    
    for(j=0; j<par_buf.num_smear;j++){
      /* S-wave part of shell source wave function for the propagator */
      IF_OK status += get_w_quark_source( stdin, prompt, 
					&par_buf.source_wqs[j] );
      IF_OK status += get_s(stdin, prompt,"wave_func_label", 
			    par_buf.source_wf_label[j]);
      for(i=0;i<MAXDIR;i++){
	IF_OK status += ask_starting_wprop( stdin, prompt, 
					    &par_buf.startflag_w[j][i],
					    par_buf.startfile_w[j][i]);
	IF_OK status += ask_ending_wprop( stdin, prompt, 
					  &par_buf.saveflag_w[j][i],
					  par_buf.savefile_w[j][i]);
      }
      /* S-wave part of relative smearing wave function for sink */
      IF_OK status += get_w_quark_sink( stdin, prompt, &par_buf.sink_wqs[j]);
      IF_OK status += get_s(stdin, prompt,"wave_func_label", 
			    par_buf.sink_wf_label[j]);
    }

    IF_OK status += get_s( stdin, prompt, "save_a0",    par_buf.a0_file    );
    IF_OK status += get_s( stdin, prompt, "save_b1",    par_buf.b1_file    );
    IF_OK status += get_s( stdin, prompt, "save_a1",    par_buf.a1_file    );
    IF_OK status += get_s( stdin, prompt, "save_a2_t2", par_buf.a2_t2_file );
    IF_OK status += get_s( stdin, prompt, "save_a2_e",  par_buf.a2_e_file  );
    IF_OK status += get_s( stdin, prompt, "save_a2",    par_buf.a2_file    );
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )
    normal_exit(0);

  sequence = par_buf.sequence;
  startflag = par_buf.startflag;
  fixflag = par_buf.fixflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  
  qic = par_buf.qic;
  dcp = par_buf.dcp;
  
  num_smear = par_buf.num_smear;
  for(j=0;j<num_smear;j++){
    for(i=0;i<MAXDIR;i++){
      startflag_w[j][i] = par_buf.startflag_w[j][i];
      strcpy(startfile_w[j][i],par_buf.startfile_w[j][i]);
      saveflag_w[j][i] = par_buf.saveflag_w[j][i];
      strcpy(savefile_w[j][i],par_buf.savefile_w[j][i]);
    }
  }
  for(i=0;i<par_buf.num_smear;i++){
    source_wqs[i] = par_buf.source_wqs[i];
    init_wqs(&source_wqs[i]);
    source_wqs[i].type = par_buf.source_wqs[i].type;
    sink_wqs[i] = par_buf.sink_wqs[i];
    init_wqs(&sink_wqs[i]);
    sink_wqs[i].type = par_buf.sink_wqs[i].type;
    strcpy(source_wf_label[i],par_buf.source_wf_label[i]);
    strcpy(sink_wf_label[i],par_buf.sink_wf_label[i]);
  }

  a_startflag_w = par_buf.a_startflag_w;
  strcpy(a_startfile_w, par_buf.a_startfile_w); 
  a_saveflag_w = par_buf.a_saveflag_w;
  strcpy(a_savefile_w, par_buf.a_savefile_w); 
  a_wqs = par_buf.a_wqs;
  init_wqs(&a_wqs);
  a_wqs.type = par_buf.a_wqs.type;

  strcpy(a0_file, par_buf.a0_file);
  strcpy(b1_file, par_buf.b1_file);
  strcpy(a1_file, par_buf.a1_file);
  strcpy(a2_t2_file, par_buf.a2_t2_file);
  strcpy(a2_e_file, par_buf.a2_e_file);
  strcpy(a2_file, par_buf.a2_file);
  
  /* Do whatever is needed to get lattice */
  if( startflag != CONTINUE ){
    startlat_p = reload_lattice( startflag, startfile );
    invalidate_this_clov(gen_clov);
  }
  
  return(0);
}

