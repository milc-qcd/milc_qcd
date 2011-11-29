/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "cl_hyb_includes.h"
#include <string.h>

/** A useful macro for copying the parameters **/
#define OUT_OF_PARAM(xx)  xx = par_buf.xx 

int  setup_operator_choice(int prompt) ;
int initial_set();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int prompt;

	/* print banner, get volume, seed */
    prompt=initial_set();
   	/* initialize the node random number generator */
    initialize_prn(&node_prn,iseed,volume+mynode());
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
	printf("SU3 Hybrid Meson Spectrum with clover fermions\n");
	printf("Microcanonical simulation with refreshing\n");
#ifdef FORWARDS_ONLY
 printf(">>>> Forward time evolution (Dong-Liu valence approximation)<<<<\n"); 
#endif
	printf("MIMD version 6\n");
	printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
	time_stamp("start");

	status=get_prompt(stdin, &prompt);

	IF_OK status += get_i(stdin,prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(stdin,prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(stdin,prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(stdin,prompt,"nt", &par_buf.nt );
	IF_OK status += get_i(stdin,prompt,"iseed", &par_buf.iseed );

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
    iseed=par_buf.iseed;
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=nx*ny*nz*nt;
    total_iters=0;
    return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  Real x;
  char savebuf[128] ;


    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

	printf("\n\n");
	status=0;
    
	/* get couplings and broadcast to nodes	*/
	/* beta, kappa */
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );
	IF_OK status += get_i(stdin, prompt,"verbose_flag", 
			      &par_buf.verbose_flag );
	IF_OK status += get_f(stdin, prompt,"kappa", &par_buf.kappa );

	/* Clover coefficient, u0 */
	IF_OK status += get_f(stdin, prompt,"clov_c", &par_buf.clov_c );
	IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );


	/* source time slice and increment */
	IF_OK status += get_i(stdin, prompt,"source_start",
			      &par_buf.source_start);
	IF_OK status += get_i(stdin, prompt,"source_inc",&par_buf.source_inc);
	IF_OK status += get_i(stdin, prompt,"n_sources",&par_buf.n_sources);


	/*** parameters controlling the smearing of the links ***/
#ifdef SMEAR
	IF_OK status += get_f(stdin, prompt,"space_simple_weight",&par_buf.space_simple_weight);
	IF_OK status += get_f(stdin, prompt,"space_norm_factor",&par_buf.space_norm_factor);
	IF_OK status += get_f(stdin, prompt,"time_simple_weight",&par_buf.time_simple_weight);
	IF_OK status += get_f(stdin, prompt,"time_norm_factor",&par_buf.time_norm_factor);
	IF_OK status += get_i(stdin, prompt,"smearing_level",
			      &par_buf.smearing_level);
#else
      smearing_level = 0   ;
#endif
    
	/* maximum no. of conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );
    
	/* error for propagator conjugate gradient */
	IF_OK status += get_f(stdin, prompt,"error_for_propagator", &x );
	IF_OK par_buf.rsqprop = x*x;
    
        /* find out what kind of starting lattice to use */
	IF_OK ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

	/*** gauge fixing yes or no choice *******/
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



	IF_OK if (prompt==1) 
	  printf("enter 'local_quark_src', or 'wall_quark_src'\n");
	IF_OK scanf("%s",savebuf);
	IF_OK printf("%s\n",savebuf);
	IF_OK {
	  if(strcmp("local_quark_src",savebuf) == 0 ){
	    par_buf.wot_src = LOCAL_QUARK_SRC  ;
	  }
	  else if(strcmp("wall_quark_src",savebuf) == 0 ) {
	    par_buf.wot_src = WALL_QUARK_SRC  ;
	  }
	  else{
	    printf("error in input: fixing_command is invalid\n"); status++;
	  }
	}


	/*** boundary conditions yes or no  *******/
	IF_OK if (prompt==1) 
	  printf("enter 'periodic_everywhere', or 'anti_periodic_in_time'\n");
	IF_OK scanf("%s",savebuf);
	IF_OK printf("%s\n",savebuf);
	IF_OK {
	  if(strcmp("periodic_everywhere",savebuf) == 0 ){
	    par_buf.boundary_flag = PERIODIC_EVERYWHERE;
	  }
	  else if(strcmp("anti_periodic_in_time",savebuf) == 0 ) {
	    par_buf.boundary_flag   =  ANTI_PERIODIC_IN_TIME;
	  }
	  else{
	    printf("error in input: boundary_command (= %s) is invalid\n",savebuf ); status++;
	  }
	}



        IF_OK status += setup_operator_choice(prompt) ; 

	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */
    
    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));
    
    if( par_buf.stopflag != 0 )
      normal_exit(0);

    startflag = par_buf.startflag;
    niter = par_buf.niter;
    rsqprop = par_buf.rsqprop;
    beta = par_buf.beta;
    kappa = par_buf.kappa;
    source_start = par_buf.source_start;
    source_inc = par_buf.source_inc;
    n_sources = par_buf.n_sources;
    strcpy(startfile,par_buf.startfile);

    boundary_flag = par_buf.boundary_flag  ;
    u0 = par_buf.u0             ;
    clov_c  = par_buf.clov_c    ;
    fixflag  = par_buf.fixflag  ;
    wot_src = par_buf.wot_src   ;
    verbose_flag = par_buf.verbose_flag   ;

  /** copy the operator flags **/
  OUT_OF_PARAM(oper_PION_SOURCE);
  OUT_OF_PARAM(oper_PION2_SOURCE);
  OUT_OF_PARAM(oper_RHO_SOURCE); 
  OUT_OF_PARAM(oper_RHO2_SOURCE);
  OUT_OF_PARAM(oper_A1P_SOURCE);
  OUT_OF_PARAM(oper_A1_SOURCE);
  OUT_OF_PARAM(oper_ZEROMP_SOURCE);
  OUT_OF_PARAM(oper_ZEROPM_SOURCE);
  OUT_OF_PARAM(oper_ZEROPMP_SOURCE);
  OUT_OF_PARAM(oper_ZEROPMB_SOURCE);
  OUT_OF_PARAM(oper_ZEROMM_SOURCE); 
  OUT_OF_PARAM(oper_ZEROMMP_SOURCE);
  OUT_OF_PARAM(oper_ONEMP_SOURCE);
  OUT_OF_PARAM(oper_ONEMP2_SOURCE);
  OUT_OF_PARAM(oper_ONEMM_SOURCE);
  OUT_OF_PARAM(oper_ONEPP_SOURCE);
  OUT_OF_PARAM(oper_QQQQ_SOURCE);


  OUT_OF_PARAM(space_simple_weight);
  OUT_OF_PARAM(space_norm_factor);
  OUT_OF_PARAM(time_simple_weight);
  OUT_OF_PARAM(time_norm_factor);
  OUT_OF_PARAM(smearing_level);



    /* Do whatever is needed to get lattice */
  if(startflag != CONTINUE){
    reload_lattice( startflag, startfile );
    invalidate_this_clov(gen_clov);
  }

    return(0);
}


/**
 **   Load in thew choice of operators to RUN
 **/


int  setup_operator_choice(int prompt)
{
  int status = 0 ; 
  char savebuf[128] ;
  int what_to_do ; 
  enum what_to_do_choices { read_operators = 50 ,  stop_reading_operators  } ; 

  /***------------------------------------------------------------****/

  /*** set all flags so that all the operatosre are NOT calculated ***/

  par_buf.oper_PION_SOURCE = DO_NOT_CALCULATE ; 
  par_buf.oper_PION2_SOURCE = DO_NOT_CALCULATE ;
  par_buf.oper_RHO_SOURCE = DO_NOT_CALCULATE  ; 
  par_buf.oper_RHO2_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_A1P_SOURCE  = DO_NOT_CALCULATE ;
  par_buf.oper_A1_SOURCE  = DO_NOT_CALCULATE ;
  par_buf.oper_ZEROMP_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_ZEROPM_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_ZEROPMP_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_ZEROPMB_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_ZEROMM_SOURCE = DO_NOT_CALCULATE; 
  par_buf.oper_ZEROMMP_SOURCE = DO_NOT_CALCULATE;
  par_buf.oper_ONEMP_SOURCE  = DO_NOT_CALCULATE;
  par_buf.oper_ONEMP2_SOURCE = DO_NOT_CALCULATE ;
  par_buf.oper_ONEMM_SOURCE = DO_NOT_CALCULATE ;
  par_buf.oper_ONEPP_SOURCE = DO_NOT_CALCULATE ;
  par_buf.oper_QQQQ_SOURCE   = DO_NOT_CALCULATE;


  /*** look at all the operators that the user wants to calculate ***/

  IF_OK if (prompt==1) 
    printf("enter 'start_of_operators' \n");
  IF_OK scanf("%s",savebuf);
  IF_OK printf("%s\n",savebuf);
	IF_OK {
	  if(strcmp("start_of_operators",savebuf) == 0 )
	  {
	    what_to_do = read_operators ;
	  }
	  else
	  {
	    printf("error in input: start_of_operators required, but %s got\n",savebuf);
	    status++ ;
	    what_to_do = stop_reading_operators ; 
	  }
	}


  /*** useful macro to save typing space ***/
 /* The following macro doesn't work with the Solaris 2.6 preprocessor */
	/*#define OPER_OPTION(a)    else if(strcmp(#a,savebuf) == 0 )  par_buf.a = CALCULATE */
 /* The following macro doesn't work on other systems, such as IBM AIX.
    We need a better one! */
#define OPER_OPTION(a)    else if(strcmp("a",savebuf) == 0 )  par_buf.a = CALCULATE

  /*** now read the operators to compute ***/
  while( what_to_do == read_operators   ) 
  {

    IF_OK if (prompt==1) 
      printf("enter an operator name to compute the spectrum of\n");
    IF_OK scanf("%s",savebuf);
    IF_OK printf("%s\n",savebuf);


    if(strcmp("end_operators",savebuf) == 0 )
    {
      what_to_do = stop_reading_operators ; 
    }
    else if(strcmp("oper_PION_SOURCE",savebuf) == 0) par_buf.oper_PION_SOURCE = CALCULATE;
    else if(strcmp("oper_PION2_SOURCE",savebuf) == 0) par_buf.oper_PION2_SOURCE = CALCULATE;
    else if(strcmp("oper_RHO_SOURCE",savebuf) == 0) par_buf.oper_RHO_SOURCE = CALCULATE;
    else if(strcmp("oper_RHO2_SOURCE",savebuf) == 0) par_buf.oper_RHO2_SOURCE = CALCULATE;
    else if(strcmp("oper_A1P_SOURCE",savebuf) == 0) par_buf.oper_A1P_SOURCE = CALCULATE;
    else if(strcmp("oper_A1_SOURCE",savebuf) == 0) par_buf.oper_A1_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROMP_SOURCE",savebuf) == 0) par_buf.oper_ZEROMP_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROPM_SOURCE",savebuf) == 0) par_buf.oper_ZEROPM_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROPMP_SOURCE",savebuf) == 0) par_buf.oper_ZEROPMP_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROPMB_SOURCE",savebuf) == 0) par_buf.oper_ZEROPMB_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROMM_SOURCE",savebuf) == 0) par_buf.oper_ZEROMM_SOURCE = CALCULATE;
    else if(strcmp("oper_ZEROMMP_SOURCE",savebuf) == 0) par_buf.oper_ZEROMMP_SOURCE = CALCULATE;
    else if(strcmp("oper_ONEMP_SOURCE",savebuf) == 0) par_buf.oper_ONEMP_SOURCE = CALCULATE;
    else if(strcmp("oper_ONEMP2_SOURCE",savebuf) == 0) par_buf.oper_ONEMP2_SOURCE = CALCULATE;
    else if(strcmp("oper_ONEMM_SOURCE",savebuf) == 0) par_buf.oper_ONEMM_SOURCE = CALCULATE;
    else if(strcmp("oper_ONEPP_SOURCE",savebuf) == 0) par_buf.oper_ONEPP_SOURCE = CALCULATE;
    else if(strcmp("oper_QQQQ_SOURCE",savebuf) == 0) par_buf.oper_QQQQ_SOURCE = CALCULATE;
    else 
    {
      printf("operator = %s is not understood\n",savebuf ) ; 
      ++status ; 
      what_to_do = stop_reading_operators ; 
    }

      


  }  /*** end the loop over the operators to calculate ***/






  return status ; 
}


