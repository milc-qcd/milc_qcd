/******** setup_p_cl.c *********/
/* MIMD code version 7 */
#define IF_OK if(status==0)

#include "arb_ov_includes.h"
int initial_set();
void make_fields();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"


int  setup_p()   {
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
#ifdef RANDOM
        /* initialize the node random number generator */
    initialize_prn(&node_prn,iseed,volume+mynode());
#endif
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* allocate space for fields */
  make_fields();

  setup_offset();
  
  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
/* stringification kludge from gnu preprocessor manual
   http://gcc.gnu.org/onlinedocs/cpp/Stringification.html */
#define XSTR(s) STR(s)
#define STR(s) #s
/* end kluge */
#ifdef NHYP
    printf("nHYP links with %d smearing level(s)\n", SMEAR_LEVEL);
#ifdef HARD_CODE_SMEAR
    printf("HARD-CODED alpha_smear={%g, %g, %g}\n",
           alpha_smear[0], alpha_smear[1], alpha_smear[2] );
#else
    printf("Reading alpha_smear parameters from infile\n");
#endif
#if (SMEAR_LEVEL<3)
    printf("REMEMBER: last %d smearing parameter(s) not used in this run\n",
           3-SMEAR_LEVEL);
#endif
    printf("IR_STAB=%e,  EPS_SQ=%e \n", (Real)IR_STAB, (Real)EPS_SQ );
#endif /* NHYP */

#ifndef TPERIODIC
    printf("Antiperiodic boundary conditions\n");
#endif
    printf("MIMD version 7 $Name:  $\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    time_stamp("start");

    
    status = get_prompt(stdin, &prompt );
    
    IF_OK status += get_i( stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i( stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i( stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i( stdin, prompt,"nt", &par_buf.nt );
#ifdef DOMAINX
    IF_OK status += get_i( stdin, prompt,"cut_x", &par_buf.cut_x );
    IF_OK status += get_i( stdin, prompt,"cut_y", &par_buf.cut_y );
    IF_OK status += get_i( stdin, prompt,"cut_z", &par_buf.cut_z );
    IF_OK status += get_i( stdin, prompt,"cut_t", &par_buf.cut_t );
#endif
#ifdef RANDOM
        IF_OK status += get_i( stdin, prompt,"iseed", &par_buf.iseed );
#endif
    
    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));

  if( par_buf.stopflag != 0 ){
    normal_exit(0);
  }
  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
#ifdef RANDOM
    iseed=par_buf.iseed;
#endif
#ifdef DOMAINX
    cut_x=par_buf.cut_x;
    cut_y=par_buf.cut_y;
    cut_z=par_buf.cut_z;
    cut_t=par_buf.cut_t;
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
  char tempstring[128];
  char save_w[128];
  int i;



  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;


        /* get couplings and broadcast to nodes */
        /* beta, mass */
        IF_OK status += get_f( stdin, prompt,"beta", &par_buf.beta );
        IF_OK status += get_f( stdin, prompt,"u0", &par_buf.u0 );



    
    /* Number of masses */
    IF_OK status += get_i( stdin, prompt,"number_of_masses", &par_buf.num_masses );
    if( par_buf.num_masses>MAX_MASSES ){
      printf("num_masses = %d must be <= %d!\n", par_buf.num_masses, MAX_MASSES);
      status++;
    }
    if( par_buf.num_masses==0) status++;
    
    /* To be safe initialize the following to zero */
    for(i=0;i<MAX_MASSES;i++){
      mass[i] = 0.0;
      resid[i] = 1.0;
    }
    
    for(i=0;i<par_buf.num_masses;i++){
      IF_OK status += get_f( stdin, prompt,"m0", &par_buf.mass[i] );
    }

      IF_OK status += get_f( stdin, prompt,"R0", &par_buf.R0 );
      IF_OK status += get_f( stdin, prompt,"scalez", &par_buf.scalez );
/* order of inner expansion */
    IF_OK status += get_f( stdin, prompt,"prec_sign", &par_buf.prec_sign );
    IF_OK status += get_f( stdin, prompt,"zolo_min", &par_buf.zolo_min );
    IF_OK status += get_f( stdin, prompt,"zolo_max", &par_buf.zolo_max );
    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i( stdin, prompt,"inner_cg_iterations", &par_buf.maxcg_inner );
     /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_f( stdin, prompt,"inner_residue", &par_buf.resid_inner );
    IF_OK status += get_f( stdin, prompt,"inner_residue_h", &par_buf.resid_inner_h );

#ifdef EIG
        IF_OK status += get_i( stdin, prompt,"Number_of_inner_eigenvals", &par_buf.Nvecs_h0r0 );
        IF_OK status += get_i( stdin, prompt,"Number_of_h0_eigenvals", &par_buf.Nvecs_h0 );
        IF_OK status += get_i( stdin, prompt,"Number_of_hov_eigenvals", &par_buf.Nvecs_hov );
	if(par_buf.Nvecs_hov>par_buf.Nvecs_h0){
node0_printf("Whoops! You need Number_of_h0_eigenvals >= Number_of_hov_eigenvals\n");
exit(1);}
        IF_OK status += get_i( stdin, prompt,"Max_Rayleigh_iters", &par_buf.MaxIter );
        IF_OK status += get_i( stdin, prompt,"Max_r0_iters", &par_buf.Maxr0Iter );
        IF_OK status += get_i( stdin, prompt,"Restart_Rayleigh", &par_buf.Restart );
        IF_OK status += get_i( stdin, prompt,"Kalkreuter_iters", &par_buf.Kiters );
         IF_OK status += get_f( stdin, prompt,"eigenvec_quality", &par_buf.eigenvec_qual );
       IF_OK status += get_f( stdin, prompt,"eigenval_tol_low",
                              &par_buf.eigenval_tol );
        IF_OK status += get_f( stdin, prompt,"error_decr_low", &par_buf.error_decr);
        IF_OK status += get_f( stdin, prompt,"eigenval_tol_high",
                              &par_buf.eigenval_tol_acc );
        IF_OK status += get_f( stdin, prompt,"error_decr_high", &par_buf.error_decr_acc);
#endif
    
    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i( stdin, prompt,"max_cg_iterations", &par_buf.niter );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i( stdin, prompt,"max_cg_restarts", &par_buf.nrestart );
    
    /* error for propagator conjugate gradient */
      for(i=0;i<par_buf.num_masses;i++){
      IF_OK status += get_f( stdin, prompt,"error_for_fermforce", &par_buf.resid[i] );
      }
      for(i=0;i<par_buf.num_masses;i++){
      IF_OK status += get_f( stdin, prompt,"error_for_fermaction", &par_buf.resid2[i] );
      }
    /* Get source type */
    /* (Same source type for each spectator) */
    init_qs(&par_buf.wqs);
    IF_OK status += get_wv_quark_source( stdin,prompt,&par_buf.wqs );
//    IF_OK status += ask_quark_source(stdin,prompt,&par_buf.wqs.type,
//				     par_buf.wqs.descrp);
//    /* width: psi=exp(-(r/r0)^2) */
//    IF_OK if (prompt==1)
//      printf("enter width(s) r0 as in: source=exp(-(r/r0)^2)\n");
//    IF_OK status += get_f( stdin, prompt,"r0", &par_buf.wqs.r0 );
//    /* (Hardwired source location for each spectator) */
//    IF_OK {
//	par_buf.wqs.x0 = source_loc[0];
//	par_buf.wqs.y0 = source_loc[1];
//	par_buf.wqs.z0 = source_loc[2];
//	par_buf.wqs.t0 = source_loc[3];

//#ifdef PROPOFF
///* code to do different source plane */
//	par_buf.wqs.t0 = source_loc[3] + nt/2 -2;
//	if(mynode()==0) printf("source plane for propagator = %d\n",
//           par_buf.wqs.t0);
//#endif
//
//    }


    /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
	par_buf.startfile );

    IF_OK status += get_i( stdin, prompt,"topology", &par_buf.topology);
    fflush(stdout);
    
       IF_OK if (prompt==1) printf(
           "enter 'no_gauge_fix', 'landau_gauge_fix', or 'coulomb_gauge_fix'\n");
       IF_OK status2=scanf("%s",savebuf);
       IF_OK {
         if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
           par_buf.fixflag = COULOMB_GAUGE_FIX;
           if(this_node==0)printf("fixing to coulomb gauge\n");
         }
         else if(strcmp("landau_gauge_fix",savebuf) == 0 ) {
           par_buf.fixflag = LANDAU_GAUGE_FIX;
           if(this_node==0)printf("fixing to landau gauge\n");
         }
         else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
           par_buf.fixflag = NO_GAUGE_FIX;
           if(this_node==0)printf("NOT fixing the gauge\n");
         }
         else{
           printf("error in input: fixing_command %s is invalid\n",savebuf);
           status++;
         }
       }

    /* find out what to do with lattice at end */
       IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
			     par_buf.savefile );
       IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				    par_buf.stringLFN );

    /* find out starting propagator(s) */

       IF_OK status += ask_starting_wprop( stdin, prompt,
			&par_buf.startflag_w[0],par_buf.startfile_w[0]);

    
    /* what to do with computed propagator */
       IF_OK status += ask_ending_wprop( stdin, prompt,
			&par_buf.saveflag_w[0],par_buf.savefile_w[0]);

      IF_OK for(i=1;i<par_buf.num_masses;i++){
	par_buf.startflag_w[i]= par_buf.startflag_w[0];
	par_buf.saveflag_w[i]= par_buf.saveflag_w[0];
      }

      /* append the mass to the startfile and savefile */
      if(par_buf.saveflag_w[0] != FORGET){
               strcpy(tempstring,par_buf.savefile_w[0]);
      IF_OK for(i=0;i<par_buf.num_masses;i++){
               sprintf(par_buf.savefile_w[i],"%s.n%f",tempstring,par_buf.mass[i]);
      }
      }
      if(par_buf.startflag_w[0] != FRESH && par_buf.startflag_w[0] != CONTINUE){
strcpy(tempstring,par_buf.startfile_w[0]);
      IF_OK for(i=0;i<par_buf.num_masses;i++){
sprintf(par_buf.startfile_w[i],"%s.n%f",tempstring,par_buf.mass[i]);
      }
      }
#ifdef H0INV
          IF_OK 
	    node0_printf("input files for writing inverses of H(0**2 + s**2 \n");
/*
      IF_OK status += ask_starting_prop(stdin, prompt,&par_buf.startflag_w3[0],
			par_buf.startfile_w3[0]);
*/

    
    /* what to do with computed propagator */
	  IF_OK status += ask_ending_prop(stdin, prompt,&par_buf.saveflag_w3[0],
		      par_buf.savefile_w3[0]);

      IF_OK for(i=1;i<par_buf.num_masses;i++){
/*
	par_buf.startflag_w3[i]= par_buf.startflag_w3[0];
*/
	par_buf.saveflag_w3[i]= par_buf.saveflag_w3[0];
      }

      /* append the mass to the startfile and savefile */
      if(par_buf.saveflag_w3[0] != FORGET){
strcpy(tempstring,par_buf.savefile_w3[0]);
      IF_OK for(i=0;i<par_buf.num_masses;i++){
sprintf(par_buf.savefile_w3[i],"%s.n%f",tempstring,par_buf.mass[i]);
      }
      }
/*
      if(par_buf.startflag_w3[0] != FORGET){
strcpy(tempstring,par_buf.startfile_w3[0]);
      IF_OK for(i=0;i<par_buf.num_masses;i++){
sprintf(par_buf.startfile_w3[i],"%s.n%f",tempstring,par_buf.mass[i]);
      }
      }
*/
#endif

      /* input hr0 */
    IF_OK if(prompt==1)
      printf("eigenmode h0(-r0) outfile:\n enter 'fresh_hr0_modes' 'iserial_hr0_modes'\n");
    IF_OK status2=scanf("%s",save_w);
    IF_OK printf("%s\n",save_w);
    IF_OK
      {
        if(strcmp("fresh_hr0_modes",save_w) == 0 )
          par_buf.in_hr0_flag = FRESH;
        else if(strcmp("iserial_hr0_modes",save_w) == 0 )
          par_buf.in_hr0_flag = RELOAD_SERIAL;
	/*
        else if(strcmp("iparallel_h0_prop",save_w) == 0 )
          par_buf.in_hr0_flag = RELOAD_PARALLEL;
	  */
        else
          {
            printf("error in input: %s is not an hr0 input command\n",save_w);
            status++;
          }
        IF_OK
          {
      if(par_buf.in_hr0_flag != FRESH){
            /*read name of file and load it */
            if(prompt==1)printf("enter name of file stem--hr0 in\n");
            status2=scanf("%s",par_buf.in_hr0);
            if(status2 !=1) {
              printf("error in input: hr0 in file  name\n"); status++;
            }
            printf("%s\n",par_buf.in_hr0);
      }
          }
      }

      /* output hr0 */
    IF_OK if(prompt==1)
      printf("h(-r0) mode outfile:\n enter 'forget_hr0_modes, 'serial_hr0_modes',  'parallel_hr0_modes'\n");
    IF_OK status2=scanf("%s",save_w);
    IF_OK printf("%s \n",save_w);
    IF_OK
      {
        if(strcmp("forget_hr0_modes",save_w) == 0 )
          par_buf.out_hr0_flag = FORGET;
        else  if(strcmp("serial_hr0_modes",save_w) == 0 )
          par_buf.out_hr0_flag = SAVE_SERIAL_FM_SC;
	/*
        else if(strcmp("parallel_hr0_modes",save_w) == 0 )
          par_buf.out_hov_flag = SAVE_PARALLEL;
	  */
        else
          {
            printf("error in input: %s is not a hr0 out file command\n",save_w);
            status++;
          }
        IF_OK
          {
      if(par_buf.out_hr0_flag != FORGET){
            /*read name of file and load it */
            if(prompt==1)printf("enter name of scratch file stem--hr0 out\n");
            status2=scanf("%s",par_buf.out_hr0);
            if(status2 !=1) {
              printf("error in input: hr0 file stem name\n"); status++;
            }
            printf("%s\n",par_buf.out_hr0);
      }
          }
      }



      /* input hov */
    IF_OK if(prompt==1)
      printf("propagator hov infile:\n enter 'fresh_hov_modes, 'iserial_hov_modes',\n");
    IF_OK status2=scanf("%s",save_w);
    IF_OK printf("%s\n ",save_w);
    IF_OK
      {
        if(strcmp("fresh_hov_modes",save_w) == 0 )
          par_buf.in_hov_flag = FRESH;
        else  if(strcmp("iserial_hov_modes",save_w) == 0 )
          par_buf.in_hov_flag = RELOAD_SERIAL;
	/*
        else if(strcmp("iparallel_hov_modes",save_w) == 0 )
          par_buf.in_hov_flag = RELOAD_PARALLEL;
	  */
        else
          {
            printf("error in input: %s is not a hov in command\n",save_w);
            status++;
          }
        IF_OK
          {
      if(par_buf.in_hov_flag != FRESH){
            /*read name of file and load it */
            if(prompt==1)printf("enter name of  file --hov in\n");
            status2=scanf("%s",par_buf.in_hov);
            if(status2 !=1) {
              printf("error in input: hov file stem name\n"); status++;
            }
            printf("%s\n",par_buf.in_hov);
      }
          }
      }


      /* output hov */
    IF_OK if(prompt==1)
      printf("H(0) mode outfile:\n enter 'forget_hov_modes, 'serial_hov_modes',  'parallel_hov_modes'\n");
    IF_OK status2=scanf("%s",save_w);
    IF_OK printf("%s\n ",save_w);
    IF_OK
      {
        if(strcmp("forget_hov_modes",save_w) == 0 )
          par_buf.out_hov_flag = FORGET;
        else  if(strcmp("serial_hov_modes",save_w) == 0 )
          par_buf.out_hov_flag = SAVE_SERIAL_FM_SC;
	/*
        else if(strcmp("parallel_hov_modes",save_w) == 0 )
          par_buf.out_hov_flag = SAVE_PARALLEL;
	  */
        else
          {
            printf("error in input: %s is not a hov out file command\n",save_w);
            status++;
          }

        IF_OK
          {
      if(par_buf.out_hov_flag != FORGET){
            /*read name of file and load it */
            if(prompt==1)printf("enter name of scratch file stem--hov out\n");
            status2=scanf("%s",par_buf.out_hov);
            if(status2 !=1) {
              printf("error in input: hov file stem name\n"); status++;
            }
            printf("%s\n",par_buf.out_hov);
      }
          }
      }
#ifdef IMAGISO
         IF_OK status += get_f(prompt,"imag_isospin", &par_buf.delta_iso );
#endif


    /* Number of (random) sources already done 
    IF_OK status += get_i(prompt,"number_sources_done", &par_buf.ndone );
    */
    if( status > 0){ /*printf("enter number_sources_done\n"); */
       par_buf.stopflag=1;} else par_buf.stopflag=0;

  } /* end if(this_node==0) */
  
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if( par_buf.stopflag != 0 ){
    normal_exit(0);
  }
  startflag = par_buf.startflag;
  fixflag = par_buf.fixflag;
  saveflag = par_buf.saveflag;
  for(i=0;i<par_buf.num_masses;i++){
    startflag_w[i] = par_buf.startflag_w[i];
    saveflag_w[i] = par_buf.saveflag_w[i];
#ifdef H0INV
/*
    startflag_w3[i] = par_buf.startflag_w3[i];
*/
    saveflag_w3[i] = par_buf.saveflag_w3[i];
#endif
  }

R0=par_buf.R0;
prec_sign=par_buf.prec_sign;

zolo_min=par_buf.zolo_min;
zolo_min_save=par_buf.zolo_min;
zolo_max=par_buf.zolo_max;
zolo_max_save=par_buf.zolo_max;

scalez=par_buf.scalez;
maxcg_inner=par_buf.maxcg_inner;
resid_inner=par_buf.resid_inner;
resid_inner_save=resid_inner;
resid_inner_h=par_buf.resid_inner_h;
#ifdef MINN
resid_inner_run=resid_inner;
#endif
#ifdef EIG
    Nvecs_h0r0 = par_buf.Nvecs_h0r0 ;
    Nvecs_h0 = par_buf.Nvecs_h0 ;
    Nvecs_hov = par_buf.Nvecs_hov ;
    MaxIter = par_buf.MaxIter ;
    Maxr0Iter = par_buf.Maxr0Iter ;
    Restart = par_buf.Restart ;
    Kiters = par_buf.Kiters ;
    eigenvec_quality=par_buf.eigenvec_qual;
    eigenval_tol_low = par_buf.eigenval_tol ;
    error_decr_low = par_buf.error_decr ;
    eigenval_tol_high = par_buf.eigenval_tol_acc ;
    error_decr_high = par_buf.error_decr_acc ;
#endif
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  num_masses = par_buf.num_masses;
  current_topology=par_buf.topology;
  for(i=0;i<par_buf.num_masses;i++){
    mass[i] = par_buf.mass[i];
    resid[i] = par_buf.resid[i];
    resid_acc[i] = par_buf.resid2[i];
  }
    wqs = par_buf.wqs;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  for(i=0;i<par_buf.num_masses;i++){
    strcpy(startfile_w[i],par_buf.startfile_w[i]);
    strcpy(savefile_w[i],par_buf.savefile_w[i]);
#ifdef H0INV
/*
    strcpy(startfile_w3[i],par_buf.startfile_w3[i]);
*/
    strcpy(savefile_w3[i],par_buf.savefile_w3[i]);
#endif
  }
  strcpy(in_hr0,par_buf.in_hr0);
  in_hr0_flag = par_buf.in_hr0_flag;
  strcpy(out_hr0,par_buf.out_hr0);
  out_hr0_flag = par_buf.out_hr0_flag;
  strcpy(out_hov,par_buf.out_hov);
  out_hov_flag = par_buf.out_hov_flag;
   strcpy(in_hov,par_buf.in_hov);
  in_hov_flag = par_buf.in_hov_flag;
  ndone = par_buf.ndone;

  beta=par_buf.beta;
  u0=par_buf.u0;
  nsmear=par_buf.nsmear;


#ifdef IMAGISO
    delta_iso=par_buf.delta_iso;
#endif


  
  /* Do whatever is needed to get lattice */
  startlat_p = reload_lattice( startflag, startfile );



	  

  return(0);
}
/* allocate all space for fields */
void make_fields() {

int memfield;

/* move here alloc for clov? */

#ifdef NHYP

    FIELD_ALLOC_VEC(gauge_field,su3_matrix,4);
    FIELD_ALLOC_VEC(gauge_field_thin,su3_matrix,4);
    FIELD_ALLOC_VEC(Staple3,su3_matrix,4);
    memfield=13;

#if (SMEAR_LEVEL>1)
    FIELD_ALLOC_MAT_OFFDIAG(Staple2,su3_matrix,4);
    FIELD_ALLOC_MAT_OFFDIAG(hyplink2,su3_matrix,4);
    memfield=37;
#endif

#if (SMEAR_LEVEL==3)
    FIELD_ALLOC_MAT_OFFDIAG(hyplink1,su3_matrix,4);
    FIELD_ALLOC_MAT_OFFDIAG(Staple1,su3_matrix,4);
    memfield=61;
#endif

    FIELD_ALLOC(tempmat_nhyp1,su3_matrix);

    if(this_node==0)printf("Mallocing %.1f MBytes per node for fields\n",
            (double)sites_on_node * memfield * sizeof(su3_matrix)/1e6);

#endif

}
