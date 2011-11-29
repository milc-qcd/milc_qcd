/******** setup_mr.c *********/
/*  set tabstop=2   for easy reading of this file */
/* $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/setup_mr.c,v 1.11 2011/11/29 19:07:09 detar Exp $  ***/
/* MIMD version 7 */
#define IF_OK if(status==0)

#include "w_static_includes.h"
#include "lattice_qdp.h"
#include <string.h>


/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int initial_set();


int setup_h()
{
  int prompt;

  /* print banner, get volume */
  prompt = initial_set();
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor gathers */
  make_nn_gathers();
  /* Create clover structure */
  gen_clov = create_clov();

  node0_printf("Finished setup\n"); fflush(stdout);
  return (prompt);
}


/* SETUP ROUTINES */
int initial_set()
{
  int prompt, status;

  /* On node zero, read lattice size and send to others */
  if (mynode() == 0)
  {
    /* print banner */
    printf("SU3 Wilson valence fermions;  heavy-light static\n");
    printf("Static variational smearing and B_B parameters correlators\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");

    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
    


    if (nt % 2 != 0)
    {
      printf("nt must be even!! \n");
      ++status ; 
    }

    if(status>0) 
      par_buf.stopflag=1; 
    else 
      par_buf.stopflag=0;

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

   
  return (prompt);
}

/* read in parameters and coupling constants	 */
int readin(int prompt)
{
  /* read in parameters for su3 monte carlo	 */
  /* argument "prompt" is 1 if prompts are to be given for input	 */


  register int i,j;		/******loop variable*****/
  int status = 0 ;
  int count;
  Real x = 0.;
  char savebuf[80];
  Real width;
  
  /* On node zero, read parameters and send to all other nodes */
  if (this_node == 0)
  {

    printf("\n\n");

    /* get number of values of kappa to be run */

    IF_OK status +=get_i(stdin, prompt, "nkap",&par_buf.nkap );
    if (par_buf.nkap > MAX_NKAP)
    {
      printf("nkap cannot be larger than MAX_NKAP!!! \n");
     ++status ; 
    }
    /*
     * get the start and end values of spin, color and kappa in their loops 
     */

     IF_OK status += get_i(stdin, prompt, "start_spin",&par_buf.start_spin);
     IF_OK status += get_i(stdin, prompt, "start_color",&par_buf.start_color);
     IF_OK status += get_i(stdin, prompt, "start_kap",&par_buf.start_kap);
     IF_OK status += get_i(stdin, prompt, "end_spin",&par_buf.end_spin);
     IF_OK status += get_i(stdin, prompt, "end_color",&par_buf.end_color);
     IF_OK status += get_i(stdin, prompt, "end_kap",&par_buf.end_kap);


    /* get couplings and broadcast to nodes	 */
    /* beta, kappa */
    IF_OK status +=get_f(stdin, prompt, "beta",&par_buf.beta);

    for (i = 0; i < par_buf.nkap; i++)
    {
      IF_OK status +=get_f(stdin, prompt, "kappa",&par_buf.cappa[i] ) ;
    }
    IF_OK status +=get_f(stdin, prompt, "approx_kappa_c",&par_buf.kappa_c );


    /* maximum no. of conjugate gradient iterations */
    IF_OK status +=get_i(stdin, prompt, "max_iterations",&par_buf.niter );

    /* maximum no. of conjugate gradient iterations */
    IF_OK status +=get_i(stdin, prompt, "max_restarts",&par_buf.nrestart );

    /* error for propagator minimal residue */
    IF_OK status +=get_f(stdin, prompt, "error_for_propagator",&x);
    par_buf.rsqprop = x * x;

    /* flag is 1 to start with psi !=0, 0 if psi=0 */
    IF_OK printf("flag = 1 to start with non-zero light quark prop, 0 if psi=0\n");
    IF_OK status += get_i(stdin, prompt, "flag",&par_buf.flag );

    /* number of hopping parameter steps */
    IF_OK status += get_i(stdin, prompt, "number_hopping_steps",&par_buf.nhop );

    /* Source type */
    IF_OK 
      status += ask_w_quark_source(stdin, prompt, 
				 &par_buf.wqs.type,par_buf.wqs.descrp);

    IF_OK if(par_buf.wqs.type != POINT_WEYL 
       && par_buf.wqs.type != CUTOFF_GAUSSIAN_WEYL)
      {
	printf("setup_mr: ERROR. %s is not a valid source for this application\n",
	       par_buf.wqs.descrp);
	printf("Use only Weyl sources.\n");
	status++;
      }
       

    /* width: psi=exp(-width^(-2)*r*r) */
    IF_OK  printf("width^(-2): source=exp(-width^(-2)*r*r)\n");
    IF_OK status += get_f(stdin, prompt, "width^(-2)",&width );
    IF_OK {
      par_buf.wqs.r0 = 0;
      if(width != 0)par_buf.wqs.r0 = sqrt(1./width);
      else if(par_buf.wqs.type != POINT_WEYL){
	printf("ERROR: If you really want zero width, specify a point_weyl source!\n");
	status++;
      }
    }

    /* get source points (if wall, these are center of wall) and check range */
    IF_OK status += get_i(stdin, prompt, "source_x",&par_buf.wqs.x0 );
    IF_OK status += get_i(stdin, prompt, "source_y",&par_buf.wqs.y0 );
    IF_OK status += get_i(stdin, prompt, "source_z",&par_buf.wqs.z0 );
    IF_OK status += get_i(stdin, prompt, "source_t",&par_buf.wqs.t0 );

    IF_OK 
    {
      if (par_buf.wqs.x0 < 0 || par_buf.wqs.x0 >= nx)
      {
	printf("wrong range for source_x\n");
	++status ;
      }
      if (par_buf.wqs.y0 < 0 || par_buf.wqs.y0 >= ny)
      {
	printf("wrong range for source_y\n");
	++status ;
      }
      if (par_buf.wqs.z0 < 0 || par_buf.wqs.z0 >= nz)
      {
	printf("wrong range for source_z\n");
	++status ;
      }
      if (par_buf.wqs.t0 < 0 || par_buf.wqs.t0 >= nt)
      {
	printf("wrong range for source_t\n");
	++status ;
      }
    } /*** end of IF_OK ****/


    /* find parity of source */
    if ((par_buf.wqs.x0 + par_buf.wqs.y0 + par_buf.wqs.z0 + par_buf.wqs.t0) % 2 == 0)
      par_buf.source_parity = EVEN;
    else
      par_buf.source_parity = ODD;
 
    IF_OK printf("source_parity= %0x\n", par_buf.wqs.parity);

    /*
     * get wall_cutoff and wall_separation and check that latter is
     * reasonable 
     */
    IF_OK status += get_i(stdin, prompt, "wall_cutoff",&par_buf.wqs.wall_cutoff );

    IF_OK status += get_i(stdin, prompt, "wall_separation",&par_buf.wall_separation );
    if (((par_buf.wqs.type == CUTOFF_GAUSSIAN_WEYL) ||
	 (par_buf.wqs.type == CUTOFF_GAUSSIAN))
	&& (nx % par_buf.wall_separation != 0 ||
	    ny % par_buf.wall_separation != 0 ||
	    nz % par_buf.wall_separation != 0))
      {
	IF_OK   printf("Warning: normally, wall_separation should\
 divide nx, ny, and nz\n");
      }

    IF_OK {
      if (prompt != 0)
	printf("enter 'no_extra_sink', or 'extra_sink'\n");
      
      if( scanf("%s", savebuf) != 1 ) 
      {
	printf("Error reading the SINK commands\n") ; 
	++status ; 
      }
    }
    IF_OK {
      if (strcmp("extra_sink", savebuf) == 0)
	{
	  if ((par_buf.wqs.type == POINT) ||
	      (par_buf.wqs.type == POINT_WEYL))
	    {
	      printf("error in input: can't have extra_sink with pt. source\n");
	      ++status ; 
	    }
	  if (nx % 4 != 0 || ny % 4 != 0 || nz % 4 != 0)
	    {
	      printf("nx,ny, or nz not divisible by 4; extra_sink not allowed\n");
	      ++status ; 
	    }
	  par_buf.nchannels = NCHANNELS + 2;
	  printf("computing extra sinks\n");
      } else
	if (strcmp("no_extra_sink", savebuf) == 0)
	  {
	    par_buf.nchannels = NCHANNELS;
	    printf("no extra sinks will be computed\n");
	} else
	  {
	    printf("error in input: sink command %s is invalid\n",savebuf);
	    ++status ; 
	  }


    } /*** end of the IF_OK condition ****/

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

    /*** load in information for static variational calculation *****/
    /* The number of smearing functions to use ***/
    IF_OK get_i(stdin, prompt, "no_smearing_func",&par_buf.nosmear );
    
    IF_OK {
      if( nosmear > MAX_SMEAR ) 
	{
	  printf("Not enough space for the smearing functions %d > %d\n",
		 par_buf.nosmear , MAX_SMEAR);
	  terminate(1);
	}
    }
    
    
    /** 
      load in the name of the file to save the  variational matrix 
      ***/
    
    IF_OK {
      if (prompt==1) 
	{
	  
	  printf("enter the name of the file to save the variational matrix to\n");
	  count=scanf("%s",par_buf.vary_out);
	  
	  if(count !=1) 
	    {
	      printf("error in input: file name read\n"); status++;
	    }
	  else printf("The variational matrix will be saved to %s\n",par_buf.vary_out);
	}
      else 
	{
	  if( scanf("%s",savebuf) != 1 )
	    {
	      printf("error in input: vary_file required\n"); 
	      status++;
	    }
	  else if(strcmp("vary_file",savebuf) != 0 ) 
	    {
	      printf("error in input: vary_file required\n"); 
	      status++;
	    }
	  IF_OK {
	    count=scanf("%s",par_buf.vary_out);
	    if(count !=1) 
	      {
		printf("error in input: file name read\n"); 
		status++;
	      }
	    else printf("The variational matrix will be saved to %s\n",par_buf.vary_out);
	  }
	  
        }
    }
    
    /** end of the read of the file name to save the smearing matrix ***/
    
    
    /***** LOAD the smearing functions *******/
    IF_OK status += get_smearing_funcs_code( savebuf, &par_buf ) ; 
    
    IF_OK {
      if (prompt==1) 
	{
	  printf("enter the name of the file to save the smeared meson operators to\n");
	  count=scanf("%s",par_buf.smear_meson_out);
	  
	  if(count !=1) 
	    {
	      printf("error in input: file name read\n"); status++;
	    }
	}
      else 
	{
	  if( scanf("%s",savebuf) != 1 )
	    {
	      printf("error in input: smear_meson_out required\n"); 
	      status++;
	    }
	  else if(strcmp("smear_meson_out",savebuf) != 0 ) 
	    {
	      printf("error in input: smear_meson_out required\n"); 
	      status++;
	    }
	  IF_OK {
	    count=scanf("%s",par_buf.smear_meson_out);
	    if(count !=1) 
	      {
		printf("error in input: file name read\n"); 
		status++;
	      }
	    else printf("The smeared meson correlators will be to saved to %s\n",
			par_buf.smear_meson_out);
	  }
	  
        }

    }
    
    /** end of the read of the file name to save the smeared  meson correlators ***/
    
    
    for (i = 0; i < par_buf.nkap; i++)
    {

      IF_OK status += ask_starting_wprop( stdin, prompt,
					  &par_buf.startflag_w[i],
					  par_buf.startfile_w[i]);

      IF_OK status += ask_ending_wprop( stdin, prompt,&par_buf.saveflag_w[i],
					par_buf.savefile_w[i]);

      IF_OK 
      {
	if (prompt != 0)
	  printf("save meson propagators:\n enter 'forget', 'save' or 'save_binary\n");
	if( scanf("%s", savebuf) != 1 )
	{
	  IF_OK printf("Error in  reading the meson output file\n") ;
	++status ;
	}
	
	if (strcmp("forget", savebuf) == 0)
	{
	  par_buf.saveflag_m = FORGET_MESON ;
	} else
	  if (strcmp("save", savebuf) == 0)
	  {
	    par_buf.saveflag_m = SAVE_MESON_ASCII;
	    /* read name of file and load it */
	    if (prompt != 0)
	      printf("enter name of file to contain mesons\n");
	    if ( scanf("%s", par_buf.savefile_m[i]) != 1)
	    {
	      printf("error in input: file name read\n");
	      ++status ; 
	    }
	  } else
	    if (strcmp("save_binary", savebuf) == 0)
	    {
	      par_buf.saveflag_m = SAVE_MESON_BINARY;
	      /* read name of file and load it */
	      if (prompt != 0)
		printf("enter name of file to contain mesons\n");
	      if( scanf("%s", par_buf.savefile_m[i])   != 1)
	      {
		printf("error in input: savefile_m file name read\n");
	  ++status ; 
	      }
	    } else
	    {
	      printf("error in input: save_MESON_command is invalid\n");
	      ++status ; 
	    }
      }				/****** end of loop over nkap ******/
      
    } /*** end of IF_OK for meson input ***/

    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    
  }  /* end if(this_node==0) */ 
  
  

  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if( par_buf.stopflag != 0 )
    normal_exit(0);

   

    nkap = par_buf.nkap;
    start_kap = par_buf.start_kap;
    start_spin = par_buf.start_spin;
    start_color = par_buf.start_color;
    end_kap = par_buf.end_kap;
    end_spin = par_buf.end_spin;
    end_color = par_buf.end_color;
    startflag = par_buf.startflag;
    saveflag_m = par_buf.saveflag_m;
    niter = par_buf.niter;
    nrestart = par_buf.nrestart;
    nhop = par_buf.nhop;
    flag = par_buf.flag;
    wqs = par_buf.wqs;
    init_qs(&wqs);
    wqs.type = par_buf.wqs.type;
    source_parity = par_buf.source_parity;
    rsqmin = par_buf.rsqmin;
    rsqprop = par_buf.rsqprop;
    beta = par_buf.beta;
    kappa_c = par_buf.kappa_c;
    wall_separation = par_buf.wall_separation;
    nchannels = par_buf.nchannels;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);
  fixflag = par_buf.fixflag ; 
  saveflag = par_buf.saveflag  ;

  /*** static variational calculation ********/
  strcpy(vary_out,par_buf.vary_out);
  strcpy(smear_meson_out,par_buf.smear_meson_out);
  nosmear = par_buf.nosmear;
  for(i=0 ; i < nosmear ; ++i )
    {
      strcpy(smearfile_in[i],par_buf.smearfile_in[i]);
      for(j=0 ; j < 5 ; ++j)
	smear_code[i][j] = par_buf.smear_code[i][j];
    }
  /********************************************/

    for (i = 0; i < nkap; i++)
    {
      cappa[i] = par_buf.cappa[i];
      strcpy(startfile_w[i], par_buf.startfile_w[i]);
      strcpy(savefile_w[i], par_buf.savefile_w[i]);
      strcpy(savefile_m[i], par_buf.savefile_m[i]);

      startflag_w[i] = par_buf.startflag_w[i];
      saveflag_w[i] = par_buf.saveflag_w[i];

    }






  /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE ){
      startlat_p = (gauge_file *) reload_lattice( startflag, startfile );  
      invalidate_this_clov(gen_clov);
    }

  return 0 ;
}

