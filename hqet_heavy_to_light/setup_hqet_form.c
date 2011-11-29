/******** setup_hqet_form.c *********/
/* MIMD version 6 */
/*  set tabstop=2   for easy reading of this file */
/* $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/setup_hqet_form.c,v 1.8 2011/11/29 18:04:41 detar Exp $   ****/
/* MIMD code version 4 */

#include "hqet_light_includes.h"
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
int  setup_hqet_form(void)   
{
int initial_set();
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


/*
 *  This function reads the some basic input
 *  parameters (the lattice dimensions), for standard input. 
 *  The code sends the parameters to all the other nodes,
 *
 *  The top of titles for the program is written out here.
 */

int initial_set()
{
  int prompt,status ;

  /* On node zero, read lattice size and send to others */
  if(mynode()==0)
  {
    /* print banner */
    printf("SU3 Wilson valence fermions;  heavy--->light form factors\n");
    printf("HQET used for the bottom quarks\n") ; 
    printf("MIMD version 6\n");
#ifdef BICG_CLOVER  
   printf("Using CLOVER fermions with the BI-CG algorithm\n") ; 
#else
   printf("Using WILSON fermions with the minimal residual algorithm\n") ; 
#endif
    time_stamp("start");


    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    
    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );

    if(par_buf.nt%2 !=0) 
    {
      printf("nt must be even!! \n"); 
      status++;
    }
    /** The final meson is always kept at the mid point of the lattice ****/
/*    par_buf.tf = par_buf.nt/2 - 1  ;   */
    par_buf.tf = par_buf.nt/2  ;
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
  tf=par_buf.tf;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  
  if(mynode()==0)
    {
      printf("Meson fixed at time slices 0 and %d\n",tf);
    }
  
  return(prompt);

}

/* 
   Node 0 reads in all parameters controlling the
   code, except the dimensions of the lattice (which are
   read in using the function int initial_set()).

   The parameters are then send to all the other nodes.


   argument "prompt" is 1 if prompts are to be given for input	

*/

int readin(int prompt)  
{
  register int i;  
  int j ; 
  int status;
  Real x;
  int ismear ;
  double ts,te ;
  int ikappa , imom ; 
  char descrp[30];
  
  char velfile[80 ] ;
  char momfile[80 ] ; 
  /*** ---------- ..........----------..........  ***/
  
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0)
    {
      printf("\n\n");
      status=0;
      
      IF_OK status += get_i(stdin, prompt,"verbose_flag",&par_buf.verbose_flag);
#ifdef BICG_CLOVER  
      /* Clover coefficient, u0 */
      IF_OK status += get_f(stdin, prompt,"clov_c",&par_buf.clov_c);
      IF_OK status += get_f(stdin, prompt,"u0",&par_buf.u0);
#endif
      
      /* find out what kind of starting lattice to use */
      IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
					    par_buf.startfile );
      
      IF_OK status += get_i(stdin, prompt,"nkap_spectator",&par_buf.no_spectator);
      if(par_buf.no_spectator  >MAX_KAPPA) 
	{
	  printf("no_spectator = %d  cannot be larger than MAX_KAPPA =%d!!! \n",
		 par_buf.no_spectator,MAX_KAPPA );
	  status++;
	}
      
      /* maximum no. of spectator conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter_spectator );
      
      /* maximum no. of spectator conjugate gradient restarts */
      IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart_spectator );
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid_spectator );
      
      /* Get source type */
      IF_OK status += ask_w_quark_source( stdin, prompt, &wallflag, descrp);
      
      /*** Load in the SPECTATOR kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_spectator  ;i++) 
	  { 
	    IF_OK status += get_f(stdin, prompt,"kappa_spectator",
				  &par_buf.kappa_spectator[i]);
	  }
	
      }
      
      /* Load in the source widths and names of the files for the
       *   spectator propagators */
      /* It would be better to code this in w_source. CD */
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
      
      IF_OK {
	for(i=0;i< par_buf.no_spectator ;i++)
	  {
	    IF_OK status += ask_starting_wprop(stdin, prompt,
					       &par_buf.startflag_spectator[i],
					       par_buf.qfile_spectator[i] );
	  }
      }
      
      
      IF_OK status += get_i(stdin, prompt,"nkap_light_zonked",
			    &par_buf.no_zonked_light);
      if(par_buf.no_zonked_light  >MAX_KAPPA) 
	{
	  printf("no_light_zonked = %d  cannot be larger than MAX_KAPPA =%d!\n",
		 par_buf.no_zonked_light,MAX_KAPPA );
	  status++;
	}
      
      /* maximum no. of zonked conjugate gradient iterations */
      IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter_zonked );
      
      /* maximum no. of zonked conjugate gradient restarts */
      IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart_zonked );
      
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			    &par_buf.resid_zonked );
      
      /* Get source type */
      IF_OK status += ask_w_quark_source( stdin, prompt, &wallflag, descrp);

      /*** Load in the zonked light kappa values ****/
      /**** MORE WORK :: this loop should count the kappa values ***/
      IF_OK {
	for(i=0;i< par_buf.no_zonked_light  ;i++) 
	  { 
	    IF_OK status += get_f(stdin,  prompt,"kappa_zonked_light",
				  &par_buf.kappa_zonked_light[i] );
	  }
      }
      
      /* Load in the source widths and names of the files for the
       *   zonked propagators */
      
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK if (prompt==1) 
	printf("enter width(s) r0 as in: source=exp(-(r/r0)^2)\n");

      for(i=0;i<par_buf.no_zonked_light;i++){
	IF_OK status += get_f(stdin, prompt,"r0", &par_buf.wqs_zonked_light[i].r0 );
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
      
      IF_OK {
	for(i=0;i< par_buf.no_zonked_light ;i++)
	  {    	
	    IF_OK status += ask_starting_wprop(stdin, prompt,
					      &par_buf.startflag_zonked[i],
					      par_buf.qfile_zonked[i]);
	  }
      }

      /*** load in the velocites   ****/
      /** 
	load in the name of the file to read the velocites from
	(( The filename is local to this routine ))
	***/
      
      IF_OK status += get_s(stdin, prompt,"vel_file",velfile);
      
      
      /** end of the read of the file name to read the velocities ***/
      
      /**** Load the velocities from disk ******/
      IF_OK
	par_buf.novel = load_velocity_from_disk(&par_buf, 
						velfile,  MAXVEL);
      
      
      /*** load in the operator momentum   ****/
      /** 
	load in the name of the file to read the momentum from
	(( The filename is local to this routine ))
	***/

      IF_OK status += get_s(stdin, prompt,"mom_file",momfile);
      
      /** end of the read of the file name to read the momentum ***/
      IF_OK
	par_buf.no_q_values = load_momentum_from_disk(par_buf.q_momstore, 
						      momfile,  MAXMOM) ;
      
      
      /** 
	load in the name of the file to save the heavy --> light form
	factors to
	***/
      
      IF_OK status += get_s(stdin, prompt,"heavy_light_out",par_buf.heavy_light_out);
      
      /** 
	load in the name of the file to save the  two point functions to
	***/
    
      IF_OK status += get_s(stdin, prompt,"twopt_out",par_buf.twopt_out);
      
      /** 
	load in the name of the file to save the sequential two point
	functions to
	***/
    
      IF_OK status += get_s(stdin, prompt,"seq_out",par_buf.seq_out);

      if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */
  
  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  

  if( par_buf.stopflag != 0 )
    normal_exit(0);


  /**** unpack the parameters *****/
  
  clov_c = par_buf.clov_c ;
  u0 = par_buf.u0 ;
  verbose_flag  = par_buf.verbose_flag  ;

  startflag = par_buf.startflag ;
  strcpy(startfile,par_buf.startfile);
  
  no_spectator  = par_buf.no_spectator  ;
  niter_spectator = par_buf.niter_spectator ;
  nrestart_spectator = par_buf.nrestart_spectator ;
  resid_spectator = par_buf.resid_spectator ;
  
  no_zonked_light = par_buf.no_zonked_light ;
  niter_zonked = par_buf.niter_zonked ;
  nrestart_zonked = par_buf.nrestart_zonked ;
  resid_zonked = par_buf.resid_zonked ;

  for(ikappa =0 ; ikappa < MAX_KAPPA ; ++ikappa )
    {
      kappa_zonked_light[ ikappa ] = par_buf.kappa_zonked_light[ ikappa ] ;
      kappa_spectator[ ikappa ]  =   par_buf.kappa_spectator[ ikappa ] ; 
      wqs_zonked_light[ ikappa ] = par_buf.wqs_zonked_light[ ikappa ];
      init_qs(&wqs_zonked_light[ ikappa ]);
      wqs_zonked_light[ ikappa ].type = par_buf.wqs_zonked_light[ ikappa ].type;
      wqs_spectator[ ikappa ] = par_buf.wqs_spectator [ ikappa ];
      init_qs(&wqs_spectator_light[ ikappa ]);
      wqs_spectator[ ikappa ].type = par_buf.wqs_spectator [ ikappa ].type;
      
      strcpy(qfile_spectator[ikappa] , par_buf.qfile_spectator[ikappa]);
      strcpy(qfile_zonked[ikappa] , par_buf.qfile_zonked[ikappa]);
      startflag_spectator[ikappa] = par_buf.startflag_spectator[ikappa] ;
      startflag_zonked[ikappa] = par_buf.startflag_zonked[ikappa]   ; 
    }
  
  
  for(imom = 0 ; imom < MAXVEL ; ++imom)
    strcpy(hqet_smear_file[imom] , par_buf.hqet_smear_file[imom] );
  
  
  strcpy(heavy_light_out, par_buf.heavy_light_out);
  strcpy(twopt_out, par_buf.twopt_out );
  strcpy(seq_out, par_buf.seq_out );
  
  no_q_values  = par_buf.no_q_values   ;
  quark_type = par_buf.quark_type  ;
  
  novel = par_buf.novel  ;
  
  for(imom = 0 ; imom < MAXMOM ; ++imom)
    {
      q_momstore[imom][0] = par_buf.q_momstore[imom][0] ;  
      q_momstore[imom][1] = par_buf.q_momstore[imom][1] ;  
      q_momstore[imom][1] = par_buf.q_momstore[imom][1] ;  
    }
  
  for(i=0 ; i< novel ;++i)
    for(j=0 ; j < 4 ;++j)
      velocity[ i ][ j ] = par_buf.velocity[ i ][ j ] ;
  
  /*** create/load the gauge configuration ****/
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile )  ;

  return(0);
}

