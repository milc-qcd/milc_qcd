/******** setup.c *********/
/* MIMD version 7 */

#define IF_OK if(status==0)

#include "smooth_inst_includes.h"
#include <string.h>

int ask_ending_topo( int prompt, int *flag, char *filename );

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int setup(void)
{
   int prompt;

   /* print banner, get volume, nflavors, seed */
   prompt=initial_set();

   /* initialize the node random number generator */
/*
   initialize_prn(&node_prn,iseed,volume+mynode());
*/

   /* Initialize the layout functions, which decide where sites live */
   setup_layout();

   /* allocate space for lattice, set up coordinate fields */
   make_lattice();

   /* set up neighbor pointers and comlink structures */
   make_nn_gathers();

   return prompt;
}


/* SETUP ROUTINES */
int initial_set(void)
{
   int prompt,status;
   /* On node zero, read lattice size, seed, nflavors and send to others */
   if(mynode()==0)
   {
      /* print banner */
      printf("SU3 Smoothing for Instantons\n");
      printf("MIMD version 7\n");
      printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
      time_stamp("start");
      status=get_prompt(stdin, &prompt);
      IF_OK status += get_i(stdin,  prompt,"nx", &par_buf.nx );
      IF_OK status += get_i(stdin,  prompt,"ny", &par_buf.ny );
      IF_OK status += get_i(stdin,  prompt,"nz", &par_buf.nz );
      IF_OK status += get_i(stdin,  prompt,"nt", &par_buf.nt );

      if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
   } /* end if(mynode()==0) */

   /* Node 0 broadcasts parameter buffer to all other nodes */
   broadcast_bytes((char *)&par_buf, sizeof(par_buf));

   if( par_buf.stopflag != 0 )
       normal_exit(0);

   nx=par_buf.nx;
   ny=par_buf.ny;
   nz=par_buf.nz;
   nt=par_buf.nt;

   this_node = mynode();
   number_of_nodes = numnodes();
   volume=nx*ny*nz*nt;
   return prompt;
}


/* read in parameters and coupling constants */
int readin(const int prompt)
{
   /* read in parameters for su3 monte carlo */
   /* argument "prompt" is 1 if prompts are to be given for input */

   int status;
   char savebuf[128];

   /* On node zero, read parameters and send to all other nodes */
   if(this_node==0)
   {
      printf("\n\n");
      status=0;

#ifdef HYP
      /* Weights for hypercube blocking */
      IF_OK status += get_f(stdin,  prompt, "alpha", &par_buf.alpha );
      IF_OK status += get_f(stdin,  prompt, "alpha2", &par_buf.alpha2 );
      IF_OK status += get_f(stdin,  prompt, "alpha3", &par_buf.alpha3 );
#else      
      /* Weight for simple APE blocking */
      IF_OK status += get_f(stdin,  prompt, "ape_weight", &par_buf.ape_weight );
#endif

      /* sweeps */
      IF_OK status += get_i(stdin,  prompt, "sweeps", &par_buf.sweeps );

      /* trajectories between propagator measurements */
      IF_OK status +=
        get_i(stdin,  prompt, "sweeps_between_meas", &par_buf.measinterval );

      /* number of hits in SU(3) projection */
      IF_OK status += get_i(stdin,  prompt, "hits_per_sweep", &par_buf.hits );

      /* find out what kind of starting lattice to use */
      IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
                                           par_buf.startfile );

      /* what kind of gauge fixing */
      IF_OK if (prompt==1)
      {
         printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
      }
      IF_OK scanf("%s",savebuf);
      IF_OK { 
	if(strcmp( "coulomb_gauge_fix",savebuf) == 0 )
	  {
	    par_buf.fixflag = COULOMB_GAUGE_FIX;
	    printf("fixing to coulomb gauge\n");
	  }
	else if(strcmp("no_gauge_fix",savebuf) == 0 )
	  {
	    par_buf.fixflag = NO_GAUGE_FIX;
	    printf("NOT fixing the gauge\n");
	  }
	else
	  {
	    printf("error in input: fixing_command is invalid\n"); status++;
	  }
      }

      /* where to store topological density */
      IF_OK status += ask_ending_topo( prompt, &(par_buf.savetopoflag),
				       par_buf.topofile);

      /* find out what to do with lattice at end */
      IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
                                         par_buf.savefile );
      IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				    par_buf.stringLFN );

      /* send parameter structure */
      if( status > 0) par_buf.stopflag=1; else par_buf.stopflag=0;
   } /* end if(this_node==0) */

   /* Node 0 broadcasts parameter buffer to all other nodes */
   broadcast_bytes((char *)&par_buf, sizeof(par_buf));

   if( par_buf.stopflag != 0 )
     normal_exit(0);


#ifdef HYP
   alpha = par_buf.alpha;
   alpha2 = par_buf.alpha2;
   alpha3 = par_buf.alpha3;
#else
   ape_weight = par_buf.ape_weight;
#endif
   sweeps = par_buf.sweeps;
   hits = par_buf.hits;
   measinterval = par_buf.measinterval;
   startflag = par_buf.startflag;
   fixflag = par_buf.fixflag;
   saveflag = par_buf.saveflag;
   savetopoflag = par_buf.savetopoflag;
   strcpy(startfile,par_buf.startfile);
   strcpy(savefile,par_buf.savefile);
   strcpy(stringLFN, par_buf.stringLFN);
   strcpy(topofile,par_buf.topofile);

   /* Do whatever is needed to get lattice */
   if( startflag != CONTINUE )
     startlat_p = reload_lattice( startflag, startfile );

   return 0;

} /*readin()*/



/* find out what do to with topo file at end, and file name if
   necessary.  This routine is only called by node 0.
*/
int ask_ending_topo( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt==1) printf(
        "'forget_topo' topo at end or 'save_topo',\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_ending_topo: ERROR IN INPUT: ending topo command\n");
        return(1);
    }
    if(strcmp("save_topo",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("forget_topo",savebuf) == 0 ) {
        *flag=FORGET;
    }
    else {
      printf("ask_ending_topo: ERROR IN INPUT: %s is not a save topo command\n",savebuf);
      return(1);
    }


    if( *flag != FORGET ){
        if(prompt==1)printf("enter filename\n");
        status=scanf("%s",filename);
        if(status !=1){
    	    printf("ask_ending_topo: ERROR IN INPUT: save filename\n"); return(1);
        }
        printf("topo density to be saved in %s\n",filename);
    }
    return(0);
}
