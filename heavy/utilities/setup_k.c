/******** setup_k.c *********/
/* MIMD version 6 */
#define IF_OK if(status==0)

#include "w_sum_includes.h"
#include <string.h>

int setup_k()
{
  int prompt;

  /* print banner, get all parameters except kappa_h */
  prompt = initial_set();

  return (prompt);
}

/* SETUP ROUTINES */
int initial_set()
{
  int prompt, status = 0 ;
  char savebuf[20];

  /* On node zero, read parameters */
  if (mynode() == 0)
  {
    /* print banner */
    printf("SU3 Wilson valence fermions;  summing heavy-light hopping\n");
    printf("MIMD version 6\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");

    printf("type 0 for no prompts, 1 for prompts or 2 for list of prompts\n");
	 
    status = get_prompt(stdin, &prompt);
    if (status != 0)
    {
      printf("error in input: initial prompt\n");
      return (-1);
    }
    IF_OK 
    {
      status += get_i(stdin, prompt, "nt",&nt);
      printf("number of time slices = %d\n", nt);
    }

    IF_OK 
    {
      status += get_i(stdin, prompt, "source_t",&source_t );
      printf("source time slice = %d\n", source_t);
    }

    /* get couplings */
    /* beta, kappa */
    IF_OK status += get_f(stdin, prompt, "beta",&beta );
    IF_OK status += get_f(stdin, prompt, "kappa",&kappa );
    IF_OK status += get_f(stdin, prompt, "approx_kappa_c",&kappa_c );

    /* number of hopping parameter steps */
    IF_OK status += get_i(stdin, prompt, "number_hopping_steps",&nhop);


    IF_OK 
    {
      if (prompt != 0)       printf("enter 'no_extra_sink', or 'extra_sink'\n");
      if( scanf("%s", savebuf)  != 1 )
      {
	printf("There was an error reading the sink parameters\n"); 
	++status ; 
      }
    }
      
    IF_OK {
      
      if (strcmp("extra_sink", savebuf) == 0)
	{
	  nchannels = NCHANNELS + 2;
	  if (this_node == 0)
	    printf("computing extra sinks\n");
      } else
	if (strcmp("no_extra_sink", savebuf) == 0)
	  {
	    nchannels = NCHANNELS;
	    if (this_node == 0)
	      printf("no extra sinks will be computed\n");
	} else
	  {
	    printf("error in input: sink command is invalid\n");
	    ++status ; 
	  }
      
    } /*** end of IF_OK ***/

    if (prompt != 0)
      printf("loading hopping expansion:\n enter 'reload', or 'reload_binary'\n");
    if( scanf("%s", savebuf) != 1 ) 
    {
      printf("ERROR during the IO option for the input file\n"); 
      ++status ; 
    }

    if (strcmp("reload", savebuf) == 0)
    {
      startflag_m = RELOAD_ASCII;
      /* read name of file and load it */
      if (prompt != 0)
	printf("enter name of file containing hopping expan.\n");
      if ( scanf("%s", startfile_m) != 1)
      {
	printf("error in input: file name read\n");
	++status ; 
      }
    } else
    if (strcmp("reload_binary", savebuf) == 0)
    {
      startflag_m = RELOAD_BINARY;
      /* read name of file and load it */
      if (prompt != 0)
	printf("enter name of file containing hopping expan.\n");
      if ( scanf("%s", startfile_m) != 1)
      {
	printf("error in input: file name read\n");
	++status ; 
      }
    } else
    {
      printf("error in input: meson file command is invalid\n");
      ++status ; 
    }


    IF_OK 
    {
      if (prompt != 0)
	printf("save results for given kappa_h\n enter 'save' or 'save_binary\n");
      if(  scanf("%s", savebuf) != 1 ) 
      {
	printf("Error reading IO option for the results\n") ; 
	++status ; 
      }

      if (strcmp("save", savebuf) == 0)
      {
	saveflag_k = SAVE_ASCII;
	/* read name of file and load it */
	if (prompt != 0)
	  printf("enter generic name of file to contain kappa_h results\n");
	if (scanf("%s", savefile)   != 1)
	{
	  printf("error in input: file name read\n");
	  ++status ;
	}
      } else
	if (strcmp("reload_binary", savebuf) == 0)
	{
	  printf(" binary output of kappa_h results not implemented\n");
	  ++status ;
	} else
	{
	  printf("error in input: save_command is invalid\n");
	  ++status ;
	}

    } /*** end of IF_OK ***/

    /*
     * writeflag=0 for just writing out sum , 1 for writing out all partial
     * sums 
     */
    IF_OK printf("writeflag=0 for final result only, 1 for partial sums\n");
    IF_OK status += get_i(stdin, prompt, "writeflag",&writeflag);
    IF_OK
    {
      if (writeflag == WRITELAST)
	printf(" writing result only \n");
      if (writeflag == WRITEALL)
	printf(" writing partial results\n");
    }


  }
   /* end if(mynode()==0) */ 
  else
  {				/* nodes other than 0 just exit */
    time_stamp("exit");
    exit(0);
  }

  this_node = mynode();
  return (status);
}




int readin(int prompt)
{
  /* read in parameters for su3 monte carlo	 */
  /* argument "prompt" is 1 if prompts are to be given for input	 */
  int status ; 
  int kappa_h_int;

  /* On node zero, read parameters */
  if (mynode() == 0)
  {

    printf("\n\n");

    /* get next kappa_h 	 */
    status = get_f(stdin, prompt, "kappa_h",&kappa_h);
    if (kappa_h < -0.0001 || status != 0)
      return 1;			/* stop for negative kappas */
    if (fabs(kappa_h - kappa) < EPS)
    {
      printf("heavy_kappa = light_kappa: using light quark\n");
      printf("\tpropagator twice; hopping expansion ignored.\n");
      printf("\tTo override use any kappa_h > 1.0 \n");
    }
    /* make output file name for this kappa_h */
    kappa_h_int = (int) (1000 * (kappa_h + EPS));
    if (kappa_h_int != 0 && kappa_h <= 1.0)
    {
      printf("heavy_kappa= %f\n", (double) kappa_h);
      if (1000.* kappa_h - (Real) kappa_h_int < 1000.* EPS)
      {				/* 3 digit kappa */
	sprintf(savefile_k, "%s.%03d", savefile, kappa_h_int);
      } else
      {
	kappa_h_int = (int) (10000 * (kappa_h + EPS));
	if (10000.* kappa_h - (Real) kappa_h_int < 10000.* EPS)
	{			/* 4 digit kappa */
	  sprintf(savefile_k, "%s.%04d", savefile, kappa_h_int);
	} else
	{
	  kappa_h_int = (int) (100000 * (kappa_h + EPS));
	  if (100000.* kappa_h - (Real) kappa_h_int < 100000.* EPS)
	  {			/* 5 digit kappa */
	    sprintf(savefile_k, "%s.%05d", savefile, kappa_h_int);
	  } else
	  {
	    printf("too many digits in kappa_h\n");
	    terminate(1);
	  }
	}
      }
    } else
    if (kappa_h_int == 0)
    {
      sprintf(savefile_k, "%s.static", savefile);
      printf("heavy_kappa=0; static assumed \n");
      printf("static result will include a heavy particle \n");
      printf("\tpropagating forward in time AND a heavy \n");
      printf("\tantiparticle propagating backward in time\n");
    } else
    {
      kappa_h_int = (int) (1000 * kappa);
      sprintf(savefile_k, "%s.%03d_hop", savefile, kappa_h_int);
      printf("degenerate meson but 2nd quark done by hopping expansion\n");
    }



  }				/* end if(mynode()==0) */
  return (0);
}
