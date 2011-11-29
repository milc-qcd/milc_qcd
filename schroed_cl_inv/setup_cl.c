/******** setup_cl.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...
   01/30/97 ANSI prototyping U.M.H.
   */

#include "schroed_cl_includes.h"
#include <string.h>

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup_cl()   {
  int initial_set();
  int prompt;

  /* print banner, get volume */
  prompt=initial_set();
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
  /* On node zero, read lattice size and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 Wilson valence fermion Schroedinger functional\n");
    printf("MIMD version 4\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    
    status=get_prompt(stdin, &prompt);
  
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
  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/

  int status;
  int i;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){

    printf("\n\n");
    status=0;

    /* Number of kappas */
    IF_OK status += get_i(stdin, prompt,"number_of_kappas", &par_buf.num_kap );
    if( par_buf.num_kap>MAX_KAP ){
      printf("num_kap = %d must be <= %d!\n", par_buf.num_kap, MAX_KAP);
      status++;
    }

    /* boundary condition flag */
    IF_OK status += get_i(stdin, prompt, "bc_flag", &par_buf.bc_flag);

    /* Number of APE smearings */
    IF_OK status += get_i(stdin, prompt, "num_smear", &par_buf.num_smear);

    /* APE smearing parameter (Boulder convention) */
    IF_OK status += get_f(stdin, prompt, "alpha", &par_buf.alpha);

    /* To be save initialize the following to zero */
    for(i=0;i<MAX_KAP;i++){
      kap[i] = 0.0;
      resid[i] = 0.0;
    }

    for(i=0;i<par_buf.num_kap;i++){
      IF_OK status += get_f(stdin, prompt,"kappa", &par_buf.kap[i] );
    }

    /* Clover coefficient */
    IF_OK status += get_f(stdin, prompt,"clov_c", &par_buf.clov_c );

    /* fermion phase factors */
    IF_OK status += get_f(stdin, prompt,"ferm_phases[0]", &par_buf.ferm_phas[0] );
    IF_OK status += get_f(stdin, prompt,"ferm_phases[1]", &par_buf.ferm_phas[1] );
    IF_OK status += get_f(stdin, prompt,"ferm_phases[2]", &par_buf.ferm_phas[2] );

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin, prompt,"max_cg_iterations", &par_buf.niter );

    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin, prompt,"max_cg_restarts", &par_buf.nrestart );

    /* error for propagator conjugate gradient */
    for(i=0;i<par_buf.num_kap;i++){
      IF_OK status += get_f(stdin, prompt,"error_for_propagator", &par_buf.resid[i] );
    }

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
	par_buf.startfile );

    /* send parameter structure */
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));

  if( par_buf.stopflag != 0 )
    normal_exit(0);

  startflag = par_buf.startflag;
  bc_flag = par_buf.bc_flag;
  num_kap = par_buf.num_kap;
  clov_c = par_buf.clov_c;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  num_smear = par_buf.num_smear;
  alpha = par_buf.alpha;
  for(i=0;i<par_buf.num_kap;i++){
    kap[i] = par_buf.kap[i];
    resid[i] = par_buf.resid[i];
  }
  for(i=0;i<3;i++){
    ferm_phases[i] = par_buf.ferm_phas[i];
  }
  strcpy(startfile,par_buf.startfile);

  for(i=0;i<par_buf.num_kap;i++){
//    wqs[i].c_src = NULL;
//    wqs[i].wv_src = NULL;
    wqs[i].type = 0;
//    wqs[i].x0 = 0;
//    wqs[i].y0 = 0;
//    wqs[i].z0 = 0;
//    wqs[i].t0 = 0;
    strcpy(wqs[i].descrp,"Schroedinger wall source");
  }

  beta = 1e20;	/* Only needed in io_helpers for setting boundary fields */
  c_t11 = 0;	/* Only needed in io_helpers for setting boundary fields */
		/* These boundary fields are actually never used here */

  /* Do whatever is needed to get lattice */
  if( startflag != CONTINUE ){
    startlat_p = reload_lattice( startflag, startfile );
    invalidate_this_clov(gen_clov);
  }

  /* put in fermion phases */
  if( startflag != CONTINUE) do_phases();
  return(0);
}

void do_phases()  {
/* Multiply gauge and boundary fields with fermionic phase factors */

  register int i,j,k,dir;
  register site *sit;
  register Real phr,phi,ddr=0.0,ddi;

  for(dir=XUP;dir<=ZUP;dir++){
    if(ferm_phases[dir] != 0.0){
      switch(dir){
	case(XUP):	ddr = ferm_phases[dir] / (Real)nx;	break;
	case(YUP):	ddr = ferm_phases[dir] / (Real)ny;	break;
	case(ZUP):	ddr = ferm_phases[dir] / (Real)nz;	break;
      }
      phr = cos((double)ddr);
      phi = sin((double)ddr);

      FORALLSITES(i,sit){
	for(j=0; j<3; j++) for(k=0; k<3; k++)  {
	  ddr = sit->link[dir].e[j][k].real;
	  ddi = sit->link[dir].e[j][k].imag;
	  sit->link[dir].e[j][k].real = phr*ddr - phi*ddi;
	  sit->link[dir].e[j][k].imag = phr*ddi + phi*ddr;
	  ddr = sit->boundary[dir].e[j][k].real;
	  ddi = sit->boundary[dir].e[j][k].imag;
	  sit->boundary[dir].e[j][k].real = phr*ddr - phi*ddi;
	  sit->boundary[dir].e[j][k].imag = phr*ddi + phi*ddr;
	}
      }
    }
  }
  node0_printf("Fermion phases set\n");
}

