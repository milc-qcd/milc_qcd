/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)

/* Modifications ...

* 5/30/07 Created from setup_cl.c */

//  $Log: setup.c,v $
//  Revision 1.6  2009/05/31 02:00:56  detar
//  Fix "continue" and NULL startlat_p bug in clover_info.c and setup*.c
//
//  Revision 1.5  2009/04/05 16:43:07  detar
//  Utilities for open meson propagators
//
//  Revision 1.4  2008/04/18 15:12:11  detar
//  Revise metadata format for correlator files.
//
//  Revision 1.3  2008/03/28 15:42:04  detar
//  Switch to indexing by momentum-operator pair.  Add another test case.
//
//  Revision 1.2  2007/11/09 16:15:58  detar
//  Support changes in propagator I/O
//
//  Revision 1.1  2007/10/07 20:02:32  detar
//  Add new application.  Generarlizes clover_invert.
//
//


#include "cl_inv_includes.h"
#include "lattice_qdp.h"
#include <string.h>
int initial_set();
void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
void make_3n_gathers();

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
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* set up K-S phase vectors, boundary conditions */
  phaseset();

#ifdef HAVE_QDP
  {
    int i;
    for(i=0; i<4; ++i) {
      shiftdirs[i] = QDP_neighbor[i];
      shiftdirs[i+4] = neighbor3[i];
    }
    for(i=0; i<8; ++i) {
      shiftfwd[i] = QDP_forward;
      shiftbck[i] = QDP_backward;
    }
  }
#endif
  
  return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
  int prompt,status;
#ifdef FIX_NODE_GEOM
  int i;
#endif
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
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   param.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_s(stdin, prompt,"job_id",param.job_id);
    
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
  
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif

  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  return(prompt);
}

/* Forward declarations */
static int hash_corr_label(char meson_label[MAX_CORR][MAX_MESON_LABEL], 
			   char mom_label[MAX_CORR][MAX_MOM_LABEL],
			   char *meson_label_in, char *mom_label_in, int *n);
char decode_parity(char *parity_label_in);
double decode_factor(char *factor_op, double factor);


/* read in parameters and coupling constants	*/
int readin(int prompt) {
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  char savebuf[128];
  int i;
  int ipair;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){
    
    printf("\n\n");
    status=0;

    /*------------------------------------------------------------*/
    /* Gauge configuration section                                */
    /*------------------------------------------------------------*/

    IF_OK status += ask_starting_lattice(stdin,  prompt, &param.startflag,
	param.startfile );
    IF_OK status += get_f(stdin, prompt,"u0", &param.u0 );

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

    /* APE smearing parameters (if needed) */
    /* Zero suppresses APE smearing */
    IF_OK status += get_f(stdin, prompt, "staple_weight", 
			  &param.staple_weight);
    IF_OK status += get_i(stdin, prompt, "ape_iter",
			  &param.ape_iter);
    
    /*------------------------------------------------------------*/
    /* Propagator inversion control                               */
    /*------------------------------------------------------------*/

    /* maximum no. of conjugate gradient iterations */
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &param.qic.max );
    
    /* maximum no. of conjugate gradient restarts */
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &param.qic.nrestart );
    
    /* error for clover propagator conjugate gradient */
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			  &param.qic.resid );
    IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			  &param.qic.relresid );
    /* Precision fixed to prevailing precision for now */
    param.qic.prec = PRECISION;
    param.qic.parity = EVENANDODD;
    

    /*------------------------------------------------------------*/
    /* Propagators and their sources                              */
    /*------------------------------------------------------------*/

    /* Number of propagators */
    IF_OK status += get_i(stdin,prompt,"number_of_propagators", 
			  &param.num_prop );
    if( param.num_prop>MAX_PROP ){
      printf("num_prop = %d must be <= %d!\n", param.num_prop, MAX_PROP);
      status++;
    }

    IF_OK for(i = 0; i < param.num_prop; i++){
    
      /* Propagator parameters */
      IF_OK status += get_s(stdin, prompt,"propagator_type", savebuf );
      IF_OK {
	if(strcmp(savebuf,"clover") == 0)param.prop_type[i] = CLOVER_TYPE;
	else if(strcmp(savebuf,"KS") == 0)param.prop_type[i] = KS_TYPE;
	else {
	  printf("Unknown quark type %s\n",savebuf);
	  status++;
	}
      }
      if(param.prop_type[i] == CLOVER_TYPE){
	IF_OK status += get_s(stdin, prompt,"kappa", param.kappa_label[i]);
	IF_OK param.dcp[i].Kappa = atof(param.kappa_label[i]);
	IF_OK status += get_f(stdin, prompt,"clov_c", &param.dcp[i].Clov_c );
	param.dcp[i].U0 = param.u0;
	IF_OK status += get_s(stdin, prompt,"check", savebuf);
	IF_OK {
	  /* Should we be checking the propagator by running the solver? */
	  if(strcmp(savebuf,"no") == 0)param.check[i] = 0;
	  else param.check[i] = 1;
	}
	IF_OK status += ask_starting_wprop( stdin, prompt, 
					    &param.startflag_w[i],
					    param.startfile_w[i]);
	
	IF_OK status += ask_ending_wprop( stdin, prompt, &param.saveflag_w[i],
					  param.savefile_w[i]);
	
	/* Get source type */
	IF_OK init_wqs(&param.src_wqs[i]);
	IF_OK status += get_w_quark_source( stdin, prompt, &param.src_wqs[i]);
	
      } else {  /* KS_TYPE */
	IF_OK status += get_s(stdin, prompt,"mass", param.mass_label[i] );
	IF_OK param.ksp[i].mass = atof(param.mass_label[i]);
	param.ksp[i].u0 = u0;
	IF_OK status += get_s(stdin, prompt,"check", savebuf);
	IF_OK {
	  /* Should we be checking the propagator by running the solver? */
	  if(strcmp(savebuf,"no") == 0)param.check[i] = 0;
	  else param.check[i] = 1;
	}
	IF_OK status += ask_starting_ksprop( stdin, prompt, 
					     &param.startflag_ks[i],
					     param.startfile_ks[i]);
	
	IF_OK status += ask_ending_ksprop( stdin, prompt, &param.saveflag_ks[i],
					   param.savefile_ks[i]);
	
	/* Get source type */
	IF_OK init_ksqs(&param.src_ksqs[i]);
	IF_OK status += get_ks_quark_source( stdin, prompt, &param.src_ksqs[i]);
      }
    }

    /*------------------------------------------------------------*/
    /* Quarks                                                     */
    /*------------------------------------------------------------*/

    /* Number of quarks */
    IF_OK status += get_i(stdin,prompt,"number_of_quarks", 
			  &param.num_qk );
    if( param.num_qk>MAX_QK ){
      printf("num_qk = %d must be <= %d!\n", param.num_qk, MAX_QK);
      status++;
    }

    IF_OK for(i = 0; i < param.num_qk; i++){
      char *check_tag;
      /* Get the propagator that we act on with the sink operator to
	 form the "quark" field used in the correlator.  It might be a
	 raw "propagator" or it might be a previously constructed
	 "quark" field */
      /* First we look for the type */
      IF_OK {
	if(prompt!=0)printf("enter 'propagator' or 'quark'\n");
	check_tag = get_next_tag(stdin, "propagator or quark", "readin");
	printf("%s ",check_tag);
	if(strcmp(check_tag,"propagator") == 0)
	  param.parent_type[i] = PROP_TYPE;
	else if(strcmp(check_tag,"quark") == 0)
	  param.parent_type[i] = QUARK_TYPE;
	else{
	  printf("\nError: expected 'propagator' or 'quark'\n");
	  status++;
	}
      }
      /* Next we get its index */
      IF_OK {
	if(prompt!=0)printf("enter the index\n");
	if(scanf("%d",&param.prop_for_qk[i]) != 1){
	  printf("\nFormat error reading index\n");
	  status++;
	}
	else{
	  printf("%d\n",param.prop_for_qk[i]);
	  if(param.parent_type[i] == PROP_TYPE && 
	     param.prop_for_qk[i] >= param.num_prop){
	    printf("Propagator index must be less than %d\n",
		   param.num_prop);
	    status++;
	  }
	  else if(param.parent_type[i] == QUARK_TYPE && 
		  param.prop_for_qk[i] >= i){
	    printf("Quark index must be less than %d here\n",i);
	    status++;
	  }
	}
      }
      /* Get sink operator attributes */
      IF_OK init_wqs(&param.snk_wqs[i]);
      IF_OK status += get_w_quark_sink( stdin, prompt, &param.snk_wqs[i]);
      IF_OK status += ask_ending_wprop( stdin, prompt, &param.saveflag_q[i],
					param.savefile_q[i]);
	
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
    
    IF_OK for(ipair = 0; ipair < param.num_pair; ipair++){
      char request_buf[MAX_SPECTRUM_REQUEST];
      
      /* Which quarks in the pair? */
      
      IF_OK status += get_vi(stdin, prompt, "pair", param.qkpair[ipair], 2);
      
      IF_OK {
	int j;
	for(j = 0; j < 2; j++){
	  if(param.qkpair[ipair][j] < 0 || param.qkpair[ipair][j] >= param.num_qk){
	    printf("Quark index %d must be in [0,%d]\n",
		   param.qkpair[ipair][j],param.num_qk-1);
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
	/* Mesons */
	if(strstr(request_buf,",meson,") != NULL)
	  param.do_meson_spect[ipair] = 1;
	else param.do_meson_spect[ipair] = 0;
	
	/* Baryons */
	if(strstr(request_buf,",baryon,") != NULL)
	  param.do_baryon_spect[ipair] = 1;
	else param.do_baryon_spect[ipair] = 0;

	/* Any meson or baryon */
	param.do_closed_hadron[ipair]  = 
	  param.do_meson_spect[ipair] |
	  param.do_baryon_spect[ipair];

	/* Open meson */
	if(strstr(request_buf,",open_meson,") != NULL)
	  param.do_open_meson[ipair] = 1;
	else param.do_open_meson[ipair] = 0;

	/* We need some correlator specification */
	if(!param.do_open_meson[ipair] && !param.do_closed_hadron[ipair]){
	  printf("Unrecognized spectrum request\n");
	  status++;
	}

	/* But we can't do both closed and open hadrons in the same pairing */
	if(param.do_open_meson[ipair] && param.do_closed_hadron[ipair]){
	  printf("Can't combine open_meson with other options ");
	  printf("in the same pairing\n");
	  status++;
	}
      }
      
      /* What file for the resulting correlators? */
      
      IF_OK status += ask_corr_file( stdin, prompt, &param.saveflag_c[ipair],
				     param.savefile_c[ipair]);
      /* Correlator values are printed with t relative to the t_offset */
      /* FT phases are computed with x,y,z relative to r_offset */
      IF_OK {
	int r[4];
	status += get_vi(stdin,prompt, "r_offset", r, 4);
	param.r_offset[ipair][0] = r[0];
	param.r_offset[ipair][1] = r[1];
	param.r_offset[ipair][2] = r[2];
	param.t_offset[ipair]    = r[3];
      }
      
      /*------------------------------------------------------------*/
      /* Table of correlator combinations actually needed           */
      /*------------------------------------------------------------*/
      
      /* Number of correlators */
      IF_OK status += get_i(stdin,prompt,"number_of_correlators", 
			    &param.num_corr[ipair] );
      IF_OK {
	if( param.num_corr[ipair]>MAX_CORR ){
	  printf("num_corr = %d must be <= %d!\n", param.num_corr[ipair], 
		 MAX_CORR);
	  status++;
	}
      }
      
      /* Sample format for correlator line:
	 correlator  A1_P5 p200 -i * 1 G5 G5X 2 0 0 E E E */
      
      param.num_corr_report[ipair] = 0;
      IF_OK for(i = 0; i < param.num_corr[ipair]; i++){
	int ok,m;
	char meson_label_in[MAX_MESON_LABEL], mom_label_in[MAX_MOM_LABEL],
	  gam_src_lab[MAXGAMMA], gam_snk_lab[MAXGAMMA], phase_lab[4],
	  factor_op[2], parity_x_in[3], parity_y_in[3], parity_z_in[3];
	double factor;
	
	IF_OK status += get_sn(stdin, prompt, "correlator", meson_label_in);
	
	/* Read the momentum label next */
	IF_OK {
	  ok = scanf("%s", mom_label_in);
	  if(ok != 1){
	    printf("Error reading momentum label\n");
	    status++;
	  }
	}
	IF_OK printf(" %s", mom_label_in);
	
	/* Add to hash table */
	m = hash_corr_label(param.meson_label[ipair], param.mom_label[ipair],
			    meson_label_in, mom_label_in,
			    &param.num_corr_report[ipair]);
	param.corr_index[ipair][i] = m;
	if(m < 0)status++;
	
	/* phase, op, factor, and gamma matrix specification */
	IF_OK {
	  ok = scanf("%s %s %lf %s %s\n",phase_lab,factor_op,&factor,
		     gam_src_lab,gam_snk_lab);
	  if(ok != 5){
	    printf("\nError reading phase, factor, and gammas\n");
	    status++;
	  }
	  else {
	    printf(" %3s %1s %6f %3s %3s",phase_lab,factor_op,factor,
		   gam_src_lab,gam_snk_lab);
	  }
	}
	
	/* decode phase for correlator */
	IF_OK {
	  param.corr_phase[ipair][i] = decode_phase(phase_lab);
	  if(param.corr_phase[ipair][i] == -1){
	    printf("\n%s is not a valid phase label\n",phase_lab);
	    status ++;
	  }
	}
	
	/* decode real factor for correlator */
	/* Permitted syntax is /ddd or *ddd where ddd is a real number */
	IF_OK {
	  param.corr_factor[ipair][i] = decode_factor(factor_op,factor);
	  if(param.corr_factor[ipair][i] == 0 ){
	    printf("\n: Decoding factor %s %g results in zero.\n",
		   factor_op,factor);
	    status++;
	  }
	}
	
	/* decode gamma matrix labels for source and sink */
	IF_OK {
	  param.gam_src[ipair][i] = gamma_index(gam_src_lab);
	  if(param.gam_src[ipair][i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",gam_src_lab);
	    status += 1;
	  }
	  param.gam_snk[ipair][i] = gamma_index(gam_snk_lab);
	  if(param.gam_snk[ipair][i] < 0){
	    printf("\n%s is not a valid gamma matrix label\n",gam_snk_lab);
	    status ++;
	  }
	}
	
	/* momentum indices */
	IF_OK {
	  if(scanf("%d %d %d",&param.corr_mom[ipair][i][0],
		   &param.corr_mom[ipair][i][1],
		   &param.corr_mom[ipair][i][2]) != 3){
	    printf("\nFormat error on momentum line\n");
	    status++;
	  }
	  
	  IF_OK printf(" %2d %2d %2d",param.corr_mom[ipair][i][0],
		       param.corr_mom[ipair][i][1],
		       param.corr_mom[ipair][i][2]);
	}
	
	/* parity for each momentum component */
	
	IF_OK {
	  ok = scanf("%s %s %s", parity_x_in, parity_y_in, parity_z_in);
	  if(ok != 3){
	    printf("\nError reading parity label\n");
	    status++;
	  }
	  IF_OK printf(" %2s %2s %2s\n", parity_x_in, parity_y_in, parity_z_in);
	}
	
	/* decode parity labels */
	IF_OK {
	  param.corr_parity[ipair][i][0] = decode_parity(parity_x_in);
	  param.corr_parity[ipair][i][1] = decode_parity(parity_y_in);
	  param.corr_parity[ipair][i][2] = decode_parity(parity_z_in);
	  if(param.corr_parity[ipair][i][0] == 0xf)status++;
	  if(param.corr_parity[ipair][i][1] == 0xf)status++;
	  if(param.corr_parity[ipair][i][2] == 0xf)status++;
	}
      } /* correlators for this pair */
    } /* pairs */
    
    /* End of input fields */
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */
    

  broadcast_bytes((char *)&param,sizeof(param));
  u0 = param.u0;
  if( param.stopflag != 0 )
    normal_exit(0);

  /* Do whatever is needed to get lattice */
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( param.startflag, param.startfile );
  /* Construct APE smeared links */
  ape_smear_3D( param.staple_weight, param.ape_iter );

  /* Initialization for KS operations if needed */
#ifdef FN
  invalidate_all_ferm_links(&fn_links);
#endif
  phases_in = OFF;

  /* make table of coefficients and permutations of paths in quark action */
  init_path_table(&ks_act_paths);
  make_path_table(&ks_act_paths, NULL);

  return(0);
}

/* decode parity label */
char decode_parity(char *parity_label_in){
  if(strcmp("E",parity_label_in) == 0)return EVEN;
  else if(strcmp("O",parity_label_in) == 0)return ODD;
  else if(strcmp("EO",parity_label_in) == 0)return EVENANDODD;

  printf("\nUnrecognized parity label %s\n",parity_label_in);
  return 0xf;
}
    
/* decode real factor label */
/* Permitted syntax is *ddd or /ddd where ddd is a real number */
double decode_factor(char *factor_op, double factor){

  if(factor == 0)
    return factor;
  else if(factor_op[0] == '*')
    return factor;
  else if(factor_op[0] == '/')
    return 1/factor;
  else
    printf("Incorrect operator for normalization factor\n");
  return 0.;
}
    

/* Manage hash table for correlator labels */

/* Look the labels up in the table.  If found, return its
   index.  Otherwise, add it to the table and assign it the new index */
static int 
hash_corr_label(char meson_label[MAX_CORR][MAX_MESON_LABEL], 
		char mom_label[MAX_CORR][MAX_MOM_LABEL],
		char *meson_label_in, char *mom_label_in, int *n){
  int i;

  for(i = 0; i < *n; i++){
    if(strcmp(meson_label_in, meson_label[i]) == 0 &&
       strcmp(mom_label_in,   mom_label[i]  ) == 0)
      return i;
  }
  if(*n >= MAX_CORR){
    printf("Too many correlators labels %d\n",*n);
    return -1;
  }
  strcpy(meson_label[*n],meson_label_in);
  strcpy(mom_label[*n],mom_label_in);
  (*n)++;
  return *n-1;
}


/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
void 
make_3n_gathers()
{
  int i;
#ifdef HAVE_QDP
  int disp[4]={0,0,0,0};
#endif
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X3UP);

#ifdef HAVE_QDP
  for(i=0; i<4; i++) {
    disp[i] = 3;
    neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
  }
#endif
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

void 
third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
	"forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+3)%nx; break;
  case XDOWN: *xp = (x+4*nx-3)%nx; break;
  case YUP: *yp = (y+3)%ny; break;
  case YDOWN: *yp = (y+4*ny-3)%ny; break;
  case ZUP: *zp = (z+3)%nz; break;
  case ZDOWN: *zp = (z+4*nz-3)%nz; break;
  case TUP: *tp = (t+3)%nt; break;
  case TDOWN: *tp = (t+4*nt-3)%nt; break;
  default: printf("third_neighb: bad direction\n"); exit(1);
  }
}
#if 0
    /*------------------------------------------------------------*/
    /* Gamma matrix pairings                                      */
    /*------------------------------------------------------------*/
    
    /* Number of gamma matrix pairings */
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
#endif
