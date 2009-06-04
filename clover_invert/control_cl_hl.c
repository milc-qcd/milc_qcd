/***************** control_cl_hl.c *****************************************/

/* Main procedure for quenched SU3 clover fermions 			*/
/* MIMD version 7 */

/* This version computes propagators for clover fermions on a
   supplied background field config */

/* Modifications ...
   
   8/15/96 Made scratch file name a variable C.D. 
   8/10/96 Installed new propagator IO and added timing C.D. */

#define CONTROL
#include "cl_inv_includes.h"
#include <string.h>
#ifdef HAVE_QDP
#include <qdp.h>
#endif

#ifndef HAVE_QIO
MUST COMPILE WITH QIO FOR THE SCRATCH FILE
#endif

/* Comment these out if you want to suppress detailed timing */
/*#define IOTIME*/
/*#define PRTIME*/

int main(int argc, char *argv[])
{
  int meascount;
  int prompt;
  Real avm_iters,avs_iters;
  
  double starttime,endtime;
#ifdef IOTIME
  double dtime;
  int iotime = 1;
#else
  int iotime = 0;
#endif
  
  int MinCG,MaxCG;
  Real RsdCG, RRsdCG;
  
  int spin,color,j,k;
  int flag;
  int status;
  
  w_prop_file *fp_in_w[MAX_KAP];        /* For reading binary propagator files */
  w_prop_file *fp_out_w[MAX_KAP];       /* For writing binary propagator files */
  w_prop_file *fp_scr[MAX_KAP];
  wilson_quark_source wqs_scr;  /* scratch file */
  char scratch_file[MAX_KAP][MAXFILENAME];
  
  wilson_vector *psi = NULL;
  wilson_prop_field quark_propagator = NULL;
  wilson_prop_field quark_prop2 = NULL;
  int cg_cl = CL_CG;
  int source_type;
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup_cl();
  /* loop over input sets */
  
  psi = create_wv_field();
  quark_propagator = create_wp_field();
  quark_prop2 = create_wp_field();

  while( readin(prompt) == 0)
    {
      
      starttime=dclock();
      MaxCG=niter;
      
      avm_iters=0.0;
      meascount=0;
      
      spectrum_cl_hl_init();
  
      if( fixflag == COULOMB_GAUGE_FIX)
	{
	  if(this_node == 0) 
	    printf("Fixing to Coulomb gauge\n");
#ifdef IOTIME
	  dtime = -dclock();
#endif
	  gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
#ifdef IOTIME
	  dtime += dclock();
	  if(this_node==0)printf("Time to gauge fix = %e\n",dtime);
#endif
	  invalidate_this_clov(gen_clov);
	}
      else
	if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
      
      /* save lattice if requested */
      if( saveflag != FORGET ){
	savelat_p = save_lattice( saveflag, savefile, stringLFN );
      }
      
      if(this_node==0)printf("END OF HEADER\n");
      
      /*	if(this_node==0) printf("num_kap = %d\n", num_kap); */
      /* Loop over kappas to compute and store quark propagator */
      for(k=0;k<num_kap;k++){
	
	kappa=kap[k];
	RsdCG=resid[k];
	RRsdCG=relresid[k];
	if(this_node==0)printf("Kappa= %g r0= %g residue= %g rel= %g\n",
			       (double)kappa,(double)wqs.r0,(double)RsdCG,
			       (double)RRsdCG);
	
	/* open files for wilson propagators */
	
#ifdef IOTIME
	dtime = -dclock();
#endif
	wqstmp = wqs;  /* For clover_info.c */
	fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
	fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				   wqs.type);
#ifdef IOTIME
	dtime += dclock();
	if(startflag_w[k] != FRESH)
	  node0_printf("Time to open prop = %e\n",dtime);
#endif
	
	/* Open scratch file and write header */
	sprintf(scratch_file[k],"%s_%02d",scratchstem_w,k);
	source_type = UNKNOWN;
	fp_scr[k] = w_open_wprop(scratchflag, scratch_file[k], source_type);
	init_wqs(&wqs_scr);

	  
	/* Loop over source spins */
	for(spin=0;spin<4;spin++){
	  /* Loop over source colors */
	  for(color=0;color<3;color++){
	    meascount ++;
	    /*if(this_node==0)printf("color=%d spin=%d\n",color,spin);*/
	    
	    if(startflag_w[k] == CONTINUE)
	      {
		if(k == 0)
		  {
		    node0_printf("Can not continue propagator here! Zeroing it instead\n");
		    startflag_w[k] = FRESH;
		  }
		else
		  {
		    copy_wv_from_wp(psi, quark_propagator, color, spin);
		  }
	      }
	    
	    /* Saves one multiplication by zero in cgilu */
	    if(startflag_w[k] == FRESH)flag = 0;
	    else 
	      flag = 1;      
	    
	    /* load psi if requested */
	    status = reload_wprop_sc_to_field( startflag_w[k], fp_in_w[k], 
					       &wqs, spin, color, psi, iotime);
	    if(status != 0)
	      {
		node0_printf("control_cl_hl: Recovering from error by resetting initial guess to zero\n");
		reload_wprop_sc_to_field( FRESH, fp_in_w[k], &wqs,
					  spin, color, psi,0);
		flag = 0;
	      }
	    
	    /* Complete the source structure */
	    wqs.color = color;
	    wqs.spin = spin;
	    
	    /* If we are starting afresh, we set a minimum number
	       of iterations */
	    if(startflag_w[k] == FRESH || status != 0)MinCG = nt/2; 
	    else MinCG = 0;
	    
	    /* Load inversion control structure */
	    qic.prec = PRECISION;
	    qic.max = MaxCG;
	    qic.nrestart = nrestart;
	    qic.resid = RsdCG;
	    qic.relresid = RRsdCG;
	    qic.start_flag = flag;
	    
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa;
	    dcp.Clov_c = clov_c;
	    dcp.U0 = u0;
	    
	    switch (cg_cl) {
	    case BICG:
	      avs_iters =
		(Real)wilson_invert_field_wqs(&wqs, w_source_field, psi,
					      bicgilu_cl_field,
					      &qic,(void *)&dcp);
	      break;
	    case HOP:
	      avs_iters = 
		(Real)wilson_invert_field_wqs(&wqs, w_source_field, psi,
					      hopilu_cl_field,
					      &qic,(void *)&dcp);
	      break;
	    case MR:
	      avs_iters = 
		(Real)wilson_invert_field_wqs(&wqs, w_source_field, psi,
					      mrilu_cl_field,
					      &qic,(void *)&dcp);
		break;
	    case CG:
	      avs_iters = 
		(Real)wilson_invert_field_wqs(&wqs, w_source_field, psi,
					      cgilu_cl_field,
					      &qic,(void *)&dcp);
	      break;
	    default:
	      node0_printf("main(%d): Inverter choice %d not supported\n",
			   cg_cl, this_node);
	    }
	    
	    avm_iters += avs_iters;
	    
	    copy_wp_from_wv(quark_propagator, psi, color, spin);
	    
	    /* Write psi to scratch disk */
#ifdef IOTIME
	    dtime = -dclock();
#endif
	    save_wprop_sc_from_field(scratchflag, fp_scr[k], &wqs_scr,
				     spin, color, psi, "Scratch record", iotime);
#ifdef IOTIME
	    dtime += dclock();
	    if(this_node==0) 
	      printf("Time to dump prop spin %d color %d %e\n",
		     spin,color,dtime);
#endif
	    
	    /* save psi if requested */
	    save_wprop_sc_from_field( saveflag_w[k],fp_out_w[k], &wqs,
				      spin,color,psi,"", iotime);
	    
	  } /* source colors */
	} /* source spins */
	
	/* Close and release scratch file */
	w_close_wprop(scratchflag, fp_scr[k]);
	
	/*if(this_node==0)printf("Dumped prop to file  %s\n",
	  scratch_file[k]); */
	
	/* close files for wilson propagators */
#ifdef IOTIME
	dtime = -dclock();
#endif
	r_close_wprop(startflag_w[k],fp_in_w[k]);
	w_close_wprop(saveflag_w[k],fp_out_w[k]);
#ifdef IOTIME
	dtime += dclock();
	if(saveflag_w[k] != FORGET)
	  node0_printf("Time to close prop = %e\n",dtime);
#endif
	
      } /* kappas */
      
      
      /* Loop over heavy kappas for the point sink spectrum */
      for(k=0;k<num_kap;k++){
	
	/* Read the propagator from the scratch file */

#ifdef IOTIME
	dtime = -dclock();
#endif
	kappa=kap[k];
	init_wqs(&wqs_scr);
	reload_wprop_to_wp_field(scratchflag, scratch_file[k], &wqs_scr,
				 quark_propagator, iotime);
#ifdef IOTIME
	dtime += dclock();
	if(this_node==0) 
	  {
	    printf("Time to read 12 spin,color combinations %e\n",dtime);
	    fflush(stdout);
	  }
#endif
	
	/*if(this_node==0)
	  printf("Closed scratch file %s\n",scratch_file[k]);
	  fflush(stdout); */
	
	/* Diagonal spectroscopy */
	spectrum_cl_hl_diag_baryon(quark_propagator, k);
	spectrum_cl_hl_diag_meson(quark_propagator, k);
	spectrum_cl_hl_diag_rot_meson(quark_propagator, k);
	if(strstr(spectrum_request,",sink_smear,") != NULL){
	  spectrum_cl_hl_diag_smeared_meson(quark_propagator, k);
	}
	
	/* Heavy-light spectroscopy */
	/* Loop over light kappas for the point sink spectrum */
	for(j=k+1;j<num_kap;j++){

#ifdef IOTIME
	  dtime = -dclock();
#endif
	  /* Read the propagator from the scratch file */
	  kappa=kap[j];
	  init_wqs(&wqs_scr);
	  reload_wprop_to_wp_field(scratchflag,  scratch_file[j], &wqs_scr,
				   quark_prop2, iotime);
#ifdef IOTIME
	  dtime += dclock();
	  if(this_node==0) 
	    {
	      printf("Time to read 12 spin,color combinations %e\n",dtime);
	      fflush(stdout);
	    }
#endif
#ifdef PRTIME
	  dtime = -dclock();
#endif
	  spectrum_cl_hl_offdiag_baryon( quark_propagator, quark_prop2, 
					 j, k);
	  spectrum_cl_hl_offdiag_meson( quark_propagator, quark_prop2, 
					j, k);
	  spectrum_cl_hl_offdiag_rot_meson( quark_propagator, quark_prop2, 
					    j, k);
	  
#ifdef PRTIME
	  dtime = -dclock();
#endif
	} /* light kappas */
	
	/* Smear the heavy propagator in place */
	sink_smear_prop( quark_propagator );
	
	/* Write the smeared propagator to the scratch file (overwriting)*/
	
	kappa=kap[k];

#ifdef IOTIME
	dtime = -dclock();
#endif
	save_wprop_from_wp_field(scratchflag, scratch_file[k], &wqs_scr,
				 quark_propagator, "Scratch propagator", 
				 iotime);
	
#ifdef IOTIME
	dtime += dclock();
	if(this_node==0) 
	  {
	    printf("Time to dump convolution %d %e\n",k,dtime);
	    fflush(stdout);
	  }
#endif
      } /* heavy kappas */
      
      /* Loop over heavy kappas for the shell sink spectrum */
      if(strstr(spectrum_request,",sink_smear,") != NULL)
	for(k=0;k<num_kap;k++){
	  
#ifdef IOTIME
	  dtime = -dclock();
#endif
	  /* Read the propagator from the scratch file */
	  kappa=kap[k];
	  init_wqs(&wqs_scr);
	  reload_wprop_to_wp_field(scratchflag,  scratch_file[k], &wqs_scr,
				   quark_propagator, iotime);
#ifdef IOTIME
	  dtime += dclock();
	  if(this_node==0) 
	    {
	      printf("Time to read convolution %d %e\n",k,dtime);
	      fflush(stdout);
	    }
#endif
	  
	  /* Diagonal spectroscopy */
	  spectrum_cl_hl_diag_smeared_meson(quark_propagator, k);
	  
	  /* Heavy-light spectroscopy */
	  /* Loop over light kappas for the shell sink spectrum */
	  for(j=k+1;j<num_kap;j++){
#ifdef PRTIME
	    dtime = -dclock();
#endif
	    /* Read the propagator from the scratch file */
	    kappa=kap[j];
	    init_wqs(&wqs_scr);
	    reload_wprop_to_wp_field(scratchflag,  scratch_file[j], &wqs_scr,
				     quark_prop2, iotime);
	      
	    /* Compute the spectrum */
	    spectrum_cl_hl_offdiag_smeared_meson( quark_propagator,
						  quark_prop2, j, k);
	    
#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0) 
	      {
		printf("Time to read and do off diagonal mesons %d %d %e\n",
		       j,k,dtime);
		fflush(stdout);
	      }
#endif
	} /* light kappas */
	
      } /* heavy kappas */
      
      spectrum_cl_hl_print(wqs.t0);
      spectrum_cl_hl_cleanup();

      if(this_node==0)printf("RUNNING COMPLETED\n");
      if(meascount>0){
	if(this_node==0)printf("total cg iters for measurement= %e\n",
			       (double)avm_iters);
	if(this_node==0)printf("cg iters for measurement= %e\n",
			       (double)avm_iters/(double)meascount);
      }
      
      endtime=dclock();
      if(this_node==0){
	printf("Time = %e seconds\n",(double)(endtime-starttime));
	printf("total_iters = %d\n",total_iters);
      }
      fflush(stdout);
    }
      
  destroy_wv_field(psi);
  destroy_wp_field(quark_propagator);
  
  return 0;
}

