/***************** control_cl.c *****************************************/

/* Main procedure for quenched SU3 clover fermions 			*/
/* MIMD version 7 */

/* This version computes propagators for clover fermions on a
   supplied background field config */

/* Modifications ...
   
   8/8/98  Rearranged loop so hadron propagators are written as soon
           as they are calculated. C.D.
   8/10/96 Installed new propagator IO and added timing C.D.
 */

#define CONTROL
#include "cl_inv_includes.h"
#include <string.h>
#ifdef HAVE_QDP
#include <qdp.h>
#endif

/* Comment these out if you want to suppress detailed timing */
/*#define IOTIME*/

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
  
  int MaxCG;
  Real RsdCG, RRsdCG;
  
  int spin,color,k;
  int flag;

  int status;

  int cl_cg = CL_CG;

  w_prop_file *fp_in_w[MAX_KAP];        /* For propagator files */
  w_prop_file *fp_out_w[MAX_KAP];       /* For propagator files */

  wilson_vector *psi = NULL;
  wilson_prop_field quark_propagator = NULL;
  
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
  
  while( readin(prompt) == 0)
    {
      
      starttime=dclock();
      MaxCG=niter;
      wqstmp = wqs;  /* For clover_info.c */
      
      avm_iters=0.0;
      meascount=0;
      
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
      /* Loop over kappas */
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
	fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
#ifdef IOTIME
	dtime += dclock();
	if(startflag_w[k] != FRESH)
	node0_printf("Time to open prop = %e\n",dtime);
#endif
	fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				   wqs.type);
	
	
	/* Loop over source spins */
	for(spin=0;spin<4;spin++){
	  /* Loop over source colors */
	  for(color=0;color<3;color++){
	  
	    meascount ++;
	    /*if(this_node==0)printf("color=%d spin=%d\n",color,spin); */
	    if(startflag_w[k] == CONTINUE)
	      {
		node0_printf("Can not continue propagator here! Zeroing it instead\n");
		startflag_w[k] = FRESH;
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
		node0_printf("control_cl: Recovering from error by resetting initial guess to zero\n");
		reload_wprop_sc_to_field( FRESH, fp_in_w[k], &wqs,
			       spin, color, psi, 0);
		flag = 0;
	      }

	    /* Complete the source structure */
	    wqs.color = color;
	    wqs.spin = spin;

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


	    /* compute the propagator.  Result in psi. */
	    
	    switch (cl_cg) {
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
			     this_node,cl_cg);
	      }

	    avm_iters += avs_iters;

	    copy_wp_from_wv(quark_propagator, psi, color, spin);
	    
	    /* save psi if requested */
	    save_wprop_sc_from_field( saveflag_w[k],fp_out_w[k], &wqs,
			     spin,color,psi,"Fill in record info here",iotime);
	  } /* source spins */
	} /* source colors */
	
	/* close files for wilson propagators */
	r_close_wprop(startflag_w[k],fp_in_w[k]);
#ifdef IOTIME
	dtime = -dclock();
#endif
	w_close_wprop(saveflag_w[k],fp_out_w[k]);
#ifdef IOTIME
      dtime += dclock();
      if(saveflag_w[k] != FORGET)
	node0_printf("Time to close prop = %e\n",dtime);
#endif
	
	/* spectrum, if requested */

	if(strstr(spectrum_request,",spectrum,") != NULL)
	  spectrum_cl(quark_propagator, wqs.t0, k);

      } /* kappas */


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
