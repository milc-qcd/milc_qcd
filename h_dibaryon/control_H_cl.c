/***************** control_H_cl.c ***************************************/
/* MIMD version 6 */

/* Calculate H-dibaryon multichannel propagators and
   associated norelativistic baryon propagators

   C. DeTar 7/6/98 Ported to version 5
   
   This version for clover propagators 
   with PERIODIC boundary conditions */

/* Modifications ...
   4/6/00 Generalized to any source time CBM CD
 */

#define CONTROL
#include "h_dibaryon_includes.h"
#include <string.h>

/**************************/
/* Timing macros */

/* Comment out either one or both of the first two defines if you want
   to suppress detailed timing */
#define IOTIME   /* For I/O related processing */
#define PRTIME   /* For propagator and correlator processing */

#ifdef IOTIME
#define STARTIOTIME dtime = -dclock();
#define STOPIOTIME(STRING) 	  dtime += dclock(); \
	  if(this_node==0){printf("Time to %s = %e\n",STRING,dtime); \
          fflush(stdout);}
#else
#define STARTIOTIME
#define STOPIOTIME(STRING)
#endif

#ifdef PRTIME
#define STARTPRTIME dtime = -dclock();
#define STOPPRTIME(STRING) 	  dtime += dclock(); \
	  if(this_node==0)printf("Time to %s = %e\n",STRING,dtime);
#else
#define STARTPRTIME
#define STOPPRTIME(STRING)
#endif
/**************************/

#define NRPROPS 2
#define HDIPROPS 4

int main(int argc,char *argv[])
{
  int meascount;
  int prompt;
  Real avm_iters,avs_iters;
  
  double starttime,endtime,dclock();
  double dtime;
  
  int MinCG,MaxCG;
  Real RsdCG;
  
  register int i;
  register site *s;
  
  int spinindex,spin,color,j,k,t,t_off;
  int kh,kl;
  int nr_fb;
  char nr_fb_label[3][2] = { "0", "F", "B" };
  int flag;
  int kprop;
  int num_prop;
  Real space_vol;

  int status;

  propagator hdibar_prop[MAX_KAP][MAX_KAP][HDIPROPS];
  propagator nrbar_prop[MAX_KAP][MAX_KAP][NRPROPS];
  
  char scratch_file[MAX_KAP][MAXFILENAME];
  
  Real norm_fac[10];
  static char *mes_kind[10] = {"PION","PS505","PS055","PS0505",
			       "RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};

  complex *pmes_prop[MAX_KAP][MAX_KAP][10];
  int pmes_prop_done[MAX_KAP][MAX_KAP];

  w_prop_file *fp_in_w[MAX_KAP];  /* For reading binary propagator files */
  w_prop_file *fp_out_w[MAX_KAP]; /* For writing binary propagator files */
  w_prop_file *fp_scr[MAX_KAP];
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup_H_cl();
  

  /* loop over input sets */
  
  while( readin(prompt) == 0)
    {

      MaxCG = niter;
      starttime=dclock();

      avm_iters=0.0;
      meascount=0;
      
      /* Allocate space for relativistic meson propagator */
      for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++)for(j=0;j<=i;j++){
	  pmes_prop[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++){
	    pmes_prop[i][j][num_prop][t] = cmplx(0.0,0.0); 
	  }
	  pmes_prop_done[i][j] = 0;
	}

      /* Allocate space for non relativistic baryon propagators */
      for(kprop=0;kprop<NRPROPS;kprop++)
	for(i=0;i<num_kap;i++)for(j=0;j<num_kap;j++){
	  nrbar_prop[i][j][kprop].c
	    = (complex *)malloc(nt*sizeof(complex));
	  if(nrbar_prop[i][j][kprop].c == NULL)
	    {
	      printf("control_H_cl: Can't malloc nrbar prop %d %d %d\n",
		     i,j,kprop);
	      terminate(1);
	    }
	  for(t=0;t<nt;t++)nrbar_prop[i][j][kprop].c[t] 
	    = cmplx(0.0,0.0); 
	  nrbar_prop[i][j][kprop].label
	    = (char *)malloc(10*sizeof(char));
	  if(nrbar_prop[i][j][kprop].c == NULL)
	    {
	      printf("control_H_cl: Can't malloc nrbar prop label %d %d %d\n",
		     i,j,kprop);
	      terminate(1);
	    }
	}
      
      /* Allocate space for H-dibaryon channel propagators */
      for(kprop=0;kprop<HDIPROPS;kprop++)
	for(kh=0;kh<num_kap_heavy;kh++)for(kl=0;kl<num_kap_light;kl++){
	  /* kappa indexing scheme is consistent with baryon propagator
	     even though we compute only the propagators with
	     one heavy (s) quark and two light (u,d) quarks */
	  i = kh; j = kl + num_kap_heavy;
	  hdibar_prop[i][j][kprop].c
	    = (complex *)malloc(nt*sizeof(complex));
	  if(hdibar_prop[i][j][kprop].c == NULL)
	    {
	      printf("control_H_cl: Can't malloc baryon prop %d %d %d\n",
		     i,j,kprop);
	      terminate(1);
	    }
	  for(t=0;t<nt;t++)hdibar_prop[i][j][kprop].c[t] 
	    = cmplx(0.0,0.0); 
	  hdibar_prop[i][j][kprop].label
	    = (char *)malloc(10*sizeof(char));
	  if(hdibar_prop[i][j][kprop].label == NULL)
	    {
	      printf("control_H_cl: Can't malloc baryon prop label %d %d %d\n",
		     i,j,kprop);
	      terminate(1);
	    }
	}
      
      if( fixflag == COULOMB_GAUGE_FIX)
	{
	  if(this_node == 0) 
	    printf("Fixing to Coulomb gauge\n");
	  STARTIOTIME;
	  gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
	  STOPIOTIME("gauge fix");
	  invalidate_this_clov(gen_clov);
	}
      else
	if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
      
      /* save lattice if requested */
      if( saveflag != FORGET ){
	/* Note: beta, kappa are kept only for save_old_binary */
	STARTIOTIME;
	savelat_p = save_lattice( saveflag, savefile, stringLFN );
	STOPIOTIME("save lattice");
      }

      if(this_node==0)printf("END OF HEADER\n");
      
      /* Loop over all kappas to compute and store quark propagator */
      for(k=0;k<num_kap;k++){
	
	kappa = kap[k];
	source_r0=wqs[k].r0;
	RsdCG=resid[k];
	if(this_node==0)printf("Kappa=%e r0=%e residue=%e\n",
			       (double)kappa,(double)source_r0,(double)RsdCG);
	
	/* open file for kth wilson propagator */
	
	fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
	fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				   wqs[k].type);
	
	/* Open scratch file and write header */
	sprintf(scratch_file[k],"%s_%02d",scratchstem_w,k);
	if(scratchflag == SAVE_CHECKPOINT)
	  {
	    fp_scr[k] = w_checkpoint_w_i(scratch_file[k]);
	    /* Close, temporarily */
	    w_checkpoint_w_c(fp_scr[k]);
	  }
	else
	  /* If serial, write header and leave it open */
	  fp_scr[k] = w_serial_w_i(scratch_file[k]);
	
	/* Loop over source colors */
	for(color=0;color<3;color++){
	  
	  for(spinindex=0;spinindex<n_spins;spinindex++){
	    spin = spins[spinindex];
	    
	    meascount ++;
	    if(this_node==0)printf("color=%d spin=%d\n",color,spin);

	    if(startflag_w[k] == CONTINUE)
	      {
		if(k == 0)
		  {
		    node0_printf("Can not continue propagator here! Zeroing it instead\n");
		    startflag_w[k] = FRESH;
		  }
		else
		  {
		    FORALLSITES(i,s)
		      copy_wvec(&(s->quark_propagator.c[color].d[spin]),
				&(s->psi));
		  }
	      }

	    /* Saves one multiplication by zero in cgilu */
	    if(startflag_w[k] == FRESH)flag = 0;
	    else 
	      flag = 1;      
	    
	    /* load psi if requested */
#ifdef IOTIME
	    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
			       spin, color, F_OFFSET(psi),1);
#else
	    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
			       spin, color, F_OFFSET(psi),0);
#endif	    
	    if(status != 0)
	      {
		node0_printf("control_H_cl: Recovering from error by resetting initial guess to zero\n");
		reload_wprop_sc_to_site( FRESH, fp_in_w[k], 
			       spin, color, F_OFFSET(psi),0);
		flag = 0;
	      }

	    
	    /* Invert to find propagator */

	    /* Complete the source structure */
	    wqs[k].color = color;
	    wqs[k].spin = spin;

	    /* For clover_info */
	    wqstmp = wqs[k];

	   /* If we are starting afresh, we set a minimum number
	      of iterations */
	   if(startflag_w[k] == FRESH || status != 0)MinCG = nt; 
	   else MinCG = 0;

	    /* Load inversion control structure */
	    qic.prec = PRECISION;
	    qic.min = MinCG;
	    qic.max = MaxCG;
	    qic.nrestart = nrestart;
	    qic.start_flag = flag;
	    qic.nsrc = 1;
	    qic.resid = RsdCG;
	    qic.relresid = 0;
	    
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa;
	    dcp.Clov_c = clov_c;
	    dcp.U0 = u0;

#ifdef BI
	    /* compute the propagator.  Result in psi. */
	    avs_iters 
	      = (Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
					  w_source,&wqs[k],
					  bicgilu_cl_site,&qic,(void *)&dcp);
#else
	    /* compute the propagator.  Result in psi. */
	    avs_iters = 
	      (Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
					w_source,&wqs[k],
					cgilu_cl_site,&qic,(void *)&dcp);
#endif
	    avm_iters += avs_iters;
	    
	    FORALLSITES(i,s)
	      copy_wvec(&(s->psi),
			&(s->quark_propagator.c[color].d[spin]));
	    
	    STARTIOTIME;
	    /* Write psi to scratch disk */
	    if(scratchflag == SAVE_CHECKPOINT)
	      {
		w_checkpoint_w_o(fp_scr[k]);
		w_checkpoint_w(fp_scr[k],spin,color,F_OFFSET(psi));
		w_checkpoint_w_c(fp_scr[k]);
	      }
	    else
	      w_serial_w(fp_scr[k],spin,color,F_OFFSET(psi));
	    STOPIOTIME("do fast quark dump");
	    /* save psi if requested */
#ifdef IOTIME
	    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
			     spin,color,F_OFFSET(psi),1);
#else
	    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
			     spin,color,F_OFFSET(psi),0);
#endif
	  } /* source spins */
	} /* source colors */
	
	/* Close and release scratch file */
	if(scratchflag == SAVE_CHECKPOINT)
	  w_checkpoint_w_f(fp_scr[k]);
	else
	  w_serial_w_f(fp_scr[k]);

	if(this_node==0)printf("Saved binary wilson_vector in file  %s\n",
			       scratch_file[k]);
	
	/* close files for wilson propagators */
	r_close_wprop(startflag_w[k],fp_in_w[k]);
	w_close_wprop(saveflag_w[k],fp_out_w[k]);
	
      } /* kappas */
      
      
      /* Loop over choice forward - backward for NR source and sink */

      for(nr_fb = 1; nr_fb <= 2; nr_fb++)if(nr_fb & nr_forw_back)
	{
	  
	  /* Reset completion flags */
	    for(i=0;i<num_kap;i++)for(j=0;j<num_kap;j++){
	      for(kprop=0;kprop<NRPROPS;kprop++)
		nrbar_prop[i][j][kprop].done = 0;
	      for(kprop=0;kprop<HDIPROPS;kprop++)
		hdibar_prop[i][j][kprop].done = 0;
	    }

	  /* Loop over heavy kappas for the point sink spectrum */
	  for(k=0;k<num_kap_heavy;k++){
	    
	    /* Read the kth heavy kappa propagator from the scratch file */
	    kappa = kappa_heavy = kap[k];
	    if(scratchflag == SAVE_CHECKPOINT)
	      fp_scr[k] = r_parallel_w_i(scratch_file[k]);
	    else
	      fp_scr[k] = r_serial_w_i(scratch_file[k]);
	    
	    STARTIOTIME;
	    for(color=0;color<3;color++) for(spin=0;spin<4;spin++){
	      if(scratchflag == SAVE_CHECKPOINT)
		r_parallel_w(fp_scr[k], spin, color,
			     F_OFFSET(quark_propagator.c[color].d[spin])); 
	      else
		r_serial_w(fp_scr[k], spin, color,
			   F_OFFSET(quark_propagator.c[color].d[spin])); 
	    }
	    STOPIOTIME("to read 12 spin-color combinations");

	    if(scratchflag == SAVE_CHECKPOINT)
	      r_parallel_w_f(fp_scr[k]); 
	    else
	      r_serial_w_f(fp_scr[k]); 
	    
	    /* Convert to NR propagator */
	    
	    STARTPRTIME;
	    nr_propagator(F_OFFSET(quark_propagator),
			  F_OFFSET(nr_prop1), nr_fb);
	    diquarkprop(F_OFFSET(nr_prop1),
			F_OFFSET(diquark_prop1));
	    STOPPRTIME("make nr and diquark");
	    
	    /* Diagonal spectroscopy - not needed */
	    
/**	    w_nrbaryon(F_OFFSET(nr_prop1), F_OFFSET(nr_prop1),
		       F_OFFSET(diquark_prop1), nrbar_prop[k][k]); **/
	    
/**	    w_hdibaryon(F_OFFSET(diquark_prop1),
			F_OFFSET(diquark_prop1), hdibar_prop[k][k]); **/
	    
	    /* Heavy-light spectroscopy */
	    /* Loop over light kappas for the point sink spectrum */
	    for(j=num_kap_heavy;j<num_kap;j++){

	      /* Read the propagator from the scratch file */
	      kappa = kappa_light = kap[j];
	      if(scratchflag == SAVE_CHECKPOINT)
		fp_scr[j] = r_parallel_w_i(scratch_file[j]);
	      else
		fp_scr[j] = r_serial_w_i(scratch_file[j]);
	      
	      STARTIOTIME;
	      for(color=0;color<3;color++) for(spin=0;spin<4;spin++){
		if(scratchflag == SAVE_CHECKPOINT)
		  r_parallel_w(fp_scr[j], spin, color,
			       F_OFFSET(quark_prop2.c[color].d[spin])); 
		else
		  r_serial_w(fp_scr[j], spin, color,
			     F_OFFSET(quark_prop2.c[color].d[spin])); 
	      }
	      STOPIOTIME("do fast quark read");
	      if(scratchflag == SAVE_CHECKPOINT)
		r_parallel_w_f(fp_scr[j]);
	      else
		r_serial_w_f(fp_scr[j]);
	      
	      /* Convert to NR propagator */
	      
	      STARTPRTIME;
	      nr_propagator(F_OFFSET(quark_prop2),
			    F_OFFSET(nr_prop2),nr_fb);
	      diquarkprop(F_OFFSET(nr_prop2),
			  F_OFFSET(diquark_prop2));
	      STOPPRTIME("make nr and diquark propagators");
	      
	      /* Diagonal spectroscopy - baryons only - done if
		 any of them was not previously done */
	    
	      for(kprop=0;kprop<NRPROPS;kprop++)
		{
		  if(nrbar_prop[j][j][kprop].done == 0)
		    {
		      STARTPRTIME;
		      w_nrbaryon(F_OFFSET(nr_prop2), F_OFFSET(nr_prop2),
				 F_OFFSET(diquark_prop2), nrbar_prop[j][j]);
		      STOPPRTIME("do diagonal baryons");
		      break;
		    }
		}
	    
	      /* Heavy-light spectroscopy - baryons and H */

	      /* We don't do baryon heavy-light if the kappa values
		 are the same, since the result is the same as the
		 diagonal light propagator */
	      
	      if(kappa_heavy != kappa_light)
		{
		  /* Relativistic meson propagator: Do only once */		  
		  if(pmes_prop_done[j][k] == 0) {
		    STARTPRTIME;
		    for(color=0;color<3;color++){
		      w_meson_site(F_OFFSET(quark_propagator.c[color]),
			      F_OFFSET(quark_prop2.c[color]), pmes_prop[j][k]);
		    }
		    pmes_prop_done[j][k] = 1;
		    STOPPRTIME("do off-diagonal relativistic meson");
		  }

		  STARTPRTIME;
		  w_nrbaryon(F_OFFSET(nr_prop2),
			     F_OFFSET(nr_prop1),F_OFFSET(diquark_prop1),  
			     nrbar_prop[j][k]);
		  
		  w_nrbaryon(F_OFFSET(nr_prop1),
			     F_OFFSET(nr_prop2),F_OFFSET(diquark_prop2),  
			     nrbar_prop[k][j]);
		  STOPPRTIME("do two sets of hl baryons");
		}
	      
	      /* For H we do only the case prop2 = u (light) index j
		 and prop1 = s (heavy) index k */

	      STARTPRTIME;
	      w_hdibaryon(F_OFFSET(diquark_prop2),
			  F_OFFSET(diquark_prop1), hdibar_prop[k][j]);
	      STOPPRTIME("do one set of hl H dibaryons");
	      
	    } /* light kappas */
	  } /* heavy kappas */

	  /* Stick with same convention as clover_invert/control_cl_hl.c */
	  space_vol = (Real)(nx*ny*nz);
	  for(num_prop=0;num_prop<10;num_prop++) norm_fac[num_prop] = space_vol;
	  norm_fac[4] *= 3.0;
	  norm_fac[5] *= 3.0;
	  norm_fac[8] *= 3.0;
	  norm_fac[9] *= 3.0;

	  /* print relativistic meson propagators */
	  for(num_prop=0;num_prop<10;num_prop++)
	    for(i=0;i<num_kap;i++)
	      for(j=0;j<=i;j++)
		if(pmes_prop_done[i][j] == 1){
		  for(t = 0; t < nt; t++){
		    t_off = (t + source_time)%nt;
		    g_floatsum( &pmes_prop[i][j][num_prop][t_off].real );
		    pmes_prop[i][j][num_prop][t_off].real  /= norm_fac[num_prop];
		    g_floatsum( &pmes_prop[i][j][num_prop][t_off].imag );
		    pmes_prop[i][j][num_prop][t_off].imag  /= norm_fac[num_prop];
		    if(this_node == 0)
		      printf("POINT%s %d %d %d  %e %e\n",
			     mes_kind[num_prop],i,j,t,
			     (double)pmes_prop[i][j][num_prop][t_off].real,
			     (double)pmes_prop[i][j][num_prop][t_off].imag);
		  }
		}

	  /* Once printed, this propagator should be neither
             calculated nor printed again */
	  for(i=0;i<num_kap;i++)
	    for(j=0;j<=i;j++)
	      if(pmes_prop_done[i][j] == 1)
		pmes_prop_done[i][j] = 2;
	  

	  /* print non-relativistic baryon propagators */
	  if(this_node == 0)
	    for(kprop=0;kprop<NRPROPS;kprop++)
	      for(i=0;i<num_kap;i++){
		for(j=0;j<i;j++)
		  if(nrbar_prop[i][j][kprop].done==1){
		    for(t = 0; t < nt; t++){
		      t_off = (t + source_time)%nt;
		      /* Periodic boundary conditions - no wraparound sign */
		      printf("%s_NR%s %d %d %d %d  %e %e\n",
			     nr_fb_label[nr_fb],
			     nrbar_prop[i][j][kprop].label,i,j,j,t,
			     (double)nrbar_prop[i][j][kprop].c[t_off].real,
			     (double)nrbar_prop[i][j][kprop].c[t_off].imag);
		    }
		  }
		
		if(nrbar_prop[i][i][kprop].done==1)
		  for(t = 0; t < nt; t++){
		    t_off = (t + source_time)%nt;
		    printf("%s_NR%s %d %d %d %d  %e %e\n",
			   nr_fb_label[nr_fb],
			   nrbar_prop[i][j][kprop].label,i,i,i,t,
			   (double)nrbar_prop[i][i][kprop].c[t_off].real,
			   (double)nrbar_prop[i][i][kprop].c[t_off].imag);
		  }
		
		for(j=i+1;j<num_kap;j++)
		  if(nrbar_prop[i][j][kprop].done==1)
		    for(t = 0; t < nt; t++){
		      t_off = (t + source_time)%nt;
		      printf("%s_NR%s %d %d %d %d  %e %e\n",
			     nr_fb_label[nr_fb],
			     nrbar_prop[i][j][kprop].label,j,j,i,t,
			     (double)nrbar_prop[i][j][kprop].c[t_off].real,
			     (double)nrbar_prop[i][j][kprop].c[t_off].imag);
		    }
	      }

	  
	  /* print H-dibaryon mixed channel propagators */
	  if(this_node == 0)
	    for(kprop=0;kprop<HDIPROPS;kprop++)
	      for(i=0;i<num_kap;i++){
		for(j=0;j<num_kap;j++)if(hdibar_prop[i][j][kprop].done==1){
		  for(t = 0; t < nt; t++){
		    t_off = (t + source_time)%nt;
		    printf("%s_%s %d %d %d %d %e %e\n",
			   nr_fb_label[nr_fb],
			   hdibar_prop[i][j][kprop].label,i,j,j,t,
			   (double)hdibar_prop[i][j][kprop].c[t_off].real,
			   (double)hdibar_prop[i][j][kprop].c[t_off].imag);
		  }
		}
	      }
	} /* Loop over nr forward - backward */

      /* Cleanup */
      for(kprop=0;kprop<NRPROPS;kprop++)
	for(i=0;i<num_kap;i++)for(j=0;j<num_kap;j++){
	  free(nrbar_prop[i][j][kprop].c);
	  free(nrbar_prop[i][j][kprop].label);
	}
      
      for(kprop=0;kprop<HDIPROPS;kprop++)
	for(kh=0;kh<num_kap_heavy;kh++)for(kl=0;kl<num_kap_light;kl++){
	  i = kh; j = kl + num_kap_heavy;
	  free(hdibar_prop[i][j][kprop].c);
	  free(hdibar_prop[i][j][kprop].label);
	}
      
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
    return 0;
} /* control_H_cl */
