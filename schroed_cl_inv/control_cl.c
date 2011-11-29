/***************** control_cl.c *****************************************/

/* Main procedure for quenched SU3 clover fermions 			*/
/* MIMD version 7 */

/* This version computes propagators for clover fermions on a
 supplied background field config with Schroedinger functional boundary
 conditions */

#define CONTROL
#include "schroed_cl_includes.h"

/* Comment these out if you want to suppress detailed timing */
#define PRTIME

#define MAX_PROP 2

int main(int argc, char *argv[])  {
  int meascount,todo;
  int prompt;
  Real avm_iters,avs_iters;
  
  double starttime,endtime;
  
  int MaxCG;
  Real RsdCG;
  
  int cl_cg = CL_CG;
  
  register int i;
  register site *s;
  
  int spin,color,k,t;
  int flag;
  int num_prop;
  Real space_vol;
  
  double dssplaq, dstplaq;
  
  /**
     Real norm_fac[5];
     static char *mes_kind[5] = {"PION","PS505","RHO33","SCALAR","PV35"};
     complex *pmes_prop_f[MAX_KAP][5];
     complex *pmes_prop_b[MAX_KAP][5];
  **/
  Real norm_fac[MAX_PROP];
  static char *mes_kind[MAX_PROP] = {"PION","PS505"};
  complex *pmes_prop_f[MAX_KAP][MAX_PROP];
  complex *pmes_prop_b[MAX_KAP][MAX_PROP];
  complex *f_V[MAX_KAP];
  Real f_1[MAX_KAP];
  
  initialize_machine(&argc,&argv);
  
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup_cl();
  /* loop over input sets */
  
  while( readin(prompt) == 0)
    {
      
      starttime=dclock();
      MaxCG=niter;
      
      avm_iters=0.0;
      meascount=0;
      node0_printf("END OF HEADER\n");
      
      if( startflag != CONTINUE && num_smear > 0){
	for(todo = num_smear; todo > 0; --todo ){
	  ape_smear_SF();
	  d_plaquette(&dssplaq, &dstplaq);
	  i = num_smear - todo + 1;
	  node0_printf("At smearing %d plaq: %e %e\n",
		       i, dssplaq, dstplaq);
	}
	node0_printf("SMEARING COMPLETED\n");
      }
      
      for(num_prop=0;num_prop<MAX_PROP;num_prop++)
	for(i=0;i<num_kap;i++){
	  pmes_prop_f[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	  pmes_prop_b[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++){
	    pmes_prop_f[i][num_prop][t] = cmplx(0.0,0.0); 
	    pmes_prop_b[i][num_prop][t] = cmplx(0.0,0.0); 
	  }
	}
      
      for(i=0;i<num_kap;i++){
	f_V[i] = (complex *)malloc(nt*sizeof(complex));
	f_1[i] = 0.0;
	for(t=0;t<nt;t++){
	  f_V[i][t] = cmplx(0.0,0.0);
	}
      }
      
      space_vol = (Real)(nx*ny*nz);
      for(num_prop=0;num_prop<MAX_PROP;num_prop++)
	norm_fac[num_prop] = 2.0*space_vol;
      
      norm_fac[1] *= -1.0;
      
      i = 2;
      if(MAX_PROP > i) norm_fac[i] *= 3.0;
      i = 4;
      if(MAX_PROP > i) norm_fac[i] *= 3.0;
      
      node0_printf("num_kap = %d\n", num_kap);
      /* Loop over kappas */
      for(k=0;k<num_kap;k++){
	
	kappa=kap[k];
	RsdCG=resid[k];
	node0_printf("Kappa=%e residue=%e\n", (double)kappa,(double)RsdCG);
	
	/* Loop over source colors */
	for(color=0;color<3;color++){
	  
	  /* Loop over source spins for "forward" propagation */
	  /* We need only two source spins, since the source is
	     proportional to a projector */
	  for(spin=0;spin<2;spin++){
	    
	    meascount ++;
	    node0_printf("color=%d spin=%d\n",color,spin);
	    wqs[k].color = color;
	    wqs[k].spin = spin;
	    wqs[k].type = PLUS;
	    wqs[k].kappa = kappa;
	    
	    if(k==0){
	      flag = 0;      /* Saves one multiplication in cgilu */
	      FORALLSITES(i,s)clear_wvec( &(s->psi));
	    }
	    else{
	      flag = 1;
	      FORALLSITES(i,s)
		copy_wvec(&(s->forw_quark_prop.c[color].d[spin]),
			  &(s->psi));
	    }
	    
	    /* Load inversion control structure */
	    qic.prec = PRECISION;
	    qic.min = 0;
	    qic.max = MaxCG;
	    qic.nrestart = nrestart;
	    qic.parity = EVENANDODD;
	    qic.start_flag = flag;
	    qic.nsrc = 1;
	    qic.resid = RsdCG;
	    qic.relresid = 0;
	    
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa;
	    dcp.Clov_c = clov_c;
	    dcp.U0 = 1.0;

	    /* Generate the source in chi */
	    w_source_sf_site(F_OFFSET(chi),&wqs[k]);
	    
	    /* compute the propagator.  Result in psi. */
	    switch (cl_cg) {
	    case BICG:
	      avs_iters = 
		(Real)bicgilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     bicgilu_cl_site,&qic,(void *)&dcp);
	      break;
	    case MR:
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     mrilu_cl_site,&qic,(void *)&dcp);
	      avs_iters = 
		(Real)mrilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
	      break;
	    case CG:
	      avs_iters = 
		(Real)cgilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     cgilu_cl_site,&qic,(void *)&dcp);
	      break;
	    default:
	      node0_printf("main(%d): Inverter choice %d not supported\n",
			   this_node,cl_cg);
	    }
	    
	    avm_iters += avs_iters;
	    node0_printf("size_r= %e, iters= %e\n",
			 (double)(qic.size_r), (double)avs_iters);
	    
	    /* In MILC conventions PLUS stands for P_- */
	    FORALLSITES(i,s){
	      copy_wvec(&(s->psi),
			&(s->forw_quark_prop.c[color].d[spin]));
	      scalar_mult_wvec(&(s->psi), -1.0,
			       &(s->forw_quark_prop.c[color].d[spin+2]));
	    }
	    
	  } /* source spins */
	  
	  /* spectrum */
	  schroed_meson(F_OFFSET(forw_quark_prop.c[color]),
			F_OFFSET(forw_quark_prop.c[color]), pmes_prop_f[k],
			MAX_PROP);
	  
	  
	  /* Loop over source spins for "backward" propagation */
	  /* We need only two source spins, since the source is
	     proportional to a projector */
	  for(spin=0;spin<2;spin++){
	    
	    meascount ++;
	    node0_printf("color=%d spin=%d\n",color,spin);
	    wqs[k].color = color;
	    wqs[k].spin = spin;
	    wqs[k].type = MINUS;
	    wqs[k].kappa = kappa;
	    
	    if(k==0){
	      flag = 0;      /* Saves one multiplication in cgilu */
	      FORALLSITES(i,s)clear_wvec( &(s->psi));
	    }
	    else{
	      flag = 1;
	      FORALLSITES(i,s)
		copy_wvec(&(s->backw_quark_prop.c[color].d[spin]),
			  &(s->psi));
	    }
	    
	    /* Load inversion control structure */
	    qic.max = MaxCG;
	    qic.nrestart = nrestart;
	    qic.resid = RsdCG;
	    qic.start_flag = flag;
	    
	    /* Load Dirac matrix parameters */
	    dcp.Kappa = kappa;
	    dcp.Clov_c = clov_c;
	    dcp.U0 = 1.0;
	    
	    /* Generate the source in chi */
	    w_source_sf_site(F_OFFSET(chi),&wqs[k]);

	    /* compute the propagator.  Result in psi. */
	    switch (cl_cg) {
	    case BICG:
	      avs_iters = 
		(Real)bicgilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
//
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     bicgilu_cl_site,&qic,(void *)&dcp);
	      break;
	    case MR:
	      avs_iters = 
		(Real)mrilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
//	      /* compute the propagator.  Result in psi. */
//	      
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     mrilu_cl_site,&qic,(void *)&dcp);
	      break;
	    case CG:
	      avs_iters = 
		(Real)cgilu_cl_site(F_OFFSET(chi),F_OFFSET(psi),&qic,(void *)&dcp);
//	      /* compute the propagator.  Result in psi. */
//	      
//	      avs_iters = 
//		(Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
//					     w_source_sf_site,&wqs[k],
//					     cgilu_cl_site,&qic,(void *)&dcp);
	      break;
	    default:
	      node0_printf("main(%d): Inverter choice %d not supported\n",
			   this_node,cl_cg);
	    }
	    avm_iters += avs_iters;
	    node0_printf("size_r= %e, iters= %e\n",
			 (double)(qic.size_r), (double)avs_iters);
	    
	    /* In MILC conventions MINUS stands for P_+ */
	    FORALLSITES(i,s){
	      copy_wvec(&(s->psi),
			&(s->backw_quark_prop.c[color].d[spin]));
	      copy_wvec(&(s->psi),
			&(s->backw_quark_prop.c[color].d[spin+2]));
	    }
	    
	  } /* source spins */
	  
	  /* spectrum */
	  schroed_meson(F_OFFSET(backw_quark_prop.c[color]),
			F_OFFSET(backw_quark_prop.c[color]), pmes_prop_b[k],
			MAX_PROP);
	  
	} /* source colors */
	
	
	/* Measurements for Z_V */
	zv_meas(F_OFFSET(forw_quark_prop), F_OFFSET(backw_quark_prop),
		&f_1[k], f_V[k], kappa);
	
      } /* kappas */
      
	/* print forward meson propagators */
      for(num_prop=0;num_prop<MAX_PROP;num_prop++)
	for(i=0;i<num_kap;i++){
	  for(t=0; t<nt; t++){
	    g_floatsum( &pmes_prop_f[i][num_prop][t].real );
	    pmes_prop_f[i][num_prop][t].real  /= norm_fac[num_prop];
	    g_floatsum( &pmes_prop_f[i][num_prop][t].imag );
	    pmes_prop_f[i][num_prop][t].imag  /= norm_fac[num_prop];
	    node0_printf("FORW_%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			 (double)pmes_prop_f[i][num_prop][t].real,
			 (double)pmes_prop_f[i][num_prop][t].imag);
	  }
	}
      
      norm_fac[1] *= -1.0;
      
      /* print backward meson propagators */
      for(num_prop=0;num_prop<MAX_PROP;num_prop++)
	for(i=0;i<num_kap;i++){
	  for(t=0; t<nt; t++){
	    if(t==0) k=t; else k=nt-t;
	    g_floatsum( &pmes_prop_b[i][num_prop][k].real );
	    pmes_prop_b[i][num_prop][k].real  /= norm_fac[num_prop];
	    g_floatsum( &pmes_prop_b[i][num_prop][k].imag );
	    pmes_prop_b[i][num_prop][k].imag  /= norm_fac[num_prop];
	    node0_printf("BACKW_%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			 (double)pmes_prop_b[i][num_prop][k].real,
			 (double)pmes_prop_b[i][num_prop][k].imag);
	  }
	}
      
      /* print Z_V measurements stuff */
      for(i=0;i<num_kap;i++){
	node0_printf("ZV_f_1 %d  %e\n", i, (double)f_1[i]);
	for(t=0; t<nt; t++){
	  g_floatsum( &f_V[i][t].real );
	  f_V[i][t].real /= space_vol;
	  g_floatsum( &f_V[i][t].imag );
	  f_V[i][t].imag /= space_vol;
	  node0_printf("ZV_f_V %d %d  %e %e\n", i, t,
		       (double)f_V[i][t].real, (double)f_V[i][t].imag);
	}
      }
      
      for(num_prop=0;num_prop<MAX_PROP;num_prop++)
	for(i=0;i<num_kap;i++){
	  free(pmes_prop_f[i][num_prop]);
	  free(pmes_prop_b[i][num_prop]);
	}
      
      for(i=0;i<num_kap;i++){
	free(f_V[i]);
      }
      
      node0_printf("RUNNING COMPLETED\n");
      if(meascount>0){
	node0_printf("total cg iters for measurement= %e\n",
		     (double)avm_iters);
	node0_printf("cg iters for measurement= %e\n",
		     (double)avm_iters/(double)meascount);
      }
      
      endtime=dclock();
      node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
      node0_printf("total_iters = %d\n",total_iters);
      fflush(stdout);
      
    }
  return 0;
}
