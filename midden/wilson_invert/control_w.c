/***************** control_w.c *****************************************/

/* Main procedure for quenched SU3 Wilson fermions 			*/
/* MIMD version 6 */

/* This version computes propagators for Wilson fermions on a
 supplied background field config */

/* Modifications ...

   8/15/96 Made scratch file name a variable C.D. 
   8/10/96 Installed new propagator IO and added timing C.D. */

#define CONTROL
#include "wi_inv_includes.h"
#include <string.h>

/* Comment these out if you want to suppress detailed timing */
#define IOTIME
#define PRTIME

int main(int argc, char *argv[])
{
int meascount;
int prompt;
Real avm_iters,avs_iters;

double starttime,endtime;
double dtime;

int MinCG,MaxCG;
Real RsdCG;

register int i;
register site *s;

int spinindex,spin,color,k,t;
int flag;
int ci,si,sf,cf;
int num_prop;
Real space_vol;

int status;

int key[4];
#define restrict rstrict /* C-90 T3D cludge */
int restrict[4];

Real norm_fac[10];

static char *mes_kind[10] = {"PION","PS505","PS055","PS0505",
		"RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};
static char *bar_kind[4] = {"PROTON","PROTON0","DELTA","DELTA0"};

complex *pmes_prop[MAX_KAP][10];
complex *smes_prop[MAX_KAP][10];
complex *bar_prop[MAX_KAP][4];

w_prop_file *fp_in_w[MAX_KAP];        /* For propagator files */
w_prop_file *fp_out_w[MAX_KAP];       /* For propagator files */

    initialize_machine(argc,argv);

    g_sync();
    /* set up */
    prompt = setup_w();
    /* loop over input sets */

    /* Set up Fourier transform */
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 0;
    setup_restrict_fourier(key, restrict);

    while( readin(prompt) == 0)
    {

	starttime=dclock();
	MaxCG=niter;

	avm_iters=0.0;
	meascount=0;

	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    pmes_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    smes_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    for(t=0;t<nt;t++){
		pmes_prop[i][num_prop][t] = cmplx(0.0,0.0); 
		smes_prop[i][num_prop][t] = cmplx(0.0,0.0); 
	    }
	}

	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    bar_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    for(t=0;t<nt;t++)bar_prop[i][num_prop][t] = cmplx(0.0,0.0); 
	}

	if( fixflag == COULOMB_GAUGE_FIX)
	  {
	    if(this_node == 0) 
	      printf("Fixing to Coulomb gauge\n");
#ifdef IOTIME
	    dtime = -dclock();
#endif
	    gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL,
		     F_OFFSET(mp),F_OFFSET(chi),0,NULL,NULL,0,NULL,NULL);
#ifdef IOTIME
	    dtime += dclock();
	    if(this_node==0)printf("Time to gauge fix = %e\n",dtime);
#endif
	  }
	else
	  if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");

        /* save lattice if requested */
        if( saveflag != FORGET ){
          save_lattice( saveflag, savefile );
        }

	if(this_node==0)printf("END OF HEADER\n");

/*	if(this_node==0) printf("num_kap = %d\n", num_kap); */
	/* Loop over kappas */
	for(k=0;k<num_kap;k++){

	  kappa=kap[k];
	  RsdCG=resid[k];
	  if(this_node==0)printf("Kappa= %g r0= %g residue= %g\n",
		(double)kappa,(double)wqs[k].r0,(double)RsdCG);

	  /* open files for wilson propagators */
	  
	  fp_in_w[k]  = r_open_prop(startflag_w[k], startfile_w[k]);
	  fp_out_w[k] = w_open_prop(saveflag_w[k],  savefile_w[k] );


	    /* Loop over source colors */
	    for(color=0;color<3;color++){

		/* Loop over source spins */
		for(spinindex=0;spinindex<n_spins;spinindex++){
		    spin = spins[spinindex];

		    meascount ++;
/*		    if(this_node==0)printf("color=%d spin=%d\n",color,spin); */
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
#ifdef IOTIME
		    status = reload_propagator( startflag_w[k], fp_in_w[k], 
				      spin, color, F_OFFSET(psi),1);
#else
		    status = reload_propagator( startflag_w[k], fp_in_w[k], 
				      spin, color, F_OFFSET(psi),0);
#endif
	    if(status != 0)
	      {
		node0_printf("control_w: Recovering from error by resetting initial guess to zero\n");
		reload_propagator( FRESH, fp_in_w[k], 
			       spin, color, F_OFFSET(psi),0);
		flag = 0;
	      }


	    /* Complete the source structure */
	    wqs[k].color = color;
	    wqs[k].spin = spin;

	    /* For wilson_info */
	    wqstmp = wqs[k];

	   /* If we are starting afresh, we set a minimum number
	      of iterations */
	   if(startflag_w[k] == FRESH || status != 0)MinCG = nt/2; 
	   else MinCG = 0;

	    /* Load inversion control structure */
	    qic.min = MinCG;
	    qic.max = MaxCG;
	    qic.nrestart = nrestart;
	    qic.resid = RsdCG;
	    qic.start_flag = flag;
	    
	    /* Load Dirac matrix parameters */
	    dwp.Kappa = kappa;

#ifdef BI
	    /* Load temporaries specific to inverter */
	    qic.wv2 = F_OFFSET(mp);
	    qic.wv3 = F_OFFSET(tmp);  /* Called rv in bicg */
	    qic.wv4 = F_OFFSET(sss);

	    /* compute the propagator.  Result in psi. */
	    avs_iters = 
	      (Real)wilson_invert_lean(F_OFFSET(chi),F_OFFSET(psi),
					w_source,&wqs[k],
					bicgilu_w,&qic,(void *)&dwp);
#else
	    /* Load temporaries specific to inverter */
	    qic.wv2 = F_OFFSET(mp);
	    
	    /* compute the propagator.  Result in psi. */
	    avs_iters = 
	      (Real)wilson_invert_lean(F_OFFSET(chi),F_OFFSET(psi),
					w_source,&wqs[k],
					cgilu_w,&qic,(void *)&dwp);
#endif
	    avm_iters += avs_iters;
		    FORALLSITES(i,s)
			copy_wvec(&(s->psi),
			    &(s->quark_propagator.c[color].d[spin]));

		    /* save psi if requested */
#ifdef IOTIME
		    save_propagator( saveflag_w[k],fp_out_w[k],
				    spin,color,F_OFFSET(psi),1);
#else
		    save_propagator( saveflag_w[k],fp_out_w[k],
				    spin,color,F_OFFSET(psi),0);
#endif
		} /* source spins */
	    } /* source colors */

	  /* close files for wilson propagators */
	  r_close_prop(startflag_w[k],fp_in_w[k]);
	  w_close_prop(saveflag_w[k],fp_out_w[k]);
	  
	  /* spectrum */
#ifdef PRTIME
	    dtime = -dclock();
#endif
	    w_baryon(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator),
		F_OFFSET(quark_propagator), bar_prop[k]);
#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for diagonal baryons %e\n",dtime);
		fflush(stdout);
	      }
#endif

#ifdef PRTIME
	    dtime = -dclock();
#endif
	    for(color=0;color<3;color++){
		w_meson(F_OFFSET(quark_propagator.c[color]),
		    F_OFFSET(quark_propagator.c[color]), pmes_prop[k]);
	    }
#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for diagonal mesons %e\n",dtime);
		fflush(stdout);
	      }
#endif
#ifdef PRTIME
	    dtime = -dclock();
#endif

	    /* Now convolute the quark propagator with a Gaussian for
	       the smeared mesons. This is done with FFT's */

	    /* fft quark_propagator (in place)--use mp as working space */
	    for(color=0;color<3;color++)
	      for(spinindex=0;spinindex<n_spins;spinindex++){
		spin = spins[spinindex];
		/* Use chi as a second working space */
		restrict_fourier(F_OFFSET(quark_propagator.c[color].d[spin]),
		    F_OFFSET(mp), F_OFFSET(chi),
		    sizeof(wilson_vector), FORWARDS);
	      }

#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for 12 FFT's of quark %d %e\n",k,dtime);
		fflush(stdout);
	      }
#endif
	    /* Use chi, for spin=color=0 for the sink wave function */
	    FORALLSITES(i,s) clear_wvec( &(s->chi));

	  wqstmp = wqs[k];
	  spin=0;color=0;
	  wqstmp.spin=spin; wqstmp.color=color;
	  wqstmp.type=GAUSSIAN;
	  w_sink(F_OFFSET(chi),&wqstmp);
	  
	    /* We want chi(-k)* -- the complex conjugate of FFT of the
	       complex conjugate of the quark sink. */
	    FORALLSITES(i,s){
		CONJG(s->chi.d[spin].c[color], s->chi.d[spin].c[color]);
	    }
#ifdef PRTIME
	    dtime = -dclock();
#endif
	    /* Use htmp[0].h[0].c[0] for scratch space */
	    restrict_fourier(F_OFFSET(chi.d[spin].c[color]),
		F_OFFSET(mp.d[spin].c[color]),
		F_OFFSET(htmp[0].h[0].c[0]), sizeof(complex), FORWARDS);
	    FORALLSITES(i,s){
		CONJG(s->chi.d[spin].c[color], s->chi.d[spin].c[color]);
	    }

	    /* Now multiply quark by sink wave function */
	    FORALLSITES(i,s)
	    for(ci=0;ci<3;ci++)for(si=0;si<4;si++)
	    for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
		CMUL(s->quark_propagator.c[ci].d[si].d[sf].c[cf],
		    s->chi.d[spin].c[color],
		    s->quark_propagator.c[ci].d[si].d[sf].c[cf]);
	    }

#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for FFT of chi and multiply %d %e\n",k,dtime);
		fflush(stdout);
	      }
#endif
	    /* Note: we do not have to FFT the quark propagator back,
	       as long as we don't compute smeared baryons,
	       since the meson construction works the same in momentum space!
	       We do need to divide by an additional (space) volume
	       factor, though! */

	    for(color=0;color<3;color++){
		w_meson(F_OFFSET(quark_propagator.c[color]),
		    F_OFFSET(quark_propagator.c[color]), smes_prop[k]);
	    }

	} /* kappas */

	space_vol = (Real)(nx*ny*nz);
	for(num_prop=0;num_prop<10;num_prop++) norm_fac[num_prop] = space_vol;
	norm_fac[4] *= 3.0;
	norm_fac[5] *= 3.0;
	norm_fac[8] *= 3.0;
	norm_fac[9] *= 3.0;

	/* print meson propagators */
	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    for(t=0; t<nt; t++){
		g_floatsum( &pmes_prop[i][num_prop][t].real );
		pmes_prop[i][num_prop][t].real  /= norm_fac[num_prop];
		g_floatsum( &pmes_prop[i][num_prop][t].imag );
		pmes_prop[i][num_prop][t].imag  /= norm_fac[num_prop];
		if(this_node == 0)
		    printf("POINT%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			(double)pmes_prop[i][num_prop][t].real,
			(double)pmes_prop[i][num_prop][t].imag);
	    }
	    for(t=0; t<nt; t++){
		g_floatsum( &smes_prop[i][num_prop][t].real );
		smes_prop[i][num_prop][t].real  /=
		    (space_vol*norm_fac[num_prop]);
		g_floatsum( &smes_prop[i][num_prop][t].imag );
		smes_prop[i][num_prop][t].imag  /=
		    (space_vol*norm_fac[num_prop]);
		if(this_node == 0)
		    printf("SMEAR%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			(double)smes_prop[i][num_prop][t].real,
			(double)smes_prop[i][num_prop][t].imag);
	    }
	}

	/* print baryon propagators */
	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    for(t=0; t<nt; t++){
		g_floatsum( &bar_prop[i][num_prop][t].real );
		bar_prop[i][num_prop][t].real  /= space_vol;
		g_floatsum( &bar_prop[i][num_prop][t].imag );
		bar_prop[i][num_prop][t].imag  /= space_vol;
		if(this_node == 0)
		    printf("POINT%s %d %d  %e %e\n",bar_kind[num_prop],i,t,
			(double)bar_prop[i][num_prop][t].real,
			(double)bar_prop[i][num_prop][t].imag);
	    }
	}

	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    free(pmes_prop[i][num_prop]);
	    free(smes_prop[i][num_prop]);
	}

	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    free(bar_prop[i][num_prop]);
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
}

