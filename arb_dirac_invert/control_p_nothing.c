/***************** control_p_nothing.c *****************************************/
/* MIMD version 6 */

/* Main procedure for quenched SU3 Wilson fermions 			*/
/* MIMD version 6 */

/* This version computes propagators for hypercubic fermions on a
 supplied background field config */

/* Modifications ... */


#define CONTROL
#include "arb_dirac_inv_includes.h"
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
Real size_r,RsdCG;

register int i;
register site *s;

int spinindex,spin,color,j,k,t;
int flag;
int ci,si,sf,cf;
int num_prop;
Real space_vol, inv_space_vol;

int status;



int key[4];
#define restrict rstrict /* C-90 T3D cludge */
int restrict[4];

Real norm_fac[10];



w_prop_file *fp_in_w[MAX_MASSES];        /* For propagator files */
w_prop_file *fp_out_w[MAX_MASSES];       /* For propagator files */

    initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

    g_sync();
    /* set up */
    prompt = setup_p();
    /* loop over input sets */


    while( readin(prompt) == 0)
    {



	starttime=dclock();
	MaxCG=niter;

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
	  }
	else
	  if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");

        /* save lattice if requested */
        if( saveflag != FORGET ){
	  /* Note: beta, kappa are kept only for save_old_binary */
	  save_lattice( saveflag, savefile, stringLFN );
        }

	if(this_node==0)printf("END OF HEADER\n");
	setup_offset();
/*
printf("warning: no fattening of links\n");
*/
	monte_block_ape_b(1);


	setup_links(SIMPLE);

/*	if(this_node==0) printf("num_masses = %d\n", num_masses); */
	/* Loop over mass */
	for(k=0;k<num_masses;k++){

	  m0=mass[k];
/* code to exit... */
	if(m0 <= -10.0){
        fflush(stdout);
	 exit(1);}

	  RsdCG=resid[k];
	  if(this_node==0)printf("mass= %g r0= %g residue= %g\n",
		(double)m0,(double)wqs[k].r0,(double)RsdCG);
	  build_params(m0);
	  make_clov1();



	  /* open files for wilson propagators */
	  
	  fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
	  fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				     wqs[k].type);


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
		    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
				      spin, color, F_OFFSET(psi),1);
#else
		    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
				      spin, color, F_OFFSET(psi),0);
#endif
	    if(status != 0)
	      {
		node0_printf("control_w: Recovering from error by resetting initial guess to zero\n");
		reload_wprop_sc_to_site( FRESH, fp_in_w[k], 
			       spin, color, F_OFFSET(psi),0);
		flag = 0;
	      }


	    /* Conjugate gradient inversion uses site structure
	       temporaries "tmpb" and "chi" */
	    /* Complete the source structure */
	    wqs[k].color = color;
	    wqs[k].spin = spin;

	    /* For wilson_info */
	    wqstmp = wqs[k];

	   /* If we are starting afresh, we set a minimum number
	      of iterations */
	   if(startflag_w[k] == FRESH || status != 0)MinCG = nt/2; 
	   else MinCG = 0;
      /* Make source */
      w_source(F_OFFSET(chi),&wqs[k]);

	    /* compute the propagator.  Result in psi. */
	    avs_iters = (Real)bicgstab(F_OFFSET(chi),F_OFFSET(psi),
				      MaxCG,RsdCG,&size_r,flag);
	    avm_iters += avs_iters;
		    FORALLSITES(i,s)
			copy_wvec(&(s->psi),
			    &(s->quark_propagator.c[color].d[spin]));

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

	  /* close files for wilson propagators */
	  r_close_wprop(startflag_w[k],fp_in_w[k]);
	  w_close_wprop(saveflag_w[k],fp_out_w[k]);
	  
	  /* spectrum */

	} /* masses (k) */


	if(this_node==0)printf("BEGIN PROPAGATORS\n");
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

void fpoly(Real x,Real *p, int np)
{
        int j;

        p[0]=1.0;
        for (j=1;j<np;j++) p[j]=p[j-1]*x;
}
