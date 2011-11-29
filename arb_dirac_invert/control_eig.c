/***************** control_eig.c *****************************************/
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

complex *pmes_prop[MAX_MASSES][10];
complex *smes_prop[MAX_MASSES][10];
complex *bar_prop[MAX_MASSES][4];

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



	if(this_node==0)printf("END OF HEADER\n");
	setup_offset();

/*
if(this_node==0)printf("warning--no fat link\n");
*/
	monte_block_ape_b(1);


/* flip the time oriented fat links */
boundary_flip(MINUS);



#ifdef PAULI
	setup_pauli();
#endif

	setup_links(SIMPLE);

/*	if(this_node==0) printf("num_masses = %d\n", num_masses); */
	/* Loop over mass */
	for(k=0;k<num_masses;k++){

	  m0=mass[k];
	if(m0 <= -10.0) exit(1);
	  RsdCG=resid[k];
	  if(this_node==0)printf("mass= %g r0= %g residue= %g\n",
		(double)m0,(double)wqs[k].r0,(double)RsdCG);
	  build_params(m0);
	  make_clov1();



	  /* open files for wilson propagators */
	  
	  fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
	  fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				     wqs[k].type);

spin=color=0;
 
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
		    node0_printf("BEGIN\n");
                meascount = f_measure2(flag);

		    /* save psi if requested */
#ifdef IOTIME
		    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
				    spin,color,F_OFFSET(psi),1);
#else
		    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
				    spin,color,F_OFFSET(psi),0);
#endif
	}

	  /* close files for wilson propagators */
	  r_close_wprop(startflag_w[k],fp_in_w[k]);
	  w_close_wprop(saveflag_w[k],fp_out_w[k]);
	  

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
