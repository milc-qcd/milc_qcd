/***************** control_eig.c *****************************************/
/* MIMD version 6 */

/* Main procedure for quenched SU3 Wilson fermions 			*/
/* MIMD version 6 */

/* This version computes propagators for hypercubic fermions on a
 supplied background field config */

/* Modifications ... */


#define CONTROL
#include "arb_dirac_eig_includes.h"
#include <string.h>

/* Comment these out if you want to suppress detailed timing */
#define IOTIME
#define PRTIME
#define JACOBI_TOL 1.110223e-16

int main(int argc, char *argv[])
{
int meascount;
int prompt;
Real avm_iters,avs_iters;

double ssplaq,stplaq;


double starttime,endtime;
double dtime;

int MinCG,MaxCG;
Real size_r,RsdCG;

register int i,j,l;
register site *s;

int spinindex,spin,color,k,kk,t;
int flag;
int ci,si,sf,cf;
int num_prop;
Real space_vol;

int status;

int source_chirality;

    wilson_vector **eigVec ;
    double *eigVal ;
    int total_R_iters ;
    double norm;
    Real re,im,re5,im5;
    complex cc;
    char label[20] ;

    double *grad, *err, max_error;
  Matrix Array,V ;

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
                /* call plaquette measuring process */
                d_plaquette(&ssplaq,&stplaq);
                if(this_node==0)printf("FATPLAQ  %e %e\n",
                    (double)ssplaq,(double)stplaq);




/* flip the time oriented fat links 
if(this_node==0) printf("Periodic time BC\n");
*/
if(this_node==0) printf("AP time BC\n");
boundary_flip(MINUS);





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


                eigVal = (double *)malloc(Nvecs*sizeof(double));
                eigVec = (wilson_vector **)malloc(Nvecs*sizeof(wilson_vector*));
                for(i=0;i<Nvecs;i++)
                  eigVec[i]=
                    (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));


          /* open files for wilson propagators */
          fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
          fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k],
				     wqs[k].type);


                    if(startflag_w[k] == FRESH)flag = 0;
                    else
                      flag = 1;
spin=color=0; /* needed by wilson writing routines */




		/* initialize the CG vectors */
		    if(flag==0){
if(this_node==0) printf("random (but chiral) initial vectors\n");
  /* Initiallize all the eigenvectors to a random vector */
  for(j=0;j<Nvecs;j++)
    {
      if(j< Nvecs/2){ source_chirality=1;}
      else{source_chirality= -1;}
      printf("source chirality %d\n",source_chirality);
      grsource_w();
      FORALLSITES(i,s){
        copy_wvec(&(s->g_rand),&(eigVec[j][i]));
	if(source_chirality==1){
	  for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	    eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);
	}

	if(source_chirality== -1){
	  for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	    eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);
	}
      }
      eigVal[j]=1.0e+16;
    }
		    }
		    else{
if(this_node==0) printf("reading in %d wilson_vectors--must be <= 12\n",Nvecs);
                    /* load psi if requested */
for(j=0;j<Nvecs;j++){
printf("reading %d %d %d\n",j,spin,color);
#ifdef IOTIME
                    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
                                      spin, color, F_OFFSET(psi),1);
#else
                    status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k], 
                                      spin, color, F_OFFSET(psi),0);
#endif

		    /* compute eigenvalue */
		    herm_delt(F_OFFSET(psi),F_OFFSET(chi));

		  re=im=0.0;
		  FORALLSITES(i,s){
		    cc = wvec_dot( &(s->chi), &(s->psi) );
		    re += cc.real ;
		  }
		  g_floatsum(&re);
		  eigVal[j]=re;
printf("trial eigenvalue of state %d %e\n",j,eigVal[j]);
		  FORALLSITES(i,s){eigVec[j][i]=s->psi;}
spin++;
if((spin %4) == 0){spin=0;color++;}
}

		    }



                total_R_iters=Kalkreuter(eigVec, eigVal, eigenval_tol, 
                                         error_decr, Nvecs, MaxIter, Restart, 
                                         Kiters, EVENANDODD) ;

		/* save the eigenvectors if requested */
spin=color=0;
                for(j=0;j<Nvecs;j++){
		  FORALLSITES(i,s){s->psi=eigVec[j][i];}
#ifdef IOTIME
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    spin,color,F_OFFSET(psi),1);
#else
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    spin,color,F_OFFSET(psi),0);
#endif
spin++;
if((spin %4) == 0){spin=0;color++;}
		}
          /* close files for wilson propagators */
          r_close_wprop(startflag_w[k],fp_in_w[k]);
          w_close_wprop(saveflag_w[k],fp_out_w[k]);





 
		/* now we do tests on eigenvectors */
                for(j=0;j<Nvecs;j++){

		  FORALLSITES(i,s){s->chi=eigVec[j][i];}

		  /* expectation value of action r=g5 dslash chi*/
		  delta0(F_OFFSET(chi),F_OFFSET(psi),PLUS);
		  FORALLSITES(i,s){
		    mult_by_gamma(&(s->psi),&(s->r),GAMMAFIVE);
		  }
		  re=im=0.0;
		  FORALLSITES(i,s){
		    cc = wvec_dot( &(s->chi), &(s->r) );
		    re += cc.real ;
		    im += cc.imag ;
		  }
		  g_floatsum(&re);
		  g_floatsum(&im);
		  
		  /* expectation value of gamma-5 */
		  re5=im5=0.0;
		  FORALLSITES(i,s){
		    mult_by_gamma(&(s->chi),&(s->r),GAMMAFIVE);
		    cc = wvec_dot( &(s->chi), &(s->r) );
		    re5 += cc.real ;
		    im5 += cc.imag ;
		  }
		  g_floatsum(&re5);
		  g_floatsum(&im5);

		  if(this_node==0) printf("F3MES %d %e %e %e %e\n",
					  j,re,im,re5,im5);
		}
		/* re-diagonalize H(m) */

  /** Allocate the array **/
  Array = AllocateMatrix(Nvecs) ;
  /** Allocate the Eigenvector matrix **/
  V = AllocateMatrix(Nvecs) ;



  for(j=0;j<Nvecs;j++){
    FORALLSITES(i,s){s->chi=eigVec[j][i];}
    delta0(F_OFFSET(chi),F_OFFSET(psi),PLUS);
    FORALLSITES(i,s) mult_by_gamma(&(s->psi),&(s->r),GAMMAFIVE);

    for(kk=0;kk<=j;kk++){
      FORALLSITES(i,s){s->psi=eigVec[kk][i];}
		  re=im=0.0;
		  FORALLSITES(i,s){
		    cc = wvec_dot( &(s->psi), &(s->r) );
		    re += cc.real ;
		    im += cc.imag ;
		  }
		  g_floatsum(&re);
		  g_floatsum(&im);
		  Array.M[j][kk].real = re ; Array.M[j][kk].imag = -im ;
		  Array.M[kk][j].real = re ; Array.M[kk][j].imag = im ;
    }
  }

      Jacobi(&Array, &V, JACOBI_TOL) ;
      RotateBasis(eigVec,&V,EVENANDODD) ;


		/* now we do tests on eigenvectors */
                for(j=0;j<Nvecs;j++){
		  eigVal[j] = Array.M[j][j].real ;
		  
		  /* expectation value of gamma-5 */
		  re5=im5=0.0;
		  FORALLSITES(i,s){
		    mult_by_gamma(&(eigVec[j][i]),&(s->r),GAMMAFIVE);
		    cc = wvec_dot( &(eigVec[j][i]), &(s->r) );
		    re5 += cc.real ;
		    im5 += cc.imag ;
		  }
		  g_floatsum(&re5);
		  g_floatsum(&im5);

		  if(this_node==0) printf("F3D5MES %d %e %e %e\n",
					  j,eigVal[j],re5,im5);
		}
  /** Deallocate the arrays **/
  deAllocate(&V) ;
  deAllocate(&Array) ;

                for(i=0;i<Nvecs;i++)
                  free(eigVec[i]) ;
                free(eigVec) ;
                free(eigVal) ;

		/* end of measurements */


      if(this_node==0)printf("total Rayleigh iters = %d\n",total_R_iters);
	  

	} /* masses */
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
