/******* congrad_multi_field.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Calls hdelta0_field.c for fermion matrix */

/* version of 26 May 02 FIELDWISE version  */

/* This is an implementation of a multi-mass conjugate gradient
algorithm a la B. Jegerlehner. The number of masses, Norder,
will be set in lattice.h, and so will the array of mass values,
shift[Norder].  */

/* The source vector is in "src", and the answer (initial guess is zero)
   in "dest". 
   "psim[jmass]" are working vectors for the conjugate gradient.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.

In this version of the code, all the scalars are real because M=D^\dagger D
*/

/* we check all vectors for convergence and quit doing the converged ones.
  */


/* ZOLOTAROV-ONLY VERSION */


#include "arb_ov_includes.h"

#ifdef PREFETCH
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#endif


#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif

void step(
    field_offset src,  /* type wilson_vector (where source is to be created)*/
    field_offset dest  /* type wilson_vector (answer and initial guess) */
)
{
  register int i;
  register site *s;
  int j, avs_iters;
  int MaxCG;
  Real RsdCG;
  Real size_r,r0a;

  wilson_vector **psim;
  wilson_vector *chi_re, *mpm,*pm0;
/* step fn testing */
#ifdef MINN
  double source_norm,dest_norm;
  Real test_epsilon;
#endif


#ifdef EIG
  double_complex cd;
  complex ctmp,*cproj = NULL;
  Real ener;
#endif

#ifdef MINN
  RsdCG=resid_inner_run;
#else
  RsdCG=resid_inner;
#endif
  MaxCG=maxcg_inner;
    
  r0a=1.0/scalez;


  chi_re=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  FORALLSITES(i,s) copy_wvec((wilson_vector *)F_PT(s,src),&(chi_re[i]));
  /* added 6 Sep */
  FORALLSITES(i,s)clear_wvec((wilson_vector *)F_PT(s,dest));

    /* test step fn */
#ifdef MINN
  source_norm=0.0;
  FORALLSITES(i,s){
    source_norm += (double)magsq_wvec((wilson_vector *)F_PT(s,src));
  }
  g_doublesum( &source_norm );
  /*
    node0_printf("step source norm %e\n",source_norm);
  */
#endif

    
#ifdef EIG
  if(Nvecs_h0r0 !=0){
    cproj = (complex *)malloc(Nvecs_h0r0*sizeof(complex));

    /* project out the eigenvectors from the source */
    for(j=0;j<Nvecs_h0r0;j++){
      cd=dcmplx((double)0.0,(double)0.0);
      FORALLSITES(i,s){
	ctmp =  wvec_dot( &(eigVec0[j][i]),(wilson_vector *)F_PT(s,src));
	CSUM(cd,ctmp);
      }
      g_dcomplexsum(&cd);
      cproj[j].real=cd.real;
      cproj[j].imag=cd.imag;

      CMULREAL(cd,-1.0,ctmp);

      FORALLSITES(i,s){
#ifdef PREFETCH
	if(i< sites_on_node - FETCH_UP){
	  prefetch_W(&(chi_re[i+FETCH_UP]));
	  prefetch_W(&(eigVec0[j][i+FETCH_UP]));
	}
#endif
	c_scalar_mult_add_wvec(&(chi_re[i]), &(eigVec0[j][i]),
			       &ctmp, &(chi_re[i]) );
       
      } /* sites */
    
    } /*Nvecs*/
    
  } /*Nvecs_h0r0 !=0 */
#endif




  psim = (wilson_vector **)malloc(Norder*sizeof(wilson_vector*));
  for(i=0;i<Norder;i++)
    psim[i]=
      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  avs_iters = congrad_multi_field( chi_re,psim,
				   MaxCG,RsdCG,&size_r,0); 

  mpm=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));



  FORALLSITES(i,s){
#ifdef PREFETCH
    if(i< sites_on_node - FETCH_UP){
      prefetch_W(&(psim[0][i+FETCH_UP]));
      prefetch_W(&(mpm[i+FETCH_UP]));
      for(j=1;j<Norder;j++)prefetch_W(&(psim[j][i+FETCH_UP]));
    }
#endif
    scalar_mult_wvec(&(psim[0][i]),coeff[0],&(mpm[i]));
    for(j=1;j<Norder;j++){
      scalar_mult_add_wvec(&(mpm[i]),&(psim[j][i]),
			   coeff[j],&(mpm[i]));
    }
  }
  free(chi_re);

  pm0=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  /* now for the numerator */
  delta0_field(mpm,pm0,PLUS);

  free(mpm);

  FORALLSITES(i,s){
#ifdef PREFETCH
    if(i< sites_on_node - FETCH_UP){
      prefetch_W(&(pm0[i+FETCH_UP]));
      prefetch_W((wilson_vector *)F_PT(s+FETCH_UP,dest));
    }
#endif
    mult_by_gamma(&(pm0[i]),(wilson_vector *)F_PT(s,dest),GAMMAFIVE);
    scalar_mult_wvec((wilson_vector *)F_PT(s,dest),
		     r0a,(wilson_vector *)F_PT(s,dest));
  }
  free(pm0);
/*
#ifdef DEBUG
if(this_node==0)printf("inner cg iters %d\n",avs_iters);
if(this_node==0)printf("inner cg res %e\n",size_r);
fflush(stdout);
#endif
*/
  for(i=0;i<Norder;i++) free(psim[i]);
  free(psim) ;


#ifdef EIG
      /* add the ''projector term'' back into dest...*/
  if(Nvecs_h0r0 !=0){
    for(j=0;j<Nvecs_h0r0;j++){
      ener=1.0;
      if(eigVal0[j]<0.0) ener=-1.0;

      CMULREAL(cproj[j],ener,ctmp); 
      FORALLSITES(i,s){
#ifdef PREFETCH
	if(i< sites_on_node - FETCH_UP){
	  prefetch_W(&(eigVec0[j][i+FETCH_UP]));
	  prefetch_W((wilson_vector *)F_PT(s+FETCH_UP,dest));
	}
#endif
	c_scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),&(eigVec0[j][i]),
                             &ctmp,(wilson_vector *)F_PT(s,dest) ); 
      }
    }
    free(cproj);
  }
#endif
    
  /* check the accuracy of the step function and raise or lower the CG accuracy accordingly */
#ifdef MINN
  if(do_minn ==1){
    dest_norm=0.0;
    FORALLSITES(i,s){
      dest_norm += (double)magsq_wvec((wilson_vector *)F_PT(s,dest));
    }
    g_doublesum( &dest_norm );
    /*
    node0_printf("step results %e %e %e",source_norm,dest_norm,dest_norm/source_norm-1.0); 
    */  
    test_epsilon = 0.5*(fabs((dest_norm/source_norm) -1.0));
    /*
      node0_printf("test_epsilon %e\n",test_epsilon);
    */
    /* you can always increase resid_inner--up to a point! */
    if(test_epsilon < resid_inner && resid_inner_run < 0.02 ) resid_inner_run *= 1.2;
    /* but it should not shrink too small */
    if((resid_inner_run >= resid_inner) &&(test_epsilon > resid_inner))
      resid_inner_run /= 1.2;
    /* 
    node0_printf("  resid_inner_run %e\n",resid_inner_run);
    */ 
 /* 
  node0_printf("congrad_multi_field: resid_inner_run %e for resid_inner %e\n",resid_inner_run,resid_inner);
 */ 
  }
#endif
}



void step_field(
    wilson_vector* src,  /* type wilson_vector (where source is to be created)*/
    wilson_vector* dest  /* type wilson_vector (answer and initial guess) */
)
{
  register int i;
  register site *s;
  int j, avs_iters;
  int MaxCG;
  Real RsdCG;
  Real size_r,r0a;

  wilson_vector **psim;
  wilson_vector *chi_re, *mpm,*pm0;
/* step fn testing */
#ifdef MINN
  double source_norm,dest_norm;
  Real test_epsilon;
#endif


#ifdef EIG
  double_complex cd;
  complex ctmp,*cproj = NULL;
  Real ener;
#endif

#ifdef MINN
  RsdCG=resid_inner_run;
#else
  RsdCG=resid_inner;
#endif
  MaxCG=maxcg_inner;
    
  r0a=1.0/scalez;


  chi_re=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  copy_Vector(src,chi_re);

  /* added 6 Sep */
  FORALLSITES(i,s)clear_wvec(&dest[i]);

    /* test step fn */
#ifdef MINN
  source_norm=0.0;
  FORALLSITES(i,s) source_norm += (double)magsq_wvec(&src[i]);
  g_doublesum( &source_norm );
  /*
    node0_printf("step source norm %e\n",source_norm);
  */
#endif

    
#ifdef EIG
  if(Nvecs_h0r0 !=0){
    cproj = (complex *)malloc(Nvecs_h0r0*sizeof(complex));

    /* project out the eigenvectors from the source */
    for(j=0;j<Nvecs_h0r0;j++){
      cd=dcmplx((double)0.0,(double)0.0);
      FORALLSITES(i,s){
	ctmp =  wvec_dot( &(eigVec0[j][i]),&src[i]);
	CSUM(cd,ctmp);
      }
      g_dcomplexsum(&cd);
      cproj[j].real=cd.real;
      cproj[j].imag=cd.imag;

      CMULREAL(cd,-1.0,ctmp);

      FORALLSITES(i,s){
#ifdef PREFETCH
	if(i< sites_on_node - FETCH_UP){
	  prefetch_W(&(chi_re[i+FETCH_UP]));
	  prefetch_W(&(eigVec0[j][i+FETCH_UP]));
	}
#endif
	c_scalar_mult_add_wvec(&(chi_re[i]), &(eigVec0[j][i]),
			       &ctmp, &(chi_re[i]) );
       
      } /* sites */
    
    } /*Nvecs*/
    
  } /*Nvecs_h0r0 !=0 */
#endif




  psim = (wilson_vector **)malloc(Norder*sizeof(wilson_vector*));
  for(i=0;i<Norder;i++)
    psim[i]=
      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  avs_iters = congrad_multi_field( chi_re,psim,
				   MaxCG,RsdCG,&size_r,0); 

  mpm=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));



  FORALLSITES(i,s){
#ifdef PREFETCH
    if(i< sites_on_node - FETCH_UP){
      prefetch_W(&(psim[0][i+FETCH_UP]));
      prefetch_W(&(mpm[i+FETCH_UP]));
      for(j=1;j<Norder;j++)prefetch_W(&(psim[j][i+FETCH_UP]));
    }
#endif
    scalar_mult_wvec(&(psim[0][i]),coeff[0],&(mpm[i]));
    for(j=1;j<Norder;j++){
      scalar_mult_add_wvec(&(mpm[i]),&(psim[j][i]),
			   coeff[j],&(mpm[i]));
    }
  }
  free(chi_re);

  pm0=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  /* now for the numerator */
  delta0_field(mpm,pm0,PLUS);

  free(mpm);

  FORALLSITES(i,s){
#ifdef PREFETCH
    if(i< sites_on_node - FETCH_UP){
      prefetch_W(&(pm0[i+FETCH_UP]));
      prefetch_W(&dest[i+FETCH_UP]);
    }
#endif
    mult_by_gamma(&(pm0[i]),&(dest[i]),GAMMAFIVE);
    scalar_mult_wvec(&(dest[i]), r0a,&(dest[i]));
  }
  free(pm0);
/*
#ifdef DEBUG
if(this_node==0)printf("inner cg iters %d\n",avs_iters);
if(this_node==0)printf("inner cg res %e\n",size_r);
fflush(stdout);
#endif
*/
  for(i=0;i<Norder;i++) free(psim[i]);
  free(psim) ;


#ifdef EIG
      /* add the ''projector term'' back into dest...*/
  if(Nvecs_h0r0 !=0){
    for(j=0;j<Nvecs_h0r0;j++){
      ener=1.0;
      if(eigVal0[j]<0.0) ener=-1.0;

      CMULREAL(cproj[j],ener,ctmp); 
      FORALLSITES(i,s){
#ifdef PREFETCH
	if(i< sites_on_node - FETCH_UP){
	  prefetch_W(&(eigVec0[j][i+FETCH_UP]));
	  prefetch_W(&(dest[i]));
	}
#endif
	c_scalar_mult_add_wvec(&(dest[i]),&(eigVec0[j][i]),
                             &ctmp,&(dest[i]) ); 
      }
    }
    free(cproj);
  }
#endif
    
  /* check the accuracy of the step function and raise or lower the CG accuracy accordingly */
#ifdef MINN
  if(do_minn ==1){
    dest_norm=0.0;
    FORALLSITES(i,s){
      dest_norm += (double)magsq_wvec(&(dest[i]));
    }
    g_doublesum( &dest_norm );
    /*
    node0_printf("step results %e %e %e",source_norm,dest_norm,dest_norm/source_norm-1.0); 
    */  
    test_epsilon = 0.5*(fabs((dest_norm/source_norm) -1.0));
    /*
      node0_printf("test_epsilon %e\n",test_epsilon);
    */
    /* you can always increase resid_inner--up to a point! */
    if(test_epsilon < resid_inner && resid_inner_run < 0.02 ) resid_inner_run *= 1.2;
    /* but it should not shrink too small */
    if((resid_inner_run >= resid_inner) &&(test_epsilon > resid_inner))
      resid_inner_run /= 1.2;
    /* 
    node0_printf("  resid_inner_run %e\n",resid_inner_run);
    */ 
    /*
  node0_printf("congrad_multi_field: resid_inner_run %e for resid_inner %e\n",resid_inner_run,resid_inner);
    */
  }
#endif
}


int congrad_multi_field( /* Return value is number of iterations taken */
    wilson_vector  *src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector **psim, /* elements of the multimass inverter */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src*src) */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    )
{
    int N_iter;
    register int i,j;
    register site *s;

  int iteration; /* counter for iterations */
  double rsq,rsqnew,source_norm,rsqmin,rsqstop;
  complex ctmp;
  double c1,c2,cd;
  double *zeta_i,*zeta_im1,*zeta_ip1;
  double *beta_i,*beta_im1,*alpha;
  int *converged;
  double rsqj;
  Real floatvar, floatvar2, *floatvarj, *floatvark; /* SSE kluge */



        /* SSE KLUGE */
void scalar_mult_add_wvecT(wilson_vector *src1,wilson_vector *src2,Real
        s, wilson_vector *dest);
void mult_mat_wilson_vecT(  su3_matrix *mat, wilson_vector *src,
        wilson_vector *dest );
void mult_adj_mat_wilson_vecT(  su3_matrix *mat, wilson_vector *src,
       wilson_vector *dest );




wilson_vector *mpm, *pm0;  /* malloc space for the vectors involved in gathers */

wilson_vector **pm;  /* malloc space for the vectors not involved in gathers */
wilson_vector *rm;  /* malloc space for the vectors not involved in gathers */


zeta_i=(double *)malloc(Norder*sizeof(double));
zeta_im1=(double *)malloc(Norder*sizeof(double));
zeta_ip1=(double *)malloc(Norder*sizeof(double));
beta_i=(double *)malloc(Norder*sizeof(double));
beta_im1=(double *)malloc(Norder*sizeof(double));
alpha=(double *)malloc(Norder*sizeof(double));

floatvarj=(Real *)malloc(Norder*sizeof(Real));
floatvark=(Real *)malloc(Norder*sizeof(Real));

converged=(int *)malloc(Norder*sizeof(int));
for(i=0;i<Norder;i++) converged[i]=0;


pm = (wilson_vector **)malloc(Norder*sizeof(wilson_vector*));
for(i=1;i<Norder;i++)
  pm[i]=
    (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

rm=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));


mpm=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
pm0=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));


  rsq = source_norm = 0.0;
  rsqmin=RsdCG*RsdCG;
  iteration=0;



    /* This algorithm only works with  psim=0... */
    if(start_flag != 0) {
        if(this_node==0)printf("need to take psi_0=0\n");
	exit(1);
   }


    /* initialize dest[j]=0,r=source, pm[j]=r; */
    FORALLSITES(i,s) {
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(src[i+FETCH_UP]));
             prefetch_W(&(rm[i+FETCH_UP]));
             prefetch_W(&(pm0[i+FETCH_UP]));
	     for(j=1;j<Norder;j++)prefetch_W(&(pm[j][i+FETCH_UP]));
	     for(j=1;j<Norder;j++)prefetch_W(&(psim[j][i+FETCH_UP]));
           }
#endif
      copy_wvec(&(src[i]),&(rm[i]) );
		  copy_wvec(&(rm[i]), &(pm0[i]));
		  clear_wvec(&(psim[0][i]));
		for(j=1;j<Norder;j++){
		  clear_wvec(&(psim[j][i]));
		  copy_wvec(&(rm[i]), &(pm[j][i]));
		}
    }

  FORALLSITES(i,s){
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(src[i+FETCH_UP]));
           }
#endif
    source_norm += (double)magsq_wvec(&(src[i]));
  }

  g_doublesum( &source_norm );
  rsq=source_norm;

  /*  
  if(this_node==0){
    printf("congrad: source_norm = %le\n",source_norm);
    fflush(stdout);
  } 
  */

  rsqstop = rsqmin * source_norm;


	for(j=0;j<Norder;j++){
	  zeta_im1[j]=zeta_i[j]=1.0;
	  alpha[j]=0.0;
	  beta_im1[j]=1.0;
	}

    for( N_iter = 0; N_iter < MaxCG && rsq>rsqstop; 
        N_iter = N_iter + 1) {

        /*   mp = (M(u) + shift[0])*pm */
        hdelta0_field(pm0,mpm,HZERO);
	
    iteration++;
    total_iters++;
    FORALLSITES(i,s){
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(mpm[i+FETCH_UP]));
             prefetch_W(&(pm0[i+FETCH_UP]));
           }
#endif
      scalar_mult_add_wvec(&(mpm[i]) ,
                            &(pm0[i]),shift[0] , &(mpm[i]));
    }

        /* beta_i[0]= - (r,r)/(pm,Mpm)  */
    cd = 0.0;
    FORALLSITES(i,s){
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(mpm[i+FETCH_UP]));
             prefetch_W(&(pm0[i+FETCH_UP]));
           }
#endif
      ctmp =  wvec_dot( &(pm0[i]), &(mpm[i]));
      cd += ctmp.real;
    }

    g_doublesum(&cd);

    beta_i[0]= -rsq/cd;

    /* beta_i(sigma), zeta_ip1(sigma) */

    zeta_ip1[0]=1.0;
	for(j=1;j<Norder;j++){
	  if(converged[j]==0){
	    zeta_ip1[j] = zeta_i[j]*zeta_im1[j]*beta_im1[0];
	    c1=beta_i[0]*alpha[0]*(zeta_im1[j]-zeta_i[j]);
	    c2=zeta_im1[j]*beta_im1[0]*(1.0-(shift[j]-shift[0])*beta_i[0]);
	    zeta_ip1[j] /= c1+c2;

	    beta_i[j]=beta_i[0]*zeta_ip1[j]/zeta_i[j];
	  }
	}



       /* psim[j] = psim[j] - beta[j]*pm[j]  */
	floatvar= -(Real)beta_i[0];
	for(j=1;j<Norder;j++)if(converged[j]==0)floatvarj[j]= -(Real)beta_i[j];


        FORALLSITES(i,s) {
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(psim[0][i+FETCH_UP]));
             prefetch_W(&(pm0[i+FETCH_UP]));
	     for(j=1;j<Norder;j++)if(converged[j]==0)prefetch_W(&(psim[j][i+FETCH_UP]));
	     for(j=1;j<Norder;j++)if(converged[j]==0)prefetch_W(&(pm[j][i+FETCH_UP]));
           }
#endif
	   scalar_mult_add_wvec(  &(psim[0][i]),
                &(pm0[i]), floatvar,  &(psim[0][i])  );
	   for(j=1;j<Norder;j++){
	     if(converged[j]==0){
	     scalar_mult_add_wvec(  &(psim[j][i]),
				    &(pm[j][i]),floatvarj[j],  &(psim[j][i])  );
	     }
	   }
	}

        /* r = r + beta[0]*mp */
	floatvar= (Real)beta_i[0];
        FORALLSITES(i,s) {
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(mpm[i+FETCH_UP]));
             prefetch_W(&(rm[i+FETCH_UP]));
           }
#endif
            scalar_mult_add_wvec( &(rm[i]),
                &(mpm[i]),floatvar, &(rm[i]) );
        }

	/* alpha_ip1[j] */
	rsqnew=0.0;
        FORALLSITES(i,s) {
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(rm[i+FETCH_UP]));
           }
#endif
            rsqnew += (double)magsq_wvec( &(rm[i]) );
        }
        g_doublesum(&rsqnew);
	alpha[0]=rsqnew/rsq;

	/*alpha_ip11--note shifted indices wrt eqn 2.43! */

	for(j=1;j<Norder;j++){
	  if(converged[j]==0){
	    alpha[j]=alpha[0]*zeta_ip1[j]*beta_i[j]/zeta_i[j]/beta_i[0];
	  }
	}

	/* pm[j]=zeta_ip1[j]r +alpha[j]pm[j] */

	floatvar= (Real)zeta_ip1[0];
	floatvar2= (Real)alpha[0];
	for(j=1;j<Norder;j++){
	  floatvarj[j]= (Real)zeta_ip1[j];
	  floatvark[j]= (Real)alpha[j];
	}
	  FORALLSITES(i,s) {
#ifdef PREFETCH
           if(i< sites_on_node - FETCH_UP){
             prefetch_W(&(rm[i+FETCH_UP]));
             prefetch_W(&(mpm[i+FETCH_UP]));
             prefetch_W(&(pm0[i+FETCH_UP]));
	     for(j=1;j<Norder;j++)if(converged[j]==0)prefetch_W(&(pm[j][i+FETCH_UP]));
           }
#endif
	      scalar_mult_wvec( &(rm[i]),floatvar, &(mpm[i]) );
	      scalar_mult_add_wvec( &(mpm[i]),
				    &(pm0[i]),floatvar2, &(pm0[i]) ); 
	      for(j=1;j<Norder;j++){
		if(converged[j]==0){
		  scalar_mult_wvec( &(rm[i]),floatvarj[j], &(mpm[i]) );
		  scalar_mult_add_wvec( &(mpm[i]),
				    &(pm[j][i]),floatvark[j], &(pm[j][i]) );
				    }
	      }
	  }



	/* test for convergence */
        rsq=rsqnew;
	for(j=1;j<Norder;j++){
	  if(converged[j]==0){
	    rsqj=rsq*zeta_ip1[j]*zeta_ip1[j];
	    if(rsqj <= rsqstop) {
	      converged[j]=1;
	      /*
	      node0_printf(" vector %d converged in %d steps %e\n",
			   j,N_iter,rsqj);
	      */
	    }
	  }
	}
	/*
        if(this_node==0 && ((N_iter / 1)*1==N_iter) ){
	printf("iter %d residue %e\n",N_iter,
	    (double)(rsq));
	fflush(stdout);}
	*/
	/* and scroll scalars */
	for(j=0;j<Norder;j++){
	  if(converged[j]==0){
	    beta_im1[j]=beta_i[j];
	    zeta_im1[j]=zeta_i[j];
	    zeta_i[j]=zeta_ip1[j];
	  }
	}

    }
    if( rsq > rsqstop ) {
      if(this_node==0)printf(" inner CONGRAD Not Converged\n rsq = %e\n",rsq);
    }
	*size_r=rsq;

	for(i=1;i<Norder;i++) free(pm[i]);
	free(pm) ;
	free(rm) ;

	free(mpm) ;
	free(pm0) ;

	free(zeta_i);
	free(zeta_ip1);
	free(zeta_im1);
	free(beta_im1);
	free(beta_i);
	free(alpha);
	free(converged);

			  free(floatvarj);
			  free(floatvark);


    return iteration ;
}


void re_setup_inner(double emin, double emax)
{
	free(shift); free(coeff);
	node0_printf("re_setup  : ");
	zolo_min=(Real)emin;
	zolo_max=(Real)emax;
	
        setup_inner();
}

#define NMAX 200
void setup_inner()
{
  double lmin,lmax,c[NMAX],b[NMAX];
  int i,l;
  double x, x1, x2, y1,y2;
  void sncndn(double, double, double*, double*, double*);
  double cel();

  double sn,cn,dn,kappa,kappap,K, Kp,uu;
  double num,den,norm;
  double  prec;
  double fzolo2(int n, double* c, double* b, double x);
  static int called=0;
  
  if (called==0)
  {
  node0_printf("shift multi mass R0 %e\n",R0);
  node0_printf("Zolotarov inner CG order %d\n",Norder);
  node0_printf("inner CG steps %d  resid %e\n",maxcg_inner,resid_inner);
  }
  if(Norder < 0){node0_printf("Can't have negative order\n");exit(1);}
  lmin=(double)zolo_min; lmax=(double)zolo_max;


    /* initialize the parameters */
  kappa=lmin/lmax;
  kappap=sqrt(1.-kappa*kappa);
  K=cel(kappap,1.,1.,1.);
  Kp=cel(kappa,1.,1.,1.);

  for (Norder=4;Norder<NMAX/2-2;Norder+=2)
  {
    /* compute the shifts c_l for approximation
     * on [1;1/kappa]  */
    for(l=1;l<2*Norder+1;l++){
      uu=((double)l)*Kp/((double)(2*Norder));
      sncndn(uu,kappa*kappa,&sn,&cn,&dn);
      num=sn*sn;
      den=cn*cn;
      c[l]=num/den;
    }


    /* rescale shifts for approximation [lmin,lmax]*/
    for(l=1;l<2*Norder+1;l++) c[l]*=lmin*lmin;

    /* compute the residues */ 
    for(l=1;l<=Norder;l++){
      b[l]=1.0;
      for(i=1;i<=Norder-1;i++) b[l] *= (c[2*i]-c[2*l-1]);
      for(i=1;i<=Norder;i++)if(i != l) b[l] /= (c[2*i-1]-c[2*l-1]);
    }
 
   
    /*compute normalization, 
     * the extrema are at x1 and x2 
     * choose norm such that the deviation is
     * equal in both directions */
    
    x=Kp/((double)(2*Norder));
    sncndn(x,kappa*kappa,&sn,&cn,&dn);

    x1=lmin;
    y1=0.0;
    y1=fzolo2(Norder,c,b,x1); 

    x2=lmin/fabs(dn);
    y2=0.0;
    y2=fzolo2(Norder,c,b,x2);

    norm=0.5*fabs(y1+y2);
    prec=fabs(y1-y2)/fabs(y1+y2);
    if (prec<prec_sign) break;
  }
    
  node0_printf("lmin %e lmax %e prec %e Norder %i\n",lmin,lmax,prec,Norder);
    
    /* rescale for ``working'' code 
     *  and print parameters*/
    for(l=1;l<=Norder;l++) b[l]/=norm;

    shift=(Real *)malloc(Norder*sizeof(Real));
    coeff=(Real *)malloc(Norder*sizeof(Real));

    for(i=0;i<Norder;i++){
      coeff[i]=b[i+1];
      shift[i]=c[2*(i+1)-1];
    }

    if (called==0)
    {

      for(i=0;i<Norder;i++)
        node0_printf("shift[ %d ]=%e coeff= %e\n",i,shift[i],coeff[i]);
      node0_printf("scale dslash by %e\n",scalez);

    }
  called=1;

}


double fzolo2(int n, double* c, double* b, double x)
{
  int l;
  double y1;	  
  y1=0.0;
  for(l=1;l<=n;l++) y1 += b[l]/(x*x + c[2*l-1]);
  y1*=x;
  return y1;
}	  

#define CA 1.e-8
#define PIO2 1.57079632679490

double cel(qqc,pp,aa,bb)
double qqc,pp,aa,bb;
{
	double a,b,e,f,g,em,p,q,qc;

	if (qqc == 0.0) {printf("Bad qqc in routine CEL");exit(1);}
	qc=fabs(qqc);
	a=aa;
	b=bb;
	p=pp;
	e=qc;
	em=1.0;
	if (p > 0.0) {
		p=sqrt(p);
		b /= p;
	} else {
		f=qc*qc;
		q=1.0-f;
		g=1.0-p;
		f -= p;
		q *= (b-a*p);
		p=sqrt(f/g);
		a=(a-b)/g;
		b = -q/(g*g*p)+a*p;
	}
	for (;;) {
		f=a;
		a += (b/p);
		g=e/p;
		b += (f*g);
		b += b;
		p=g+p;
		g=em;
		em += qc;
		if (fabs(g-qc) <= g*CA) break;
		qc=sqrt(e);
		qc += qc;
		e=qc*em;
	}
	return PIO2*(b+a*em)/(em*(em+p));
}


#define CA 1.e-8

void sncndn(uu,emmc,sn,cn,dn)
double uu,emmc;
double *sn,*cn,*dn;
{
	double a,b,c,d,emc,u;
	double em[14],en[14];
	int i,ii,l,bo;
	d=0.;

	emc=emmc;
	u=uu;
	if (emc) {
		bo=(emc < 0.0);
		if (bo) {
			d=1.0-emc;
			emc /= -1.0/d;
			u *= (d=sqrt(d));
		}
		a=1.0;
		*dn=1.0;
		for (i=1;i<=13;i++) {
			l=i;
			em[i]=a;
			en[i]=(emc=sqrt(emc));
			c=0.5*(a+emc);
			if (fabs(a-emc) <= CA*a) break;
			emc *= a;
			a=c;
		}
		u *= c;
		*sn=sin(u);
		*cn=cos(u);
		if (*sn) {
			a=(*cn)/(*sn);
			c *= a;
			for (ii=l;ii>=1;ii--) {
				b=em[ii];
				a *= c;
				c *= (*dn);
				*dn=(en[ii]+a)/(b+a);
				a=c/b;
			}
			a=1.0/sqrt(c*c+1.0);
			*sn=(*sn >= 0.0 ? a : -a);
			*cn=c*(*sn);
		}
		if (bo) {
			a=(*dn);
			*dn=(*cn);
			*cn=a;
			*sn /= d;
		}
	} else {
		*cn=1.0/cosh(u);
		*dn=(*cn);
		*sn=tanh(u);
	}
}

#undef CA
