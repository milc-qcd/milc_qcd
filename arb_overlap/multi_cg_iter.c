/******* multi_cg_iter.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Calls hdelta0.c for fermion matrix */

/* version of 18 Oct 01  */

/* This is an implementation of a multi-mass conjugate gradient
algorithm a la B. Jegerlehner. The number of masses, num_masses,
will be set in lattice.h, and so will the array of mass values,
shift0[num_masses].  */

/* The source vector is in "src", and the answer (initial guess is zero)
   in "dest". 
   "psim[jmass]" are working vectors for the conjugate gradient.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.

In this version of the code, all the scalars are real because M=D^\dagger D

This is the actual CG, by itself, to invert H(0)^2 +shift^2.
 All the ``projection'' etc is in a
driver routine for the particular problem.


Conventions:
massless case
H=R0(gamma-5 + epsion(h(-R0)))
H_0^2 = R)^2(2+gamma-5 eps + eps gamma-5)

Massive:
D(m)=(1-m/(2R0))D(0) +m
H(m)^2 = C(H^2+S)
C=(1-m*m/(4R0*R0))) S=m^2/C

Multimass inverter finds psi = (H^2+S)^{-1} chi
Before subtraction
D^{-1} = D(m)^\dagger (H(m))^{-2} = [(1-m/(2R0))H(0) gamma-5 + m] psi/C


Then subtraction for psibar-psi, etc...

Dtilde chi = (D^{-1} chi - chi/(2R0))/(1-m/(2R0))


NOTE: This code uses trial vectors of H(0)^2 in the projection, with
the two chirality states packed in to the upper and lower halves of
eigVecH0 (since eigenvectors of H(0)^2 can be taken to be chiral...)
*/

/* we check all vectors for convergence and quit doing the converged ones.
  */

#include "arb_ov_includes.h"



int congrad_multi_o( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector **psim, /* elements of the multimass inverter */
    int MaxCG,          /* maximum number of iterations per restart */
    Real *RsdCG,        /* desired residual -one for each mass 
                           normalized as sqrt(r*r)/sqrt(src*src) */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    )
{
    int N_iter;
    register int i,j;
    register site *s;

  int iteration; /* counter for iterations */
  double rsq,rsqnew,source_norm,*rsqmin,*rsqstop;
  complex ctmp;
  double c1,c2,cd;
  double *zeta_i,*zeta_im1,*zeta_ip1;
  double *beta_i,*beta_im1,*alpha;
  int *converged;
  double rsqj;

  wilson_vector *pm0,*mpm0;
  wilson_vector **pm;  /* malloc space for the vectors not involved in gathers */
  wilson_vector *rm;  /* malloc space for the vectors not involved in gathers */

  zeta_i=(double *)malloc(num_masses*sizeof(double));
  zeta_im1=(double *)malloc(num_masses*sizeof(double));
  zeta_ip1=(double *)malloc(num_masses*sizeof(double));
  beta_i=(double *)malloc(num_masses*sizeof(double));
  beta_im1=(double *)malloc(num_masses*sizeof(double));
  alpha=(double *)malloc(num_masses*sizeof(double));


  rsqmin=(double *)malloc(num_masses*sizeof(double));
  rsqstop=(double *)malloc(num_masses*sizeof(double));

  converged=(int *)malloc(num_masses*sizeof(int));
  for(i=0;i<num_masses;i++) converged[i]=0;


  pm = (wilson_vector **)malloc(num_masses*sizeof(wilson_vector*));
  for(i=1;i<num_masses;i++)
    pm[i]=
      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  rm=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  pm0=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  mpm0=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  rsq = source_norm = 0.0;
  for(i=0;i<num_masses;i++) rsqmin[i]=RsdCG[i]*RsdCG[i];
  iteration=0;



    /* If num_masses != 1, this algorithm only works with  psim=0...
   we can clearly set up the CG so it reloads one mass. */
    if(start_flag != 0) {
        if(this_node==0)printf("need to take psi_0=0\n");
	exit(1);
   }


    /* initialize dest[j]=0,r=source, pm[j]=r; */
    FORALLSITES(i,s) {
      copy_wvec(((wilson_vector *)F_PT(s,src)),&(rm[i]) );
		  copy_wvec(&(rm[i]), &(pm0[i]));
		  /* clear_wvec(&(psim[0][i])); */
		for(j=1;j<num_masses;j++){
		  /* clear_wvec(&(psim[j][i])); */
		  copy_wvec(&(rm[i]), &(pm[j][i]));
		}
    }

  FORALLSITES(i,s){
    source_norm += (double)magsq_wvec(((wilson_vector *)F_PT(s,src))  );
  }

  g_doublesum( &source_norm );
  rsq=source_norm;

  if(this_node==0){
    printf("congrad: source_norm = %e\n",source_norm);
    fflush(stdout);
  } 

  for(i=0;i<num_masses;i++)rsqstop[i] = rsqmin[i] * source_norm;

	for(j=0;j<num_masses;j++){
	  zeta_im1[j]=zeta_i[j]=1.0;
	  alpha[j]=0.0;
	  beta_im1[j]=1.0;
	}

    for( N_iter = 0; N_iter < MaxCG && rsq>rsqstop[0]; 
        N_iter = N_iter + 1) {


        /*   mp = (M(u) + shift0[0])*pm */
        hdelta0_field(pm0,mpm0,HOVERLAP);
    iteration++;
    total_iters++;
    FORALLSITES(i,s){scalar_mult_add_wvec(&(mpm0[i]) ,
                            &(pm0[i]),shift0[0] , &(mpm0[i]));
    }

        /* beta_i[0]= - (r,r)/(pm0,Mpm0)  */
    cd = 0.0;
    FORALLSITES(i,s){
      ctmp =  wvec_dot( &(pm0[i]), &(mpm0[i]));
      cd += ctmp.real;
    }

    g_doublesum(&cd);

    beta_i[0]= -rsq/cd;

    /* beta_i(sigma), zeta_ip1(sigma) */

    zeta_ip1[0]=1.0;
	for(j=1;j<num_masses;j++){
          if(converged[j]==0){
	    zeta_ip1[j] = zeta_i[j]*zeta_im1[j]*beta_im1[0];
	    c1=beta_i[0]*alpha[0]*(zeta_im1[j]-zeta_i[j]);
	    c2=zeta_im1[j]*beta_im1[0]*(1.0-(shift0[j]-shift0[0])*beta_i[0]);
	    zeta_ip1[j] /= c1+c2;

	    beta_i[j]=beta_i[0]*zeta_ip1[j]/zeta_i[j];
	  }
	}

	/*
	for(j=0;j<num_masses;j++){
	  printf("beta[%d] %e \n",j,beta_i[j]);
	  printf("zeta_ip1[%d] %e \n",j,zeta_ip1[j]);
	} 
*/ 

       /* psim[j] = psim[j] - beta[j]*pm[j]  */
        FORALLSITES(i,s) {
	    scalar_mult_add_wvec(  &(psim[0][i]),
                &(pm0[i]), -(Real)beta_i[0],  &(psim[0][i])  );
	  for(j=1;j<num_masses;j++){
	    if(converged[j]==0){
	      scalar_mult_add_wvec(  &(psim[j][i]),
				     &(pm[j][i]), -(Real)beta_i[j],  &(psim[j][i])  );
	    }
	  }
	}

        /* r = r + beta[0]*mp */
        FORALLSITES(i,s) {
            scalar_mult_add_wvec( &(rm[i]),
                &(mpm0[i]),(Real)beta_i[0], &(rm[i]) );
        }

	/* alpha_ip1[j] */
	rsqnew=0.0;
        FORALLSITES(i,s) {
            rsqnew += (double)magsq_wvec( &(rm[i]) );
        }
        g_doublesum(&rsqnew);
	alpha[0]=rsqnew/rsq;

	/*alpha_ip11--note shifted indices wrt eqn 2.43! */

	for(j=1;j<num_masses;j++){
          if(converged[j]==0){
	    alpha[j]=alpha[0]*zeta_ip1[j]*beta_i[j]/zeta_i[j]/beta_i[0];
	  }
	}
	/*
	for(j=0;j<num_masses;j++){ 
	  printf("alpha %d %e\n",j,alpha[j]);
	}  */

	/* pm[j]=zeta_ip1[j]r +alpha[j]pm[j] */


	  FORALLSITES(i,s) {
	      scalar_mult_wvec( &(rm[i]),(Real)zeta_ip1[0], &(mpm0[i]) );
	      scalar_mult_add_wvec( &(mpm0[i]),
				    &(pm0[i]),(Real)alpha[0], &(pm0[i]) ); 
	      for(j=1;j<num_masses;j++){
                if(converged[j]==0){
		  scalar_mult_wvec( &(rm[i]),(Real)zeta_ip1[j], &(mpm0[i]) );
		  scalar_mult_add_wvec( &(mpm0[i]),
				    &(pm[j][i]),(Real)alpha[j], &(pm[j][i]) );
				    }
	      }
	  }



	/* test higher masses for convergence */
        rsq=rsqnew;
        for(j=1;j<num_masses;j++){
          if(converged[j]==0){
            rsqj=rsq*zeta_ip1[j]*zeta_ip1[j];
            if(rsqj <= rsqstop[j]) {
              converged[j]=1;
              
              node0_printf(" vector %d converged in %d steps %e\n",
                           j,N_iter,rsqj);
                           
            }
          }
        }

	 
        if(this_node==0 && ((N_iter / 5)*5==N_iter) ){
	printf("iter %d residue %e",N_iter,
	    (double)(rsq));
#ifdef MINN
	printf("  %e",resid_inner_run);
#endif
	printf("\n");
	fflush(stdout);}
       
	/* and scroll scalars */
	for(j=0;j<num_masses;j++){
          if(converged[j]==0){
	    beta_im1[j]=beta_i[j];
	    zeta_im1[j]=zeta_i[j];
	    zeta_i[j]=zeta_ip1[j];
	  }
	}

    }
    if( rsq > rsqstop[0] ) {
        if(this_node==0)printf(" outer CONGRAD Not Converged\n");
    }

	*size_r=rsq;

	for(i=1;i<num_masses;i++) free(pm[i]);
	free(pm) ;
	free(rm) ;
	free(pm0) ;
	free(mpm0) ;

	free(zeta_i);
	free(zeta_ip1);
	free(zeta_im1);
	free(beta_im1);
	free(beta_i);
	free(alpha);
        free(converged);
        free(rsqstop);
        free(rsqmin);


    return iteration ;
}
