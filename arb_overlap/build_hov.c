/************* build_hov.c ******************************/
/* MIMD version 7 */

/*
 This routine computes eigenvectors and eigenvalues of H(0)^2 in one
chirality sector (set by the parameter trial_chirality) and returns them. It computes
the eigenvalues and eigenvectors of H(0) by diagonalizing the nonchiral modes
two-by-two: i.e, take each nonchiral mode psi,
 compute |chi> = H(0)|psi> - |psi><psi<H(0)|psi>, normalize |chi>, and
diagonalize H(0) in the 2x2 basis. ``Chiral'' modes are assumend not to
 be degenerate, so the state is just the eigenstate of H(0), the
eigenvalue is just the square root of the eigenvalue of H(0)^2 (both had
better be close to zero).

This routine returns ``packed'' eigenmodes of H(0)^2 in a Wilson vector;
the upper 2 components are the chirality +1 eigenmode of H(0)^2, the
lower 2 components are the chirality -1 eigenmode.



*/



#include "arb_ov_includes.h"


int build_hov(int *trial_chirality, int* jcount)
{


  wilson_vector **eigVec1;
  double *eigVal1;
  register int i,j;
  register site *s;

  int kk,jj,l,j0max;
  int total_R_iters ;
  Real norm,ener;
  Real re,im,re5,im5;

  Real test5[2];

  complex cc;
  Matrix Array,V;

  /* for guessing eigenvectors */
  half_wilson_vector hwvec;
  complex ctmp;
  int source_chirality = 0;
  int ndel;

  eigenval_tol=eigenval_tol_low;
  error_decr=error_decr_high;
  ndel=ndelta0;


#ifdef OPP
  node0_printf("flipping selected chirality sector\n");
  *trial_chirality *= -1;
#endif
  node0_printf("\n-----------------eigenvalues of H_ov-----------------------\n");
  node0_printf("pass 2, want %d eigenmodes in chirality sector %d\n",Nvecs_hov,
	       *trial_chirality);
  build_params(-R0);

  chirality_flag = *trial_chirality;


  /* reload only the desired chirality of the eigenmodes */
  for(j=0;j<Nvecs_hov;j++){
    if( *trial_chirality == 1){
      FORALLSITES(i,s){
        w_to_hw(&eigVec[j][i],0,&hwvec);
        hw_to_w(&hwvec,0,&eigVec[j][i]);
      }
    }
    if( *trial_chirality == -1){
      FORALLSITES(i,s){
        w_to_hw(&eigVec[j][i],2,&hwvec);
        hw_to_w(&hwvec,2,&eigVec[j][i]);
      }
    }
  }


  /*  we are ready to begin finding eigenvectors of
      the overlap hamiltonian. We reload the pole-crossing (negative) mass */	

  build_params(-R0);
  make_clov1();
  kind_of_h0=HOVERLAP;
  total_R_iters=Kalkreuter(eigVec, eigVal, eigenval_tol, 
			   error_decr, Nvecs_hov, MaxIter,
			   Restart, 
			   Kiters, EVENANDODD) ;
  if(this_node==0)printf("total Rayleigh iters = %d\n",total_R_iters);


  /* now we presumably have a set of eigenvectors of H_0^2 */
  node0_printf("calling last analysis package\n");
  node0_printf("\n");
  node0_printf("START_EIGVALS of H_ov and other properties\n");
  fflush(stdout) ;


  for(j=0;j<Nvecs_hov;j++){

  FORALLSITES(i,s){s->chi=eigVec[j][i];}
  /* expectation value of gamma-5 */
  re5=im5=0.0;
  FORALLSITES(i,s){
    mult_by_gamma(&(s->chi),&(s->r),GAMMAFIVE);
    cc = wvec_dot( &(s->chi), &(s->r) );
    re5 += cc.real ;
  }
  g_floatsum(&re5);
  if(j==0){if(re5>0.0) source_chirality=1;
  if(re5<0.0) source_chirality= -1;
  node0_printf("source chirality %d\n",
	       source_chirality);}
  if(this_node==0) printf("F3HO2V %d %e %e\n",
			  j,(double)eigVal[j],(double)re5);
  }

  /* now we want to diagonalize our eigenmodes. The
     chiral modes assumed to be unmixed among themselves and the nonchiral modes
     are mixed two by two (and we have one chirality) */

  /* recall, eigValcut = variational bound on lowest paired eigenmode. 
     If eigVal[j]-eigValcut > delta, we certainly have a paired mode, for some delta.
     If eigVal[j]-eigValcut < -delta
  */
  j0max=0;
  for(j=0;j<Nvecs_hov;j++){
    if(eigVal[j]-(eigValcut+eigenval_tol) < 0.0 ) j0max++;
  }
  node0_printf("%d  modes below threshold %e\n",j0max,eigValcut);
 /* there are at most j0_max true zero modes */

  *jcount=0;
  for(j=0;j<Nvecs_hov;j++){
      node0_printf("Eigenvalue_hov(%i)= %e\n",j,eigVal[j]);

    if(eigVal[j] > 1.e-6) /* doubly degenerate, non-topological modes, this is our cut which assumes
			 we have a zero mode without further processing */
      {

	eigVal1 = (double *)malloc(2*sizeof(double));
	eigVec1 = (wilson_vector **)malloc(2*sizeof(wilson_vector*));
	for(i=0;i<2;i++)
	  eigVec1[i]=
	    (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));


	FORALLSITES(i,s){ eigVec1[0][i]  =eigVec[j][i]; s->chi=eigVec[j][i];}
	hoverlap(F_OFFSET(chi),F_OFFSET(r));
	/* orthogonalize wrt eigVec1[0] */
	ctmp =cmplx(0.0,0.0);
	FORALLSITES(i,s){
	  cc =  wvec_dot( &(s->chi),&(s->r));
	  CSUM(ctmp,cc);
	}
	g_complexsum(&ctmp);
	CMULREAL(ctmp,-1.0,ctmp);
	FORALLSITES(i,s)
	  c_scalar_mult_add_wvec(&(s->r), &(s->chi),
                             &ctmp, &(s->r) );


	/* and renormalize */
	norm=0.0;
	FORALLSITES(i,s){
	  cc = wvec_dot( &(s->r), &(s->r) );
	  norm+= cc.real;
	}
	g_floatsum(&norm);
	norm=1.0/sqrt(norm);
	FORALLSITES(i,s)scalar_mult_wvec(&(s->r),norm,&(s->r));
	FORALLSITES(i,s){eigVec1[1][i]=s->r;}



	/* Now re-pack the states eigVec1 into eigVec, in the ``empty''
	   chirality slot */
	FORALLSITES(i,s){
	  if(source_chirality== 1){
	    for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	      eigVec[j][i].d[kk].c[l]=eigVec1[1][i].d[kk].c[l];
	  }
	  if(source_chirality== -1){
	    for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	      eigVec[j][i].d[kk].c[l]=eigVec1[1][i].d[kk].c[l];
	  }
	}




	/** Allocate the array **/
	Array = AllocateMatrix(2) ;

	/** Allocate the Eigenvector matrix **/
	V = AllocateMatrix(2) ;


	for(jj=0;jj<2;jj++){
	  FORALLSITES(i,s){s->r=eigVec1[jj][i];}
	  node0_printf("ARRAY ");   

	  for(kk=0;kk<=jj;kk++){
	    FORALLSITES(i,s){s->chi=eigVec1[kk][i];}
	    hoverlap(F_OFFSET(chi),F_OFFSET(psi));

	    re=im=0.0;
	    FORALLSITES(i,s){

	      cc = wvec_dot( &(s->r), &(s->psi) );
	      re += cc.real ;
	      im += cc.imag ;
	    }
	    g_floatsum(&re);
	    g_floatsum(&im);
	    Array.M[jj][kk].real = (double)re ; Array.M[jj][kk].imag = (double)im ;
	    Array.M[kk][jj].real = (double)re ; Array.M[kk][jj].imag = -(double)im ;
	    node0_printf("%e %e       ",Array.M[jj][kk].real,Array.M[jj][kk].imag);   
	  }
	  node0_printf("\n");
	}
	Jacobi(&Array, &V, JACOBI_TOL) ;
	sort_eigenvectors(&Array,&V);
	RotateBasis(eigVec1,&V) ;
	/* now we do tests on eigenvectors */
	for(jj=0;jj<2;jj++){
	  eigVal1[jj] = Array.M[jj][jj].real ;
	  FORALLSITES(i,s){s->chi=eigVec1[jj][i];}
                  
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

	  if(this_node==0) printf("F3OGH02X2 %d %e  %e %e  %e\n",
				jj,(double)eigVal1[jj],(double)re5,(double)im5,
				  (double)(eigVal1[jj]/2.0/R0));
	test5[jj]= fabs(re5 - eigVal1[jj]/2.0/R0);
	}
      /* final check: if eigenvalue is large enough, don't split, no matter what. Only
	 j0max(at most) zero modes can exist, so once you've found them all, quit */
	if(eigVal[j] <eigValcut)if(*jcount < j0max )if( (fabs(eigVal1[0]+eigVal1[1]) > 0.0008 )   )
	{
	  (*jcount)++;
	  node0_printf(" NOT keeping the wrong chirality eigenvector!\n");
          FORALLSITES(i,s){
	    if(source_chirality== 1){
	      for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
		eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);;
	    }
	    if(source_chirality== -1){
	      for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
		eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);         ;
	    }
	  }
	}
	/** Deallocate the arrays **/
	deAllocate(&V) ;
	deAllocate(&Array) ;

	for(i=0;i<2;i++)
	  free(eigVec1[i]) ;
	free(eigVec1) ;
	free(eigVal1) ;
      }

    else{ 
      /* near zero mode--can be taken from H_ov^2 ...
	 Note: the state is assumed to be unmixed, the energy is the square root*/
      FORALLSITES(i,s){ s->chi=eigVec[j][i];}
      ener=sqrt(fabs(eigVal[j]));


                  /* expectation value of gamma-5 */
      re5=im5=0.0;
      FORALLSITES(i,s){
	mult_by_gamma(&(s->chi),&(s->r),GAMMAFIVE);
	cc = wvec_dot( &(s->chi), &(s->r) );

	re5 += cc.real ;
      }
      g_floatsum(&re5);

      if(this_node==0) printf("F3OGH01X1 %d %e  %e   %e\n",
				j,(double)ener,(double)re5,
				(double)(ener/2.0/R0));

      (*jcount)++;

	/* load zeros into wrong chirality (should be redundant) */
      FORALLSITES(i,s){
	if( *trial_chirality==1){
	  for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	    eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);
	}
	if( *trial_chirality== -1){
	  for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	    eigVec[j][i].d[kk].c[l]=cmplx(0.0,0.0);
	}
      }

    }

  }
 
  node0_printf("ZEROMODES %d CHIRALITY %d\n",*jcount, *trial_chirality);
  node0_printf("END_EIGVALS of H_ov and other properties\n");
  node0_printf("\n");
  node0_printf("mxv operations in build_hov: %i\n",(int)(ndelta0-ndel));
  (*jcount) *= *trial_chirality;
  return 0;
}
