/************* build_h0.c ******************************/
/* MIMD version 7 */

#include "arb_ov_includes.h"



/* Comment these out if you want to suppress detailed timing */
#define IOTIME
#define PRTIME

int build_h0()
{

  register int i,j;
  register site *s;

  int kk,l;

  int total_R_iters ;

  wilson_vector **eigVec1;
  double *eigVal1;


  /* for guessing eigenvectors */
  wilson_vector wtmp;
  complex ctmp,cn;
  Real chirality,cd;
  int jplus,jminus,source_chirality = 0;
  int *ov_chirality;

  void normalize(wilson_vector *vec) ;
  void project_out(wilson_vector *vec, wilson_vector **vector, int Num);


  eigenval_tol=eigenval_tol_low;
  error_decr=error_decr_low;

	  /* begin eigenmodes and eigenvals of h(0) */

  node0_printf("\n-------------------eigenvalues of H(0)---------------------\n");

      eigVal1 = (double *)malloc(Nvecs_h0*sizeof(double));
      eigVec1 = (wilson_vector **)malloc(Nvecs_h0*sizeof(wilson_vector*));
      for(i=0;i<Nvecs_h0;i++)
        eigVec1[i]=
          (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  /*  we find some trial eigenvalues of te massless approximate GW action */
  build_params(0.0);
  make_clov1();
  kind_of_h0=HZERO;

 

  for(j=0;j<Nvecs_h0;j++){
    if(j< Nvecs_h0/2){ source_chirality=1;}
    else{source_chirality= -1;}
/*
    if(this_node==0)printf("source chirality %d\n",source_chirality);
*/
    grsource_dirac(EVENANDODD);
    FORALLSITES(i,s){
      copy_wvec(&(s->g_rand),&(eigVec1[j][i]));
      if(source_chirality==1){
	for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	  eigVec1[j][i].d[kk].c[l]=cmplx(0.0,0.0);
      }
      if(source_chirality== -1){
	for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	  eigVec1[j][i].d[kk].c[l]=cmplx(0.0,0.0);
      }		  
    }
    eigVal1[j]=1.0e+16;

  }


  /* next we find some trial eigenvalues of our
		   approximate GW action */
		/* we'll live at m=0 for now */
  total_R_iters=Kalkreuter(eigVec1, eigVal1, eigenval_tol, 
                                         error_decr, Nvecs_h0, MaxIter, Restart, 
                                         Kiters, EVENANDODD) ;
  if(this_node==0)printf("total Rayleigh iters = %d\n",total_R_iters);

  /** Print out various properties of h(0) **/
  /*
  for(j=0;j<Nvecs_h0;j++){

    FORALLSITES(i,s){s->chi=eigVec1[j][i];}

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
			    j,(double)re,(double)im,
			    (double)re5,(double)im5);
  }
  */


  node0_printf("\n");
  node0_printf("START_EIGVALS of h(0)\n");
  for(j=0;j<Nvecs_h0;j++){
    node0_printf("Eigval %d %e\n", j, (double)eigVal1[j]);
  }
  node0_printf("END_EIGVALS of h(0)\n");
  node0_printf("\n");
  fflush(stdout);


      /* our overlap iteration for H_0^2 wants chiral eigenvectors. Measure
the chirality of our vectors and chop the appropriate part off... */
  ov_chirality=malloc(Nvecs_h0*sizeof(int));
  for(j=0;j<Nvecs_h0;j++){
    cn=cmplx(0.0,0.0);
    cd=0.0;
    FORALLSITES(i,s){
      mult_by_gamma(&(eigVec1[j][i]),&wtmp,GAMMAFIVE);
      cd += magsq_wvec(&(eigVec1[j][i]));
      ctmp =  wvec_dot(&(eigVec1[j][i]),&wtmp);
      CSUM(cn,ctmp);
    }
    g_complexsum(&cn);
    g_floatsum(&cd);
    chirality=cn.real/cd;

    if(chirality >0.0) {source_chirality=1;ov_chirality[j]=1;}
    if(chirality <0.0) {source_chirality= -1;ov_chirality[j]= -1;}
/*
    node0_printf("chirality %d %e\n",source_chirality,(double)chirality);
*/
    if(source_chirality==1){
      FORALLSITES(i,s){
	for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	  eigVec1[j][i].d[kk].c[l]=cmplx(0.0,0.0);
      }
    }

    if(source_chirality== -1){
      FORALLSITES(i,s){
	for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	  eigVec1[j][i].d[kk].c[l]=cmplx(0.0,0.0);
      }
    }
  }



  /* now pack vectors into chiral pairs, to treat ''like'' H(0) eigenmodes */
  jplus=jminus=0;
  for(j=0;j<Nvecs_h0;j++){
    if(ov_chirality[j] == 1 && jplus<Nvecs_hov){
      FORALLSITES(i,s){
	for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	  eigVec[jplus][i].d[kk].c[l]=eigVec1[j][i].d[kk].c[l];
      }
      jplus++;
    }
    if(ov_chirality[j] == -1 && jminus<Nvecs_hov){
      FORALLSITES(i,s){
	for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	  eigVec[jminus][i].d[kk].c[l]=eigVec1[j][i].d[kk].c[l];
      }
      jminus++;
    }
  }
  /* now top off with random vectors (and warn the user) */
  if(jplus<Nvecs_hov){
    node0_printf("Need %d random  positive chirality vectors\n",Nvecs_hov-jplus);
    for(j=jplus;j<Nvecs_hov;j++){
      grsource_dirac(EVENANDODD);
      FORALLSITES(i,s){
	for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
	  eigVec[j][i].d[kk].c[l]=s->g_rand.d[kk].c[l];
      }
    }
  }
  if(jminus<Nvecs_hov){
    node0_printf("Need %d random negative chirality vectors\n",Nvecs_hov-jminus);
    for(j=jminus;j<Nvecs_hov;j++){
      grsource_dirac(EVENANDODD);
      FORALLSITES(i,s){
	for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
	  eigVec[j][i].d[kk].c[l]=s->g_rand.d[kk].c[l];
      }
    }
  }


  for(j=0;j<Nvecs_hov;j++){
    eigVal[j]=1.e16;
  }


  for(i=0;i<Nvecs_h0;i++)
    free(eigVec1[i]) ;
  free(eigVec1) ;
  free(eigVal1) ;
  free(ov_chirality);
  build_params(-R0);

  return total_R_iters;
}
