/************* build_lowest_chi.c ******************************/
/* MIMD version 7 */

/*
This routine takes a set of trial vectors (assumed to be ``packed'' into
chirality pairs), takes one of each chirality, and diagonalizes the overlap
operator H(0)^2  in the 2-dimensional vector space. It repacks the two resulting
chiral vectors as the eigVal[0] of the array and returns the lowest chiraliy.


*/



#include "arb_ov_includes.h"



int build_lowest_chi(int *trial_chirality)
{


  wilson_vector **eigVec1;
  double *eigVal1;


  register int i,j;
  register site *s;

  int kk,l;

  int total_R_iters ;

  half_wilson_vector hwvec;





  /* for guessing eigenvectors */
  wilson_vector wtmp;
  complex ctmp,cn;
  Real chirality,cd;

  int Nvecs_hove;

  int lowest;
/* for splitting eigenmodes */
  int j0max;



  node0_printf("\n----------------find sector of lowest chirality----------------------\n");

  eigVal1 = (double *)malloc(Nvecs_hov*sizeof(double));
  eigVec1 = (wilson_vector **)malloc(Nvecs_hov*sizeof(wilson_vector*));

  for(i=0;i<Nvecs_hov;i++)
    eigVec1[i]=
      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  for(j=0;j<Nvecs_hov;j++){
    FORALLSITES(i,s){eigVec1[j][i]=eigVec[j][i];}
    eigVal1[j]=eigVal[j];
  }


  Nvecs_hove=2;

  node0_printf("loading the lowest %d eigenvecs  of H_0^2 in the chirality sector of the minimum \n",Nvecs_hove);
  /* eigVecs are stored ''packed'' so must only check the norm of the packed vector */

  /*source_chirality = 1 */
  j0max=0;

  for(j=0;j<Nvecs_hov;j++)if(j0max<1){

    cd=0.0;
    FORALLSITES(i,s){
      w_to_hw(&eigVec1[j][i],0,&hwvec);
      cd += magsq_hwvec(&hwvec);
    }
    g_floatsum(&cd);
    if(cd != 0.0){
      j0max=1;
      node0_printf("taking positive chirality part of vector %d\n",j);
      FORALLSITES(i,s){
	w_to_hw(&eigVec1[j][i],0,&hwvec);
	hw_to_w(&hwvec,0,&eigVec[0][i]);
      }
    }
  }
  if(j0max == 0){
    node0_printf("could not find a positive chirality vector in build_lowest_chi. Exiting\n");
    exit(1);
  }


  /*source_chirality = -1 */
  j0max=0;

  for(j=0;j<Nvecs_hov;j++)if(j0max<1){

    cd=0.0;
    FORALLSITES(i,s){
      w_to_hw(&eigVec1[j][i],2,&hwvec);
      cd += magsq_hwvec(&hwvec);
    }
    g_floatsum(&cd);
    if(cd != 0.0){
      j0max=1;
      node0_printf("taking negative chirality part of vector %d\n",j);
      FORALLSITES(i,s){
	w_to_hw(&eigVec1[j][i],2,&hwvec);
	hw_to_w(&hwvec,2,&eigVec[1][i]);
      }
    }
  }
  if(j0max == 0){
    node0_printf("could not find a negative chirality vector in build_lowest_chi. Exiting\n");
    exit(1);
  }




		/* finally we are ready to begin finding eigenvectors of
the overlap hamiltonian. We reload the pole-crossing (negative) mass */	

  build_params(-R0);
  /*
  make_clov1();
  refresh_links();
  */
  kind_of_h0=HOVERLAP;

  /*find lowest mode in chirality sector plus*/
  chirality_flag=1;
  total_R_iters=Kalkreuter(eigVec, eigVal, eigenval_tol, 
                                         error_decr, 1, MaxIter,
					 Restart, 
                                         Kiters, EVENANDODD) ;

  /*find lowest mode in chirality sector minus*/
  chirality_flag=-1;
  total_R_iters=Kalkreuter(&eigVec[1], &eigVal[1], eigenval_tol, 
                                         error_decr, 1, MaxIter,
					 Restart, 
                                         Kiters, EVENANDODD) ;
  if(this_node==0)printf("total Rayleigh iters = %d\n",total_R_iters);

  node0_printf("\n");
  node0_printf("START_EIGVALS for finding chirality of lowest eigenvalue\n");
  for(i=0;i<2;i++){
    node0_printf("eigval %d %e\n",i,eigVal[i]);
  }
  /* pick the lowest one */
  if(eigVal[0]<=eigVal[1]){
    lowest=0;
  }
  else{
    lowest=1;
  }

  
  cn=cmplx(0.0,0.0);
  cd=0.0;
  FORALLSITES(i,s){
    mult_by_gamma(&(eigVec[lowest][i]),&wtmp,GAMMAFIVE);
    cd += magsq_wvec(&(eigVec[lowest][i]));
    ctmp =  wvec_dot(&(eigVec[lowest][i]),&wtmp);
    CSUM(cn,ctmp);
  }
  g_complexsum(&cn);
  g_floatsum(&cd);
  chirality=cn.real/cd;
  if(chirality >0.0) *trial_chirality=1;
  if(chirality <0.0) *trial_chirality= -1;


  node0_printf("2 state min energy H0^2 chirality %e %e  %d\n",
	       (double)eigVal[lowest],(double)chirality,*trial_chirality);

	  /* the upper eigenmode we have found is (within fluctuations) a
variational bound on the energy of all doubled modes. Keep it to help us split
the modes after the diagonalization in the lower subspace */
  eigValcut=eigVal[1-lowest];
  node0_printf("eigValcut %e\n",eigValcut);


  /* repack states for minimization --lowest two in eigVec[0]: Note that 
    eigVec[0] is the lowest chirality +1 state and eigVec[1] is the lowest chi=-1 state */
  FORALLSITES(i,s){
    for(kk=0;kk<2;kk++)for(l=0;l<3;l++)
      eigVec1[0][i].d[kk].c[l]=eigVec[0][i].d[kk].c[l];
    for(kk=2;kk<4;kk++)for(l=0;l<3;l++)
      eigVec1[0][i].d[kk].c[l]=eigVec[1][i].d[kk].c[l];
  }
  


  /* copy the original vectors back */
  for(j=0;j<Nvecs_hov;j++){
    FORALLSITES(i,s){eigVec[j][i]=eigVec1[j][i];}
    eigVal[j]=eigVal1[j];
  }




  for(i=0;i<Nvecs_hov;i++)free(eigVec1[i]) ;
  free(eigVec1) ;
  free(eigVal1) ;


  node0_printf("TRIAL CHIRALITY %d\n", *trial_chirality);
  node0_printf("END_EIGVALS for finding chirality of lowest eigenvalue\n");
  return 0;
}
