/************************** integrals.c ************************/
/* MIMD version 6 */
/* test integrals over gauge group */
/* DT 10/1/97
  complile:
  cc -g -o integrals -I../include \
    integrals.c ../libraries/su3.a ../libraries/complex.a -lm
*/
#include "su3_dense_includes.h"

Real make_change(su3_matrix *mat, Real scale);
double drand48(); void srand48( int );

main(int argc, char**argv){
    su3_matrix A,B,C;
    complex cc1,cc2,cc3;
    Real beta,sum;
    int i,j,k,nsamples;
    Real Delta,theta,z1,z2,z3,phi1,phi2,phi3;


    srand48(12345);
    beta=1.0;
    nsamples = 1000;

    /* Make a "source matrix" in C */
    make_change( &C, 1.0 );
    for(i=0;i<50;i++){
	make_change( &A, 1.5 );
	add_su3_matrix( &C, &A, &C );
    }
    scalar_mult_su3_matrix( &C, 1/50.0, &C );
    dumpmat( &C );

    /* integrate exp(beta*Tr( A_adjoint*C + C_adjoint*A ) )  numerically */
    sum = 0.0;
    make_change( &A, 1.5 );
    for(i=0;i<nsamples;i++){
	random_su3_matrix( &A );
	/**mult_su3_na( &A, &A, &B );
	printf("CHECK\n");
	dumpmat( &A ); dumpmat( &B );**/
	sum += exp( beta*2.0*realtrace_su3( &A, &C ) );
    }
    sum /= nsamples;
    printf("beta = %e, nsamples = %d, average = %e\n",beta,nsamples,sum);

    /* integrate exp(beta*Tr( A*C + C_adjoint*A_adjoint ) )  analytically */
    /* find eigenvalues of C_adjoint*C */
    mult_su3_na( &C, &C, &A );

}

/* Make a reasonably random SU3 matrix - take the current matrix and
make a bunch of random changes */
random_su3_matrix( su3_matrix *M ){
    int i;
    su3_matrix M2,M3;
    for(i=0;i<10;i++){
	make_change( &M2, 1.5 );
	mult_su3_nn( M, &M2, &M3 );
	*M = M3;
    }
    reunit_su3( M );
}


/* Make change matrix for Metropolis step.  Choose random generator
  and rotation angle */
Real make_change(su3_matrix *mat, Real scale){
  register int ia,ib,gen;
  register Real theta,c,s;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  /* random (symmetric) change */
  /* angles range from +- theta/4 to theta/2 */
  do{ theta = scale * (drand48() - 0.5); }
	while( fabs(theta) < 0.25*scale );
  c = cos(theta); s = sin(theta);
  /* random generator, range 1 through 8 */
  gen = ((int)(8.0*drand48()))+1;

  switch(gen){
    case 1:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][1].imag = mat->e[1][0].imag = s;
      break;
    case 2:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][1].real = s;
      mat->e[1][0].real = -s;
      break;
    case 3:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][0].imag = s;
      mat->e[1][1].imag = -s;
      break;
    case 4:
      mat->e[0][0].real = mat->e[2][2].real = c;
      mat->e[0][2].imag = mat->e[2][0].imag = s;
      break;
    case 5:
      mat->e[0][0].real = mat->e[2][2].real = c;
      mat->e[0][2].real = s;
      mat->e[2][0].real = -s;
      break;
    case 6:
      mat->e[1][1].real = mat->e[2][2].real = c;
      mat->e[1][2].imag = mat->e[2][1].imag = s;
      break;
    case 7:
      mat->e[1][1].real = mat->e[2][2].real = c;
      mat->e[1][2].real = s;
      mat->e[2][1].real = -s;
      break;
    case 8:
      mat->e[0][0].real = mat->e[1][1].real = cos(theta/sqrt(3.));
      mat->e[0][0].imag = mat->e[1][1].imag = sin(theta/sqrt(3.));
      mat->e[2][2].real = cos(-2.*theta/sqrt(3.));
      mat->e[2][2].imag = sin(-2.*theta/sqrt(3.));
      break;
  }
/*printf("make_change: theta = %e\n",theta);*/
  return(theta);
}
