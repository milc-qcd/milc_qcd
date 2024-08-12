/********************** check_reunit_deriv.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the fermion force calculation */

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#else
#error "Requires QIO compilation"
#endif

#ifdef HAVE_GRID
#include "../include/mGrid/mGrid.h"  /* For GRID_info_t */
extern GRID_4Dgrid *grid_full;
#endif

static void ranmat(su3_matrix *mat)
{
  int i,j,k,dir;
  site *s;

  FORALLSITES(i,s) {
    for(dir=XUP;dir<=TUP;dir++) {
      for(j=0; j<3; j++)  {
	for(k=0; k<3; k++)  {
	  Real x = 0.7*gaussian_rand_no(&(s->site_prn));
	  Real y = 0.7*gaussian_rand_no(&(s->site_prn));
	  if (j != k)  {
	    mat[4*i+dir].e[j][k] = cmplx(x,y);
	  }
	  else  {
	    mat[4*i+dir].e[j][k] = cmplx(1.0+x,y);
	  }
	}
      }
    }
  }
}

static void check_answer(su3_matrix *dW, su3_matrix *ansdW, Real tol){
  /* Unpack the answer and compare if possible */
  Real diff, reldiff;
  Real maxdiff = 0;
  Real maxnorm = 0;
  Real norm = 0;
  int i, dir;
  su3_matrix diffmat;
  FORALLFIELDSITES(i){
    FORALLUPDIR(dir){
      su3_matrix *tmat = ansdW + 4*i + dir;
      sub_su3_matrix( dW + 4*i + dir, tmat, &diffmat);
      diff = sqrt(realtrace_su3( &diffmat, &diffmat ));
      norm = sqrt(realtrace_su3( tmat, tmat));
      if(diff > tol * norm){
	printf("Intolerable relative difference %e (norm %e) node %d site %d dir %d\n",
	       diff/norm,norm,this_node,i,dir);
	dumpmat(dW + 4*i + dir);
	dumpmat(tmat);
      }
      if(maxdiff < diff)maxdiff = diff;
      if(maxnorm < norm)maxnorm = norm;
      /* In any case, copy the new result to the answer matrix */
      ansdW[4*i + dir] = *tmat;
    }
  }
  g_floatmax(&maxdiff);
  g_floatmax(&maxnorm);
  if(maxnorm > 0){
    reldiff = maxdiff/maxnorm;
    node0_printf("Relative difference %e\n",reldiff);
  }
  else
    node0_printf("Absolute difference %e but norm is 0???\n",maxdiff);
}

/* Pretend that Q contains the old accumulated force, and V contains the link matrices
   to be unitarized. Then W contains the unitarized link matrices.
   Calculate the derivatives of W with respect to respect to V:
   dW/dV and d(W^+)/dV (at fixed V^+ !), where W=V(V^+V)^-1/2 
   Return dW = Tr[Q dW/dV) + Tr(Q^+ dW+/dV) */

static void milc_u3_reunit_deriv( info_t info, su3_matrix *V, su3_matrix *dW, su3_matrix *Q ){
  su3_tensor4 dwdv, dwdagdv;
  int i;
  complex ftmp;
  su3_matrix tmat;
  FORALLFIELDSITES(i) {
    for(int dir=XUP;dir<=TUP;dir++) {
      /* Calculate derivative for a single link */
      u3_unit_der_analytic( &info, V+4*i+dir, &dwdv, &dwdagdv );
      clear_su3mat( dW + 4*i+dir );
      su3_adjoint( Q+4*i+dir, &tmat );
      /* Take trace of Q times derivatives */
      for( int m=0; m<3; m++) {
        for( int n=0; n<3; n++) {
	  for(int k=0; k<3; k++) {
	    for(int l=0; l<3; l++) {
	      CMUL( dwdv.t4[k][m][n][l], Q[4*i+dir].e[l][k], ftmp );
	      CSUM( dW[4*i+dir].e[n][m], ftmp );
	      /* CAREFUL with the adjoint part here! */
	      CMUL( dwdagdv.t4[k][m][n][l], tmat.e[l][k], ftmp );
	      CSUM( dW[4*i+dir].e[n][m], ftmp );
	      
	    }
	  }
	}
      }
    }
  }
}

#ifdef HAVE_GRID
static void
extract_grid_info( info_t *info, GRID_info_t *grid_info ){
  info->final_sec = grid_info->final_sec;
  info->final_flop = grid_info->final_flop;
  info->status = grid_info->status;
  info->count1 = grid_info->count1;
  info->count2 = grid_info->count2;
}
#endif

void check_reunitarization_derivative( char *ansfilein, int ansflagin,
				       char *ansfileout, int ansflagout )
{
  int i, dir;
  char *filexml;
  char recxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test fermion unitarization derivative</title>";
#if (MILC_PRECISION == 1)
  Real tol = 1e-3;
#else
  Real tol = 1e-8;
#endif
  int ff_prec = MILC_PRECISION;  /* Just use prevailing precision for now */

  /* Set up dummy inputs to check reunitarization */

  su3_matrix *V    = create_G();
  su3_matrix *Q    = create_G();
  su3_matrix *dW   = create_G();
  su3_matrix *dWans = create_G();
    
  if(V == NULL || Q == NULL || dW == NULL || dWans == NULL){
    node0_printf("No room for test matrices\n");
    terminate(1);
  }

  node0_printf("Creating random link matrices for U(3) unitarization\n"); fflush(stdout);
  ranmat(V);
  ranmat(Q);

  info_t info = INFO_ZERO;

#ifdef HAVE_GRID
  {
    GRID_info_t grid_info = GRID_INFO_ZERO;
    
    node0_printf("Calling for unitarization derivative via Grid\n"); fflush(stdout);
    if(MILC_PRECISION == 1)
      GRID_F3_reunit_deriv( &grid_info, V, dW, Q, grid_full);
    else
      GRID_D3_reunit_deriv( &grid_info, V, dW, Q, grid_full);

    extract_grid_info( &info, &grid_info );
  }
#else
  {
    node0_printf("Calling for unitarization derivative via MILC\n"); fflush(stdout);
    milc_u3_reunit_deriv( info, V, dW, Q );
  }
#endif

#ifdef HISQ_SVD_COUNTER
  printf("hisq_svd_counter = %d\n", hisq_svd_counter);
#endif
  
  /* If an answer file is given, read it for comparison */
  if(ansflagin == RELOAD_SERIAL){
    node0_printf("Reading the derivatives of the reunitarized links from %s\n", ansfilein);
    restore_color_matrix_scidac_to_field(ansfilein, dWans, 4, MILC_PRECISION);
    node0_printf("Checking the answer\n"); fflush(stdout);
    
    /* Check the answer */
    check_answer(dW, dWans, tol);
  }

  /* Save the updated "force" from the derivative of the reunitarized
     links if requested. */
  if (ansflagout != FORGET ){
#ifdef HAVE_QIO
    filexml = create_QCDML();
    node0_printf("Saving the derivatives of the reunitarized links to %s\n", ansfileout);
    save_color_matrix_scidac_from_field( ansfileout, filexml, 
					 "dWans", QIO_SINGLEFILE, dW, 4, MILC_PRECISION,
					 stringLFNfat);
    free_QCDML(filexml);
  }
#endif
}
