/************************* sphaleron.c *******************************/
/* module functions exclusive to the sphaleron rate calculation      */

#include "wilson_flow_includes.h"
// #include "../include/field_strength.h"
#include <string.h>

#ifdef SPHALERON

#if 0

  // double maxdev = 0.0;
  // int ig;
  // node0_printf("END RESULT?\n");
  // FORALLSITES(i, s)
  // IF_BLOCKED(s, block_stride) 
  // IF_ACTIVE(s) {
  // 	sub_su3_matrix( &(s->fieldstrength[FS_XY]), &(fstrength[FS_XY][i]), &tmp );
  //   for ( ig = 0; ig < 9; ig ++ ) {
  //    	if (maxdev < cabs( ((complex*)(&tmp) + ig) )) {
  //     	node0_printf("i %d: s %d %d %d %d; ig %d\n",i,s->x,s->y,s->z,s->t, ig);
  //     	dumpmat( &tmp );
  //     }
  //    	maxdev = ( maxdev > cabs( ((complex*)(&tmp) + ig) ) 
  //    					 ? maxdev : cabs( ((complex*)(&tmp) + ig) ) );
  //   }
 	// }
  // node0_printf("maxdev %.16g in FS_XY\n",maxdev);

#endif

/* Special unitarize a unitary matrix:
   B=A/det(A)^1/3
   A is a U(3) matrix and B is an SU(3) matrix.
   This function also returns det(A) */
void su3_iterative_special_unitarize( su3_matrix *a, su3_matrix *b, complex *detA ) {
  complex s;
  Real argdet, r;

  #define MAXITER 20
  #if (MILC_PRECISION==1)
  #define TOL 1.e-6
  #else
  #define TOL 1.e-14
  #endif
  int iter,ic,jc;
  double epsmax = 0;
  for (ic=0; ic<3; ic++) {
		for (jc=0; jc<3; jc++) {
			epsmax += cabs( &(a->e[ic][jc]) );
		}
	}
	double lastepsmax = epsmax + 1.;
	su3_matrix aadag, eps, aprime, temp;
	su3_matrix *anew = &aprime, *aold = a; 
	for ( iter=0; iter < MAXITER && epsmax > TOL && lastepsmax > epsmax; iter++ ) {
		aold = ( (iter % 2) ? &aprime : a);
		anew = ( (iter % 2) ? a : &aprime );
		mult_su3_na( aold, aold, &aadag );
		su3mat_copy( &aadag, &eps );
		eps.e[0][0].real -= 1.0;
		eps.e[1][1].real -= 1.0;
		eps.e[2][2].real -= 1.0;
		mult_su3_nn( &eps, aold, &temp );
		scalar_mult_sub_su3_matrix( aold, &temp, 0.5, anew );
		lastepsmax = epsmax;
		epsmax = 0.;
		for (ic=0; ic<3; ic++) {
			for (jc=0; jc<3; jc++) {
				epsmax = ( epsmax > cabs( &(eps.e[ic][jc]) ) 
					       ? epsmax : cabs( &(eps.e[ic][jc]) ) );
			}
		}
		// dumpmat( aold );
		// dumpmat( &aadag );
		// dumpmat( &eps );
		// dumpmat( anew );
		// node0_printf("iter %02d epsmax %.6e, TOL %.1e\n",iter,epsmax,TOL);
	}

  (*detA) = det_su3( anew );

  /* take third root, choose argument in [-pi/3, pi/3) branch */
  argdet = carg( detA );
  r = cabs( detA );
  r = exp( log(r)/3 );
  argdet /= 3;

  /* multiply the matrix by inverse root */
  s.real = cos( argdet ) / r;
  s.imag = -sin( argdet ) / r;
  c_scalar_mult_su3mat( anew, &s, b );
}

su3_matrix ** new_half_links( void ) {

  int i;
  su3_matrix **this = new_field( N_HALF );
  renew_half_links( this );

  return ( this );
}

void renew_half_links( su3_matrix **link_half ) {
	register int i;
	register int ic,jc;
	int dir;
	register site *s = NULL;
	su3_matrix tmp;
	msg_tag *tag=NULL;
	complex detA;

	if ( link_half == NULL ) {
		node0_printf("Tried to clear unallocated link_half\n"); 
		terminate(1);
	}
	if ( link_half[0] == NULL ) {
		node0_printf("Tried to clear unallocated link_half\n"); 
		terminate(1);
	}

	FORALLUPDIRBUT(TUP,dir) {

		tag = start_gather_site( F_OFFSET(link[dir]) , sizeof(su3_matrix),
														TUP, EVENANDODD, gen_pt[0]);

    wait_gather(tag);
		FORALLSITES(i,s) 
		IF_BLOCKED(s, block_stride)
		// IF_ACTIVE(s) 
		// {
		// 	su3mat_copy( &(s-> link[dir]), &tmp );
		IF_LOWER_BULK(s) 
		{
			add_su3_matrix( &(s-> link[dir]), (su3_matrix *)(gen_pt[0][i]), &tmp );
			for (ic=0; ic<3; ic++) 
				for (jc=0; jc<3; jc++) 
				 	CMULREAL(tmp.e[ic][jc],0.5,tmp.e[ic][jc]);
			su3_iterative_special_unitarize( &tmp, &(link_half[dir][i]), &detA );
		}
		cleanup_gather( tag );
	}
}

void update_last_flow_links( su3_matrix **link_last_flow ) {
	register int i;
	int dir;
	register site *s = NULL;

	if ( link_last_flow== NULL ) {
		node0_printf("Tried to use unallocated link_last_flow1\n"); 
		terminate(1);
	}
	if ( link_last_flow[0] == NULL ) {
		node0_printf("Tried to use unallocated link_last_flow1\n"); 
		terminate(1);
	}
	FORALLUPDIRBUT(TUP,dir) 
		FORALLSITES(i, s) 
		IF_BLOCKED(s, block_stride)
		IF_BOUNDARY(s)
			su3mat_copy( &(s->link[dir]), &(link_last_flow[dir][i]) );
}

/** set all inactive links to the identity to test 
 *  whether they contribute where they should not **/
static 
void kill_all_inactive_links( void ) {
  register int i, dir;
  register site *s = NULL;

  FORALLUPDIR(dir) {
    FORALLSITES(i,s) {
      if (!ACTIVE_COND_LINK(s, dir)) {
        set_identity( &(s -> link[TUP]) );
      }
    } 
  }
}

void prepare_bulk_links( void ) {
	register int i;
	register site *s = NULL;

	ax_gauge();
	node0_printf("Fixed to axial gauge\n");
	fflush(stdout);

	FORALLSITES(i,s) {
		if (s->t < nt -1) {
			set_identity( &(s -> link[TUP]) );
		} 
	}
  // kill_all_inactive_links();
}

void report_bulk( Real time, Real *q_bulk ) {

  /* Wilson flow output variables */
  double Et_WS[2], Es_WS[2];
  double Et_C[2], Es_C[2], charge[4];
  double q_ret;
  
  /* Print flow output column labels */
  node0_printf("#LABEL      time Clover_t Clover_s iClover_t iClover_s Plaq_t Plaq_s Rect_t Rect_s charge icharge\n");
  fflush(stdout);

  fmunu_fmunu_full( Et_C, Es_C, charge );
  gauge_action_w_s_full( Et_WS, Es_WS );
  print_observables( "REPT_FULL:", time, Et_WS, Es_WS, Et_C, Es_C, charge );

  fmunu_fmunu_bulk( Et_C, Es_C, charge );
  gauge_action_w_s_bulk( Et_WS, Es_WS );
  print_observables( "REPT_BULK:", time, Et_WS, Es_WS, Et_C, Es_C, charge );
  q_bulk[0] = charge[1];
  fmunu_fmunu_half( Et_C, Es_C, charge );  
  gauge_action_w_s_half( Et_WS, Es_WS );
  print_observables( "REPT_HALF:", time, Et_WS, Es_WS, Et_C, Es_C, charge );
  fflush(stdout);
  q_bulk[1] = charge[1];
}

void bulk_flow( Real *q_bulk ) {

	int i, ibulk,thr_bulk = 0;
	Real q_last[2] = {0.0, 0.0};
	Real delta_q2[2]; 
	Real thresh_q2 = ( qthr_bulk * qthr_bulk < qs_tol * qs_tol 
                   ? qthr_bulk * qthr_bulk : qs_tol * qs_tol );
	prepare_bulk_links();
  report_bulk( stoptime, q_last );    
  for ( ibulk = 1; 
        ibulk <= maxnflow_bulk && thr_bulk < minqthr_bulk; 
        ibulk++ ) {
    /* integrate the flow */
    run_gradient_flow( BULK );
    /* report the status */
  	report_bulk( stoptime + ibulk * stoptime_bulk, q_bulk );    
    /* check convergence */
    for ( i = 0; i < 2; i++ ) 
    	delta_q2[i] = ( q_bulk[i] - q_last[i] ) * ( q_bulk[i] - q_last[i] );
    if ( delta_q2[0] > thresh_q2 && delta_q2[1] > thresh_q2 ) {
      thr_bulk = 0;
    } else 
      // thr_bulk++;
    for ( i = 0; i < 2; i++ ) 
    	q_last[i] = q_bulk[i];
  }
  q_bulk[2] = ( q_bulk[1] - q_bulk[0] ) * ( q_bulk[1] - q_bulk[0] );
}

void bdry_flow( Real *q_bulk ) {

	int i, ibdry, thr_bdry = 0, thr_qs2 = 0;
  Real q_last[2] = {0.0, 0.0}, delta_q2[2];
  Real q_s[2], q_int[2], delta_qs2[2]; 
  Real thresh_qs2 = qs_tol * qs_tol;
  Real thresh_q2 = ( qthr_bdry * qthr_bdry < thresh_qs2 
                ? qthr_bdry * qthr_bdry : thresh_qs2 );
  thresh_qs2 = ( q_bulk[2] > thresh_qs2 ? q_bulk[2] : thresh_qs2 );
  thresh_qs2 = ( thresh_q2 > thresh_qs2 ? thresh_q2 : thresh_qs2 );
   
#ifdef BLOCKING
  Real block_time, rest_time;    
#endif
  for ( ibdry = 1, q_acc[0] = 0.0, q_acc[1] = 0.0;
        ( ibdry <= maxnflow_bdry
       &&   thr_bdry < 2 * minqthr_bdry
       && ( thr_bdry < minqthr_bdry 
         || thr_qs2  < minqthr_bdry ) ); 
        ibdry++ ) {
    #ifdef BLOCKING
      block_time = ( block_stride == 1 ? block_1to2_time : 
                     block_stride == 2 ? block_2to4_time :
                     block_stride == 4 ? block_4to8_time : -1.0 );
      rest_time = 0.0;
      if ( ( ibdry - 1 ) * stoptime_bdry <= block_time 
        && block_time < ibdry  * stoptime_bdry 
        ) {
        rest_time = ibdry * stoptime_bdry - block_time;
        block_time = stoptime_bdry - rest_time;
        stoptime_bdry = block_time;
      }
    #endif

    if ( stoptime_bdry > 0.0 ) {
      /* integrate the flow until blocking */
      run_gradient_flow( BOUNDARY );
    }
    #ifdef BLOCKING
      if ( rest_time > 0.0 ) {
        /* block one step */
        spatial_blocking();
        stoptime_bdry = rest_time;
        /* integrate the flow for the rest */
        run_gradient_flow( BOUNDARY );
        stoptime_bdry = block_time + rest_time;
      }
    #endif

    /* check convergence */
    for ( i = 0; i < 2; i++ ) {
      delta_q2[i] = ( q_acc[i] - q_last[i] ) * ( q_acc[i] - q_last[i] );
      q_s[i] = q_bulk[i] + q_acc[1];
      q_int[i] = round( q_s[i] );
      delta_qs2[i] = (q_s[i] - q_int[i] ) * ( q_s[i] - q_int[i] );
      q_last[i] = q_acc[i];
    }
    node0_printf("FULL_ACCUM @ %02d",
      ibdry);
    node0_printf(" q_acc (%.6g,%.6g) delta_q2 (%.6g,%.6g)",
      q_acc[0],q_acc[1], delta_q2[0],delta_q2[1]);
    if ( delta_q2[0] < thresh_q2 || delta_q2[1] < thresh_q2 ) {
      // thr_bdry++;
      node0_printf(" PASS");
    } else {
      thr_bdry = 0;
    }
    node0_printf(" q_s = (%.6g, %.6g) delta_qs2 = (%.6g, %.6g)",
      q_s[0],q_s[1], delta_qs2[0],delta_qs2[1]);
    if ( delta_qs2[0] < thresh_qs2 || delta_qs2[1] < thresh_qs2 ) {
      // thr_qs2++;
      node0_printf(" PASS");
    } else {
      thr_qs2 = 0;
    }
    node0_printf("\n");
  }
}


#endif