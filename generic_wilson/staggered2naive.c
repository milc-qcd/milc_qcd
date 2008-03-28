/******************* staggered2naive.c *********************************/
/* MIMD version 7 */

/* Convert a staggered propagator to a naive propagator (Wilson prop field) */
/* Good only for a staggered propagator originating from a point source */

/* 6/02/07 C. DeTar */

#include "generic_wilson_includes.h"
#include "../include/gammatypes.h"
#include <string.h>

static gamma_matrix_t omega0_dag;
gamma_matrix_t g1;
gamma_matrix_t gx;
gamma_matrix_t gy;
gamma_matrix_t gz;
gamma_matrix_t gt;

static void init_gamma(void){
  g1 = gamma_mat(G1);
  gx = gamma_mat(GX);
  gy = gamma_mat(GY);
  gz = gamma_mat(GZ);
  gt = gamma_mat(GT);
}

static void mult_omega_l(gamma_matrix_t *omega, int x, int y, int z, int t)
{
  gamma_matrix_t omtmp;
  
  /* Definition follows our KS phase conventions */

  /**
  if(z % 2 == 0){ mult_gamma_by_gamma( &gz, omega, &omtmp ); *omega = omtmp; }
  if(y % 2 == 0){ mult_gamma_by_gamma( &gy, omega, &omtmp ); *omega = omtmp; }
  if(x % 2 == 0){ mult_gamma_by_gamma( &gx, omega, &omtmp ); *omega = omtmp; }
  if(t % 2 == 0){ mult_gamma_by_gamma( &gt, omega, &omtmp ); *omega = omtmp; }
  **/
  
  if(t % 2 == 0){ mult_gamma_by_gamma( &gt, omega, &omtmp ); *omega = omtmp; }
  if(x % 2 == 0){ mult_gamma_by_gamma( &gx, omega, &omtmp ); *omega = omtmp; }
  if(y % 2 == 0){ mult_gamma_by_gamma( &gy, omega, &omtmp ); *omega = omtmp; }
  if(z % 2 == 0){ mult_gamma_by_gamma( &gz, omega, &omtmp ); *omega = omtmp; }
  
}

static void m_by_omega_swv(spin_wilson_vector *swv, su3_vector *ksp, site *s,
			   int ks_source_r[]){

  int s0;
  spin_wilson_vector swvtmp;
  gamma_matrix_t omega;

  /* Initialize the propagator */
  
  memset(swv, 0, sizeof(spin_wilson_vector));

  /* Start with WP = I_spin x KS_color where I_spin is a unit matrix
     in spin */

  for(s0 = 0; s0 < 4; s0++)
    su3vec_copy(ksp, &swv->d[s0].d[s0]);


  /* Construct (omega(sink) omega0^+(source))^* */
  omega = omega0_dag;
  mult_omega_l(&omega, s->x, s->y, s->z, s->t);
  //  gamma_conj( &omega, &omega );
  
  mult_sw_by_gamma_mat_l( swv, &swvtmp, &omega); 
  *swv = swvtmp;
}

static void m_by_omega(wilson_propagator *wp, su3_matrix *ksp, site *s,
		int ks_source_r[]){

  int s0,c0,c1;
  wilson_propagator wptmp;
  gamma_matrix_t omega;

  /* Initialize the propagator */
  
  memset(wp, 0, sizeof(wilson_propagator));

  /* Start with WP = I_spin x KS_color where I_spin is a unit matrix
     in spin */

  for(s0 = 0; s0 < 4; s0++)
    for(c0 = 0; c0 < 3; c0++)
      for(c1 = 0; c1 < 3; c1++)
	wp->c[c0].d[s0].d[s0].c[c1] = ksp[c0].e[c0][c1];


  /* Construct (omega(sink) omega0^+(source))^* */
  omega = omega0_dag;
  mult_omega_l(&omega, s->x, s->y, s->z, s->t);
  //  gamma_conj( &omega, &omega );
  
  for(c0=0;c0<3;c0++){
    mult_sw_by_gamma_mat_l( &(wp->c[c0]), &(wptmp.c[c0]), &omega); 
    wp->c[c0] = wptmp.c[c0]; 
  }
}

/* r is the KS source point */
void convert_ksprop_to_wprop(wilson_propagator *wp, su3_matrix *ksp, int r[])
{
  int i;
  site *s;
  gamma_matrix_t omega0;

  init_gamma();

  /* Construct source omega adjoint */
  omega0 = g1;   /* Unit gamma */
  mult_omega_l(&omega0, r[0], r[1], r[2], r[3]);
  gamma_adj(&omega0_dag, &omega0);

  FORALLSITES(i,s){
    m_by_omega( wp+i, ksp+3*i, s, r);
  }

}
/* r is the KS source point */
void convert_ksprop_to_wprop_swv(spin_wilson_vector *swv, 
				 su3_vector *ksp, int r[])
{
  int i;
  site *s;
  gamma_matrix_t omega0;

  init_gamma();

  /* Construct source omega adjoint */
  omega0 = g1;   /* Unit gamma */
  mult_omega_l(&omega0, r[0], r[1], r[2], r[3]);
  gamma_adj(&omega0_dag, &omega0);

  FORALLSITES(i,s){
    m_by_omega_swv( swv+i, ksp+i, s, r);
  }

}
