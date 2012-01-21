/******************* staggered2naive.c *********************************/
/* MIMD version 7 */

/* Convert a staggered propagator to a naive propagator (Wilson prop field) */
/* Good only for a staggered propagator originating from a source with support
   only on the origins of hypercubes */

/* 6/02/07 C. DeTar */

#include "generic_wilson_includes.h"
#include "../include/gammatypes.h"
#include <string.h>

gamma_matrix_t g1;
gamma_matrix_t gx;
gamma_matrix_t gy;
gamma_matrix_t gz;
gamma_matrix_t gt;
gamma_matrix_t g5;

static void init_gamma(void){
  g1 = gamma_mat(G1);
  gx = gamma_mat(GX);
  gy = gamma_mat(GY);
  gz = gamma_mat(GZ);
  gt = gamma_mat(GT);
  g5 = gamma_mat(G5);
}

#ifdef OLD_STAGGERED2NAIVE

/* Kept for backward compatibility */

static gamma_matrix_t gamma0_dag;

/* omega <-- gz^z gy^y gx^x gt^t * omega */

static void mult_omega_l(gamma_matrix_t *omega, int x, int y, int z, int t)
{
  gamma_matrix_t omtmp;
  
  /* Definition follows our KS phase conventions */

  if(t % 2 == 0){ mult_gamma_by_gamma( &gt, omega, &omtmp ); *omega = omtmp; }
  if(x % 2 == 0){ mult_gamma_by_gamma( &gx, omega, &omtmp ); *omega = omtmp; }
  if(y % 2 == 0){ mult_gamma_by_gamma( &gy, omega, &omtmp ); *omega = omtmp; }
  if(z % 2 == 0){ mult_gamma_by_gamma( &gz, omega, &omtmp ); *omega = omtmp; }
  
}

static void m_by_omega_swv(spin_wilson_vector *swv, su3_vector *ksp, site *s, int r0[]){

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
  omega = gamma0_dag;
  mult_omega_l(&omega, s->x-r0[0], s->y-r0[1], s->z-r0[2], s->t-r0[3]);
  
  mult_sw_by_gamma_mat_l( swv, &swvtmp, &omega); 
  *swv = swvtmp;
}

/* Convert a color vector field (usually from a propagator from a
   single color source) to a naive propagator as a spin-wilson vector
   field */
/* r is the KS source point */

void convert_ksprop_to_wprop_swv(spin_wilson_vector *swv, 
				 su3_vector *ksp, int r[], int r0[])
{
  int i;
  site *s;
  gamma_matrix_t gamma0;

  init_gamma();

  /* Construct source gamma adjoint */
  gamma0 = g1;   /* Unit gamma */
  mult_omega_l(&gamma0, r[0]-r0[0], r[1]-r0[1], r[2]-r0[2], r[3]-r0[3]);
  gamma_adj(&gamma0_dag, &gamma0);

  FORALLSITES(i,s){
    m_by_omega_swv( swv+i, ksp+i, s, r0);
  }
}

#else

/* Multiply omega on the right by the adjoint of the Kawamoto-Smit
   transformation matrix */

/* omega <--  gt^t gx^x gy^y gz^z * omega */

static void mult_omega_l(gamma_matrix_t *omega, int x, int y, int z, int t)
{
  gamma_matrix_t omtmp;
  
  /* Definition follows our KS phase conventions */

  if(z & 0x1){ mult_gamma_by_gamma( &gz, omega, &omtmp ); *omega = omtmp; }
  if(y & 0x1){ mult_gamma_by_gamma( &gy, omega, &omtmp ); *omega = omtmp; }
  if(x & 0x1){ mult_gamma_by_gamma( &gx, omega, &omtmp ); *omega = omtmp; }
  if(t & 0x1){ mult_gamma_by_gamma( &gt, omega, &omtmp ); *omega = omtmp; }
  
}

/* Convert a color vector field (usually from a propagator from a
   single color source) to a naive propagator as a spin-wilson vector
   field */
/* r is the KS source point */

void convert_ksprop_to_wprop_swv(spin_wilson_vector *swv, 
				 su3_vector *v, int r[], int r0[])
{
  int i;
  site *s;
  gamma_matrix_t gamma0, gamma0_dag, omega, omega5;

  init_gamma();

  /* Initialize the propagator */
  memset(swv, 0, sites_on_node*sizeof(spin_wilson_vector));
  
  /* Construct gamma0 = Gamma(source)= g3^(r[3]-r0[3]) g0^(r[0]-r0[0]) g1^(r[1]-r0[1]) g2^(r[2]-r0[2]) */

  gamma0 = g1;   /* Unit gamma */
  mult_omega_l(&gamma0, r[0]-r0[0], r[1]-r0[1], r[2]-r0[2], r[3]-r0[3]);

  /* Then gamma0_dag = g2^(r[2]-r0[2]) g1^(r[1]-r0[1]) g0^(r[0]-r0[0]) g3^(r[3]-r0[3]) */
  gamma_adj(&gamma0_dag, &gamma0);

  /* Complete the transformation. Take the direct product with the
     color vector */

  FORALLSITES(i,s){

    /* Construct omega = g5 Gamma(sink) Gamma(source)^\dagger g5 */

    /* First construct Gamma(sink) Gamma(source)^\dagger */

    omega = gamma0_dag;
    mult_omega_l(&omega, s->x-r0[0], s->y-r0[1], s->z-r0[2], s->t-r0[3]);

    /* Then conjugate by gamma5 to map the MILC staggered Dslash convention
       to the MILC Dirac Dslash convention */
    
    mult_gamma_by_gamma(&g5, &omega, &omega5);
    mult_gamma_by_gamma(&omega5, &g5, &omega);

    /* Take the direct product with the color vector */
    direct_prod_gamma_su3_vector(swv+i, v+i, &omega);
  }
}

#endif

/* Apply the inverse Kawamoto Smit transformation */
/* This is applied in place to a Dirac vector field */
/* r is the KS source point */

void convert_naive_to_staggered_wv(wilson_vector *wv, int r[], int r0[])
{
  int i;
  site *s;
  gamma_matrix_t gamma0, gamma0_dag;
  gamma_matrix_t omega, omega_dag5, omega_dag;
  wilson_vector wvtmp;

  init_gamma();

  /* Construct gamma0 = g3^(r[3]-r0[3]) g0^(r[0]-r0[0]) g1^(r[1]-r0[1]) g2^(r[2]-r0[2]) */

  gamma0 = g1;   /* Unit gamma */
  mult_omega_l(&gamma0, r[0]-r0[0], r[1]-r0[1], r[2]-r0[2], r[3]-r0[3]);

  /* Then gamma0_dag = g2^(r[2]-r0[2]) g1^(r[1]-r0[1]) g0^(r[0]-r0[0]) g3^(r[3]-r0[3]) */
  gamma_adj(&gamma0_dag, &gamma0);

  FORALLSITES(i,s){

    /* Construct omega = g5 Gamma(sink) Gamma(source)^\dagger g5 */
    
    omega = gamma0_dag;
    mult_omega_l(&omega, s->x-r0[0], s->y-r0[1], s->z-r0[2], s->t-r0[3]);
    
    /* Then omega_dag = Gamma(source) Gamma(sink)^\dagger */
    gamma_adj(&omega_dag, &omega);
    
    /* Conjugate by gamma5 to map the MILC Dirac Dslash convention
       to the MILC staggered Dslash convention */
    
    mult_gamma_by_gamma(&g5, &omega_dag, &omega_dag5);
    mult_gamma_by_gamma(&omega_dag5, &g5, &omega_dag);

    mult_w_by_gamma_mat_l( wv+i, &wvtmp, &omega_dag ); 
    wv[i] = wvtmp;
  }
}

#ifdef DEBUG_NAIVE

#ifdef OLD_STAGGERED2NAIVE

void dslash_naive(wilson_vector *src, wilson_vector *dst){
  int i;
  wilson_vector *tmp1 = create_wv_field();
  wilson_vector *tmp2 = create_wv_field();

  /* tmp1 <- (1 + 2 Dslash src) */
  dslash_w_field(src, tmp1, +1, EVENANDODD);
  /* tmp2 <- (1 - 2 Dslash src) */
  dslash_w_field(src, tmp2, -1, EVENANDODD);

  FORALLFIELDSITES(i){
    /* dst <- -2 Dslash src */
    sub_wilson_vector(tmp2+i, tmp1+i, dst+i);
    scalar_mult_wvec(dst+i, 0.5, dst+i);
  }

  destroy_wv_field(tmp1);
  destroy_wv_field(tmp2);
}

#else

void dslash_naive(wilson_vector *src, wilson_vector *dst){
  int i;
  wilson_vector *tmp1 = create_wv_field();
  wilson_vector *tmp2 = create_wv_field();

  /* tmp1 <- (1 + 2 Dslash src) */
  dslash_w_field(src, tmp1, +1, EVENANDODD);
  /* tmp2 <- (1 - 2 Dslash src) */
  dslash_w_field(src, tmp2, -1, EVENANDODD);

  FORALLFIELDSITES(i){
    /* dst <- -2 Dslash src */
    sub_wilson_vector(tmp1+i, tmp2+i, dst+i);
    scalar_mult_wvec(dst+i, 0.5, dst+i);
  }

  destroy_wv_field(tmp1);
  destroy_wv_field(tmp2);
}

#endif

void check_naive(wilson_vector *dst, wilson_vector *src, Real mass, Real tol){
  int i;
  wilson_vector *tmp = create_wv_field();
  wilson_vector diff;
  Real d,dmag = 0;
  char myname[] = "check_naive";

  dslash_naive(dst, tmp);

  FORALLFIELDSITES(i){
    /* tmp <- Dslash dst +  m dst */
    scalar_mult_add_wvec(tmp + i, dst+i, 2.*mass, tmp+i);
    sub_wilson_vector(tmp + i, src + i, &diff);
    d = magsq_wvec(&diff);
    if(d > tol){
      printf("site %d expected \n",i);
      dump_wvec(src + i);
      printf("got\n");
      dump_wvec(tmp + i);
    }
    dmag += d;
  }
  g_floatsum(&dmag);
  dmag = sqrt(dmag);
  node0_printf("%s: mass %g. Norm diff is %g\n", myname, mass, dmag);

  destroy_wv_field(tmp);
}

#endif

