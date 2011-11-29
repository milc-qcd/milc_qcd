/********************* cl_solver_utilities.c  ***************************/
/* MIMD version 7 */

#include "generic_clover_includes.h"

/* Notes for ILU preconditioning */
/* ------------------------------------------------------------
   The matrix to be inverted is written in block even-odd form as


   M = ( R_o     -K D_oe )
       ( -K D_eo  R_e    )

   where R_o and R_e are 1 - K Clov_c/u_0^3 i sigma_mu,nu F_mu,nu
   are the site-diagonal hermitian clover matrices on odd and even
   sites.

   The dslash operators D_oe src and D_eo src are given by

   SUM_dirs ( 
      ( 1 + gamma[dir] ) * U(x,dir) * src(x+dir)
    + ( 1 - gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
   )

   with gammas defined in libraries/mb_gamma.c.

   The ILU decomposition results in M = L A U

            L                 A              U

   M = ( 1            0 ) ( R_o  0   ) ( 1  -K/R_o D_oe )
       ( -K D_eo/R_o  1 ) ( 0    M_e ) ( 0     1        )

   where

         M_e = R_e - K^2 D_eo/R_o D_oe

   acts only on even sites.

   So

       M dest = src

   implies

       A (U dest) = L^(-1) src

   or

       R_o (U dest)_o  = [L^(-1) src]_o = src_o
       M_e dest_e  = [L^(-1) src]_e

   We solve the second equation and substitute into the first:

       R_o (dest_o - K/R_o D_oe dest_e = src_o

    or      dest_o = 1/R_o (src_o + K D_oe dest_e)

 ------------------------------------------------------------ */

/* The Fermilab relative residue - squared norms used here! */

Real relative_residue(wilson_vector *p, wilson_vector *q, int parity)
{
  double residue, num, den;
  int i;
  site *s;
  
  residue = 0;
  FORSOMEPARITY(i,s,parity){
    num = (double)magsq_wvec( &(p[i]) );
    den = (double)magsq_wvec( &(q[i]) );
    residue += (den==0) ? 1.0 : (num/den);
  }

  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return residue/volume;
  else
    return 2*residue/volume;
}

/* ---------  src = L^(-1)*src  ------------- */

double ilu_xfm_source( 
     wilson_vector *dest,
     wilson_vector *r,
     wilson_vector *mp,
     Real Kappa,
     int *is_startede,
     msg_tag *tage[]
     )
{
  double sr, size_src;
  int i;
  site *s;

  /* src_e = srce_e + K D_eo/R_o srce_o */
  /* (leaving src_o = src_o)   */

  /* mp_o = 1/R_o srce_o */
  mult_this_ldu_field(gen_clov, r, mp, ODD);
  /* mp_e = D_eo/R_o srce_o */
  dslash_w_field_special(mp, mp, PLUS, EVEN, tage, *is_startede);
  *is_startede = 1;
  
  /* Normalization  */
  sr = 0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(r[i]), &(mp[i]), Kappa, &(r[i]) );
    sr += (double)magsq_wvec( &(r[i]) );
  }
  g_doublesum(&sr);
  size_src = (Real)sqrt(sr);

  return size_src;
}

/* ------------------------------------------------------------ */

/* Computes two parts of the ILU linear operator:

   drd_in_e = D_eo/R_o D_oe in_e

      and

   r_in_e = R_e in_e

*/

void ilu_DRD(
  wilson_vector *in,
  wilson_vector *drd_in,
  wilson_vector *r_in,
  wilson_vector *tmpo,
  int isign,
  msg_tag* tago[],
  int *is_startedo,
  msg_tag* tage[],
  int *is_startede
  )
{
  /* r_in_e = R_e in_e */
  mult_this_ldu_field(gen_clov, in, r_in, EVEN);
  /* drd_in_o = D_oe in_e */
  dslash_w_field_special(in, tmpo, isign, ODD, tago, *is_startedo);
  *is_startedo = 1;
  /* r_in_o = 1/R_o D_oe in_e */
  mult_this_ldu_field(gen_clov, tmpo, r_in, ODD);
  /* drd_in_e = D_eo/R_o D_oe in_e */
  dslash_w_field_special(r_in, drd_in, isign, EVEN, tage, *is_startede);
  *is_startede = 1;
}

/* ------------------------------------------------------------ */
/* Reconstruct solution on odd sites */

void ilu_xfm_dest(
  wilson_vector *dest,   /* dest_o = src_o and dest_e is solution */
  wilson_vector *mp,    /* Temporary */
  Real Kappa, 
  int *is_startedo,
  msg_tag *tago[])
{
  int i;
  site *s;

  /* --------- dest_o = U^(-1)_oo/R_o dest_o + U^(-1)_oe dest_e --------- */
  /* mp_o = D_oe * dest_e */
  dslash_w_field_special(dest, mp, PLUS, ODD, tago, *is_startedo);
  *is_startedo = 1;

  /* mp_o = dest_o + K D_oe * dest_e (remember dest_o = original src_o still)*/
  FORODDSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(dest[i]), &(mp[i]), Kappa, &(mp[i]) );
  }
  /* dest_o = 1/R_o dest_o + K/R_o D_oe * dest_e */
  mult_this_ldu_field(gen_clov, mp, dest, ODD);
}

/* ------------------------------------------------------------ */
/* Convert Dirac Wilson to Dirac clover parameters */

void map_dwp_to_dcp(dirac_clover_param *dcp, dirac_wilson_param *dwp)
{
  dcp->Kappa = dwp->Kappa;     /* hopping */
  dcp->Clov_c = 0;             /* No clover term */
  dcp->U0 = 1;                 /* No tadpole factor */
}

