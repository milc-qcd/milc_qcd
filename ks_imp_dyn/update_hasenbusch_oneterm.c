/********** update_hasenbusch.c *****************************************/
/* MIMD version 7 */

/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.

 In this version the second mass mass2 is the Hasenbusch
 regulator mass for the first mass mass1.  The regulator carries
 no Naik epsilon, so the number of Naik terms must be 1.

 ----------------------------------------------------------------------
 ALGORITHM (see phi_hasenbusch_algorithm.ltx)

 The light-quark determinant is split as

     det(M1^dag M1) = det(N^dag N) det(M2^dag M2),      N = M1 M2^{-1}

 with one pseudofermion field for each factor,

     S1 = Phi1^dag (N^dag N)^{-1} Phi1
     S2 = Phi2^dag (M2^dag M2)^{-1} Phi2

 For staggered fermions M1^dag M1 = D^dag D + 4 m1^2, so

     (N^dag N)^{-1} = 1 + 4 (m2^2 - m1^2) (M1^dag M1)^{-1}
     S1             = Phi1^dag Phi1 + 4 (m2^2 - m1^2) Phi1^dag X1,
                               X1 = (M1^dag M1)^{-1} Phi1
     dS1            = -4 (m2^2 - m1^2) X1^dag [M1^dag dM1 + h.c.] X1

 The Phi1^dag Phi1 term is independent of the gauge field, so it is a
 constant of the trajectory and cancels in Delta S.

 IMPLEMENTATION NOTE.  Rather than carry the factor 4(m2^2 - m1^2)
 through d_action() and update_h(), we absorb it into the pseudofermion
 field at the time it is generated:

     phi1 = 2 sqrt(m2^2 - m1^2) Phi1     (requires m2 > m1)

 Then, with xxx1 = (M1^dag M1)^{-1} phi1,

     phi1^dag xxx1     = 4 (m2^2 - m1^2) Phi1^dag X1        (the S1 term)
     xxx1 (...) xxx1^dag = 4 (m2^2 - m1^2) X1 (...) X1^dag    (the dS1 term)

 Both are exactly what d_action.c and update_h.c already compute for an
 ordinary two-mass phi-algorithm term.  The fermion force is bilinear in
 xxx1, which is why the square root appears.  Consequently NEITHER
 d_action.c NOR update_h.c NOR lattice.h requires any modification, and
 downstream of the pseudofermion generation this routine is identical to
 update.c.

 FLAVOR WEIGHTS.  Both pseudofermions describe the SAME light quark, so
 the input file must set nflavors1 = nflavors2 = the number of light
 flavors (4 for exact HMC with this code).
 ----------------------------------------------------------------------
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"

/*---------------------------------------------------------------------*/
/* Generate the pseudofermion field for the Hasenbusch-corrected light
   quark term.  Implements, for Gaussian random R on all sites,

       Z    = M1^dag R1                       (both parities)
       Phi1  = N^dag R1 = (M2^dag M2)^{-1} M2 Z          (even sites)
       phi1 = 2 sqrt(m2^2 - m1^2) Phi1                   (even sites)

   The operator order used for Phi1 is the cheap one: applying M2 to Z
   first requires Z on both parities but leaves only an even-parity
   solve.  (Solving first and applying M2 afterwards would need the
   solution on both parities.)

   Returns the number of CG iterations used.                           */

static int hasenbusch_pseudofermion(Real m1, Real m2, imp_ferm_links_t *fn,
				    Real *fermion_action_const){
  Real final_rsq;
  int iters = 0;

  if(m2 <= m1){
    node0_printf("hasenbusch_pseudofermion: need mass2 > mass1, but got %e <= %e\n",
		 (double)m2, (double)m1);
    terminate(1);
  }

  /* phi1 <- Z = M1^dag R1 on ALL sites.  Both parities of Z are needed to
     form M2 Z on the even sites below. */
  grsource_imp( F_OFFSET(phi1), m1, EVENANDODD, fn );

  /* xxx1 <- M2 Z on the even sites.  xxx1 is used here only as scratch;
     it is cleared again before it is needed as a CG solution. */
  ks_dirac_op_site( F_OFFSET(phi1), F_OFFSET(xxx1), m2, EVEN, fn );

  /* Clear phi1: discards the odd sites of Z, which are no longer
     needed, and supplies a zero starting guess for the solve. */
  clear_latvec( F_OFFSET(phi1), EVENANDODD );

  /* phi1 <- (M2^dag M2)^{-1} M2 Z = M2^{dag,-1} Z = Phi1   (even sites) */
  iters += ks_congrad( F_OFFSET(xxx1), F_OFFSET(phi1), m2,
		       niter, nrestart, rsqmin, MILC_PRECISION, EVEN,
		       &final_rsq, fn );

  /* phi1 <- 2 sqrt(m2^2 - m^2) Phi1 */
  Real cc = (Real)(2.0*sqrt((double)m2*(double)m2
		       - (double)m1*(double)m1));
  scalar_mult_latvec( F_OFFSET(phi1), cc, F_OFFSET(phi1), EVEN );
  *fermion_action_const = 1./(cc*cc); /* For the inert constant additive term in fermion action */

  /* Zero starting guess for the mass1 solve */
  clear_latvec( F_OFFSET(xxx1), EVENANDODD );

  return iters;
}

/*---------------------------------------------------------------------*/

int update()  {
  int step, iters=0;
  int n;
  Real final_rsq;
#ifdef HMC_ALGORITHM
  double startaction=0., endaction, d_action(Real cc);
  Real xrandom;
  Real fermion_action_const = 0;
#endif
  imp_ferm_links_t *fn;

  /* refresh the momenta */
  ranmom();

  /* The Hasenbusch regulator carries no Naik epsilon, so only one set
     of fat/long links is allowed */
  n = fermion_links_get_n_naiks(fn_links);
  if(n != 1){
    node0_printf("update: Naik epsilon is not supported with the Hasenbusch term\n");
    terminate(1);
  }

  /* do "steps" microcanonical steps  */
  for(step=1; step <= steps; step++){

#ifdef PHI_ALGORITHM
    /* generate the pseudofermion fields only at the start of the
       trajectory.  Also clear xxx, since zero is our best guess for the
       solution with a new random phi field. */
    if(step==1){
      restore_fermion_links_from_site(fn_links, MILC_PRECISION);
      fn = get_fm_links(fn_links, 0);
      /* Hasenbusch-corrected light term:
	 phi1 = 2 sqrt(m2^2 - m1^2) N^dag R1 */
      clear_latvec( F_OFFSET(xxx1), EVENANDODD );
      iters += hasenbusch_pseudofermion( mass1, mass2, fn, &fermion_action_const );
      /* Regulator term: Phi2 = M2^dag R', the ordinary construction */
      clear_latvec( F_OFFSET(xxx2), EVENANDODD );
      grsource_imp( F_OFFSET(phi2), mass2, EVEN, fn);
      destroy_fn_links(fn);
    }

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (M1^dag M1)inverse * phi */
    if(step==1){
      restore_fermion_links_from_site(fn_links, MILC_PRECISION);
      fn = get_fm_links(fn_links, 0);
      iters += ks_congrad_two_src( F_OFFSET(phi1), F_OFFSET(phi2),
				   F_OFFSET(xxx1), F_OFFSET(xxx2),
				   mass1, mass2, niter, nrestart, rsqmin,
				   MILC_PRECISION, EVEN, &final_rsq, fn );
      destroy_fn_links(fn);
      startaction=d_action(fermion_action_const);
      /* copy link field to old_link */
      gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
    }
#endif

    /* update U's to middle of interval */
    update_u(0.5*epsilon);

#else /* "R" algorithm */
    /* first update the U's to special time interval */
    /* and generate a pseudofermion configuration */
    /* nflavors1 and nflavors2 are equal here, so the middle update_u()
       is a no-op; it is kept for symmetry with update.c */

    update_u(epsilon*(0.5-nflavors1/8.0));
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links, 0);
    clear_latvec( F_OFFSET(xxx1), EVENANDODD );
    iters += hasenbusch_pseudofermion( mass1, mass2, fn );
    destroy_fn_links(fn);

    update_u(epsilon*((nflavors1-nflavors2)/8.0));

    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links, 0);
    clear_latvec( F_OFFSET(xxx2), EVENANDODD );
    grsource_imp( F_OFFSET(phi2), mass2, EVEN, fn);
    destroy_fn_links(fn);

    /* update U's to middle of interval */
    update_u(epsilon*nflavors2/8.0);
#endif

    /* do conjugate gradient to get xxx1 = (M1^dag M1)^{-1} * phi1 */
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links, 0);
    iters += ks_congrad_two_src( F_OFFSET(phi1), F_OFFSET(phi2),
				 F_OFFSET(xxx1), F_OFFSET(xxx2),
				 mass1, mass2, niter, nrestart, rsqmin,
				 MILC_PRECISION, EVEN, &final_rsq, fn );
    /* Combine xxx on the even sites with Dslash xxx on the odd sites,
       as the fermion force routine expects */
    dslash_site( F_OFFSET(xxx1), F_OFFSET(xxx1), ODD, fn);
    dslash_site( F_OFFSET(xxx2), F_OFFSET(xxx2), ODD, fn);
    destroy_fn_links(fn);

    /* now update H by full time interval */
    update_h(epsilon);

    /* update U's by half time step to get to even time */
    update_u(epsilon*0.5);

    /* reunitarize the gauge field */
    reunitarize_ks();

  }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
  /* find action */
  /* do conjugate gradient to get xxx1 = (M1^dag M1)^{-1} * phi1 and
     (M2^dag M2)^{-1} * phi2*/
  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
  fn = get_fm_links(fn_links, 0);
  iters += ks_congrad_two_src( F_OFFSET(phi1), F_OFFSET(phi2),
			       F_OFFSET(xxx1), F_OFFSET(xxx2),
			       mass1, mass2, niter, nrestart, rsqmin,
			       MILC_PRECISION, EVEN, &final_rsq, fn );
  destroy_fn_links(fn);
  endaction=d_action(fermion_action_const);
  /* decide whether to accept, if not, copy old link field back */
  /* careful - must generate only one random number for whole lattice */
  if(this_node==0)xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if( exp( (double)(startaction-endaction) ) < xrandom ){
    if(steps > 0)
      gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
#ifdef FN
    invalidate_fermion_links(fn_links);
#endif
    node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
  }
  else {
    node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
  }
#endif

  if(steps > 0)return (iters/steps);
  else return(-99);
}


/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_int_alg_opt_chr( void )
{
  return "INT_ALG_NEEDS_TO_BE_FIXED";
}
