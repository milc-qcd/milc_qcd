/************************ grsource.c *****************************/
/* MIMD version 7 */
/* Clover fermions */

#include "cl_dyn_includes.h"

/* construct a gaussian random vector, g_rand, and chi=M(dagger)*g_rand  */
/* also clear psi, since zero is our best guess for the solution with a
   new random chi field. */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:

  without LU:
   M = A - kappa*( Dslash_eo + Dslash_oe )
  with LU:
   M = A_even - kappa^2 * Dslash_eo * (A_odd)^{-1} * Dslash_oe
*/
#ifdef LU
#define FORMYSITES FOREVENSITES
#else
#define FORMYSITES FORALLSITES
#endif

void grsource_w() {
register int i,j,k;
register site *s;

    FORMYSITES(i,s){
        for(k=0;k<4;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else							
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
	    s->psi.d[k].c[j] = cmplx(0.0,0.0);
        }
    }

#ifdef LU
    /*  chi <- M^dagger g_rand  */
    mult_this_ldu_site( gen_clov, F_OFFSET(g_rand), F_OFFSET(tmp), EVEN );
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(tmp), MINUS, ODD);
    mult_this_ldu_site( gen_clov, F_OFFSET(tmp), F_OFFSET(g_rand), ODD );
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi), MINUS, EVEN);
    FOREVENSITES(i,s){
        scalar_mult_add_wvec( &(s->tmp), &(s->chi), -kappa*kappa,
            &(s->chi) );
    }
#else
    mult_this_ldu_site( gen_clov, F_OFFSET(g_rand), F_OFFSET(tmp), EVENANDODD );
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi), MINUS, EVENANDODD);
    cleanup_dslash_wtemps();
    FORALLSITES(i,s){
        scalar_mult_add_wvec( &(s->tmp), &(s->chi), -kappa,
	    &(s->chi) );
    }
#endif
}/* grsource_w */


/* Check congrad by multiplying psi by Madj*M, compare result to chi */
/* Before calling checkmul() you should call grsource_cl() and congrad_w() */
void checkmul() {
register int i,j,k;
register site *s;
    printf("CHECKMUL: starting\n");

    /* multiply by M_adjoint*M */
#ifdef LU
    mult_this_ldu_site(gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVEN );
    dslash_w_site( F_OFFSET(psi), F_OFFSET(psi), PLUS, ODD );
    mult_this_ldu_site(gen_clov, F_OFFSET(psi), F_OFFSET(mp),ODD );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(mp), PLUS, EVEN);
    FOREVENSITES(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->mp), -kappa*kappa, &(s->mp) );
    mult_this_ldu_site(gen_clov,F_OFFSET(mp), F_OFFSET(tmp), EVEN);
    dslash_w_site( F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD );
    mult_this_ldu_site(gen_clov,F_OFFSET(mp), F_OFFSET(tmp), ODD);
    dslash_w_site( F_OFFSET(tmp), F_OFFSET(p) , MINUS, EVEN);
    FOREVENSITES(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa*kappa, &(s->p));
#else
    mult_this_ldu_site(gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVENANDODD );
    dslash_w_site( F_OFFSET(psi), F_OFFSET(mp), PLUS, EVENANDODD);
    FORALLSITES(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->mp), -kappa, &(s->mp) );
    mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), EVENANDODD );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(p), MINUS, EVENANDODD);
    FORALLSITES(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa, &(s->p) );
#endif

    /* Compare to source */
    FORMYSITES(i,s){
	/**printf("Site %d %d %d %d\n",s->x,s->y,s->z,s->t);**/
	for(k=0;k<4;k++)for(j=0;j<3;j++){
if( fabs((double)s->chi.d[k].c[j].real-(double)s->p.d[k].c[j].real) > 1e-4 )
	    printf("%d %d %d\t%e\t%e\t%e\n",i,k,j,
		(double)s->chi.d[k].c[j].real,
		(double)s->p.d[k].c[j].real,
		(double)s->chi.d[k].c[j].real - (double)s->p.d[k].c[j].real);
if( fabs((double)s->chi.d[k].c[j].imag-(double)s->p.d[k].c[j].imag) > 1e-4 )
	    printf("%d %d %d\t%e\t%e\t%e\n",i,k,j,
		(double)s->chi.d[k].c[j].imag,
		(double)s->p.d[k].c[j].imag,
		(double)s->chi.d[k].c[j].imag - (double)s->p.d[k].c[j].imag);
	}
	printf("\n");
    }
}/* checkmul */
