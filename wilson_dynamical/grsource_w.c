/************************ grsource_w.c *****************************/
/* MIMD version 6 */

#include "wi_dyn_includes.h"

/* construct a gaussian random vector, g_rand, and chi=M(dagger)*g_rand  */
/* also clear psi, since zero is our best guess for the solution with a
   new random chi field. */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:

  without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
  with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
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
    /* use g_rand on odd sites as temporary storage:
        g_rand(odd) = Dslash^adjoint * g_rand(even) */
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(g_rand), MINUS, ODD);
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi)   , MINUS, EVEN);
    FOREVENSITES(i,s){
        scalar_mult_add_wvec( &(s->g_rand), &(s->chi), -kappa*kappa,
            &(s->chi) );
    }
#else
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi), MINUS, EVENANDODD);
    FORALLSITES(i,s){
        scalar_mult_add_wvec( &(s->g_rand), &(s->chi), -kappa,
	    &(s->chi) );
    }
#endif
}/* grsource_w */


/* Check congrad by multiplying psi by Madj*M, compare result to chi */
/* Before calling checkmul() you should call grsource_w() and congrad_w() */
void checkmul() {
register int i,j,k;
register site *s;
    printf("CHECKMUL: starting\n");

    /* multiply by M_adjoint*M */
#ifdef LU
    dslash_w_site( F_OFFSET(psi), F_OFFSET(psi), PLUS, ODD );
    dslash_w_site( F_OFFSET(psi), F_OFFSET(mp ), PLUS, EVEN);
    FOREVENSITES(i,s)
        scalar_mult_add_wvec( &(s->psi), &(s->mp), -kappa*kappa, &(s->mp) );
   dslash_w_site( F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(p) , MINUS, EVEN);
    FOREVENSITES(i,s)
        scalar_mult_add_wvec( &(s->mp), &(s->p), -kappa*kappa, &(s->p));
#else
    dslash_w_site( F_OFFSET(psi), F_OFFSET(mp), PLUS, EVENANDODD);
    FORALLSITES(i,s)
        scalar_mult_add_wvec( &(s->psi), &(s->mp), -kappa, &(s->mp) );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(p), MINUS, EVENANDODD);
    FORALLSITES(i,s)
        scalar_mult_add_wvec( &(s->mp), &(s->p), -kappa, &(s->p) );
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
