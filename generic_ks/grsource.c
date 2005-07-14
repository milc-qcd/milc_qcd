/************************ grsource.c *****************************/
/* MIMD version 7 */

/* 1/27/00 combined with Schroedinger functional version - UMH */

#include "generic_ks_includes.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*g_rand  */
/* also clear xxx, since zero is our best guess for the solution with a
   new random phi field. */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */
void grsource(int parity) {
register int i,j;
register site *s;
    FORALLSITES(i,s){
#ifdef SCHROED_FUN
	if(s->t > 0){
#endif
	    for(j=0;j<3;j++){
#ifdef SITERAND
		s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
		s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
		s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
		s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
	    }
#ifdef SCHROED_FUN
	}
	else{	/* Set all fermion vectors to zero at t=0 */
	    for(j=0;j<3;j++){
		s->phi.c[j].real = s->phi.c[j].imag = 0.0;
		s->resid.c[j].real = s->resid.c[j].imag = 0.0;
		s->cg_p.c[j].real = s->cg_p.c[j].imag = 0.0;
		s->xxx.c[j].real = s->xxx.c[j].imag = 0.0;
		s->ttt.c[j].real = s->ttt.c[j].imag = 0.0;
		s->g_rand.c[j].real = s->g_rand.c[j].imag = 0.0;
	    }
	}
#endif
    }
    clear_latvec( F_OFFSET(xxx), EVENANDODD );
    dslash_site( F_OFFSET(g_rand), F_OFFSET(phi), parity);
    scalar_mult_latvec( F_OFFSET(phi), -1.0, F_OFFSET(phi), parity );
    scalar_mult_add_latvec( F_OFFSET(phi), F_OFFSET(g_rand), 2.0*mass,
	F_OFFSET(phi), parity );
}/* grsource */


/* Check congrad by multiplying xxx by M, compare result to g_rand */
/* Before calling checkmul() you should call grsource(EVENANDODD) and
   congrad(...,EVENANDODD) */
void checkmul() {
register int i,j;
register site *s;
    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), EVENANDODD);
    scalar_mult_add_latvec( F_OFFSET(ttt), F_OFFSET(xxx), 2.0*mass,
	F_OFFSET(ttt), EVENANDODD );
    FORALLSITESDOMAIN(i,s){
	printf("Site %d %d %d %d\n",s->x,s->y,s->z,s->t);
	for(j=0;j<3;j++){
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].real,
		(double)s->ttt.c[j].real,(double)s->g_rand.c[j].real - (double)s->ttt.c[j].real);
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].imag,
		(double)s->ttt.c[j].imag,(double)s->g_rand.c[j].imag - (double)s->ttt.c[j].imag);
	}
	printf("\n");
    }
}/* checkmul */
