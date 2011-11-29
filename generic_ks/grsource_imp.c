/************************ grsource_imp.c *****************************/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

#if 0
/* OBSOLETE VERSION */

void grsource_imp_old( field_offset dest, Real mass, int parity,
		       ferm_links_t *fn ) {
register int i,j;
register site *s;
    FORALLSITES(i,s){
        for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
    }
    dslash_site( F_OFFSET(g_rand), dest, parity, fn);
    scalar_mult_latvec( dest, -1.0, dest, parity );
    scalar_mult_add_latvec( dest, F_OFFSET(g_rand), 2.0*mass,
	dest, parity );
}/* grsource_imp_old */
#endif

/* This version uses a temporary field vector */

/* DEPRECATED VERSION */

void grsource_imp( field_offset dest, Real mass, int parity,
		   imp_ferm_links_t *fn ) {
  su3_vector *g_rand = create_v_field();
  su3_vector *tdest = create_v_field();
  
  grsource_plain_field( g_rand, EVENANDODD );
  ks_dirac_adj_op( g_rand, tdest, mass, parity, fn );

  copy_site_member_from_v_field(dest, tdest);
  
  destroy_v_field(tdest);
  destroy_v_field(g_rand);

}/* grsource_imp */


/* This version works entirely with field vectors */

void grsource_imp_field( su3_vector *dest, Real mass, int parity,
			 imp_ferm_links_t *fn ) {
  su3_vector *g_rand = create_v_field();

  grsource_plain_field( g_rand, EVENANDODD);
  ks_dirac_adj_op( g_rand, dest, mass, parity, fn );

  destroy_v_field(g_rand);

}/* grsource_imp_field */


#if 0
/* OBSOLETE VERSION */

void z2rsource_imp_old( field_offset dest, Real mass, int parity,
			imp_ferm_links_t *fn ) {
register int i,j;
register site *s;
    FORALLSITES(i,s){
        for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.c[j].real = z2_rand_no(&(s->site_prn));
            s->g_rand.c[j].imag = z2_rand_no(&(s->site_prn));
#else
            s->g_rand.c[j].real = z2_rand_no(&node_prn);
            s->g_rand.c[j].imag = z2_rand_no(&node_prn);
#endif
        }
    }
    dslash_site( F_OFFSET(g_rand), dest, parity, fn);
    scalar_mult_latvec( dest, -1.0, dest, parity );
    scalar_mult_add_latvec( dest, F_OFFSET(g_rand), 2.0*mass,
	dest, parity );
}/* z2rsource_imp_old */
#endif


/* DEPRECATED VERSION */

void z2rsource_imp( field_offset dest, Real mass, int parity,
		    imp_ferm_links_t *fn ) {
  su3_vector *z2_rand = create_v_field();
  su3_vector *tdest = create_v_field();

  z2rsource_plain_field( z2_rand, EVENANDODD );
  ks_dirac_adj_op( z2_rand, tdest, mass, parity, fn );

  copy_site_member_from_v_field(dest, tdest);

  destroy_v_field(tdest);
  destroy_v_field(z2_rand);

}/* z2rsource_imp */


/* construct a plain gaussian random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

/* DEPRECATED VERSION */

void grsource_plain( field_offset dest, int parity ) {
  int i,j;
  site *s;
  su3_vector *rand;
  FORSOMEPARITY(i,s,parity){
    for(j=0;j<3;j++){
#ifdef SITERAND
      rand = (su3_vector *)F_PT(s,dest);
            rand->c[j].real = gaussian_rand_no(&(s->site_prn));
            rand->c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            rand->c[j].real = gaussian_rand_no(&node_prn);
            rand->c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }    
}/* grsource_plain */


/* construct a plain gaussian random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

void grsource_plain_field( su3_vector *dest, int parity ) {
  int i,j;
  site *s;

  FORSOMEPARITY(i,s,parity){
    for(j=0;j<3;j++){
#ifdef SITERAND
      dest[i].c[j].real = gaussian_rand_no(&(s->site_prn));
      dest[i].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      dest[i].c[j].real = gaussian_rand_no(&node_prn);
      dest[i].c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }    
}/* grsource_plain_field */


/* construct a plain Z(2) random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

/* DEPRECATED VERSION */

void z2rsource_plain( field_offset dest, int parity ) {
  int i,j;
  site *s;
  su3_vector *rand;
  FORSOMEPARITY(i,s,parity){
    for(j=0;j<3;j++){
#ifdef SITERAND
      rand = (su3_vector *)F_PT(s,dest);
            rand->c[j].real = z2_rand_no(&(s->site_prn));
            rand->c[j].imag = z2_rand_no(&(s->site_prn));
#else
            rand->c[j].real = z2_rand_no(&node_prn);
            rand->c[j].imag = z2_rand_no(&node_prn);
#endif
    }
  }    
}/* z2rsource_plain */


/* construct a plain Z(2) random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

void z2rsource_plain_field( su3_vector *dest, int parity ) {
  int i,j;
  site *s;

  FORSOMEPARITY(i,s,parity){
    for(j=0;j<3;j++){
#ifdef SITERAND
      dest[i].c[j].real = z2_rand_no(&(s->site_prn));
      dest[i].c[j].imag = z2_rand_no(&(s->site_prn));
#else
      dest[i].c[j].real = z2_rand_no(&node_prn);
      dest[i].c[j].imag = z2_rand_no(&node_prn);
#endif
    }
  }    
}/* z2rsource_plain */


/* Check congrad by multiplying src by M, compare result to g_rand */
/* Before calling checkmul() you should call grsource(EVENANDODD) and
   congrad(...,EVENANDODD) */

/* DEPRECATED VERSION */

#if 0
void checkmul_imp( field_offset src, Real mass,
		   imp_ferm_links_t *fn )
{
register int i,j;
register site *s;
    dslash_site( src, F_OFFSET(ttt), EVENANDODD, fn);
    scalar_mult_add_latvec( F_OFFSET(ttt), src, 2.0*mass,
	F_OFFSET(ttt), EVENANDODD );
    FORALLSITES(i,s){
	printf("Site %d %d %d %d\n",s->x,s->y,s->z,s->t);
	for(j=0;j<3;j++){
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].real,
		(double)s->ttt.c[j].real,(double)s->g_rand.c[j].real - (double)s->ttt.c[j].real);
	    printf("%d %d\t%e\t%e\t%e\n",i,j,(double)s->g_rand.c[j].imag,
		(double)s->ttt.c[j].imag,(double)s->g_rand.c[j].imag - (double)s->ttt.c[j].imag);
	}
	printf("\n");
/**sleep(2);**/
    }
}/* checkmul */

#endif
