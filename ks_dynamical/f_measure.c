/**************** f_measure.c ***************************************/
/* MIMD version 6 */

/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermionic energy and pressure
    fermion action

    Before calling this routine, a Gaussian random vector should be
    placed in g_rand, M_adjoint*g_rand in phi, and M_inverse*g_rand in xxx.
    phi and xxx are defined on even sites only.
*/


#include "ks_dyn_includes.h"

void f_measure(Real *r_psi_bar_psi_even, Real *r_psi_bar_psi_odd, 
	       Real *r_ferm_energy, Real *r_ferm_pressure, 
	       Real *r_ferm_action)
{
/* local variables for accumulators */
register int i,dir;
register site *st;
msg_tag *tag0,*tag1;
Real rpbp_e,rpbp_o,rfenergy,rfpressure,rfaction;
complex cc;

    rpbp_e = rpbp_o = rfenergy = rfpressure = rfaction = 0.0;

    /* fermion action = phi.xxx */
    /* psi-bar-psi on even sites = g_rand.xxx */
    FOREVENSITES(i,st){
        cc = su3_dot( &(st->phi), &(st->xxx) );
	rfaction += cc.real;
        cc = su3_dot( &(st->g_rand), &(st->xxx) );
	rpbp_e += cc.real;
    }
    /* psi-bar-psi on odd sites (energy and pressure will be added later) */
    FORODDSITES(i,st){
	cc = su3_dot( &(st->g_rand), &(st->g_rand) );
	rpbp_o += cc.real;
    }

    /* fermion energy and pressure */
    for(dir=XUP;dir<=TUP;dir++){
        /* Bring g_rand down to even site */
        tag1 = start_gather_site( F_OFFSET(g_rand), sizeof(su3_vector), dir,
	    EVEN, gen_pt[1] );
        FORODDSITES(i,st){
	    mult_adj_su3_mat_vec( &(st->link[dir]), &(st->g_rand),
		&(st->tempvec[0]));
        }
        /* Bring result up to even site */
        tag0 = start_gather_site( F_OFFSET(tempvec[0]), sizeof(su3_vector),
	    OPP_DIR(dir), EVEN, gen_pt[0] );
        wait_gather(tag0);
        wait_gather(tag1);
	/* now gen_pt[1] points to g_rand at the site in the plus direction, 
	   and gen_pt[0] to g_rand at the site in the minus direction
	   multiplied by the link matrix.  It is OK to overwrite *gen_pt[0]
	   because it points to tempvec[0] on some site. */
        FOREVENSITES(i,st){
	    mult_su3_mat_vec_nsum( &(st->link[dir]), (su3_vector *)gen_pt[1][i],
		(su3_vector *)gen_pt[0][i] );
	    cc = su3_dot( (su3_vector *)gen_pt[0][i], &(st->xxx) );
	    if(dir==TUP)rfenergy += cc.real;
	    else      rfpressure += cc.real;

        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
    }
    g_floatsum( &rpbp_o );
    g_floatsum( &rpbp_e );
    g_floatsum( &rfenergy );
    g_floatsum( &rfpressure );
    g_floatsum( &rfaction );

    *r_psi_bar_psi_odd =  (rpbp_o-rfenergy-rfpressure) *
				(2.0/((double)volume*2.0*mass));
    *r_psi_bar_psi_even =  rpbp_e*(2.0/(double)volume) ;
    *r_ferm_energy =  rfenergy*((double)nflavors*0.5/(double)volume) ;
    *r_ferm_pressure =  rfpressure*(-(double)nflavors/(6.0*(double)volume)) ;
    *r_ferm_action =  rfaction*(1.0/(double)volume) ;
}
