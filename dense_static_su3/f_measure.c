/**************** f_measure.c ***************************************/
/* MIMD version 4 */
/* version for dense_su3 - everything is complex */

/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermionic energy and pressure
    fermion action

    Before calling this routine, a Gaussian random vector should be
    placed in g_rand, M_adjoint*g_rand in phi, and M_inverse*g_rand in xxx.
    phi and xxx are defined on even sites only.
*/

#include "su3_dense_includes.h"

void f_measure() {
/* local variables for accumulators */
register int i,dir;
register site *st;
msg_tag *tag0,*tag1;
complex pbp_e,pbp_o,fenergy,fpressure,faction;
complex cc;

    pbp_e = pbp_o = fenergy = fpressure = faction = cmplx(0.0,0.0);

    /* fermion action = phi.xxx */
    /* psi-bar-psi on even sites = g_rand.xxx */
    FOREVENSITES(i,st){
        cc = su3_dot( &(st->phi), &(st->xxx) );
	CSUM(faction,cc);
        cc = su3_dot( &(st->g_rand), &(st->xxx) );
	CSUM(pbp_e,cc);
    }
    /* psi-bar-psi on odd sites (energy and pressure will be added later) */
    FORODDSITES(i,st){
	cc = su3_dot( &(st->g_rand), &(st->g_rand) );
	CSUM(pbp_o,cc);
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
	    if(dir==TUP){CSUM(fenergy,cc)}
	    else        {CSUM(fpressure,cc)}
        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
    }
    g_complexsum( &pbp_o );
    g_complexsum( &pbp_e );
    g_complexsum( &fenergy );
    g_complexsum( &fpressure );
    g_complexsum( &faction );

    CSUB(pbp_o,fenergy,pbp_o);
    CSUB(pbp_o,fpressure,pbp_o);
    CMULREAL(pbp_o,2.0/((double)volume*2.0*mass),pbp_o);
    CMULREAL(pbp_e,2.0/(double)volume,pbp_e);
    CMULREAL(fenergy,(double)nflavors*0.5/(double)volume,fenergy);
    CMULREAL(fpressure,-(double)nflavors/(6.0*(double)volume),
	fpressure);
    CMULREAL(faction,1.0/(double)volume,faction);
    if(this_node==0)printf("ZFMES %e %e %e %e %e %e\n",
	(double)pbp_e.real, (double)pbp_e.imag,
	(double)pbp_o.real, (double)pbp_o.imag,
	(double)fenergy.real, (double)fenergy.imag );
}
