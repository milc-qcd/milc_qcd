/*********************** rand_gauge_u1.c ***************************/
/* MIMD version 7 */
/* original code by UMH */
/* 2/19/98 Version 5 port CD */
/* 11/03/16 by YL */

/* Makes a random gauge transformation on the U(1)
    gauge fields */
/* Warning KS fermion applications: Must be done with KS phases OUT! */
/* Requires field workspace G */
/* NOT TESTED YET */

#include "generic_includes.h"

void randomize_u1(complex *G, Real radius);
void gauge_trans_u1(complex *G, complex *tlink_u1);

void rand_gauge_u1(complex *G, complex *tlink_u1)
/* G holds the gauge transformation matrices */
{
    randomize_u1(G, 1.0);
    gauge_trans_u1(G, tlink_u1);
}

void randomize_u1(complex *G, Real radius)
{
    register int i;
    site *s;

    FORALLSITES(i, s) {
        /* is this really random? */
        /* G[i] = ce_itheta(radius*((Real)drand48()-0.5)); */
        G[i] = ce_itheta(radius * ((Real)rand() - 0.5));
    }
}

void gauge_trans_u1(complex *G, complex *tlink_u1)
{
    register int i, mu;
    site *s;
    complex tmp;
    msg_tag *tag[4];

    FORALLUPDIR(mu)
    tag[mu] = start_gather_field(G, sizeof(complex), mu, EVENANDODD,
                                 gen_pt[mu]);

    FORALLUPDIR(mu) {
        wait_gather(tag[mu]);
        FORALLSITES(i, s) {
            node0_printf("random G i %d 4i+dir %d"
                         "real link_u1 %e"
                         "imag link_u1 %e"
                         "realG %e imagG %e\n",
                         i, 4 * i + mu,
                         tlink_u1[i].real,
                         tlink_u1[i].imag,
                         G[i].real, G[i].imag);

            CMULJ_(G[i], tlink_u1[i], tmp);
            CMUL(tmp, *(complex *)gen_pt[mu][i], tlink_u1[i]);
        }
        cleanup_gather(tag[mu]);
    }
}
