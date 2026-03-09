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

void randomize_u1_A(Real *G, Real radius);
void gauge_trans_u1_A(Real *G, Real *tu1_A);

void rand_gauge_u1_A(Real *G, Real *tu1_A)
/* G holds the gauge transformation matrices */
{
    randomize_u1_A(G, 1.0);
    gauge_trans_u1_A(G, tu1_A);
}

void randomize_u1_A(Real *G, Real radius)
{
    register int i;
    site *s;

    FORALLSITES(i, s) {
        /* is this reall random? */
        /* G[i] = ce_itheta(radius*((Real)drand48()-0.5)); */
        G[i] = (int)rand() % 100;
    }
}

void gauge_trans_u1_A(Real *G, Real *tu1_A)
{
    register int i, mu;
    site *s;
    Real tmp;
    msg_tag *tag[4];

    FORALLUPDIR(mu)
    tag[mu] = start_gather_field(G, sizeof(Real), mu, EVENANDODD,
                                 gen_pt[mu]);

    FORALLUPDIR(mu) {
        wait_gather(tag[mu]);
        FORALLSITES(i, s) {
            node0_printf("random G i %d 4i+dir %d"
                         "A %e  G %e G+ %e\n",
                         i, 4 * i + mu,
                         tu1_A[4 * i + mu],
                         G[i],
                         *(Real *)gen_pt[mu][i]);

            tmp = tu1_A[4 * i + mu] - G[i];
            tu1_A[4 * i + mu] = tmp + *(Real *)gen_pt[mu][i];
        }
        cleanup_gather(tag[mu]);
    }
}
