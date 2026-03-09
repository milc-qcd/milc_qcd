/***************** gauge_force_u1_cpu.c *************************************/
/* Calculate u1 gauge force */
/* MIMD version 7 */
/* 11/03/16 by YL */
/* **************************************************************/
/*                                                              */
/*                                                              */
/* The non-compact U1 gauge action is                           */
/* S_A = 1/(4*e^2) \sum_{n,\mu,\nu} F(n)_{\mu\nu} F(n)^{\mu\nu} */
/*     = 1/(4*e^2) \sum_{n,\mu,\nu} F(n)_{\mu\nu} F(n)^{\mu\nu} */
/*     = 1/(4*e^2) \sum_{n,\mu,\nu} [A(n)_\mu + A(n+\mu)_\nu    */
/*                            - A(n+\nu)_\mu - A(n)_\nu]^2      */
/* 1/e^2 = beta_u1                                              */
/* dS_A/dA(n)_\mu = 1/(e^2) \sum_\nu[                           */
/*  (A(n)_\mu + A(n+\mu)_\nu - A(n+\nu)_\mu - A(n)_\nu)         */
/* -(A(n-\nu)_\mu + A(n-\nu+\mu)_\nu - A(n)_\mu - A(n-\nu)_\nu)]*/
/*   = beta_u1 * \sum_\nu [ F(n)_{\mu\nu} - F(n-\nu)_{\mu\nu} ] */
/* The conjugate momentum h is defined as h = dA/dt.            */
/* The U(1) gauge force is then                                 */
/* dh/dt = - dS_A/dA                                            */
/*       =beta_u1 * \sum_\nu[ -(d+a-b-c)-(d-e-g+f) ]            */
/*  The horizontal portial of the force at site n               */
/*  pointing to horizontal dir1 direction is                    */
/*          b                                                   */
/*      o --<-- o                                               */
/*      :       :                                               */
/*    c v       ^ a                                             */
/*      :       :                                               */
/*      o -->-- o                                               */
/* site n   d                                                   */
/*      o -->---o                                               */
/*      :       :                                               */
/*    f ^       ve                                              */
/*      :       :                                               */
/*      o --<-- o                                               */
/*          g                                                   */
/*      :                                                       */
/*      ^                                                       */
/*      :                                                       */
/*  dir2 \nu                                                    */
/*  dir1 \mu -->                                                */
/* Author: YL                                                   */
/* Last updated on 11/03/2016                                   */
/* **************************************************************/

#include "generic_u1_includes.h"

void gauge_force_u1_cpu(Real eps)
{
    register int i, dir1, dir2;
    register site *s;
    msg_tag *mtag0, *mtag1, *mtag2;
    int start;
    Real t1, t2;
    register Real eb;
    Real *Astpl;
    Real *temp1;

#ifdef GFTIME
    /* For Symanzik1 action */
    /* int nflop = 153004;  */
    double dtime = -dclock();
#endif

    Astpl = create_r_field();
    temp1 = create_r_field();
    eb = eps * beta_u1;
    eb = (-1.0 / 3.0 / pseudo_charges[0]) * eps * beta_u1;

    for (dir1 = XUP; dir1 <= TUP; dir1++) {
        start = 1;
        for (dir2 = XUP; dir2 <= TUP; dir2++)if (dir2 != dir1) {
            /*  lower staple */
            mtag0 = declare_strided_gather(u1_A + dir2,
                                           4 * sizeof(Real), sizeof(Real),
                                           dir1, EVENANDODD, gen_pt[0]);
            prepare_gather(mtag0);
            do_gather(mtag0);
            /* = a */
            /* = A(n+\nu)_\mu */
            mtag2 = declare_strided_gather(u1_A + dir1,
                                           4 * sizeof(Real), sizeof(Real),
                                           dir2, EVENANDODD, gen_pt[2]);
            prepare_gather(mtag2);
            do_gather(mtag2);
            /* = b */
            /* = A(n+\mu)_\nu */
            wait_gather(mtag0);
            FORALLFIELDSITES(i) {
                t1 = u1_A[4 * i + dir2] - u1_A[4 * i + dir1];
                temp1[i] = t1 - (*(Real *)gen_pt[0][i]);
            }        /* = c(d^)(a^) */
            /* Gather lower staple "up to home site" */
            mtag1 = start_gather_field(temp1,
                                       sizeof(Real), OPP_DIR(dir2),
                                       EVENANDODD, gen_pt[1]);
            /* c(d^)(a^) becomes f(g^)(e^) */

            /*  upper staple */
            /* mtag0 has already been gathered */
            wait_gather(mtag2);
            if (start) {
                FORALLFIELDSITES(i) {
                    t1 = -(*(Real *)gen_pt[2][i]) - u1_A[4 * i + dir2];
                    Astpl[i] = t1 + (*(Real *)gen_pt[0][i]); /* = (c^)(b^)a */
                }
                start = 0;
            } else {
                FORALLFIELDSITES(i) {
                    t1 = -(*(Real *)gen_pt[2][i]) - u1_A[4 * i + dir2];
                    t2 = t1 + (*(Real *)gen_pt[0][i]);  /* = (c^)(b^)a */
                    Astpl[i] += t2;
                }
            }
            wait_gather(mtag1);
            FORALLFIELDSITES(i) {
                Astpl[i] += (*(Real *)gen_pt[1][i]); /* = (c^)(b^)a+f(g^)(e^) */
            }
            cleanup_gather(mtag0);
            cleanup_gather(mtag1);
            cleanup_gather(mtag2);
        }; /* dir2 */

        /* Now complete the two plaquettes and update the momentum */

        FORALLSITES(i, s) {
            t1 = 6.0 * u1_A[4 * i + dir1] + Astpl[i];
            s->mom_u1[dir1] -= eb * t1;

#ifdef U1_DEBUG
            if (dir1 == XUP && i == 0) {
                node0_printf("../generic_u1/gauge_force_u1_cpu.c eps, "
                             "beta_u1, eb, t1 %e %e %e %e\n",
                             eps, beta_u1, eb, t1);
            }
#endif

        }
    } /* dir1 */
#ifdef GFTIME
    dtime += dclock();
    node0_printf("GFTIME:   time = %e (U1) mflops = %e\n", dtime, 0.0);
#endif

    destroy_r_field(Astpl);
    destroy_r_field(temp1);
    Astpl = NULL;
    temp1 = NULL;
} /* gauge_force_u1_cpu */
