/***************** u1Fmunu.c *************************************/
/* Calculate u1 Fmunu */
/* MIMD version 7 */
/* 11/03/16 by YL */
/* **************************************************************/
/*                                                              */
/*                                                              */
/* Code for average space-space and space-time Fmunu            */
/*                                                              */
/*                      b                                       */
/*                   o ->-- o                                   */
/*                   :      :                                   */
/*           ^     c ^      v a                                 */
/*           |       :      :                                   */
/*      dir1         o==<= o                                    */
/*                      d                                       */
/*                                                              */
/*            dir2->                                            */
/*                                                              */
/* Author: S. Basak                                             */
/* Last updated on 07.24.07                                     */
/* CD modified 5/24/12                                          */
/*                                                              */
/* ************************************************************ */

#include "generic_u1_includes.h"

void u1Fmunu(Real *sFmunu, Real *tFmunu)
{

    int i, dir1, dir2;
    double ssFmunu, stFmunu;
    double s2Fmunu, t2Fmunu;
    msg_tag *mtag0, *mtag1;
    Real pre;
    Real *Atmp = create_r_field();
    Real *Astpl = create_r_field();

    ssFmunu = stFmunu = 0.0;
    s2Fmunu = t2Fmunu = 0.0;

    /* Fmunu. in dir1-dir2 plane */
    for (dir1 = YUP; dir1 <= TUP; dir1++) {
        for (dir2 = XUP; dir2 < dir1; dir2++) {
            mtag0 = declare_strided_gather(u1_A + dir2, 4 * sizeof(Real),
                                           sizeof(Real), dir1,
                                           EVENANDODD, gen_pt[0]);
            prepare_gather(mtag0);
            do_gather(mtag0);
            /* = b */
            mtag1 = declare_strided_gather(u1_A + dir1, 4 * sizeof(Real),
                                           sizeof(Real), dir2,
                                           EVENANDODD, gen_pt[1]);
            prepare_gather(mtag1);
            do_gather(mtag1);
            /* = a */
            FORALLFIELDSITES(i) {
                Atmp[i] = u1_A[4 * i + dir1] - u1_A[4 * i + dir2];
            }  /* = (d^)c */

            wait_gather(mtag0);
            FORALLFIELDSITES(i) {
                Astpl[i] = Atmp[i] + *((Real *)gen_pt[0][i]);
            }  /* = (d^)cb */

            wait_gather(mtag1);
            FORALLFIELDSITES(i) {
                Atmp[i] = Astpl[i] - *(Real *)gen_pt[1][i];  /* = (d^)cb(a^) */
                /* check */
                /* Atmp[i] = 1.0; */
            }

            FORALLFIELDSITES(i) {
                pre = Atmp[i];
                if (dir1 == TUP) stFmunu += pre;
                else             ssFmunu += pre;
                if (dir1 == TUP) t2Fmunu += pre * pre;
                else             s2Fmunu += pre * pre;
            }

            cleanup_gather(mtag0);
            cleanup_gather(mtag1);

        } /* dir2-loop ends */
    } /* dir1-loop ends */

    destroy_r_field(Atmp);
    destroy_r_field(Astpl);

    g_doublesum(&s2Fmunu);
    g_doublesum(&t2Fmunu);
    /* sFmunu=ssFmunu/((Real)(nx*ny*nz*nt));
     * tFmunu=stFmunu/((Real)(nx*ny*nz*nt));
     */
    *sFmunu = s2Fmunu; /* this is Fmunu^2 */
    *tFmunu = t2Fmunu;

} /* end of Fmunu() */

/* ************************************************************    */
