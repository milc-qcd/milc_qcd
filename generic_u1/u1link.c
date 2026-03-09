/***************** u1link.c *****************************************/

/* Utilities inserting the U(1) field into the SU(3) gauge links */

/* MIMD version 7 */

/* Author: S. Basak ? 09.18.07 */
/* CD modified 5/24/12 */

#include "generic_u1_includes.h"

/*
* can use create_r_array_field(4) and
* create_c_array_field(4) to replace
* the following create_u1_A_field()
* and create_u1_link_field() routines.
* By YL 08/24/2016.
*/
/* ************************ */
Real *create_u1_A_field(void) {
    Real *A;
    A = (Real *)malloc(4 * sizeof(Real) * sites_on_node);
    if (A == NULL) {
        node0_printf("create_u1_A_field: No room\n");
        terminate(1);
    }
    return A;
}

void destroy_u1_A_field(Real *A) {
    if (A != NULL)
        free(A);
}

#ifndef NO_GAUGE_FIELD
static su3_matrix **tlink;

/* ************************ */
void u1phase_on(Real charge, Real *A) {

    int i, dir;
    site *s;
    complex clink;
    Real theta;

    /* save a copy of the original field */
    tlink = gauge_field_copy_site_to_field(F_OFFSET(link[0]));

    /* Insert the U(1) phase */
    FORALLSITES(i, s) {
        FORALLUPDIR(dir) {
            theta = charge * A[4 * i + dir];  /* Qe * A(mu, x) */
            clink = ce_itheta(theta);
            c_scalar_mult_su3mat( &tlink[dir][i], &clink, &s->link[dir]);
        }
    }
} /* u1_phase_on */

/* ************************ */
void u1phase_off(void) {
    int dir;

    gauge_field_copy_field_to_site(tlink, F_OFFSET(link[0]));
    FORALLUPDIR(dir) {
        free(tlink[dir]);
    }

    free(tlink);

} /* u1phase_off */

#endif

/* ************************ */
void construct_u1_links_from_A(Real charge, Real *A, complex *clink) {
    int i, dir;
    site *s;
    Real theta;

    FORALLSITES(i, s) {
        FORALLUPDIR(dir) {
            theta = charge * A[4 * i + dir];  /* Qe * A(mu, x) */
            clink[4 * i + dir] = ce_itheta(theta);
        }
    }
} /* construct_u1_links_from_A */

/* copy a gauge field - an array of four complex numbers */
void u1_link_copy_field_to_field(complex *src, complex *dest) {
    register int i, dir;
    register site *s;
    FORALLSITES(i, s) {
        for (dir = XUP; dir <= TUP; dir++) {
            dest[4 * i + dir] = src[4 * i + dir];
        }
    }
#ifdef FN
/* free_fn_links_u1(fn_links_u1); */
#endif
}

/* copy a gauge potential - an array of four real numbers */
void u1_A_copy_field_to_field(Real *src, Real *dest) {
    register int i, dir;
    register site *s;
    FORALLSITES(i, s) {
        for (dir = XUP; dir <= TUP; dir++) {
            dest[4 * i + dir] = src[4 * i + dir];
        }
    }
#ifdef FN
/*   free_fn_links_u1(&fn_links); */
#endif
}
/* ************************************************************ */
