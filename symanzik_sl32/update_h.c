/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 7 */

/* Modifications:
   2/17/98  ANSI prototyping,
            improved gauge updating taken from "ks_imp_dyn" U.M.H.
   */

#include "symanzik_sl32_includes.h"

void update_h(Real eps) {
    /* gauge field force */
    imp_gauge_force(eps,F_OFFSET(mom));
} /* update_h */
