/******************* gauge_force_imp_qphix.c ************************/
/* For the QPhiX interface */
/* MIMD version 7 */

#include "../include/generic_qphix.h"
#include "../include/generic_ks_qphix.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*! \brief call to the qphix_ks_congrad_parity.
 *
 * Choose the inversion precision
 */
void
imp_gauge_force_qphix ( Real eps, field_offset mom_off )
{

  if(MILC_PRECISION == 1)
    imp_gauge_force_qphix_F( eps, mom_off);
  else
    imp_gauge_force_qphix_D( eps, mom_off);
}
