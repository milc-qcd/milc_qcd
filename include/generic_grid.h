/* Grid-MILC Interface */
#ifndef _GENERIC_GRID_H
#define _GENERIC_GRID_H

#include "../include/mGrid/mGrid.h"
#include "../include/config.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/su3.h"
#include "../include/dirs.h"
#include "../include/fermion_links.h"
#include <stdbool.h>

/* milc_to_grid_utilities.c */
GRID_evenodd_t milc2grid_parity(int milc_parity);
int grid2milc_parity(GRID_evenodd_t grid_parity);
GRID_status_t initialize_grid(void);
void finalize_grid(void);
int *query_grid_node_mapping(int peGrid[]);

#endif // _GENERIC_GRID_H
