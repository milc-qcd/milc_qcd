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
void finalize_grid();
int grid_initialized(void);
void setup_grid_communicator(int peGrid[]);
int *query_grid_node_mapping(void);
int grid_lexicographic_to_worldrank(int lexrank);
int grid_rank_from_processor_coor(int x, int y, int z, int t);
void grid_coor_from_processor_rank(int coords[], int worldrank);

#endif // _GENERIC_GRID_H
