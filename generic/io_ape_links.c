/********************** io_ape_links.c **********************************/
/* MILC/QIO interface  -- reads and writes an APE links file */
/* Requires QIO */
/*   MIMD version 7 */
/* CD 8/2024
*/

#include "generic_includes.h"
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include <string.h>
#define LATDIM 4

gauge_file *
save_apelinks( int flag, su3_matrix *links, const char *filename){
  double dtime;
  gauge_file *gf = NULL;

  dtime = -dclock();
  switch( flag ){
  case FORGET:
    gf = NULL;
    break;
  case SAVE_SERIAL_SCIDAC:
    gf = save_serial_scidac(links, filename, 1);
    break;
  case SAVE_SERIAL_SCIDAC_DP:
    gf = save_serial_scidac(links, filename, 2);
    break;
  case SAVE_PARALLEL_SCIDAC:
    gf = save_parallel_scidac(links, filename, 1);
    break;
  case SAVE_PARALLEL_SCIDAC_DP:
    gf = save_parallel_scidac(links, filename, 2);
    break;
  case SAVE_PARTFILE_SCIDAC:
    gf = save_partfile_scidac(links, filename, 1);
    break;
  case SAVE_PARTFILE_SCIDAC_DP:
    gf = save_partfile_scidac(links, filename, 2);
    break;
  default:
    node0_printf("save_apelinks: Unrecognized save flag %d\n", flag);
    terminate(1);
  }
  
  dtime += dclock();
  if(flag != FORGET)
    node0_printf("Time to save APE links = %e\n", dtime);
  
  return gf;
}

gauge_file *
reload_apelinks( int flag, su3_matrix *links, const char *filename){

  gauge_file *gf = NULL;

  double dtime = -dclock();
  
  switch(flag){
  case FRESH:	/* Don't read anything */
    break;
  case RELOAD_SERIAL:	/* read APE links serially */
    gf = restore_serial_scidac(links, filename);
    break;
  case RELOAD_PARALLEL: /* read APE links in parallel */
    gf = restore_parallel_scidac(links, filename);
    break;
  default:
   node0_printf("reload_apelinks: Unrecognized reload flag %d\n", flag);
   terminate(1);
  }
  dtime += dclock();
  if(flag != FRESH)
    node0_printf("Time to reload gauge configuration = %e\n",dtime);
  
  return gf;
}

