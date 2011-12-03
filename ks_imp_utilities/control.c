/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the KS inversion and performs
   and/or compares with a standard file, the fermion force calculation */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#endif
#include "lattice_qdp.h"

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
  int prompt;
  char *filexml;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

  /* loop over input sets */
  while( readin(prompt) == 0){
    
    node0_printf("BEGIN\n");
#ifdef CHECK_INVERT
    
    check_ks_invert( srcfile, srcflag, F_OFFSET(phi), 
		     ansfile, ansflag, F_OFFSET(xxx), 
		     F_OFFSET(g_rand), mass);
    
#else
#ifndef HAVE_QIO
BOMB Checking the fermion force requires QIO compilation
#endif
    
    check_fermion_force( srcfile, srcflag, F_OFFSET(xxx), 
			 ansfile, ansflag, mass);
    node0_printf("Done checking fermion force\n");
#endif
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      node0_printf("Saving the lattice\n");
      save_lattice( saveflag, savefile, NULL );
      rephase( ON );
    }
    
#ifdef FN
    /* save longlinks if requested */
    if (savelongflag != FORGET ){
#ifdef HAVE_QIO
      filexml = create_QCDML();
      node0_printf("Saving the long links\n");
      save_color_matrix_scidac_from_field( savelongfile, filexml, 
         "Long links", QIO_SINGLEFILE, get_lnglinks(get_fm_links(fn_links)[0]), 4);
      free_QCDML(filexml);
#else
      printf("ERROR: Can't save the longlinks.  Recompile with QIO\n");
#endif
    }
    
    /* save fatlinks if requested */
    if (savefatflag != FORGET ){
#ifdef HAVE_QIO
      filexml = create_QCDML();
      node0_printf("Saving the fat links\n");
      save_color_matrix_scidac_from_field( savefatfile, filexml, 
	   "Fat links", QIO_SINGLEFILE, get_fatlinks(get_fm_links(fn_links)[0]), 4);
      free_QCDML(filexml);
#else
      printf("ERROR: Can't save the fatlinks.  Recompile with QIO\n");
#endif
    }
#endif
    
  }
  node0_printf("RUNNING COMPLETED\n");
  return 0;
}

