/************************* control.c *******************************/
/* MIMD version 6 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the KS inversion and performs
   and/or compares with a standard file, the fermion force calculation */

#define CONTROL
#include "ks_imp_includes_qdp.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#endif

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
  int i;
  int prompt;
  char *filexml;
  double dtime, dclock();
  
  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  initialize_machine(argc,argv);
  QDP_initialize(&argc, &argv);
  g_sync();
  /* set up */
  prompt = setup();
  /* loop over input sets */
  while( readin(prompt) == 0){
    
#ifdef CHECK_INVERT
    
    check_ks_invert( srcfile, srcflag, F_OFFSET(phi), 
		     ansfile, ansflag, F_OFFSET(xxx), 
		     F_OFFSET(g_rand), mass);
    
#else
    
    check_fermion_force( srcfile, srcflag, F_OFFSET(xxx), 
			 ansfile, ansflag, mass);
#endif
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile );
      rephase( ON );
    }
    
#ifdef FN
    /* save longlinks if requested */
    if (savelongflag != FORGET ){
      filexml = create_ks_XML();
#ifdef DSLASH_TMP_LINKS
      save_color_matrix_scidac_from_temp( savelongfile, filexml, 
			  "Long links", QIO_SINGLEFILE, t_longlink, 4);
#else
      save_color_matrix_scidac( savelongfile, filexml, 
		"Long links", QIO_SINGLEFILE, F_OFFSET(longlink[0]), 4);
#endif
      free_ks_XML(filexml);
    }
    
    /* save fatlinks if requested */
    if (savefatflag != FORGET ){
      filexml = create_ks_XML();
#ifdef DSLASH_TMP_LINKS
      save_color_matrix_scidac_from_temp( savefatfile, filexml, 
		  "Fat links", QIO_SINGLEFILE, t_fatlink, 4);
#else
      save_color_matrix_scidac( savefatfile, filexml, 
		"Fat links", QIO_SINGLEFILE, F_OFFSET(fatlink[0]), 4);
#endif
      free_ks_XML(filexml);
    }
#endif
    
  }
  return 0;
}
