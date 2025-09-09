/********************** check_fermion_force.c *******************************/
/* MIMD version 7 */

/* This code performs and/or checks the link fattening calculation */
/* It compares standard MILC link fattening with either a trial input
   file or, if none is provided, with results returned by a 
   non-MILC subroutine call. */

/* Instructions for use ...
   1. Edit the test quark action file ../generic_ks/imp_actions/hisq/hisq_u3_test_action.h
   to select the desired combination of link components to include in the test.
   2. Build the check_link_fattening code with the MILC version of link fattening, based
   on the test quark action above.
   3. Run the code, specifying that the input fat and long links are "fresh" and 
   the computed fat and long links are to be saved in files.
   4. Rebuild the code, now linking to the proposed link fattening code (e.g. by setting
   WANTGRID = true.
   5. Run the code, this time reading the above links from the file and specifying
   "forget" for the output files.
   The compares the resulting fat and long links and reports differences.
*/   
   

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#if defined(HAVE_QIO) || defined(HAVE_GRID)
#include <qio.h>
#else
#error Requires HAVE_QIO or HAVE_GRID
#endif

void check_link_fattening( char *lngansfile, int lngansflag, char *fatansfile, int fatansflag )
{
#if (MILC_PRECISION == 1)
  Real tol = 1e-6;
#else
  Real tol = 1e-14;
#endif

  ks_action_paths_hisq *ap = get_action_paths_hisq(&fn_links[0]);

  /* Get fat and long link pointers for links generated in main by the
     code specified in the Makefile. e.g. WANTGRID=true generates them
     with GRID */
  imp_ferm_links_t *fn = get_fm_links(fn_links, 0);
  su3_matrix *lng = get_lnglinks(fn);
  su3_matrix *fat = get_fatlinks(fn);
  
  /*------------------------------------------------------------*/
  /* Read fiducial result if given */

  su3_matrix *fattrust;
  su3_matrix *lngtrust;
  
  if(lngansflag != FRESH && fatansflag != FRESH){

    /* If the answer file is given, read it for comparison */
    node0_printf("Reading the fat and long links from files\n"); fflush(stdout);
    lngtrust = (su3_matrix *)malloc(4*sites_on_node*sizeof(su3_matrix));
    if(lngtrust == NULL){
      node0_printf("No room for lngtrust\n");
      terminate(1);
    }
    memset(lngtrust, '\0', 4*sites_on_node*sizeof(su3_matrix));
    
    fattrust = (su3_matrix *)malloc(4*sites_on_node*sizeof(su3_matrix));
    if(fattrust == NULL){
      node0_printf("No room for fattrust\n");
      terminate(1);
    }
    memset(fattrust, '\0', 4*sites_on_node*sizeof(su3_matrix));

    restore_color_matrix_scidac_to_field(lngansfile, lngtrust, 4, MILC_PRECISION);
    restore_color_matrix_scidac_to_field(fatansfile, fattrust, 4, MILC_PRECISION);

  } else {

    fattrust = fat;
    lngtrust = lng;

  }

  /*------------------------------------------------------------*/
  if(lngansflag != FRESH){

    node0_printf("Checking the long-link answer\n"); fflush(stdout);

    Real maxdiff = 0;
    Real maxnorm = 0;
    int i, dir;

    FORALLFIELDSITES(i){
      FORALLUPDIR(dir){
	su3_matrix diffmat;
	su3_matrix *tmat = lng + 4*i + dir;
	sub_su3_matrix( lngtrust + 4*i + dir, tmat, &diffmat);
	Real diff = sqrt(realtrace_su3( &diffmat, &diffmat ));
	Real norm = sqrt(realtrace_su3( tmat, tmat));
	if(diff > tol * norm){
	  printf("Large relative difference %e node %d coord %d %d %d %d dir %d\n",
		 diff/norm,this_node,
		 lattice[i].x,lattice[i].y,lattice[i].z,lattice[i].t,dir);
	  
	  printf("From file\n");
	  dumpmat(lngtrust + 4*i + dir);
	  printf("From link calculation\n");
	  dumpmat(tmat);
	}
	if(maxdiff < diff)maxdiff = diff;
	if(maxnorm < norm)maxnorm = norm;
      }
    }
    
    g_floatmax(&maxdiff);
    g_floatmax(&maxnorm);
    if(maxnorm > 0){
      Real reldiff = maxdiff/maxnorm;
      node0_printf("Relative difference %e\n",reldiff);
    }
  }

  /*------------------------------------------------------------*/
  if(fatansflag != FRESH){

    node0_printf("Checking the fat-link answer\n"); fflush(stdout);
    
    Real maxdiff = 0;
    Real maxnorm = 0;
    
    int i, dir;
    FORALLFIELDSITES(i){
      FORALLUPDIR(dir){
	su3_matrix diffmat;
	su3_matrix *tmat = fat + 4*i + dir;
	sub_su3_matrix( fattrust + 4*i + dir, tmat, &diffmat);
	Real diff = sqrt(realtrace_su3( &diffmat, &diffmat ));
	Real norm = sqrt(realtrace_su3( tmat, tmat));
	if(diff > tol * norm){
	  printf("Large relative difference %e node %d coord %d %d %d %d dir %d\n",
		 diff/norm,this_node,
		 lattice[i].x,lattice[i].y,lattice[i].z,lattice[i].t,dir);
	  printf("From file\n");
	  dumpmat(fattrust + 4*i + dir);
	  printf("From link calculation\n");
	  dumpmat(tmat);
	}
	if(maxdiff < diff)maxdiff = diff;
	if(maxnorm < norm)maxnorm = norm;
      }
    }
    
    g_floatmax(&maxdiff);
    g_floatmax(&maxnorm);
    if(maxnorm > 0){
      Real reldiff = maxdiff/maxnorm;
      node0_printf("Relative difference %e\n",reldiff);
    } else {
      node0_printf("Absolute difference %e but maxnorm is 0???\n",maxdiff);
    }      
  }
  
  /*------------------------------------------------------------*/
  /* save longlinks if requested */
  if (savelongflag != FORGET ){
#ifdef HAVE_QIO
    su3_matrix *lng = get_lnglinks(fn);
    char *filexml = create_QCDML();
    if(!param.withKSphases){
      rephase_field_offset( lng, OFF, NULL, param.coord_origin);
      node0_printf("Saving the long links with KS phases OUT buttime BC IN\n");
    } else {
      node0_printf("Saving the long links with KS phases IN and time BC IN\n");
    }
    node0_printf("Saving the long links with LFN\n '%s'\n", stringLFNlong);
    save_color_matrix_scidac_from_field( savelongfile, filexml, 
					 "Long links", QIO_SINGLEFILE, lng, 4, MILC_PRECISION,
					 stringLFNlong);
    free_QCDML(filexml);
    if(!param.withKSphases)
      rephase_field_offset( lng, ON, NULL, param.coord_origin);
#else
    printf("ERROR: Can't save the longlinks.  Recompile with QIO\n");
#endif
  }
  
  /*------------------------------------------------------------*/
  /* save fat links if requested */
  if (savefatflag != FORGET ){
#ifdef HAVE_QIO
    su3_matrix *fat = get_fatlinks(fn);
    char *filexml = create_QCDML();
    if(!param.withKSphases){
      rephase_field_offset( fat, OFF, NULL, param.coord_origin);
      node0_printf("Saving the fat links with KS phases OUT but time BC IN\n");
    } else {
      node0_printf("Saving the fat links with KS phases IN and time BC IN\n");
    }
    node0_printf("Saving the fat links with LFN\n '%s'\n", stringLFNfat);
    save_color_matrix_scidac_from_field( savefatfile, filexml, 
					 "Fat links", QIO_SINGLEFILE, fat, 4, MILC_PRECISION,
					 stringLFNfat);
    free_QCDML(filexml);
    if(!param.withKSphases)
      rephase_field_offset( fat, ON, NULL, param.coord_origin);
#else
    printf("ERROR: Can't save the fatlinks.  Recompile with QIO\n");
#endif

  }

  destroy_fn_links(fn);
}      
