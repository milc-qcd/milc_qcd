/********************** check_fermion_force.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the fermion force calculation */

#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#else
#error Must compile with QIO
#endif


void check_gauge_force( char *ansfile, int ansflag )
{
  Real diff, maxdiff, norm, maxnorm, reldiff;
  int i, dir;
  site *s;
  Real epsilon = 1.;
  su3_matrix tmat, diffmat;
  char *filexml;
  char recxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test gauge force field</title>";
#if (MILC_PRECISION == 1)
  Real tol = 1e-3;
#else
  Real tol = 1e-7;
#endif
  int ff_prec = MILC_PRECISION;  /* Just use prevailing precision for now */

  beta = 1.;
  n_dyn_masses = 3;
  dyn_flavors[0] = 2;
  dyn_flavors[1] = 1;
  dyn_flavors[2] = 1;

  su3_matrix *ansmom = (su3_matrix *)malloc(4*sites_on_node*sizeof(su3_matrix));
  if(ansmom == NULL){
    node0_printf("No room for ansmom\n");
    terminate(1);
  }

  /* Compute the gauge action */
  double g_action = (beta/3.0)*imp_gauge_action_ks();
  node0_printf("The gauge action is %.14e\n", g_action);

  /* Compute the gauge force */

  node0_printf("Computing the gauge force\n"); fflush(stdout);
  
  /* Just to be safe, clear the answer */
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      s->mom[dir].m00im = 0.0;
      s->mom[dir].m11im = 0.0;
      s->mom[dir].m22im = 0.0;
      s->mom[dir].m01.real = 0.0;
      s->mom[dir].m01.imag = 0.0;
      s->mom[dir].m02.real = 0.0;
      s->mom[dir].m02.imag = 0.0;
      s->mom[dir].m12.real = 0.0;
      s->mom[dir].m12.imag = 0.0;
    }
  }

  imp_gauge_force_ks(epsilon,F_OFFSET(mom));

  /* If the answer file is given, read it for comparison */
  if(ansflag == RELOAD_SERIAL){
    restore_color_matrix_scidac_to_field(ansfile, ansmom, 4, MILC_PRECISION);
    node0_printf("Checking the answer\n"); fflush(stdout);
  }

  /* Unpack the answer and compare if possible */
  maxdiff = 0;
  maxnorm = 0;
  norm = 0;
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat );
      /* If we have loaded an answer file, do the comparison */
      if(ansflag == RELOAD_SERIAL){
	sub_su3_matrix( ansmom + 4*i + dir, &tmat, &diffmat);
	diff = sqrt(realtrace_su3( &diffmat, &diffmat ));
	norm = sqrt(realtrace_su3( &tmat, &tmat));
	//printf("DIFF %g %g\n",norm,diff);
	if(diff > tol * norm){
	  printf("Intolerable relative difference %e node %d site %d dir %d\n",
		 diff/norm,this_node,i,dir);
	  dumpmat(ansmom + 4*i + dir);
	  dumpmat(&tmat);
	}
	if(maxdiff < diff)maxdiff = diff;
	if(maxnorm < norm)maxnorm = norm;
      }
      /* In any case, copy the new result to the answer matrix */
      ansmom[4*i + dir] = tmat;
    }
  }

  if(ansflag == RELOAD_SERIAL){
    g_floatmax(&maxdiff);
    g_floatmax(&maxnorm);
    if(maxnorm > 0){
      reldiff = maxdiff/maxnorm;
      node0_printf("Relative difference %e\n",reldiff);
    }
    else
      node0_printf("Absolute difference %e but norm is 0???\n",maxdiff);
  }      

  /* Save ansmom if requested */

#ifdef HAVE_QIO

  if(ansflag == SAVE_SERIAL){
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
	recxml, QIO_SINGLEFILE, ansmom, 4, MILC_PRECISION, NULL);
    free_QCDML(filexml);
  }
  else if(ansflag == SAVE_PARTFILE_SCIDAC){
    node0_printf("Saving the momentum matrix\n");
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
	recxml, QIO_PARTFILE, ansmom, 4, MILC_PRECISION, NULL);
    free_QCDML(filexml);
  }
#endif

  free(ansmom);
}      
