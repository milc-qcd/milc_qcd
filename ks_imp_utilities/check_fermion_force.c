/********************** check_fermion_force.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This code performs and/or checks the fermion force calculation */

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#else
BOMB THE COMPILE
#endif

void check_fermion_force( char srcfile[MAX_MASS][MAXFILENAME], int srcflag,
			  char *ansfile, int ansflag, int nmass, ks_param *ksp)
{
  Real diff, maxdiff, norm, maxnorm, reldiff;
  int i, dir;
  site *s;
  Real eps = 1.;
  int nflavors = 4;
  su3_matrix tmat, diffmat;
  char *filexml;
  char recxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test fermion force field</title>";
#if (PRECISION == 1)
  Real tol = 1e-3;
#else
  Real tol = 1e-5;
#endif
  int ff_prec = PRECISION;  /* Just use prevailing precision for now */
  /* Supports only asqtad at the moment */
  imp_ferm_links_t **fn = get_fm_links(fn_links);
  Real *residues = (Real *)malloc(nmass*sizeof(Real));;
  su3_vector **src = (su3_vector **)malloc(nmass*sizeof(su3_vector *));
  su3_matrix *ansmom = (su3_matrix *)malloc(4*sites_on_node*sizeof(su3_matrix));

  if(residues == NULL){
    node0_printf("No room for residues\n");
    terminate(1);
  }

  /* For testing, we just set them to a constant */
  for(i = 0; i < nmass; i++)
    residues[i] = ((Real)nflavors)/4.;
  
  if(src == NULL){
    node0_printf("No room for src\n");
    terminate(1);
  }
  
  if(ansmom == NULL){
    node0_printf("No room for ansmom\n");
    terminate(1);
  }

  for(i = 0; i < nmass; i++)
    src[i] = create_v_field();

  /* Make a random source in src if we don't reload it */

  for(i = 0; i < nmass; i++){
    if(srcflag == RELOAD_SERIAL){
      restore_ks_vector_scidac_to_field(srcfile[i], QIO_SERIAL, src[i], 1);
      fflush(stdout);
    }  else {
      /* generate g_rand random; phi = Mdagger g_rand */
      node0_printf("Generating random sources\n");
      grsource_imp_field( src[i], mass, EVENANDODD, fn[i]);
    }
  }
      
  node0_printf("Computing the fermion force\n"); fflush(stdout);
  
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

  eo_fermion_force_multi( eps, residues, src, nmass, ff_prec, fn_links );

  /* If the answer file is given, read it for comparison */
  if(ansflag == RELOAD_SERIAL){
    restore_color_matrix_scidac_to_field(ansfile, ansmom, 4, PRECISION);
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
	printf("DIFF %g %g\n",norm,diff);
	if(diff > tol * norm){
	  printf("Intolerable relative difference %e node %d site %d\n",
		 diff/norm,this_node,i);
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

  /* Save source and answer if requested */
  if(srcflag == SAVE_SERIAL || srcflag == SAVE_PARTFILE_SCIDAC)
    for(i = 0; i < nmass; i++){
#ifdef HAVE_QIO
      if(srcflag == SAVE_SERIAL)
	save_ks_vector_scidac_from_field(srcfile[i], "check fermion force",
					"source color vector field", 
					QIO_SINGLEFILE, QIO_SERIAL, src[i], 1, PRECISION);
      else if(srcflag == SAVE_PARTFILE_SCIDAC)
	save_ks_vector_scidac_from_field(srcfile[i], "check fermion force",
					"source color vector field",
					QIO_PARTFILE, QIO_SERIAL, src[i], 1, PRECISION);
    }

  if(ansflag == SAVE_SERIAL){
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
        recxml, QIO_SINGLEFILE, ansmom, 4, PRECISION);
    free_QCDML(filexml);
  }
  else if(ansflag == SAVE_PARTFILE_SCIDAC){
    node0_printf("Saving the momentum matrix\n");
    filexml = create_QCDML();
    save_color_matrix_scidac_from_field(ansfile, filexml, 
       recxml, QIO_PARTFILE, ansmom, 4, PRECISION);
    free_QCDML(filexml);
  }
#endif

}      
