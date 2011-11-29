/*********************** check_ks_invert.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* NOTE: REMOVE MdaggerM inverse CHECK */

/* This code performs and/or checks the KS inversion */

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#endif

void check_ks_invert( char *srcfile, int srcflag, field_offset src, 
		      char *ansfile, int ansflag, field_offset ans, 
		      field_offset tmp, Real mass)
{
  Real final_rsq = 0;
  /* Note: these are absolute, not relative errors. */
#if (PRECISION == 1)
  Real tol_M = 1e-2;
  Real tol_MdagM = 1e-3;
#else
  Real tol_M = 5e-6;
  Real tol_MdagM = 1e-7;
#endif
  int iters = 0;
  char srcfilexml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>KS Invert Test</title>";
  char srcrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Sample source color vector field</title>";
  char ansMrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test answer = M^-1 source</title>";
  char ansMdMrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test answer = answer = (MdaggerM)^-1 source</title>";
  imp_ferm_links_t *fn = get_fm_links(fn_links)[0];
  
  /* Make a random source in phi if we don't reload it */
  if(srcflag == RELOAD_SERIAL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_site (srcfile, QIO_SERIAL, src, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else if(srcflag == RELOAD_PARALLEL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_site (srcfile, QIO_PARALLEL, src, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else {
    /* generate g_rand random; phi = Mdagger g_rand */
    grsource_imp( src, mass, EVENANDODD, fn );
    node0_printf("Generating a random source\n");
  }
  
  /* Do the inversion if we aren't reloading the answer */
  if(ansflag == RELOAD_SERIAL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_site (ansfile, QIO_SERIAL, ans, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else if(ansflag == RELOAD_PARALLEL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_site (ansfile, QIO_PARALLEL, ans, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else {
    node0_printf("Doing the inversion\n");
    clear_latvec( ans, EVENANDODD );
   if(inverttype == INVERT_M){
      /* Compute M^-1 phi */
      iters = mat_invert_uml( src, ans, tmp, mass, PRECISION, fn );
      node0_printf("Inversion required %d iters\n",iters);
    }
    else {
      /* Compute (M^dagger M)^-1 phi */
      iters = ks_congrad( src, ans, mass,
			  niter, nrestart, rsqprop, PRECISION, 
			  EVENANDODD, &final_rsq, fn );
      node0_printf("Inversion required %d iters resid %e\n",iters,final_rsq);
    }
  }
  
  /* Check the inversion */
  node0_printf("Checking the inversion\n");
  if(inverttype == INVERT_M)
    /* Is M xxx = phi ? */
    check_invert( ans, src, mass, tol_M, fn );
  else
    /* Is MdaggerM xxx = phi ? */
    check_invert2( ans, src, tmp, mass, tol_MdagM, EVENANDODD, fn );
  
  /* Save source and answer if requested */
#ifdef HAVE_QIO
  if(srcflag == SAVE_SERIAL)
    save_ks_vector_scidac_from_site(srcfile, srcfilexml, srcrecxml, 
				    QIO_SINGLEFILE, QIO_SERIAL, src, 1);
  
  if(ansflag == SAVE_SERIAL){
    if(inverttype == INVERT_M)
      save_ks_vector_scidac_from_site(ansfile, srcfilexml, ansMrecxml, 
				      QIO_SINGLEFILE, QIO_SERIAL, ans, 1);
    else
      save_ks_vector_scidac_from_site(ansfile, srcfilexml, ansMdMrecxml,
				      QIO_SINGLEFILE, QIO_SERIAL, ans, 1);
  }
#endif
}      
