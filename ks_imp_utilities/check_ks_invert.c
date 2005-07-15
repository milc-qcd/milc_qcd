/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

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
  Real tol = 1e-5;
  int iters = 0;
  char *filexml;
  
  /* Make a random source in phi if we don't reload it */
  if(srcflag == RELOAD_SERIAL){
    restore_ks_vector_scidac_to_site (srcfile, src, QIO_SERIAL, 1);
  }
  else {
    /* generate g_rand random; phi = Mdagger g_rand */
    grsource_imp( src, mass, EVENANDODD );
    node0_printf("Generating a random source\n");
  }
  
  /* Do the inversion if we aren't reloading the answer */
  if(ansflag == RELOAD_SERIAL){
    restore_ks_vector_scidac_to_site (ansfile, ans, QIO_SERIAL, 1);
  }
  else {
    node0_printf("Doing the inversion\n");
    clear_latvec( ans, EVENANDODD );
    if(inverttype == INVERT_M){
      /* Compute M^-1 phi */
      iters = mat_invert_uml( src, 
			      ans, tmp, mass);
      node0_printf("Inversion required %d iters\n",iters);
    }
    else {
      /* Compute (M^dagger M)^-1 phi */
      iters = ks_congrad( src, ans, mass,
			  niter, rsqprop, EVENANDODD, &final_rsq);
      node0_printf("Inversion required %d iters resid %e\n",iters,final_rsq);
    }
  }
  
  /* Check the inversion */
  node0_printf("Checking the inversion\n");
  if(inverttype == INVERT_M)
    /* Is M xxx = phi ? */
    check_invert( ans, src, mass, tol);
  else
    /* Is MdaggerM xxx = phi ? */
    check_invert2( ans, src, tmp, mass, tol, EVENANDODD);
  
  /* Save source and answer if requested */
#ifdef HAVE_QIO
  if(srcflag == SAVE_SERIAL)
    save_ks_vector_scidac_from_site(srcfile, "source color vector field", 
				    QIO_SINGLEFILE, QIO_SERIAL, src, 1);
  
  if(ansflag == SAVE_SERIAL){
    if(inverttype == INVERT_M)
      save_ks_vector_scidac_from_site(ansfile, "answer = M^-1 source", 
				      QIO_SINGLEFILE, QIO_SERIAL, ans, 1);
    else
      save_ks_vector_scidac_from_site(ansfile, "answer = (MdaggerM)^-1 source",
				      QIO_SINGLEFILE, QIO_SERIAL, ans, 1);
  }
#endif
}      
