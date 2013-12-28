/*********************** check_ks_invert.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* NOTE: REMOVE MdaggerM inverse CHECK */

/* This code performs and/or checks the KS inversion */

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "params.h"
#ifdef HAVE_QIO
#include <qio.h>
#endif

void check_ks_invert( char *srcfile, int srcflag, 
		      char ansfile[MAX_MASS][MAXFILENAME],
		      int ansflag[MAX_MASS],
		      int nmass, ks_param ksp[], 
		      quark_invert_control qic[])
{
  char myname[] = "check_ks_invert";
  /* Note: these are absolute, not relative errors. */
#if (PRECISION == 1)
  Real tol_M = 1e-2;
  Real tol_MdagM = 1e-3;
#else
  Real tol_M = 5e-6;
  Real tol_MdagM = 1e-7;
#endif
  int iters = 0;
  int i;
  char srcfilexml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>KS Invert Test</title>";
  char srcrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Sample source color vector field</title>";
  char ansMrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test answer = M^-1 source</title>";
  char ansMdMrecxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test answer = answer = (MdaggerM)^-1 source</title>";

  imp_ferm_links_t **fn = (imp_ferm_links_t **)malloc(nmass*sizeof(imp_ferm_links_t *));
  su3_vector *src = create_v_field();
  su3_vector **ans = (su3_vector **)malloc(nmass*sizeof(su3_vector *));

  /* Set fn links for the inversion (all the same here) */
  if(fn == NULL){
    printf("%s(%d): No room for fn\n", myname, this_node);
    terminate(1);
  }
  
  for(i = 0; i < nmass; i++)
    fn[i] = get_fm_links(fn_links)[ksp[i].naik_term_epsilon_index];

  if(ans == NULL){
    node0_printf("%s(%d)No room for ans\n", myname, this_node);
    terminate(1);
  }
  
  for(i = 0; i < nmass; i++)
    ans[i] = create_v_field();

  /* Convert masses to offsets for ks_multicg_mass_field */
  for(i = 0; i < nmass; i++){
    ksp[i].offset = 4.0*ksp[i].mass*ksp[i].mass;
  }

  /* Make a random source in phi if we don't reload it */
  if(srcflag == RELOAD_SERIAL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_field (srcfile, QIO_SERIAL, src, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else if(srcflag == RELOAD_PARALLEL){
#ifdef HAVE_QIO
    restore_ks_vector_scidac_to_field (srcfile, QIO_PARALLEL, src, 1);
#else
    printf("QIO compilation is required for loading source or answer\n");
    terminate(1);
#endif
  }
  else {
    /* generate g_rand random; phi = Mdagger g_rand */
    /* This is just for testing, so it doesn't matter that we arbitrarily pick fn[0] */
    grsource_imp_field( src, mass, EVENANDODD, fn[0] );
    node0_printf("Generating a random source\n");
  }
  
  /* Do the inversion if we aren't reloading the answer */
  if(ansflag[0] == RELOAD_SERIAL || ansflag[0] == RELOAD_PARALLEL){
    for(i = 0; i < nmass; i++){
      if(ansflag[0] == RELOAD_SERIAL){
#ifdef HAVE_QIO
	restore_ks_vector_scidac_to_field (ansfile[i], QIO_SERIAL, ans[i], 1);
#else
	printf("QIO compilation is required for loading source or answer\n");
	terminate(1);
#endif
      }
      else if(ansflag[0] == RELOAD_PARALLEL){
#ifdef HAVE_QIO
	restore_ks_vector_scidac_to_field (ansfile[i], QIO_PARALLEL, ans[i], 1);
#else
	printf("QIO compilation is required for loading source or answer\n");
	terminate(1);
#endif
      }
    }

  } else {
    node0_printf("Doing the inversion\n");
    if(inverttype == INVERT_M){
      /* Compute M^-1 phi */
      iters += mat_invert_multi( src, ans, ksp, nmass, qic, fn );
    } else {
      /* Compute (M^dagger M)^-1 phi */
      if(nmass == 1){
	iters += ks_congrad_field( src, ans[0], qic+0, ksp[0].mass, fn[0]);
      }
      else{
	int save_parity = qic->parity;
	qic->parity = EVEN;
	iters += ks_multicg_field( src, ans, ksp, nmass, qic, fn );
	qic->parity = ODD;
	iters += ks_multicg_field( src, ans, ksp, nmass, qic, fn );
	qic->parity = save_parity;
      }
    }
  }
  
  /* Save source if requested */
#ifdef HAVE_QIO
  if(srcflag == SAVE_SERIAL)
    save_ks_vector_scidac_from_field(srcfile, srcfilexml, srcrecxml, 
				     QIO_SINGLEFILE, QIO_SERIAL, src, 1, PRECISION);
#endif

  /* Check the inversion */
  for(i = 0; i < nmass; i++){
    node0_printf("Checking inversion %d\n",i);
    if(inverttype == INVERT_M)
      /* Is M xxx = phi ? */
      check_invert_field( ans[i], src, ksp[i].mass, tol_M, fn[i] );
    else
      /* Is MdaggerM xxx = phi ? */
      check_invert2( ans[i], src, ksp[i].mass, tol_MdagM, EVENANDODD, fn[i] );
  
    /* Save answer if requested */
#ifdef HAVE_QIO
    if(ansflag[0] == SAVE_SERIAL){
      if(inverttype == INVERT_M)
	save_ks_vector_scidac_from_field(ansfile[i], srcfilexml, ansMrecxml, 
					QIO_SINGLEFILE, QIO_SERIAL, ans[i], 1, PRECISION);
      else
	save_ks_vector_scidac_from_field(ansfile[i], srcfilexml, ansMdMrecxml,
					QIO_SINGLEFILE, QIO_SERIAL, ans[i], 1, PRECISION);
    }
#endif
  }

  node0_printf("total_iters = %d\n", iters);

  for(i = 0; i < nmass; i++)
    destroy_v_field(ans[i]);
  free(ans);
  free(fn);

  destroy_v_field(src);
}      
