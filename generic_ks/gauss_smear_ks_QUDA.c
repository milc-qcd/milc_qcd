/************************** gauss_smear_ks_QUDA.c ********************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source using QUDA
 */

#include "generic_ks_includes.h"

#if defined(HAVE_QUDA) && defined(USE_GSMEAR_GPU)

#include <string.h>
#include <assert.h>

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

/* #define GS_TIME */
/* #define GS_DEBUG */

static int twolink = 0; /* Presumed QUDA state of twolink field. 0 = not computed */
static int recompute_2link = 1; 

void
gauss_smear_reuse_2link_QUDA( int flag )
{
  recompute_2link = flag? 0 : 1;
}

/* Delete saved two-link.
 */
void
gauss_smear_delete_2link_QUDA()
{
  if( initialize_quda() )
  {
    node0_printf( "%s: FATAL. QUDA has not been initialized.\n", __func__ );
    terminate(1);
  }

  qudaFreeTwoLink();
  twolink = 0;
  gauss_smear_reuse_2link_QUDA(0);

  return ;
}

/* Perform gauss_smear_v_field (in gauss_smear_ks.c) on GPU
   using performTwoLinkGaussianSmearNStep() in QUDA.
 */
void
gauss_smear_v_field_QUDA(su3_vector *src, su3_matrix *t_links,
                         Real width, int iters, int t0)
{
  char myname[] = "gauss_smear_v_field_QUDA";
  static su3_matrix *t_links_last = NULL;

#ifdef GS_DEBUG
  node0_printf( "%s: Start\n", myname );
#endif

  if(t_links == NULL){
    printf("%s(%d): NULL t_links\n", __func__, this_node);
    terminate(1);
  }

  if ( t_links != t_links_last){
    /* If the t_link pointer changed, refresh the two-link */
    node0_printf("%s: [Warning] Input field for two-links changed.  Will recompute them\n", __func__);
    gauss_smear_delete_2link_QUDA();
    t_links_last = t_links;
  }

#ifdef GS_TIME
  double dtimec;
#endif
  
  /* Initialize QUDA */
  if( initialize_quda() )
  {
    node0_printf( "%s: FATAL. QUDA has not been initialized.\n", myname );
    terminate(1);
  }
  
  /* Input parameters ***************************/
  int laplaceDim = 3;
  /**********************************************/

  if( laplaceDim > 3 && t0 != ALL_T_SLICES )
  {
    node0_printf( "%s: [Warning] t0 is ignored for d>3 dimensional Laplacian.\n", myname );
    t0 = ALL_T_SLICES;
  }

  if( ! recompute_2link && twolink == 0 )
  {
    node0_printf( "%s: [Warning] recompute_2link is false but there is no saved two-link. Two-link will be calculated.\n", myname );
    gauss_smear_reuse_2link_QUDA(0);
  }

#ifdef GS_TIME
  dtimec = -dclock();
#endif
#ifdef GS_DEBUG
  node0_printf( "%s: Gaussian smearing starts.\n", myname );
#endif

  /* Quark smearing parameters */
  QudaTwoLinkQuarkSmearArgs_t qsmear_args;
  qsmear_args.n_steps = iters;
  qsmear_args.width = width;
  qsmear_args.compute_2link = recompute_2link ? 1 : 0;
  qsmear_args.delete_2link = 0;
  qsmear_args.t0 = t0;
  qsmear_args.laplaceDim = laplaceDim;

  /* Run gaussian smearing */
  qudaTwoLinkGaussianSmear( MILC_PRECISION, MILC_PRECISION, (void*) t_links, (void*) src, qsmear_args );

  /* two-link is saved. */
  twolink = 1;
  gauss_smear_reuse_2link_QUDA(1);


#ifdef GS_DEBUG
  node0_printf( "%s: Gaussian smearing ends.\n", myname );
#endif  
#ifdef GS_TIME
  dtimec += dclock();
  node0_printf( "[GS_TIME] QUDA two-link Gaussian smearing: time = %g s, iters = %d\n", dtimec, iters );
#endif

  return ;
}

#else /* #ifdef USE_GSMEAR_GPU */

void
gauss_smear_v_field_QUDA(su3_vector *src, su3_matrix *t_links,
                         Real width, int iters, int t0)
{
  char myname[] = "gauss_smear_v_field_QUDA";
  node0_printf( "%s: Requires compilation with the QUDA library\n", myname );
  terminate(1);

  return ;  
}

#endif
