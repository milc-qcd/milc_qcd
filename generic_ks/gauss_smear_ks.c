/************************** gauss_smear_ks.c ********************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source, Chroma style
 * 
 * CD 4/07 Stolen from Chroma. MILC version.
 */

#include "generic_ks_includes.h"

/*------------------------------------------------------------*/

/* Compute
   src <- cov_gauss src
   
   using t_links
*/

#if defined(HAVE_QUDA) && defined(USE_GSMEAR_GPU)

void 
gauss_smear_v_field(su3_vector *src, su3_matrix *t_links,
		    Real width, int iters, int t0)
{
  gauss_smear_reuse_2link_QUDA( 1 );
  gauss_smear_v_field_QUDA( src, t_links, width, iters, t0 );
}

#else

void 
gauss_smear_v_field(su3_vector *src, su3_matrix *t_links,
		    Real width, int iters, int t0)
{
#ifdef GAUSS_SMEAR_KS_TWOLINK
  gauss_smear_v_field_cpu_twolink( src, t_links, width, iters, t0 );
#else
  gauss_smear_v_field_cpu( src, t_links, width, iters, t0 );
#endif
}

#endif

/*------------------------------------------------------------*/

/* Computes 
   src <- Lapl_3d src 
*/

void 
laplacian_v_field(su3_vector *src, su3_matrix *t_links, int t0)
{
  su3_vector *tmp;

  if(t_links == NULL){
    printf("laplacian_v_field(%d): NULL t_links\n",this_node);
    terminate(1);
  }

  tmp = create_v_field();
  copy_v_field(tmp, src);

#ifdef GAUSS_SMEAR_KS_TWOLINK
  klein_gord_field_twolink(tmp, src, t_links, 0., t0);
#else
  klein_gord_field(tmp, src, t_links, 0., t0);
#endif

  destroy_v_field(tmp);
}

void 
gauss_smear_ks_prop_field(ks_prop_field *src, su3_matrix *t_links,
			  Real width, int iters, int t0)
{
  int color;

  for(color = 0; color < src->nc; color++)
    gauss_smear_v_field(src->v[color], t_links, width, iters, t0);

}
