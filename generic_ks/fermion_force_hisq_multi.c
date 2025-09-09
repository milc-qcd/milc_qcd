/****** fermion_force_hisq_multi.c  -- ******************/

/* MIMD version 7 */
/* Multisource fermion force.  */

/* May 5, 2024 C. DeTar: Most code moved to fermion_force_hisq_multi.c 
   and fermion_force_hisq_cpi.c*/

/* External entry points in this file

   eo_fermion_force_multi
   ks_multiff_opt_chr

 */

#include "generic_ks_includes.h"
#include "../include/fermion_links.h"
#include <string.h>

/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
#ifdef FFTIME
static const char *prec_char[2] = {"F", "D"};
#endif

void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, fermion_links_t *fl) {
  info_t info = INFO_ZERO;

  double dtime = -dclock();

  if(KS_MULTIFF != FNMAT){
    node0_printf("Warning: Bad KS_MULTIFF value: the only HISQ choice is FNMAT\n");
    node0_printf("Reverting to FNMAT\n");
  }

#ifdef USE_FF_GPU

#if defined(HAVE_QUDA)
  fermion_force_multi_hisq_quda(&info, prec, eps, residues, xxx, nterms, fl );
#elif defined(HAVE_GRID)
  fermion_force_multi_hisq_grid(&info, prec, eps, residues, xxx, nterms, fl );
#endif

#else

  fermion_force_multi_hisq_cpu(&info, prec, eps, residues, xxx, nterms, fl );

#endif

  /* Accumulate SVD statistics */
  g_intsum(&INFO_HISQ_SVD_COUNTER(&info));
  hisq_svd_counter += INFO_HISQ_SVD_COUNTER(&info);
  
  g_intsum(&INFO_HISQ_FORCE_FILTER_COUNTER(&info));
  hisq_force_filter_counter += INFO_HISQ_FORCE_FILTER_COUNTER(&info);

  dtime += dclock();
  info.final_sec = dtime;
#ifdef FFTIME
#ifdef USE_FF_GPU

#if defined(HAVE_QUDA)
  node0_printf("FFTIME:  time = %e (HISQ QUDA %s) terms = %d flops/site = %d mflops = %e\n",
	       info.final_sec,prec_char[MILC_PRECISION-1],nterms,
	       (int)(info.final_flop*numnodes()/volume),
	       info.final_flop/(1e6*info.final_sec) );
  #elif defined(HAVE_GRID)
  node0_printf("FFTIME:  time = %e (HISQ GRID %s) terms = %d flops/site = %d mflops = %e\n",
	       info.final_sec,prec_char[MILC_PRECISION-1],nterms,
	       (int)(info.final_flop*numnodes()/volume),
	       info.final_flop/(1e6*info.final_sec) );
  #endif

#else
  node0_printf("FFTIME:  time = %e (HISQ MILC %s) terms = %d flops/site = %d mflops = %e\n",
	       info.final_sec,prec_char[MILC_PRECISION-1],nterms,
	       (int)(info.final_flop*numnodes()/volume),
	       info.final_flop/(1e6*info.final_sec) );
#endif
#endif
}

#if 0
n
/**********************************************************************/
/* These routines are left here for consistency, but unsupported
   until we decide what a HISQ one-term or two-term force means */
/**********************************************************************/

void 
eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *temp_x,
			  int prec, fermion_links_t *fl )
{
  su3_vector *xxx[1] = { temp_x };
  int nterms = 1;
  Real residues[1] = { weight };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
  
}

void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
				    int prec, fermion_links_t *fl)
{
  su3_vector *temp_x = create_v_field_from_site_member(x_off);
  su3_vector *xxx[1] = { temp_x };
  int nterms = 1;
  Real residues[1] = { weight };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

  destroy_v_field(temp_x);
}

void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, fermion_links_t *fl ){
  su3_vector *xxx[2] = { x1_off, x2_off };
  int nterms = 2;
  Real residues[2] = { weight1, weight2 };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
}

void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2,
				     field_offset x1_off, field_offset x2_off,
				     int prec, fermion_links_t *fl)
{
  su3_vector *x1 = create_v_field_from_site_member(x1_off);
  su3_vector *x2 = create_v_field_from_site_member(x2_off);
  su3_vector *xxx[2] = { x1, x2 };
  int nterms = 2;
  Real residues[2] = { weight1, weight2 };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

  destroy_v_field(x1);
  destroy_v_field(x2);
}

#endif /* #if 0 */

/**********************************************************************/
/*   Accessor for string describing the multiff option                */
/**********************************************************************/
const char 
*ks_multiff_opt_chr( void )
{
  return "FNMAT";
}
