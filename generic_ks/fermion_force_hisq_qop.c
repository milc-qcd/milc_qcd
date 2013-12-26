/******* fermion_force_hisq_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 5/25/09  C. DeTar and A Bazavov */

/* External entry points in this file

   eo_fermion_force_oneterm
   eo_fermion_force_twoterms
   eo_fermion_force_multi
   ks_multiff_opt_chr

 */

/* Compile with fermion_force_hisq_qop_[FD].c */

/*
 * $Log: fermion_force_hisq_qop.c,v $
 * Revision 1.3  2013/12/26 05:21:07  detar
 * Support eo_fermion_force_oneterm, eo_fermion_force_oneterm_site,
 *   eo_fermion_force_twoterms and eo_fermion_force_twoterms_site
 *
 * Revision 1.2  2012/01/21 21:04:12  detar
 * Urs Heller's upgrades to eigen_stuff*.c
 *
 * Revision 1.1  2011/11/29 20:42:29  detar
 * Add
 *
 */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"
#include "../include/fermion_links.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

//static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_hisq_qop.c,v 1.3 2013/12/26 05:21:07 detar Exp $";

/**********************************************************************/
/*   Parallel transport nterms source vectors                        */
/**********************************************************************/

static void 
fermion_force_multi( Real eps, Real *residues, 
		     su3_vector **xxx, int nterms, int prec,
		     fermion_links_t *fl ) 
{

  if(prec == 1)
    fermion_force_multi_hisq_F( eps, residues, xxx, n_orders_naik, fl );
  else
    fermion_force_multi_hisq_D( eps, residues, xxx, n_orders_naik, fl );

}

static void set_qop_hisq_force_opts( int fnmat_src_min, int veclength ) {
  /* Note: the want_deps and want_aux options are set in hisq_links_qop.c */

  /* Set values */
  QOP_opt_t qop_hf_opt[2] = {

    //    {.tag = "fnmat_src_min",.value=fnmat_src_min},
    {.tag = "veclength",.value=veclength},


#ifdef HISQ_FORCE_FILTER
    {.tag = "force_filter",.value=HISQ_FORCE_FILTER}
#else
    {.tag = "force_filter",.value=0},
#endif

  };

  /* Set links options, overriding defaults */

  if(QOP_hisq_force_set_opts(qop_hf_opt, 2) != QOP_SUCCESS)
    node0_printf("eo_fermion_force_multi: error setting QOP options\n");
  
}

// /**********************************************************************/
// /*   Parallel transport vectors in blocks of veclength.               */
// /**********************************************************************/
// /* Requires the xxx1 and xxx2 terms in the site structure */
// 
// void 
// fermion_force_block( Real eps, Real *residues, 
// 		     su3_vector **xxx, int nterms, int veclength, 
// 		     int prec, fermion_links_t *fl) 
// {
//   if(prec == 1)
//     fermion_force_block_F( eps, residues, xxx, nterms, veclength, fl);
//   else
//     fermion_force_block_D( eps, residues, xxx, nterms, veclength, fl);
// 
// }

/**********************************************************************/
/*   Standard MILC interface for fermion force with multiple sources  */
/**********************************************************************/

void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, fermion_links_t *fl ) {

  int veclength, fnmat_src_min;

#ifdef VECLENGTH
  veclength = VECLENGTH;
#else
  veclength = 4;
#endif

 // Discontinue use as of qopqdp-0.17.5
  switch(KS_MULTIFF){
  case ASVEC:
    fnmat_src_min = nterms + 1;
    break;
  default: /* FNMAT */
    fnmat_src_min = 4;
  }

  set_qop_hisq_force_opts(fnmat_src_min, veclength);

  fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
}

/**********************************************************************/
/* Standard MILC interface for the single-species HISQ fermion force
   routine */
/**********************************************************************/ 

void eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *x_off_one,
			       int prec, fermion_links_t *fl)
{

  Real residues[1] = { weight };
  su3_vector *xxx[1] = { x_off_one };
  int nterms = 1;


  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

}

void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off_site,
				    int prec, fermion_links_t *fl)
{

  su3_vector *x_off_one = create_v_field_from_site_member(x_off_site);

  eo_fermion_force_oneterm(eps, weight, x_off_one, prec, fl);

  destroy_v_field(x_off_one);

}

/**********************************************************************/
/* Standard MILC interface for the two-species HISQ fermion force
   routine */
/**********************************************************************/
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, fermion_links_t *fl)
{

  Real residues[2] = { weight1, weight2 };
  su3_vector *xxx[2] = { x1_off, x2_off };
  int nterms = 2;

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

}


void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2, 
				     field_offset x1_off_site, field_offset x2_off_site,
				     int prec, fermion_links_t *fl)
{

  su3_vector *x1_off = create_v_field_from_site_member(x1_off_site);
  su3_vector *x2_off = create_v_field_from_site_member(x2_off_site);

  eo_fermion_force_twoterms(eps, weight1, weight2, x1_off, x2_off, prec, fl);

  destroy_v_field(x2_off);
  destroy_v_field(x1_off);

}


/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char 
*ks_multiff_opt_chr( void )
{
  switch(KS_MULTIFF){
  case ASVEC:
    return "ASVEC";
    break;
  default:
    return "FNMAT";
  }
  return NULL;
}

