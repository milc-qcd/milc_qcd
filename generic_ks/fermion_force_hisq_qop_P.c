/******* fermion_force_hisq_qop_P.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 01/31/2011 A. Bazavov, with support for the HISQ fermion force */
/* 11/2005 C. DeTar and D. Renner */
/* 5/09/07 C. DeTar Support selected precision */
/* NOTE: This is an include file for fermion_force_hisq_qop_F.c and
   fermion_force_hisq_qop_D.c.  Any changes must take this usage
   into account */

/* External entry points in this file are private to
   fermion_force_hisq_qop.c

   eo_fermion_force_oneterm_[FD]  -- DOES NOT WORK, TEMPORARILY EXCLUDED
   eo_fermion_force_twoterms_[FD] -- DOES NOT WORK, TEMPORARILY EXCLUDED
   fermion_force_multi_[FD]

*/

/*
 * $Log: fermion_force_hisq_qop_P.c,v $
 * Revision 1.4  2012/11/24 00:02:51  detar
 * Add placeholders for HYPISQ action.  Support HISQ action within ks_imp_dyn.
 *
 * Revision 1.3  2012/02/22 03:46:31  detar
 * Add precision tag to timing line
 *
 * Revision 1.2  2012/02/16 16:30:29  detar
 * Initialize QOP_info_t structure.
 *
 * Revision 1.1  2011/11/29 20:42:28  detar
 * Add
 *
 *
 */

#if (QOP_PrecisionInt==1)

#define CREATE_F_FROM_SITE4        create_F_F_from_site4
#define CREATE_V_FROM_SITE         create_F_V_from_site
#define CREATE_V_FROM_FIELD        create_F_V_from_field
#define EO_FERMION_FORCE_ONETERM   eo_fermion_force_oneterm_F
#define EO_FERMION_FORCE_TWOTERMS  eo_fermion_force_twoterms_F
#define FERMION_FORCE_MULTI        fermion_force_multi_hisq_F
#define FERMION_FORCE_BLOCK        fermion_force_block_F
#define MY_REAL                    float
#define UNLOAD_F_TO_SITE4          unload_F_F_to_site4
#define UNLOAD_G_TO_SITE4          unload_F_G_to_site4
#define GET_HISQ_LINKS             get_F_hisq_links
#define HL                         hl_F
				   
#else				   
				   
#define CREATE_F_FROM_SITE4        create_D_F_from_site4
#define CREATE_V_FROM_SITE         create_D_V_from_site
#define CREATE_V_FROM_FIELD        create_D_V_from_field
#define EO_FERMION_FORCE_ONETERM   eo_fermion_force_oneterm_D
#define EO_FERMION_FORCE_TWOTERMS  eo_fermion_force_twoterms_D
#define FERMION_FORCE_MULTI        fermion_force_multi_hisq_D
#define FERMION_FORCE_BLOCK        fermion_force_block_D
#define MY_REAL                    double
#define UNLOAD_F_TO_SITE4          unload_D_F_to_site4
#define UNLOAD_G_TO_SITE4          unload_D_G_to_site4
#define GET_HISQ_LINKS             get_D_hisq_links
#define HL                         hl_G

#endif

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

#ifdef FFTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

//static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_hisq_qop_P.c,v 1.4 2012/11/24 00:02:51 detar Exp $";

#if 0  /* Not supported until we decide what a HISQ one-term or two-term force means */
/**********************************************************************/
/* Standard MILC interface for the single-species Hisq fermion force
   routine */
/**********************************************************************/
void EO_FERMION_FORCE_ONETERM( Real eps, Real weight, su3_vector *x_off,
			       fermion_links_t *fl)
{

  /* For example weight = nflavors/4 */

  QOP_FermionLinksHisq *HL = GET_HISQ_LINKS(fl);
  QOP_hisq_coeffs_t *coeff = get_action_coeffs_hisq(fl);
  QOP_Force *mom;
  QOP_ColorVector *vecx;
  
  QOP_info_t info = {0., 0., 0, 0, 0};

  double remaptime = -dclock();

  QOP_opt_t qop_ff_opt[2] = {
    {.tag = "fnmat_src_min",.value=4},
    {.tag = "veclength",.value=4}
  };

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    node0_printf("eo_fermion_force_oneterm: Error initializing QOP\n");
    terminate(1);
  }

  /* Convert momentum to QOP format */

  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Convert color vector to QOP format */
  vecx = CREATE_V_FROM_FIELD(x_off,EVENANDODD);

  /* Compute fermion force */
  remaptime += dclock();
  QOP_hisq_force_set_opts(qop_ff_opt, 2);

#ifdef HISQ_SVD_COUNTER
  QOP_info_hisq_svd_counter(&info) = hisq_svd_counter;
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  QOP_info_hisq_force_filter_counter(&info) = hisq_force_filter_counter;
#endif

  /* The coefficients are already loaded with weight = 0.5 */
  QOP_hisq_force(&info, HL, mom, coeff, 2.*weight*eps, vecx);

#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = QOP_info_hisq_svd_counter(&info);
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = QOP_info_hisq_force_filter_counter(&info);
#endif

  remaptime -= dclock();

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_F(mom);   mom = NULL;
  QOP_destroy_V(vecx);  vecx = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop %s) terms = 1 mflops = %e\n",
	       info.final_sec, qop_prec[QOP_PrecisionInt-1],
	       (Real)info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}

/**********************************************************************/
/* Standard MILC interface for the two-species HISQ fermion force
   routine */
/**********************************************************************/
void EO_FERMION_FORCE_TWOTERMS( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				fermion_links_t *fl )
{

  /* For example weight1 = nflavor1/4; weight2 = nflavor2/4 */

  QOP_FermionLinksHisq *HL = GET_HISQ_LINKS(fl);
  QOP_hisq_coeffs_t *coeff = get_action_coeffs_hisq(fl);
  QOP_Force *mom;
  QOP_ColorVector *vecx[2];
  
  MY_REAL epsv[2];
  QOP_info_t info = {0., 0., 0, 0, 0};
  QOP_opt_t qop_ff_opt[2] = {
    {.tag = "fnmat_src_min",.value=4},
    {.tag = "veclength",.value=4}
  };

  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force_twoterms: Error initializing QOP\n");
    terminate(1);
  }

  /* Convert momentum to QOP format */

  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Loop over fermion species */
  /* Convert color vectors in site structure to QOP format */
  vecx[0] = CREATE_V_FROM_FIELD(x1_off,EVENANDODD);
  vecx[1] = CREATE_V_FROM_FIELD(x2_off,EVENANDODD);

  /* Load coefficients */
  epsv[0] = 2.*eps*weight1;  epsv[1] =2.* eps*weight2;

#ifdef HISQ_SVD_COUNTER
  QOP_info_hisq_svd_counter(&info) = hisq_svd_counter;
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  QOP_info_hisq_force_filter_counter(&info) = hisq_force_filter_counter;
#endif

  /* Compute fermion force */
  remaptime += dclock();
  QOP_hisq_force_set_opts(qop_ff_opt, 2);
  QOP_hisq_force_multi(&info, HL, mom, coeff, epsv, vecx, 2);
  remaptime -= dclock();

#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = QOP_info_hisq_svd_counter(&info);
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = QOP_info_hisq_force_filter_counter(&info);
#endif

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_F(mom);   mom = NULL;
  QOP_destroy_V(vecx[0]);  vecx[0] = NULL;
  QOP_destroy_V(vecx[1]);  vecx[1] = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop %s) terms = %d mflops = %e\n",
	       info.final_sec, qop_prec[QOP_PrecisionInt-1],
	       2, (Real)info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}

#endif  /* #if 0 */

/**********************************************************************/
/*   Version for HISQ.  Parallel transport nterms source vectors    */
/**********************************************************************/

void FERMION_FORCE_MULTI( Real eps, Real *residues, 
			  su3_vector **xxx, int my_n_orders_naik[],
			  fermion_links_t *fl ) 
{

  QOP_FermionLinksHisq *HL = GET_HISQ_LINKS(fl);
  QOP_hisq_coeffs_t *coeff = get_action_coeffs_hisq(fl);
  int n_naiks = coeff->n_naiks;
  QOP_Force *mom;
  QOP_ColorVector **vecx;
  
  int i;
  MY_REAL *epsv;
  QOP_info_t info = {0., 0., 0, 0, 0};
  int nterms;
  char myname[] = "fermion_force_multi_hisq";

  double remaptime = -dclock();

  nterms = 0;
  for(i = 0; i < n_naiks; i++)
    nterms += my_n_orders_naik[i];

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter %s in fermion_force_qop_P.c\n",myname);fflush(stdout);
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  /* Convert momentum to QOP format */

  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Make space for QOP vector pointers */
  vecx = (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nterms);
  if(vecx == NULL){
    printf("%s(%d): no room for vecx pointers\n",myname,this_node);
    terminate(1);
  }

  /* Convert color vectors to QOP format */
  for(i = 0; i < nterms; i++){
    vecx[i] = CREATE_V_FROM_FIELD(xxx[i], EVENANDODD);
  }

  /* Compute weights, following the QOP convention for factors of 2 */
  epsv = (MY_REAL *)malloc(sizeof(MY_REAL)*nterms);
  for(i = 0; i < nterms; i++) epsv[i] = eps*residues[i];

#ifdef HISQ_SVD_COUNTER
  QOP_info_hisq_svd_counter(&info) = hisq_svd_counter;
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  QOP_info_hisq_force_filter_counter(&info) = hisq_force_filter_counter;
#endif

  /* Compute fermion force */
  remaptime += dclock();
  QOP_hisq_force_multi(&info, HL, mom, coeff, epsv, vecx, my_n_orders_naik);
  remaptime -= dclock();

#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = QOP_info_hisq_svd_counter(&info);
#endif

#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = QOP_info_hisq_force_filter_counter(&info);
#endif

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_F(mom);   mom = NULL;
  for(i = 0; i < nterms; i++){
    QOP_destroy_V(vecx[i]);  vecx[i] = NULL;
  }
  free(vecx);  vecx = NULL;
  free(epsv);  epsv = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop %s) terms = %d flops/site = %d mflops = %e\n",
	       info.final_sec, qop_prec[QOP_PrecisionInt-1],
	       nterms, (int)(info.final_flop/sites_on_node),
	       info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  %s in fermion_force_qop_P.c\n",myname);fflush(stdout);
#endif

}


// /**********************************************************************/
// /*   Parallel transport vectors in blocks of veclength. Hisq only   */
// /**********************************************************************/
// 
// void FERMION_FORCE_BLOCK( Real eps, Real *residues, 
// 			  su3_vector **xxx, int nterms, int veclength,
// 			  fermion_links_t *fl ) {
// 
//   int j;
// 
//   /* First do blocks of size veclength */
//   for( j = 0;  j <= nterms-veclength; j += veclength )
//     FERMION_FORCE_MULTI( eps, &(residues[j]), xxx+j, veclength, fl );
//   
// #ifndef ONEMASS
//   /* Continue with pairs if needed */
//   for( ; j <= nterms-2 ; j+=2 ){
//     EO_FERMION_FORCE_TWOTERMS( eps, residues[j], residues[j+1],
// 			       xxx[j], xxx[j + 1], fl );
//   }
// 
//   /* Finish with a single if needed */
//   for( ; j <= nterms-1; j++ ){
//     EO_FERMION_FORCE_ONETERM( eps, residues[j], xxx[j], fl );
//   }
// #else
//   /* Thrown in just to make it compile.  Probably never used. */
//   /* Finish with a single if needed */
//   for( ; j <= nterms-1; j++ ){
//     EO_FERMION_FORCE_ONETERM( eps, residues[j], xxx[j], fl );
//   }
// #endif
// 
// }

