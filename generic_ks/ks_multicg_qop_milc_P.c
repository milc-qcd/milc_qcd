/******* ks_multicg_qop_milc_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This version implements the QOP API with the standard MILC algorithm
 * It is intended for comparing results with other QOP routines.
 * It is based only on d_congrad5_fn_qop_milc.c
 * At present it does not use the MILC multicg inverter when it might.
 */

/* 8/2006 C. DeTar Created */

/*
 * $Log: ks_multicg_qop_milc_P.c,v $
 * Revision 1.3  2012/11/24 00:02:54  detar
 * Add placeholders for HYPISQ action.  Support HISQ action within ks_imp_dyn.
 *
 * Revision 1.2  2011/12/03 22:30:09  detar
 * Cosmetic: Fix file name in top comment line
 *
 * Revision 1.1  2006/12/09 13:52:40  detar
 * Add mixed precision capability for KS inverter in QOP and QDP
 *
 * Revision 1.1  2006/08/13 15:01:50  detar
 * Realign procedures to accommodate ks_imp_rhmc code
 * Add Level 3 wrappers and MILC dummy Level 3 implementation
 *
 */

#if ( QOP_PrecisionInt == 1 )
#define MYREAL float
#else
#define MYREAL double
#endif

#include "generic_ks_includes.h"
#include <qop.h>

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/ks_multicg_qop_milc_P.c,v 1.3 2012/11/24 00:02:54 detar Exp $";

/* Just do repeated single inversions */
void QOP_asqtad_invert_multi(QOP_info_t *info,
			     QOP_FermionLinksAsqtad *links,
			     QOP_invert_arg_t *inv_arg,
			     QOP_resid_arg_t **res_arg[],
			     MYREAL *masses[], int nmass[],
			     QOP_ColorVector **dest_pt[],
			     QOP_ColorVector *src_pt[],
			     int nsrc){
  int isrc, imass;

  for(isrc = 0; isrc < nsrc; isrc++){
    for(imass = 0; imass < nmass[isrc]; imass++){

      /* Do the inversion */
      QOP_asqtad_invert(info, links, inv_arg, res_arg[isrc][imass],
			masses[isrc][imass],
			dest_pt[isrc][imass], src_pt[isrc]);
      if(info->status != QOP_SUCCESS)break;
    }
    if(info->status != QOP_SUCCESS)break;
  }
}
