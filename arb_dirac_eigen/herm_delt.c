/******* herm_delt.c  ****/
/* MIMD version 6 */

/*
dest= herm_delt*src

Put your favorite Hermetian Dirac operator here!
*/



#include "arb_dirac_eig_includes.h"




void herm_delt(field_offset src,field_offset dest)
{
	register int i;
	register site *s;

	/* D^\dagger D */
delta0(src,F_OFFSET(r),PLUS);
delta0(F_OFFSET(r),dest,MINUS);

}
/* herm_delt.c */
