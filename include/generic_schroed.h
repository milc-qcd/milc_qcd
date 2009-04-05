#ifndef _GENERIC_SCHROED_H
#define _GENERIC_SCHROED_H
/************************ generic_schroed.h *****************************
*									*
*  Macros and declarations for miscellaneous generic_schroed routines   *
*  This header is for codes that call generic_schroed routines          *
*  MIMD version 7 							*
*									*
*/

/* Most routines have their "generic" equivalent and are
   prototyped in generic.h.
*/

/* coupling.c */
void coupling(double *ds_deta, double *bd_plaq);

/* make_schroed_lattice.c */
void set_boundary_fields(void);

#endif	/* _GENERIC_SCHROED_H */

