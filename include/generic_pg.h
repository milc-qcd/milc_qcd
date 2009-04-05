#ifndef _GENERIC_PG_H
#define _GENERIC_PG_H
/************************ generic_pg.h **********************************
*									*
*  Macros and declarations for generic_pg routines                      *
*  This header is for codes that call generic_pg routines               *
*  MIMD version 7 							*
*									*
*/

#include "../include/macros.h"

int update(void);
void update_h(Real eps);
void update_u(Real eps);
void relax(int NumStp);
void monte(int NumStp);
void dsdu_qhb(int dir1, int parity);
double d_action(void);
void gauge_field_copy(field_offset src, field_offset dest);
#endif	/* _GENERIC_PG_H */
