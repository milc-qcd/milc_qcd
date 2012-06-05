/* ************************************************************ */
/*                           PARAMS.H                           */
/*                                                              */
/* Parameter buffer structure for reading in parameters		*/
/*                                                              */
/* Last Updated on 07.19.07                                     */
/*                                                              */
/* ************************************************************ */
#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"	/* for MAXFILENAME */

typedef struct {

	int stopflag;		/* if anything goes wrong */

	int nx,ny,nz,nt;	/* lattice dimensions */
	int iseed;		/* for random numbers */

	Real echarge;		/* electron charge */

	int start_u1flag;	/* what to do for beginning u(1) lattice */
	char start_u1file[MAXFILENAME];
	int save_u1flag;	/* what to do with u(1) lattice at end */
	char save_u1file[MAXFILENAME];

} params;

#endif /* _PARAMS_H */

/* ************************************************************ */

