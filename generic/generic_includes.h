/************************ generic_includes.h ****************************
*									*
*  This header is included in all codes in this directory               *
*  MIMD version 7 							*
*									*
*/
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#ifdef HAVE_U1
/* for measurements routines */
/* in gauge_stuff.c et al */
#include "../include/generic_u1.h"
#endif
#include "../include/int32type.h"
#include "../include/dirs.h"
#include <lattice.h>
#include <time.h>
#include <string.h>
#ifdef SCHROED_FUN
#include "../include/generic_schroed.h"
#endif
#include "../include/check_malloc.h"

#ifdef GB_BARYON
#ifdef blind
#include "../include/blind_data.h"
#endif
#endif
