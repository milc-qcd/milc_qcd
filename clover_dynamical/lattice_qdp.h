#ifndef _LATTICE_QDP_H
#define _LATTICE_QDP_H

#ifdef HAVE_QDP

#include <qdp.h>

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN QDP_Shift shiftdirs[4];
EXTERN QDP_ShiftDir shiftfwd[4], shiftbck[4];

#endif
#endif /* _LATTICE_QDP_H */
