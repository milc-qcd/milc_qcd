#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define AUTO_STOPTIME -1

#define WILSON 0
#define SYMANZIK 1

// low-storage schemes
// 3-stage third order
#define INTEGRATOR_LUSCHER 0
// 5-stage fourth order
#define INTEGRATOR_CK 1

// Runge-Kutta-Munthe-Kaas schemes
// 3-stage third order
#define INTEGRATOR_RKMK3 10
// 4-stage fourth order
#define INTEGRATOR_RKMK4 11
// 6-stage fifth order
#define INTEGRATOR_RKMK5 12
// 13-stage eighth order (Dormand-Prince)
#define INTEGRATOR_RKMK8 13



#endif /* _DEFINES_H */
