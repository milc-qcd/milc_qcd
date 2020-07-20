#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define AUTO_STOPTIME -1

#define WILSON 0
#define SYMANZIK 1
#define ZEUTHEN 2

// low-storage schemes
// 3-stage third order
#define INTEGRATOR_LUSCHER 0
// 5-stage fourth order
#define INTEGRATOR_CK 1
// 6-stage fourth order
#define INTEGRATOR_BBB 2

// Runge-Kutta-Munthe-Kaas schemes
// 3-stage third order
#define INTEGRATOR_RKMK3 10
// 4-stage fourth order
#define INTEGRATOR_RKMK4 11
// 6-stage fifth order
#define INTEGRATOR_RKMK5 12
// 13-stage eighth order (Dormand-Prince)
#define INTEGRATOR_RKMK8 13

// adaptive schemes
#define INTEGRATOR_ADAPT_LUSCHER 20
#define INTEGRATOR_ADAPT_BS 21

// Safety factor for adaptive schemes
// to prevent too many rejected steps
#define SAFETY 0.95


/* one step of the flow, here branching into different integrators happens */
#if GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
 || GF_INTEGRATOR==INTEGRATOR_BBB
#define flow_step integrate_RK_2N
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
#define flow_step integrate_RKMK3
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
#define flow_step integrate_RKMK_generic
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER
#define flow_step integrate_adapt_RK_2N
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#define flow_step integrate_adapt_bs
#endif

#endif /* _DEFINES_H */
