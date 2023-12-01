#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define AUTO_STOPTIME -1

#define NO_REPORT 0
#define OLDREPORT 1
#define NEWREPORT 2

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
// 3-stage third order
#define INTEGRATOR_CF3 3

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
#define INTEGRATOR_ADAPT_CF3 22

// Safety factor for adaptive schemes
// to prevent too many rejected steps
#define SAFETY 0.95

// 1-stage first order
#define INTEGRATOR_EULER 99

/* one step of the flow, here branching into different integrators happens */
#if GF_INTEGRATOR==INTEGRATOR_EULER || \
 		GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
 || GF_INTEGRATOR==INTEGRATOR_BBB || GF_INTEGRATOR==INTEGRATOR_CF3
#define flow_step integrate_RK_2N
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
#define flow_step integrate_RKMK3
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
#define flow_step integrate_RKMK_generic
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
#define flow_step integrate_adapt_RK_2N
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#define flow_step integrate_adapt_bs
#endif

// for tuning of the coefficients of generic third-order and adaptive
//#define READ_CF3_FROM_FILE
//#define READ_ADPT_CF3_FROM_FILE

// dump lattice in double precision for debugging purposes
// BE CAREFULL with this
//#define DEBUG_FIELDS

#ifdef DEBUG_FIELDS
#define REPACK_TO_DOUBLE
#ifdef REPACK_TO_DOUBLE
#define MATRIX_TYPE dsu3_matrix
#else
#define MATRIX_TYPE fsu3_matrix
#endif
#endif

#ifdef BLOCKING

#define MAX_BLOCK_STRIDE 4
/* defines for 2nd and 4th nearest neighbor (BLOCKING) stuff
 * NOTE: INCOMPATIBLE with the 3n gathers due to the NAIK term */
#define X2UP 8
#define Y2UP 9
#define Z2UP 10
#define T2UP 11
#define T2DOWN 12
#define Z2DOWN 13
#define Y2DOWN 14
#define X2DOWN 15

#define OPP_2_DIR(dir) (23-(dir))
#define DIR2(dir) ((dir)+8)
#define FORALL2UPDIR(dir) for(dir=X2UP; dir<=T2UP; dir++)

#define X4UP 16
#define Y4UP 17
#define Z4UP 18
#define T4UP 19
#define T4DOWN 20
#define Z4DOWN 21
#define Y4DOWN 22
#define X4DOWN 23

#define OPP_4_DIR(dir) (39-(dir))
#define DIR4(dir) ((dir)+16)
#define FORALL4UPDIR(dir) for(dir=X4UP; dir<=T4UP; dir++)

// blocked condition in volume for loops
#define IF_BLOCKED(s, stride) \
	if ( s-> x % stride ==0 && s-> y % stride ==0 &&s-> z % stride ==0 )
#else
#define IF_BLOCKED(s, stride)
#endif

#ifdef REGIONS

// specify regions
#define FULLVOL 0
#define BOUNDARY 1
#define BULK 2
#define LOWER_BULK 3
#define ACTIVE 4
#define LOWER_BOUNDARY 5
#define UPPER_BOUNDARY 6

// define boundaries of bulk region
#define LWR_BDRY ( 0 )
#ifdef HALF_LATTICE_TEST
// #define UPR_BDRY ( ( 1 ) )
#define UPR_BDRY ( ( nt / 2 ) )
#else
#define UPR_BDRY ( ( nt / 2 ) )
#endif
#define ACTV_VOL_LEN ( UPR_BDRY - LWR_BDRY + 1 > nt ? nt : UPR_BDRY - LWR_BDRY + 1 )
#define LWR_BULK_LEN ( UPR_BDRY - LWR_BDRY - 0 > nt ? nt : UPR_BDRY - LWR_BDRY - 0 )
#define INR_BULK_LEN ( UPR_BDRY - LWR_BDRY - 1 > nt ? nt : UPR_BDRY - LWR_BDRY - 1 )
#define BOUNDARY_COND(s) ( (s-> t == LWR_BDRY || s -> t == UPR_BDRY ) )
#define BULK_COND(s) ( (s-> t > LWR_BDRY && s -> t < UPR_BDRY ) )
#define LOWER_BULK_COND(s) ( (s-> t >= LWR_BDRY && s -> t < UPR_BDRY ) )
#define ACTIVE_COND(s) ( ( s-> t >= LWR_BDRY && s -> t <= UPR_BDRY) )
#define LWR_BDRY_COND(s) ( s-> t == LWR_BDRY )
#define UPR_BDRY_COND(s) ( s-> t == UPR_BDRY )

#define IF_BOUNDARY(s) if BOUNDARY_COND(s)
#define IF_LWR_BDRY(s) if ( s-> t == LWR_BDRY )
#define IF_UPR_BDRY(s) if ( s-> t == UPR_BDRY )
#define IF_BULK(s) if BULK_COND(s)
#define IF_LOWER_BULK(s)  if LOWER_BULK_COND(s)
#define IF_ACTIVE(s) if ( s-> t >= LWR_BDRY && s -> t <= UPR_BDRY)
#define IF_INACTIVE(s) if ( s -> t > UPR_BDRY)

// region condition in volume for loops
#define IF_REGION(s, region) \
 	if ( region == FULLVOL \
 	|| (region == BOUNDARY && BOUNDARY_COND(s) ) \
 	|| (region == BULK && BULK_COND(s) ) \
 	|| (region == LOWER_BULK && LOWER_BULK_COND(s) ) \
 	|| (region == ACTIVE && ACTIVE_COND(s) ) \
 	|| (region == LOWER_BOUNDARY && LWR_BDRY_COND(s) ) \
 	|| (region == UPPER_BOUNDARY && UPR_BDRY_COND(s) ) \
 	)

#define BOUNDARY_COND_LINK(s, dir) \
 	( ( s-> t == LWR_BDRY || s -> t == UPR_BDRY ) && dir < TUP )
#define ACTIVE_COND_LINK(s, dir) \
 	( ( s-> t >= LWR_BDRY ) \
 	&& ( s -> t < UPR_BDRY || ( s -> t == UPR_BDRY && dir < TUP ) ) ) 
#define BULK_COND_LINK(s, dir) \
 	( ( s-> t > LWR_BDRY ) \
 	&& ( s -> t < UPR_BDRY - 1 || ( s -> t == UPR_BDRY - 1 && dir < TUP ) ) )
#define LWR_BDRY_COND_LINK(s, dir) \
 	( s-> t == LWR_BDRY && dir < TUP )
#define UPR_BDRY_COND_LINK(s, dir) \
 	( s-> t == UPR_BDRY && dir < TUP )
#define IF_LINK_IN_REGION(s, dir, region) \
 	if ( region == FULLVOL \
 	|| (region == BOUNDARY && BOUNDARY_COND_LINK(s, dir) ) \
 	|| (region == BULK && BULK_COND_LINK(s, dir) ) \
 	|| (region == ACTIVE && ACTIVE_COND_LINK(s, dir) ) \
 	|| (region == LOWER_BOUNDARY && LWR_BDRY_COND_LINK(s,dir) ) \
 	|| (region == UPPER_BOUNDARY && UPR_BDRY_COND_LINK(s,dir) ) \
 	)	

#define FORALLDIRSUP(dir, region) \
 	FORALLUPDIR(dir) \
 	if ( ( region == FULLVOL || region == BULK ) \
 	|| ( region == BOUNDARY && dir < TUP ) ) 

#else
#define IF_REGION(s, region)
#define FORALLDIRSUP(dir, region) \
	FORALLUPDIR(dir)
#define IF_LINK_IN_REGION(s, dir, region)
#endif

#ifdef SPHALERON
#ifndef REGIONS
BOMB THE COMPILE
#endif

// drop or keep time links in axial gauge
// #define DROP_TIME_LINKS

#define N_HALF 3
#define N_LAST_FLOW 3
#endif

#endif /* _DEFINES_H */
