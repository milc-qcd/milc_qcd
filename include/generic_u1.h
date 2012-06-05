#ifndef _GENERIC_U1
#define _GENERIC_U1

/************************ generic_ks.h **********************************
*									*
*  Macros and declarations for generic_u1 routines                      *
*  This header is for codes that call generic_u1 routines               *
*  MIMD version 7 							*
*									*
*/

/* u1link.c */
Real *create_u1_A_field(void);
void destroy_u1_A_field(Real *A);
void u1phase_on(Real charge, Real *A);
void u1phase_off(void);

/* u1plaq.c */
complex u1ploop(void);
void u1plaq(Real *ssplq,Real *stplq);


#endif /* _GENERIC_U1 */
