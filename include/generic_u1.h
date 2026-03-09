#ifndef _GENERIC_U1
#define _GENERIC_U1

/************************ generic_u1.h **********************************
*									*
*  Macros and declarations for generic_u1 routines                      *
*  This header is for codes that call generic_u1 routines               *
*  MIMD version 7                                                       *
*                                                                       *
*/
/* gauge_stuff_u1.c */
void g_measure_nc_u1(void);

/* u1avlink.c */
void u1avlink(double *sLink, double *tLink, double charge);

/* u1link.c */
/* this can be replaced by create_r_array_field(4);
 * ../generic/field_utilities.c
 */
Real *create_u1_A_field(void);
/* this can be replaced by destroy_r_array_field(4);
 * ../generic/field_utilities.c
 */
void destroy_u1_A_field(Real *A);
void u1phase_on(Real charge, Real *A);
void u1phase_off(void);
void construct_u1_links_from_A(Real charge, Real *A, complex *clink);
/* this can be replaced bycopy_c_array_field(complex *dst, complex *src, 4);
 * ../generic/field_utilities.c
 */
void u1_link_copy_field_to_field(complex *src, complex *dest);
/* this can be replaced by copy_r_array_field(Real *dst, Real *src, 4);
 * ../generic/field_utilities.c
 */
void u1_A_copy_field_to_field(Real *src, Real *dest);

/* u1plaq.c */
complex u1ploop(Real charge);
void u1plaq(Real *ssplq,Real *stplq,Real charge);
void u1Fmunu(Real *ssFmunu, Real *stFmunu);

/* gauge_force_u1_cpu.c */
void gauge_force_u1_cpu( Real eps);
#ifdef USE_GF_GPU
#define gauge_force_u1 gauge_force_u1_cpu
#else
#define gauge_force_u1 gauge_force_u1_cpu
#endif

#endif /* _GENERIC_U1 */
