/* QPhiX-MILC Interface */
#ifndef _GENERIC_QPHIX_H
#define _GENERIC_QPHIX_H

#include <qphix.h>
#include "../include/config.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/su3.h"
#include "../include/dirs.h"
#include "../include/fermion_links.h"
#include <stdbool.h>

/* milc_to_qphix_utilities.c */
QPHIX_evenodd_t milc2qphix_parity(int milc_parity);
int qphix2milc_parity(QPHIX_evenodd_t qphix_parity);
QPHIX_status_t initialize_qphix(int precision);
void finalize_qphix(void);

/* d_congrad5_fn_qphix_P.c */

/* Generic mapping */

#if ( QPHIX_PrecisionInt == 1 )

#define create_qphix_L_from_fn_links  create_qphix_F_L_from_fn_links

#else

#define create_qphix_L_from_fn_links  create_qphix_D_L_from_fn_links

#endif

QPHIX_F3_FermionLinksAsqtad *
create_qphix_F_L_from_fn_links (fn_links_t *fn, int parity);

QPHIX_D3_FermionLinksAsqtad *
create_qphix_D_L_from_fn_links (fn_links_t *fn, int parity);

/* map_milc_to_qphix*.c */

/* Generic mapping */

#if ( QPHIX_PrecisionInt == 1 )

#define create_qphix_raw4_G  create_qphix_raw4_F_G
#define create_qphix_raw4_F  create_qphix_raw4_F_F
#define create_qphix_raw_V   create_qphix_raw_F_V
#define create_qphix_raw_D   create_qphix_raw_F_D

#define destroy_qphix_raw4_G  destroy_qphix_raw4_F_G
#define destroy_qphix_raw4_F  destroy_qphix_raw4_F_F
#define destroy_qphix_raw_V   destroy_qphix_raw_F_V
#define destroy_qphix_raw_D   destroy_qphix_raw_F_D

#define create_qphix_L_from_fields    create_qphix_F_L_from_fields

#define create_qphix_raw4_G_from_field create_qphix_raw4_F_G_from_field
#define create_qphix_raw4_F_from_field create_qphix_raw4_F_F_from_field
#define create_qphix_raw_V_from_field  create_qphix_raw_F_V_from_field
#define create_qphix_raw_D_from_field  create_qphix_raw_F_D_from_field

#define create_qphix_F_from_field4  create_qphix_F_F_from_field4
#define create_qphix_G_from_field4  create_qphix_F_G_from_field4
#define create_qphix_V_from_field   create_qphix_F_V_from_field
#define create_qphix_D_from_field   create_qphix_F_D_from_field

#define unload_qphix_raw4_G_to_field unload_qphix_raw4_F_G_to_field
#define unload_qphix_raw4_F_to_field unload_qphix_raw4_F_F_to_field
#define unload_qphix_raw_V_to_field  unload_qphix_raw_F_V_to_field
#define unload_qphix_raw_D_to_field  unload_qphix_raw_F_D_to_field

#define unload_qphix_F_to_field4  unload_qphix_F_F_to_field4
#define unload_qphix_G_to_field4  unload_qphix_F_G_to_field4

#define unload_qphix_V_to_field unload_qphix_F_V_to_field
#define unload_qphix_D_to_field unload_qphix_F_D_to_field

#define map_milc_clov_to_qphix_raw map_milc_clov_to_qphix_raw_F

#else

#define create_qphix_raw4_G  create_qphix_raw4_D_G
#define create_qphix_raw4_F  create_qphix_raw4_D_F
#define create_qphix_raw_V   create_qphix_raw_D_V
#define create_qphix_raw_D   create_qphix_raw_D_D

#define destroy_qphix_raw4_G  destroy_qphix_raw4_D_G
#define destroy_qphix_raw4_F  destroy_qphix_raw4_D_F
#define destroy_qphix_raw_V   destroy_qphix_raw_D_V
#define destroy_qphix_raw_D   destroy_qphix_raw_D_D

#define create_qphix_L_from_fields    create_qphix_D_L_from_fields

#define create_qphix_raw4_G_from_field create_qphix_raw4_D_G_from_field
#define create_qphix_raw4_F_from_field create_qphix_raw4_D_F_from_field
#define create_qphix_raw_V_from_field  create_qphix_raw_D_V_from_field
#define create_qphix_raw_D_from_field  create_qphix_raw_D_D_from_field

#define create_qphix_F_from_field4   create_qphix_D_F_from_field4
#define create_qphix_G_from_field4   create_qphix_D_G_from_field4
#define create_qphix_V_from_field   create_qphix_D_V_from_field
#define create_qphix_D_from_field   create_qphix_D_D_from_field

#define unload_qphix_raw4_G_to_field unload_qphix_raw4_D_G_to_field
#define unload_qphix_raw4_F_to_field unload_qphix_raw4_D_F_to_field
#define unload_qphix_raw_V_to_field  unload_qphix_raw_D_V_to_field
#define unload_qphix_raw_D_to_field  unload_qphix_raw_D_D_to_field


#define unload_qphix_F_to_field4  unload_qphix_D_F_to_field4
#define unload_qphix_G_to_field4  unload_qphix_D_G_to_field4

#define unload_qphix_V_to_field unload_qphix_D_V_to_field
#define unload_qphix_D_to_field unload_qphix_D_D_to_field

#define map_milc_clov_to_qphix_raw map_milc_clov_to_qphix_raw_D

#endif

QPHIX_F3_ColorVector *create_qphix_F_V_from_field(su3_vector *src, int parity);
fsu3_matrix *create_qphix_raw4_F_G_from_field(su3_matrix *links, int parity);

fsu3_vector *create_qphix_raw_F_V_from_field(su3_vector *src, int milc_parity);

void unload_qphix_F_V_to_field(su3_vector *dest, QPHIX_F3_ColorVector* src, 
			int parity);

void destroy_qphix_raw4_F_G(fsu3_matrix *raw);

QPHIX_D3_ColorVector *create_qphix_D_V_from_field(su3_vector *src, int parity);
dsu3_matrix *create_qphix_raw4_D_G_from_field(su3_matrix *links, int parity);

dsu3_vector *create_qphix_raw_D_V_from_field(su3_vector *src, int milc_parity);

void unload_qphix_D_V_to_field(su3_vector *dest, QPHIX_D3_ColorVector* src, 
			int parity);

void destroy_qphix_raw4_D_G(dsu3_matrix *raw);


#if 0
/*! \brief Gather backward long-links for a site, but do not perform adjoint.
 * 
 * Copy of load_longbacklinks (fermion_links_helpers.c), without the adjoint.
 * We do the adjoint on the backlinks in the QPhiX-codegen generated code for 
 * the dslash kernels.
 */
su3_matrix *
load_longbacklinks_without_adjoint (fn_links_t *fn);

/*! \brief Gather backward fat-links for a site, but do not perform the adjoint.
 *
 * Copy of load_fatbacklinks (fermion_links_helpers.c), without the adjoint.
 * We do adjoint on the backlinks in the QPhiX-codegen generated code for the 
 * dslash kernels.
 */
su3_matrix *
load_fatbacklinks_without_adjoint (fn_links_t *fn);

/*! \brief Extract su3vector(spinors) from the lattice.
 *
 * Convenience function to create and array of su3_vectors from the lattice 
 * site array. The argument defines if we want to copy over the source or the
 * destination spinors.
 */
su3_vector *
get_su3_vectors_from_lattice (field_offset ft);


/*! \brief Convert the su3vector array into QPhiX's SOA layout.
 *
 * Get the su3_vectors from the site major lattice and convert them into an
 * array of ks_spinors for qphix
 */
void
get_ks_spinors_from_lattice (su3_vector *ks, void* ks_spinor, int parity);

void
gather_su3vectors_from_lattice (su3_vector *ks, void *ks_spinor, int parity);

/*!
 * Convert MILC's links into QPhiX's Gauge18 type. 
 */
void
get_links_from_lattice (void* fgauge[2], void* lgauge[2], fn_links_t *fn);

/*!
 * Convert MILC's fatlinks into QPhiX's Gauge18 type. 
 */
void
get_fatlinks_from_lattice (void* gauge18s, int parity, fn_links_t *fn);

/*!
 * Convert MILC's longlinks into QPhiX's Gauge type. 
 */
void
get_longlinks_from_lattice (void* gauge18s, int parity, fn_links_t *fn);

/*! Write back ks_spinors to the lattice */
void
set_ks_spinors_into_lattice (su3_vector *ks, void* ks_spinor, int parity);

/*! \brief Public API for dslash_fn_site. */ 
void 
qphix_dslash_fn_site ( field_offset src
                     , field_offset dest
                     , int parity
                     , fn_links_t *fn 
                     );

/*! \brief Wrapper to call QPhiX's CG for MILC. */
int
qphix_ks_congrad ( field_offset src
                 , field_offset dest
                 , Real mass
                 , int niter
                 , int nrestart
                 , Real rsqmin
                 , int prec
                 , int parity
                 , Real *final_rsq
                 , fn_links_t *fn
                 );

/*! \brief Wrapper to call QPhiX's two-src CG for MILC. */
int
qphix_ks_congrad_two_src ( field_offset src1
                         , field_offset src2
                         , field_offset dest1
                         , field_offset dest2
                         , Real mass1
                         , Real mass2
                         , int niter
                         , int nrestart
                         , Real rsqmin
                         , int prec
                         , int parity
                         , Real *final_rsq
                         , fn_links_t *fn
                         );

int
qphix_ks_congrad_site ( field_offset src
                      , field_offset dest
                      , quark_invert_control *qic
                      , Real mass
                      , fn_links_t *fn
                      );                  
#endif
#endif // _GENERIC_QPHIX_H
