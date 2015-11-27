/* QPhiX-MILC Interface */
#ifndef _GENERIC_QPHIX_H
#define _GENERIC_QPHIX_H

#include "../include/config.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/su3.h"
#include "../include/dirs.h"
#include "../include/fermion_links.h"
#include <ks_long_dslash.h>      /* qphix_ks_dslash */
#include <ks_d_congrad_fn.h>     /* congrad */
#include <stdbool.h>

/*! \brief QPhiX environment. */
typedef struct qphix_env
{
    void *ks_src1;
    void *ks_src2;
    void *ks_dest1;
    void *ks_dest2;
    void *gll[2];
    void *gfl[2];
} qphix_env_t;

int
initialize_qphix(void);

/*! \brief Constructor to create the single global qphix_env_t object. */ 
void
create_qphix_env ( int nx
                 , int ny
                 , int nz
                 , int nt
                 , int ncores
                 , int threads_per_core
                 , int min_ct
                 //, int soalen
                 //, int by
                 //, int bz
                 //, bool use_compressed12
              );

/*! \brief Destroy the global mbench object */
void
destroy_qphix_env (void);

/* Declare a global instance of the qphix environment. \fixme - Do away!  */
extern qphix_env_t *mbench;
/* Global flag to indicate if the mbench instance is usable. */
extern bool is_qphix_env_setup;

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
Sfigather_su3vectors_from_lattice (su3_vector *ks, void *ks_spinor, int parity);

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

#if 0
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
