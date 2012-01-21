#ifndef _FLAVOR_OPS_H
#define _FLAVOR_OPS_H

#include "../include/gammatypes.h"

/* flavor_ops.c */
void mult_spin_pseudoscalar(field_offset src, field_offset dest ) ;
#if 0
void sym_shift(int dir, field_offset src,field_offset dest) ;
void zeta_shift(int n, int *d, field_offset src, field_offset dest ) ;
void eta_shift(int n, int *d, field_offset src, field_offset dest ) ;


void mult_flavor_vector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudoscalar(field_offset src, field_offset dest ) ;

void mult_spin_vector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest ) ;

/* flavor_ops2.c */

void mult_flavor_vector_field(int mu, int r0[],
			      su3_vector *src, su3_vector *dest ) ;
void mult_flavor_tensor_field(int mu, int nu, int r0[],
			      su3_vector *src, su3_vector *dest ) ;
void mult_flavor_pseudovector_field(int mu, int r0[],
				    su3_vector *src, su3_vector *dest ) ;
void mult_flavor_pseudoscalar_field(int r0[], 
				    su3_vector *src, su3_vector *dest ) ;

void mult_spin_vector_field(int mu, int r0[],
			    su3_vector *src, su3_vector *dest ) ;
void mult_spin_tensor_field(int mu, int nu, int r0[],
			    su3_vector *src, su3_vector *dest ) ;
void mult_spin_pseudovector_field(int mu, int r0[],
				  su3_vector *src, su3_vector *dest ) ;
void mult_spin_pseudoscalar_field(int r0[],
				  su3_vector *src, su3_vector *dest ) ;
#endif

void mult_pion5_field( int r0[], su3_vector *src, su3_vector *dest );
void mult_pion05_field( int r0[], su3_vector *src, su3_vector *dest );
void mult_pioni5_field( int fdir, int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_pionij_field( int fdir, int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_pioni_field( int fdir, int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_pioni0_field( int fdir, int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_pions_field(int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_pion0_field(int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_rhoi_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest );
void mult_rhoi0_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest );
void mult_rhos_field( int fdir,  int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );
void mult_rho0_field( int fdir,  int r0[], su3_vector *src, su3_vector *dest, su3_matrix *links );

void
general_spin_taste_op(enum gammatype spin_index, enum gammatype taste_index, int r0[],
		      su3_vector *dest, su3_vector *src, su3_matrix *links);
void spin_taste_op(int index, int r0[], su3_vector *dest, su3_vector *src);

int spin_taste_index(char *label);
char *spin_taste_label(int index);
#endif /* _FLAVOR_OPS_H */
