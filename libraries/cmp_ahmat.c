/*****************  cmp_ahmat.c  (in su3.a) *****************************
*									*
*  Make an anti_hermitmat (anti Hermitian matrix in compressed form)	*
*  from an SU3 matrix (3x3 complex matrix). 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void compress_anti_hermitian( su3_matrix *mat_su3,
	anti_hermitmat *mat_antihermit ) {
	mat_antihermit->m00im=mat_su3->e[0][0].imag;
	mat_antihermit->m11im=mat_su3->e[1][1].imag;
	mat_antihermit->m22im=mat_su3->e[2][2].imag;
	mat_antihermit->m01.real=mat_su3->e[0][1].real;
	mat_antihermit->m02.real=mat_su3->e[0][2].real;
	mat_antihermit->m12.real=mat_su3->e[1][2].real;
	mat_antihermit->m01.imag=mat_su3->e[0][1].imag;
	mat_antihermit->m02.imag=mat_su3->e[0][2].imag;
	mat_antihermit->m12.imag=mat_su3->e[1][2].imag;
}
