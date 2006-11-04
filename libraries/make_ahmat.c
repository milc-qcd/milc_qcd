/*****************  make_ahmat.c  (in su3.a) ****************************
*									*
* void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3)	*
* take the traceless and anti_hermitian part of an su3 matrix 		*
* and compress it 							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 ) {
Real temp;
	
	temp =
	    (m3->e[0][0].imag + m3->e[1][1].imag + m3->e[2][2].imag)*0.33333333333333333;
	ah3->m00im = m3->e[0][0].imag - temp;
	ah3->m11im = m3->e[1][1].imag - temp;
	ah3->m22im = m3->e[2][2].imag - temp;
	ah3->m01.real = (m3->e[0][1].real - m3->e[1][0].real)*0.5;
	ah3->m02.real = (m3->e[0][2].real - m3->e[2][0].real)*0.5;
	ah3->m12.real = (m3->e[1][2].real - m3->e[2][1].real)*0.5;
	ah3->m01.imag = (m3->e[0][1].imag + m3->e[1][0].imag)*0.5;
	ah3->m02.imag = (m3->e[0][2].imag + m3->e[2][0].imag)*0.5;
	ah3->m12.imag = (m3->e[1][2].imag + m3->e[2][1].imag)*0.5;

}/* make_anti_hermitian */

#else
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 ) {
Real temp,temp2;
	
	temp =
	    (m3->e[0][0].imag + m3->e[1][1].imag);
	temp2 = temp + m3->e[2][2].imag;
	temp = temp2*0.333333333333333333;
	ah3->m00im = m3->e[0][0].imag - temp;
	ah3->m11im = m3->e[1][1].imag - temp;
	ah3->m22im = m3->e[2][2].imag - temp;
	temp = m3->e[0][1].real - m3->e[1][0].real; ah3->m01.real = temp*0.5;
	temp = m3->e[0][2].real - m3->e[2][0].real; ah3->m02.real = temp*0.5;
	temp = m3->e[1][2].real - m3->e[2][1].real; ah3->m12.real = temp*0.5;
	temp = m3->e[0][1].imag + m3->e[1][0].imag; ah3->m01.imag = temp*0.5;
	temp = m3->e[0][2].imag + m3->e[2][0].imag; ah3->m02.imag = temp*0.5;
	temp = m3->e[1][2].imag + m3->e[2][1].imag; ah3->m12.imag = temp*0.5;

}/* make_anti_hermitian */
#endif /*end ifdef FAST */
