/*****************  su3_proj.c  (in su3.a) ******************************
*									*
* void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c )	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#define _inline_C_su3_projector( a, b, dest ){\
  int _i,_j;\
  Real _tmp,_tmp2;\
    for(_i=0;_i<3;_i++)for(_j=0;_j<3;_j++){\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].real;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].imag;\
	(dest)->e[_i][_j].real = _tmp + _tmp2;\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].imag;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].real;\
	(dest)->e[_i][_j].imag = _tmp - _tmp2;\
    }\
}
