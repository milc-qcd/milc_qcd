/*****************  su3_s_m_a_hwvec_proj.c  (in su3.a) *******************
*									 *
* void scalar_mult_add hwvec_proj( su3_matrix *a, half_wilson_vector *b, * 
*    half_wilson_vector *c, Real *s, su3_matrix *d )	                 *
*  D  <- A + B[0] * C[0]_adj * s[0] + B[1] *C[1]_adj * s[1]              *
*/
#define _inline_C_scalar_mult_add_hwvec_proj(aa,bb,cc,ss,dd) \
{\
    int _i,_j;\
    Real tmp0,tmp1;\
    for(_i=0;_i<3;_i++)for(_j=0;_j<3;_j++){\
	tmp1 = (bb)->h[0].c[_i].real * (cc)->h[0].c[_j].real;\
	tmp0 = (bb)->h[0].c[_i].imag * (cc)->h[0].c[_j].imag;\
	(dd)->e[_i][_j].real = (aa)->e[_i][_j].real + (tmp0 + tmp1)*(ss)[0];\
\
	tmp1 = (bb)->h[0].c[_i].real * (cc)->h[0].c[_j].imag;\
	tmp0 = (bb)->h[0].c[_i].imag * (cc)->h[0].c[_j].real;\
	(dd)->e[_i][_j].imag = (aa)->e[_i][_j].imag + (tmp0 - tmp1)*(ss)[0];\
\
	tmp1 = (bb)->h[1].c[_i].real * (cc)->h[1].c[_j].real;\
	tmp0 = (bb)->h[1].c[_i].imag * (cc)->h[1].c[_j].imag;\
	(dd)->e[_i][_j].real += (tmp0 + tmp1)*(ss)[1];\
\
	tmp1 = (bb)->h[1].c[_i].real * (cc)->h[1].c[_j].imag;\
	tmp0 = (bb)->h[1].c[_i].imag * (cc)->h[1].c[_j].real;\
	(dd)->e[_i][_j].imag += (tmp0 - tmp1)*(ss)[1];\
    }\
}
