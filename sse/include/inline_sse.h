#ifndef _INLINE_SSE_H_
#define _INLINE_SSE_H_

/* These assembly inline versions work only with the gnu C compiler on
   P3 or P4 processors */

#include "inline_headers.h"

#if defined SSE_INLINE

#define mult_su3_nn(...) _inline_sse_mult_su3_nn(__VA_ARGS__)
#define mult_su3_na(...) _inline_sse_mult_su3_na(__VA_ARGS__)
#define mult_su3_an(...) _inline_sse_mult_su3_an(__VA_ARGS__)
#define mult_su3_mat_vec(...) _inline_sse_mult_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec(...) _inline_sse_mult_adj_su3_mat_vec(__VA_ARGS__)
#define mult_su3_mat_vec_sum_4dir(...) _inline_sse_mult_su3_mat_vec_sum_4dir(__VA_ARGS__)
#define mult_adj_su3_mat_vec_4dir(...) _inline_sse_mult_adj_su3_mat_vec_4dir(__VA_ARGS__)
#define mult_adj_su3_mat_4vec(...) _inline_sse_mult_adj_su3_mat_4vec(__VA_ARGS__)
#define su3_projector(...) _inline_sse_su3_projector(__VA_ARGS__)
#define mult_su3_mat_hwvec(...) _inline_sse_mult_su3_mat_hwvec(__VA_ARGS__)
#define mult_adj_su3_mat_hwvec(...) _inline_sse_mult_adj_su3_mat_hwvec(__VA_ARGS__)
#define sub_four_su3_vecs(...) _inline_sse_sub_four_su3_vecs(__VA_ARGS__)
#define add_su3_vector(...) _inline_sse_add_su3_vector(__VA_ARGS__)
#define scalar_mult_add_su3_vector(a,b,c,d) {Real _temp = c; _inline_sse_scalar_mult_add_su3_vector(a,b,_temp,d);}
#define scalar_mult_add_su3_matrix(a,b,c,d) {Real _temp = c; _inline_sse_scalar_mult_add_su3_matrix(a,b,_temp,d);}

#endif

/* Corrections by A. Alexandru */
typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((__aligned__ (16)));

static sse_mask _sse_sgn13 __attribute__ ((__aligned__ (16), __unused__)) ={0x80000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn24 __attribute__ ((__aligned__ (16), __unused__)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn3 __attribute__  ((__aligned__ (16), __unused__)) ={0x00000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn4 __attribute__  ((__aligned__ (16), __unused__)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};


#endif
