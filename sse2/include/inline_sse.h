#ifndef _INLINE_SSE_H_
#define _INLINE_SSE_H_

#include "sse_m_amat_hwvec.h"
#include "sse_m_amatvec.h"
#include "sse_m_amv_4dir.h"
#include "sse_m_amv_4vec.h"
#include "sse_m_mat_an.h"
#include "sse_m_mat_hwvec.h"
#include "sse_m_mat_na.h"
#include "sse_m_mat_nn.h"
#include "sse_m_matvec.h"
#include "sse_m_mv_s_4dir.h"
#include "sse_su3_proj.h"

/** NOT AVAILABLE, YET **/
/** #include "sse_addvec.h" **/
/** #include "sse_sub4vecs.h" **/
/** #include "sse_s_m_a_vec.h" **/
/** #include "sse_s_m_a_mat.h" **/

/* Corrections by A. Alexandru */
typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((aligned (16)));

static sse_mask _sse_sgn13 __attribute__ ((unused)) ={0x80000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn24 __attribute__ ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn3 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn4 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn2 __attribute__  ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x00000000};

#endif
