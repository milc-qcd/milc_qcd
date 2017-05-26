#ifndef QPHIXJ_ARCH_H
#define QPHIXJ_ARCH_H

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 4
#endif

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
#define VECLEN_SP 16 
#define VECLEN_HP 16 
#define VECLEN_DP 8
#endif

#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)
#define VECLEN_HP 8
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_SCALAR_SOURCE) 
#define VECLEN_SP 1
#define VECLEN_DP 1
#endif

#if defined(QPHIX_QPX_SOURCE) 
#define VECLEN_SP 4
#define VECLEN_DP 4
#endif

#if defined(QPHIX_SSE_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 2
#endif

#define COMPRESS false  // Hard wired for now

#endif // QPHIXJ_ARCH_H
