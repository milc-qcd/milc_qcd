#ifndef _PREFETCH_ASM_H
#define _PREFETCH_ASM_H

/* Macros for prefetching when inline assembly coding is available */
/* So far we have this only for GNU gcc on the P3 and P4 */

#if defined P3

/***************************************************************************/
/*               GNU C Cache Manipulation Macros                           */
/*            Assumes 8 byte alignment and 32 byte cache line              */
/*                Appropriate for Intel P3 Processor                       */
/***************************************************************************/

#define _pftch_M(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+64)))

#define _pftch_V(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+16)))

#define _pftch_W(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+64)), \
                      "m" (*(((char*)(a))+88)))

#define _pftch_H(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+40)))

#define _pftch_4V(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+64)), \
                      "m" (*(((char*)(a))+88)))

#define _pftch_4M(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3  \n\t"  \
                      "prefetchnta %4  \n\t"  \
                      "prefetchnta %5  \n\t"  \
                      "prefetchnta %6  \n\t"  \
                      "prefetchnta %7  \n\t"  \
                      "prefetchnta %8  \n\t"  \
                      "prefetchnta %9" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+64)), \
                      "m" (*(((char*)(a))+96)), \
                      "m" (*(((char*)(a))+128)), \
                      "m" (*(((char*)(a))+160)), \
                      "m" (*(((char*)(a))+192)), \
                      "m" (*(((char*)(a))+224)), \
                      "m" (*(((char*)(a))+256)), \
                      "m" (*(((char*)(a))+280)))

#define _pftch_4W(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3  \n\t"  \
                      "prefetchnta %4  \n\t"  \
                      "prefetchnta %5  \n\t"  \
                      "prefetchnta %6  \n\t"  \
                      "prefetchnta %7  \n\t"  \
                      "prefetchnta %8  \n\t"  \
                      "prefetchnta %9  \n\t"  \
                      "prefetchnta %10 \n\t"  \
                      "prefetchnta %11 \n\t"  \
                      "prefetchnta %12" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+32)), \
                      "m" (*(((char*)(a))+64)), \
                      "m" (*(((char*)(a))+96)), \
                      "m" (*(((char*)(a))+128)), \
                      "m" (*(((char*)(a))+160)), \
                      "m" (*(((char*)(a))+192)), \
                      "m" (*(((char*)(a))+224)), \
                      "m" (*(((char*)(a))+256)), \
                      "m" (*(((char*)(a))+288)), \
                      "m" (*(((char*)(a))+320)), \
                      "m" (*(((char*)(a))+352)), \
                      "m" (*(((char*)(a))+376)))

#elif defined P4 

/***************************************************************************/
/*               GNU C Cache Manipulation Macros                           */
/*            Assumes 8 byte alignment and 128 byte fetch result           */
/*                Appropriate for Intel P4 Processor                       */
/***************************************************************************/

#define _pftch_M(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  " \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+64)))


#define _pftch_V(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+16)))

#define _pftch_W(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+88)))

#define _pftch_H(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+40)))

#define _pftch_4V(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+88)))

#define _pftch_4M(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+128)), \
                      "m" (*(((char*)(a))+256)), \
                      "m" (*(((char*)(a))+280)))

#define _pftch_4W(a) \
__asm__ __volatile__ ("prefetchnta %0  \n\t"  \
                      "prefetchnta %1  \n\t"  \
                      "prefetchnta %2  \n\t"  \
                      "prefetchnta %3" \
                      : \
                      : \
                      "m" (*(((char*)(a)))), \
                      "m" (*(((char*)(a))+128)), \
                      "m" (*(((char*)(a))+256)), \
                      "m" (*(((char*)(a))+376)))

#else

CREATE YOUR OWN DEFINITION HERE!

#endif  /* P3 or P4 */

#define prefetch_M(a0) \
               _pftch_M(a0)

#define prefetch_V(a0) \
               _pftch_V(a0)

#define prefetch_W(a0) \
               _pftch_W(a0)

#define prefetch_H(a0) \
               _pftch_H(a0)

#define prefetch_VV(a0,a1) \
               _pftch_V(a0); \
               _pftch_V(a1)

#define prefetch_VVV(a0,a1,a2) \
               _pftch_V(a0); \
               _pftch_V(a1); \
               _pftch_V(a2)

#define prefetch_VVVV(a0,a1,a2,a3) \
               _pftch_V(a0); \
               _pftch_V(a1); \
               _pftch_V(a2); \
               _pftch_V(a3)

#define prefetch_VVVVV(a0,a1,a2,a3,a4) \
               _pftch_V(a0); \
               _pftch_V(a1); \
               _pftch_V(a2); \
               _pftch_V(a3); \
               _pftch_V(a4);

#define prefetch_WWW(a0,a1,a2) \
               _pftch_W(a0); \
               _pftch_W(a1); \
               _pftch_W(a2)

#define prefetch_WWWW(a0,a1,a2,a3) \
               _pftch_W(a0); \
               _pftch_W(a1); \
               _pftch_W(a2); \
               _pftch_W(a3)

#define prefetch_WWWWW(a0,a1,a2,a3,a4) \
               _pftch_W(a0); \
               _pftch_W(a1); \
               _pftch_W(a2); \
               _pftch_W(a3); \
               _pftch_W(a4);

#define prefetch_4MVVVV(a0,a1,a2,a3,a4) \
               _pftch_4M(a0); \
               _pftch_V(a1); \
               _pftch_V(a2); \
               _pftch_V(a3); \
               _pftch_V(a4);

#define prefetch_4MWWWW(a0,a1,a2,a3,a4) \
               _pftch_4M(a0); \
               _pftch_W(a1); \
               _pftch_W(a2); \
               _pftch_W(a3); \
               _pftch_W(a4);

#define prefetch_4MV4V(a0,a1,a2) \
               _pftch_4M(a0); \
               _pftch_V(a1); \
               _pftch_4V(a2);

#define prefetch_4MW4W(a0,a1,a2) \
               _pftch_4M(a0); \
               _pftch_W(a1); \
               _pftch_4W(a2);

#endif /* _PREFETCH_ASM_H */
