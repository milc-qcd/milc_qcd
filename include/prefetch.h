#ifndef _PREFETCH_H
#define _PREFETCH_H

/***************************************************************************/
/*                                                                         */
/*          Cache manipulation for a variety of architectures              */
/*                                                                         */
/***************************************************************************/


/* Abbreviations for data types...
   M  SU(3) matrix
   V  SU(3) vector
   W  Wilson vector
   H  Half Wilson vector
*/

#if (defined P3 || defined P4) && defined __GNUC__

#include "../include/prefetch_asm.h"

#else

#if ! defined PREFETCH

/***************************************************************************/
/*              Ignore Cache Manipulation Macros                           */
/*                                                                         */
/***************************************************************************/

#define prefetch_M(a0)                    
#define prefetch_V(a0)                    
#define prefetch_W(a0)                    
#define prefetch_H(a0)                    
#define prefetch_VV(a0,a1)            
#define prefetch_VVV(a0,a1,a2)            
#define prefetch_VVVV(a0,a1,a2,a3)        
#define prefetch_VVVVV(a0,a1,a2,a3,a4)    
#define prefetch_WWW(a0,a1,a2)            
#define prefetch_WWWW(a0,a1,a2,a3)        
#define prefetch_WWWWW(a0,a1,a2,a3,a4)    
#define prefetch_4MVVVV(a0,a1,a2,a3,a4)   
#define prefetch_4MWWWW(a0,a1,a2,a3,a4)   
#define prefetch_4MV4V(a0,a1,a2)          
#define prefetch_4MW4W(a0,a1,a2)          

#else

/***************************************************************************/
/*                  Cache Manipulation Macros                              */
/*                Prefetch via subroutine calls                            */
/*                                                                         */
/***************************************************************************/

void _prefetch_M( su3_matrix * );
void _prefetch_V( su3_vector * );
void _prefetch_W( wilson_vector *);
void _prefetch_H( half_wilson_vector *);
void _prefetch_VV( su3_vector *, su3_vector *);
void _prefetch_VVV( su3_vector *, su3_vector *, su3_vector *);
void _prefetch_VVVV( su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_VVVVV( su3_vector *, su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_WWW( wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_WWWW( wilson_vector *, wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_WWWWW( wilson_vector *, wilson_vector *, wilson_vector *, 
		      wilson_vector *, wilson_vector *);
void _prefetch_4MVVVV( su3_matrix *, su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_4MWWWW( su3_matrix *, wilson_vector *, wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_4MV4V( su3_matrix *, su3_vector *, su3_vector *);
void _prefetch_4MW4W( su3_matrix *, wilson_vector *, wilson_vector *);

#define prefetch_M(a0)                    _prefetch_M(a0)
#define prefetch_V(a0)                    _prefetch_V(a0)
#define prefetch_W(a0)                    _prefetch_W(a0)
#define prefetch_H(a0)                    _prefetch_H(a0)
#define prefetch_VV(a0,a1)                _prefetch_VV(a0,a1)
#define prefetch_VVV(a0,a1,a2)            _prefetch_VVV(a0,a1,a2)
#define prefetch_VVVV(a0,a1,a2,a3)        _prefetch_VVVV(a0,a1,a2,a3)
#define prefetch_VVVVV(a0,a1,a2,a3,a4)    _prefetch_VVVVV(a0,a1,a2,a3,a4)
#define prefetch_WWW(a0,a1,a2)            _prefetch_WWW(a0,a1,a2)
#define prefetch_WWWW(a0,a1,a2,a3)        _prefetch_WWWW(a0,a1,a2,a3)
#define prefetch_WWWWW(a0,a1,a2,a3,a4)    _prefetch_WWWWW(a0,a1,a2,a3,a4)
#define prefetch_4MVVVV(a0,a1,a2,a3,a4)   _prefetch_4MVVVV(a0,a1,a2,a3,a4)
#define prefetch_4MWWWW(a0,a1,a2,a3,a4)   _prefetch_4MWWWW(a0,a1,a2,a3,a4)
#define prefetch_4MV4V(a0,a1,a2)          _prefetch_4MV4V(a0,a1,a2)
#define prefetch_4MW4W(a0,a1,a2)          _prefetch_4MW4W(a0,a1,a2)

#endif /* NOPREFETCH */

#endif /* P4 or P3 */

#endif /* _PREFETCH_H */

