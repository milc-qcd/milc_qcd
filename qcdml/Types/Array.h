#ifndef __Array_h__
#define __Array_h__

#include <stddef.h>
#include <Types.h>

#define _mPrototypeArray(T) \
typedef struct T ## Array_tag T ## Array; \
\
struct T ## Array_tag \
{ \
  size_t size; \
  T*     elem; \
}; \
\
T ## Array* new_ ## T ## Array ( size_t size ); \
\
void delete_ ## T ## Array ( T ## Array* obj )

#ifdef __cplusplus
extern "C" {
#endif

_mPrototypeArray(double);
_mPrototypeArray(float);
_mPrototypeArray(int);
_mPrototypeArray(long);
_mPrototypeArray(short);
_mPrototypeArray(unsigned);
_mPrototypeArray(char);
_mPrototypeArray(bool);
_mPrototypeArray(CstringRef);
_mPrototypeArray(voidPtr);
_mPrototypeArray(doubleComplex);
_mPrototypeArray(floatComplex);

#ifdef __cplusplus
}
#endif

#undef _mPrototypeArray

#endif
