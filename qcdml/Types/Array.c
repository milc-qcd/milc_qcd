#include <stdlib.h>
#include <cast.h>
#include <Array.h>

#define _mMethodsArray(T) \
\
T ## Array* new_ ## T ## Array ( size_t size ) \
{ \
  T ## Array* self = static_cast ( T ## Array*, \
                          calloc ( 1, sizeof(T ## Array) ) ); \
  if ( size ) \
     self -> elem = static_cast ( T*, \
                         calloc ( size, sizeof(T) ) ); \
  self -> size = size; \
  return self; \
} \
\
void delete_ ## T ## Array ( T ## Array* obj ) \
{ \
  if ( obj ) \
  { \
    obj -> size = 0; \
    if ( obj -> elem ) \
      { \
        free ( obj -> elem ); \
      } \
    free ( obj ); \
  } \
}

_mMethodsArray(double);
_mMethodsArray(float);
_mMethodsArray(int);
_mMethodsArray(long);
_mMethodsArray(short);
_mMethodsArray(unsigned);
_mMethodsArray(char);
_mMethodsArray(bool);
_mMethodsArray(CstringRef);
_mMethodsArray(voidPtr);
_mMethodsArray(doubleComplex);
_mMethodsArray(floatComplex);
