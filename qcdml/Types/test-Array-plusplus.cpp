#include <stdio.h>
#include <Array.h>

#define _mPrint(T,F,...) \
void T ## Print ( T const* obj ) \
{ \
  printf ( F, __VA_ARGS__ ); \
}

#define _mArrayPrint(T) \
void T ## ArrayPrint ( T ## Array const* obj ) \
{ \
  unsigned j; \
  printf ( "[ " ); \
  for ( j = 0; j < obj -> size; ++j ) \
    { \
      printf ( " " ); \
      T ## Print ( &(obj -> elem [j]) ); \
    } \
  printf ( "]\n" ); \
}


#define _mmain(T) \
int main ( ) \
{ \
  int j; \
  int size = 10; \
  T ## Array* array = new_## T ## Array ( size ); \
\
  for ( j = 0; j < size; ++j ) \
    array -> elem [j] = j; \
\
  T ## ArrayPrint ( array ); \
\
  delete_ ## T ## Array ( array ); \
\
  return 0; \
}

_mPrint(float,"%g", *obj);
_mPrint(floatComplex,"[ %g, %g ]", real(*obj), imag(*obj));

_mArrayPrint(float);
_mArrayPrint(floatComplex);

_mmain(floatComplex);

