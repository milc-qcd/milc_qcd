#ifndef __xmlWriter_h__
#define __xmlWriter_h__

#include <Types.h>
#include <Array.h>
#include <outStream.h>

#define _mPrototypeEncode(T) \
void T ## Encode ( outStream* os, char const* name, void const* obj )

#ifdef __cplusplus
extern "C" {
#endif

// scalar types
_mPrototypeEncode(double);
_mPrototypeEncode(float);
_mPrototypeEncode(int);
_mPrototypeEncode(long);
_mPrototypeEncode(short);
_mPrototypeEncode(unsigned);
_mPrototypeEncode(char);
_mPrototypeEncode(bool);
_mPrototypeEncode(CstringRef);
_mPrototypeEncode(doubleComplex);
_mPrototypeEncode(floatComplex);

// Array types
_mPrototypeEncode(doubleArray);
_mPrototypeEncode(floatArray);
_mPrototypeEncode(intArray);
_mPrototypeEncode(longArray);
_mPrototypeEncode(shortArray);
_mPrototypeEncode(unsignedArray);
_mPrototypeEncode(charArray);
_mPrototypeEncode(boolArray);
_mPrototypeEncode(CstringRefArray);
_mPrototypeEncode(doubleComplexArray);
_mPrototypeEncode(floatComplexArray);

#ifdef __cplusplus
}
#endif
#endif
