#include <xmlWriter.h>
#include <cast.h>

// template
#define _mDefineEncode(T,FMT) \
void T ## Encode ( outStream* os, char const* name, void const* obj ) \
{ \
  os -> print ( os, "<%s>" #FMT "</%s>", \
		name, \
		*static_cast ( T *, obj ), \
		name ); \
}

#define _mDefineArrayEncode(T) \
void T ## ArrayEncode ( outStream* os, char const* name, void const* obj ) \
{ \
  T ## Array const* array = static_cast ( T ## Array const*, obj ); \
  size_t j; \
\
  os -> print ( os, "<%s>", name ); \
\
  for ( j = 0; j < array -> size; ++j ) \
    T ## Encode ( os, "elem", &(array -> elem [j] ) ); \
\
  os -> print ( os, "</%s>", name ); \
}

// scalar types
_mDefineEncode(double, %.17g);
_mDefineEncode(float, %.8g);
_mDefineEncode(int, %d);
_mDefineEncode(long, %ld);
_mDefineEncode(short, %d);
_mDefineEncode(unsigned, %u);
_mDefineEncode(char, %c);
_mDefineEncode(bool, %u);
_mDefineEncode(CstringRef, %s);

// Array types
_mDefineArrayEncode(double);
_mDefineArrayEncode(float);
_mDefineArrayEncode(int);
_mDefineArrayEncode(long);
_mDefineArrayEncode(short);
_mDefineArrayEncode(unsigned);
_mDefineArrayEncode(char);
_mDefineArrayEncode(bool);
_mDefineArrayEncode(CstringRef);
_mDefineArrayEncode(doubleComplex);
_mDefineArrayEncode(floatComplex);

// specializations

void doubleComplexEncode ( outStream* os, char const* name, void const* obj )
{
  os -> print ( os, "<%s><re>%.17g</re><im>%.17g</im></%s>",
		name,
		creal (*static_cast ( doubleComplex*, obj )),
		cimag (*static_cast ( doubleComplex*, obj )),
		name );
}

void floatComplexEncode ( outStream* os, char const* name, void const* obj )
{
  os -> print ( os, "<%s><re>%.17g</re><im>%.8g</im></%s>",
		name,
		creal (*static_cast ( floatComplex*, obj )),
		cimag (*static_cast ( floatComplex*, obj )),
		name );
}
