#include <xmlWriter.h>
#include <stdlib.h>
#include <cast.h>

// type testStruct
typedef struct testStruct_tag testStruct;

struct testStruct_tag
{
  // scalar
  double        Double;
  float         Float;
  int           Int;
  long          Long;
  short         Short;
  unsigned      Unsigned;
  char          Char;
  bool          Bool;
  CstringRef    Stringg;
  doubleComplex DoubleComplex;
  floatComplex  FloatComplex;

  // array
  doubleArray*        DoubleArray;
  floatArray*         FloatArray;
  intArray*           IntArray;
  longArray*          LongArray;
  shortArray*         ShortArray;
  unsignedArray*      UnsignedArray;
  charArray*          CharArray;
  boolArray*          BoolArray;
  CstringRefArray*    StringArray;
  doubleComplexArray* DoubleComplexArray;
  floatComplexArray*  FloatComplexArray;

};

// initialize type testStruct
void initTestStruct ( testStruct* obj )
{
  unsigned j;
  unsigned size = 3;

  // scalar
  obj -> Double = 3.14159265358979323844;
  obj -> Float  = 2.71828182845904523536;
  obj -> Int    = 123456;
  obj -> Long   = 756789123;
  obj -> Short  = 16000;
  obj -> Unsigned = 126;
  obj -> Char   = 'X';
  obj -> Bool   = true;
  obj -> Stringg = "Hello World!";
  obj -> DoubleComplex = 5 + 4*I;
  obj -> FloatComplex = 10 - 3*I;

  // array
  obj -> DoubleArray = new_doubleArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> DoubleArray -> elem [j] = j;

  obj -> FloatArray = new_floatArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> FloatArray -> elem [j] = j;

  obj -> IntArray = new_intArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> IntArray -> elem [j] = j;

  obj -> LongArray = new_longArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> LongArray -> elem [j] = j;

  obj -> ShortArray = new_shortArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> ShortArray -> elem [j] = j;

  obj -> UnsignedArray = new_unsignedArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> UnsignedArray -> elem [j] = j;

  obj -> CharArray = new_charArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> CharArray -> elem [j] = 'A' + j;

  obj -> BoolArray = new_boolArray ( size );
  for ( j = 0; j < size; ++j )
    obj -> BoolArray -> elem [j] = j % 2;

  obj -> StringArray = new_CstringRefArray ( size );
  obj -> StringArray -> elem [0] = "Array";
  obj -> StringArray -> elem [1] = "of";
  obj -> StringArray -> elem [2] = "Cstring";

  obj -> DoubleComplexArray = new_doubleComplexArray ( size );
  for ( j = 0; j < size; ++j )
    {
      obj -> DoubleComplexArray -> elem [j] = j + I * j;
    }
  obj -> FloatComplexArray = new_floatComplexArray ( size );
  for ( j = 0; j < size; ++j )
    {
      obj -> FloatComplexArray -> elem [j] = j + I * j;
    }

}

// encoder for type testStruct
void testStructEncode ( outStream* os, char const* name, void const* obj )
{
  testStruct const* ts = static_cast ( testStruct const*, obj );

  os -> print ( os, "<%s>", name );

  // invoke encoder for data members
  doubleEncode ( os, "Double", &ts -> Double );
  floatEncode  ( os, "Float",  &ts -> Float );
  intEncode ( os, "Int", &ts -> Int );
  longEncode ( os, "Long", &ts -> Long );
  shortEncode ( os, "Short", &ts -> Short );
  unsignedEncode ( os, "Unsigned", &ts -> Unsigned );
  charEncode   ( os, "Char",   &ts -> Char );
  boolEncode ( os, "Bool", &ts -> Bool );
  CstringRefEncode ( os, "Stringg", &ts -> Stringg );
  doubleComplexEncode ( os, "DoubleComplex", &ts -> DoubleComplex );
  floatComplexEncode ( os, "FloatComplex", &ts -> FloatComplex );

  os -> print ( os, "\n" );
  doubleArrayEncode ( os, "DoubleArray", ts -> DoubleArray );
  os -> print ( os, "\n" );
  floatArrayEncode  ( os, "FloatArray",  ts -> FloatArray );
  os -> print ( os, "\n" );
  intArrayEncode ( os, "IntArray", ts -> IntArray );
  os -> print ( os, "\n" );
  longArrayEncode ( os, "LongArray", ts -> LongArray );
  os -> print ( os, "\n" );
  shortArrayEncode ( os, "ShortArray", ts -> ShortArray );
  os -> print ( os, "\n" );
  unsignedArrayEncode ( os, "UnsignedArray", ts -> UnsignedArray );
  os -> print ( os, "\n" );
  charArrayEncode   ( os, "CharArray",   ts -> CharArray );
  os -> print ( os, "\n" );
  boolArrayEncode ( os, "BoolArray", ts -> BoolArray );
  os -> print ( os, "\n" );
  CstringRefArrayEncode ( os, "StringArray", ts -> StringArray );
  os -> print ( os, "\n" );
  doubleComplexArrayEncode ( os, "DoubleComplexArray", ts -> DoubleComplexArray );
  os -> print ( os, "\n" );
  floatComplexArrayEncode ( os, "FloatComplexArray", ts -> FloatComplexArray );
  os -> print ( os, "\n" );

  os -> print ( os, "</%s>", name );
}

void testFile ( testStruct const* obj )
{
  outStream* os = new_outStreamFile ( stdout );

  printf ( "\noutput to FILE*\n" );
  testStructEncode ( os, "myTestStruct", obj );
  printf ( "\n" );

  delete_outStream ( os );
}

void testString ( testStruct const* obj )
{
  size_t size = 2048;
  String* st = new_String ( size, "" );
  outStream* os = new_outStreamString ( st );

  printf ( "\noutput to String\n" );
  testStructEncode ( os, "myTestStruct", obj );
  printf ( "%s\n", st -> chars );

  delete_String ( st );
  delete_outStream ( os );
}

int main ( void )
{
  testStruct myTestStruct;

  initTestStruct ( &myTestStruct );

  testFile ( &myTestStruct );

  testString ( &myTestStruct );

  return 0;
}
