#include <outStream.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdint.h>
#include <cast.h>

static void printToFile ( outStream* self, const char *format, ... )
{
  va_list ap;
  va_start ( ap, format );
  vfprintf ( self -> _fp, format, ap );
  va_end ( ap );
}

// private prototype
int va_printfToString ( String* dest, char const* format, va_list ap );

static void printToString ( struct outStream_tag* self,
			    const char *format, ... )
{
  va_list ap;
  size_t size = 256;
  String* tmp = new_String ( size, "" );

  va_start ( ap, format );
  va_printfToString ( tmp, format, ap );
  va_end ( ap );

  StringConcat ( self -> _st, tmp );
  delete_String ( tmp );
}

outStream* new_outStreamFile ( FILE* fp )
{
  outStream* self = static_cast ( outStream*,
				  calloc ( 1, sizeof(outStream) ) );
  self -> _fp = fp;
  self -> print = printToFile;

  return self;
}

outStream* new_outStreamString ( String* str )
{
  outStream* self = static_cast ( outStream*,
				  calloc ( 1, sizeof(outStream) ) );
  self -> _st = str;
  self -> print = printToString;

  return self;
}

void delete_outStream ( outStream* obj )
{
  obj -> _fp = NULL;
  obj -> _st = NULL;

  free ( obj );
}
