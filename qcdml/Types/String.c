#include <String.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <cast.h>

void dumpString ( FILE* out, char const* name, String const* self )
{
  fprintf ( out, "dumpString: %s (%p)\n", name, self );
  if ( self )
    {
      fprintf ( out, "\tsize  = %d\n", self -> size );
      fprintf ( out, "\talloc = %d\n", self -> alloc );
      fprintf ( out, "\tchars = %p",   self -> chars );
      if ( self -> chars )
	{
	  size_t len = strlen ( self -> chars );
	  if ( self -> size != len )
	    fprintf ( out, " ERROR inconsistent size" );
	  fprintf ( out, " ( %d chars )\n", len );
	  fprintf ( out, "\t\"%s\"\n", self -> chars );
	}
      fprintf ( out, "\n" );
    }
}

size_t StringSetCapacity ( String* self, size_t newSize )
{
  if ( self -> chars )
    realloc ( self -> chars, newSize );
  else
    self -> chars = static_cast ( char*, calloc ( newSize, sizeof(char) ) );

  if ( self -> chars )
    {
      self -> alloc = newSize;
    }
  else
    {
      // realloc failed
      self -> alloc = 0;
      self -> size  = 0;
      fprintf ( stderr,
		"StringSetCapacity realloc %d failed\n",
		newSize );
    }
  return self -> size;
}


String* new_String ( size_t size, char const* s )
{
  String* self = static_cast ( String*, calloc ( 1, sizeof(String) ) );
  size_t len = strlen ( s );
  size_t alloc = 1 + size;

  if ( alloc < 15 ) // tune for performance
    alloc = 16;

  if ( len >= alloc )
    alloc = len + 1;

  self -> chars = 0;
  StringSetCapacity ( self, alloc );

  strncpy ( self -> chars, s, len + 1 ); // copy terminating null
  self -> size  = len;

  return self;
}

void delete_String ( String* self )
{
  if ( self )
    {
      if ( self -> chars )
	{
	  free ( self -> chars );
	  self -> chars = 0;
	}

      self -> alloc = 0;
      self -> size  = 0;

      free ( self );
    }
}

String* StringClone ( String const* self )
{
  String* clone = static_cast ( String*, calloc ( 1, sizeof(String) ) );

  clone -> size  = self -> size;
  clone -> chars = 0;

  StringSetCapacity ( clone, self -> alloc );

  strncpy ( clone -> chars, self -> chars, 1 + self -> size );

  return clone;
}

void StringCopy ( String* self, String const* src )
{
  if ( self -> alloc <= src -> size )
    {
      StringSetCapacity ( self, 1 + src -> size );
    }

  self -> size = src -> size;
  strncpy ( self -> chars, src -> chars, 1 + src -> size );
}

void StringCopyChars ( String* self, char const* src )
{
  size_t src_size = strlen ( src );

  if ( self -> alloc <= src_size )
    {
      StringSetCapacity ( self, 1 + src_size );
    }

  self -> size = src_size;
  strncpy ( self -> chars, src, 1 + src_size );
}

void StringConcat ( String* self, String const* src )
{
  size_t new_size = self -> size + src -> size;

  if ( self -> alloc <= new_size )
    {
      StringSetCapacity ( self, 1 + new_size );
    }

  self -> size = new_size;
  strncat ( self -> chars, src -> chars, 1 + src -> size );
}

void StringConcatChars ( String* self, char const* src )
{
  size_t src_size = strlen ( src );
  size_t new_size = self -> size + src_size;

  if ( self -> alloc <= new_size )
    {
      StringSetCapacity ( self, 1 + new_size );
    }

  self -> size = new_size;
  strncat ( self -> chars, src, 1 + src_size );
}

int va_printfToString ( String* dest, char const* format, va_list ap )
{
  int ret;
  va_list aq;

  va_copy ( aq, ap );
  ret = vsnprintf ( dest -> chars, dest -> alloc, format, aq );
  va_end ( aq );

  if ( ret < 0 )
    {
      // glibc versions through libc 2.0.6 truncates and returns negative
      ret = dest -> size = dest -> alloc - 1;
      *(dest -> chars + ret ) = 0; // null terminate
      fprintf ( stderr,
		"printfToString: Warning formatted output truncated\n" );
    }
  else
    {
      // set the size
      int size = dest -> size =
	( ret < dest -> alloc ) ? ret : -1 + dest -> alloc;

      if ( ret > size )
	{
	  // resize and retry ( >= glibc 2.1 )
	  size_t alloc = 1 + ret;
	  StringSetCapacity ( dest, alloc );

	  ret = vsnprintf ( dest -> chars, dest -> alloc, format, ap );
	  dest -> size = ret;
	}
    }
  return ret;
}

int printfToString ( String* dest, char const* format, ... )
{
  va_list ap;
  int ret;

  va_start ( ap, format );
  ret = va_printfToString ( dest, format, ap );
  va_end ( ap );

  return ret;
}
