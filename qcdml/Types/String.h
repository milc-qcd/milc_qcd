#ifndef __String_h__
#define __String_h__

#include <stdio.h>
#include <stddef.h>

typedef struct String_tag {
  char*  chars;
  size_t size;
  size_t alloc;
} String;

typedef struct {
  char* ptr;
  String* str;
} String_pos;

String* new_String ( size_t size, char const* s );

void delete_String ( String* self );

String* StringClone ( String const* self );

size_t StringSetCapacity ( String* self, size_t newSize );

void StringCopy ( String* self, String const* src );

void StringCopyChars ( String* self, char const* src );

void StringConcat ( String* self, String const* src );

void StringConcatChars ( String* self, char const* src );

int printfToString ( String* dest, char const* format, ... );

void dumpString ( FILE* out, char const* name, String const* self );

// TODO: implement functions below ==========================

// Return an iterator to the first character of string
String_pos StringFirst ( String const* self );

// Return an iterator to the last character of string
String_pos StringLast ( String const* self );

// The  strchr()  function  returns  an iterator  to the first
// occurrence of the character chr
String_pos StringIndex ( String const* self, int chr );

// The rindex function  returns  an interator to  the  last
// occurrence of the character c
String_pos StringRindex ( String const* self, int chr );

void insertString ( String_pos self, String const* src );

void replaceString ( String_pos from, String_pos to, String const* src );

#endif
