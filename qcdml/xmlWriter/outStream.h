#ifndef __outStream_h__
#define __outStream_h__
#include <stdio.h>
#include <String.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct outStream_tag outStream;

struct outStream_tag
{
  // public:

  // printf to stream
  void (*print) ( outStream* self, const char *format, ... );

  // private:

  FILE*   _fp;
  String* _st;
};

outStream* new_outStreamFile ( FILE* fp );

outStream* new_outStreamString ( String* str );

void delete_outStream ( outStream* obj );

#ifdef __cplusplus
}
#endif
#endif

