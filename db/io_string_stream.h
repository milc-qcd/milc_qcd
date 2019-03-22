#ifndef __IO_STRING_STREAM_H__
#define __IO_STRING_STREAM_H__

#include <stdlib.h>
#include <stdarg.h>
#include "milc_datatypes.h"

#ifndef static_cast
#define static_cast(t,v) ((t)(v))
#endif

typedef struct io_string_stream_t {
  char *base;     // application owns this pointer
  size_t length;
} io_string_stream;

char*  io_string_stream_alloc(io_string_stream *strm, size_t buf_size);

char*  io_string_stream_realloc(io_string_stream *strm, size_t additional);

int io_string_stream_free(io_string_stream *strm);

int io_printf_to_string_stream(io_string_stream *strm, int max_len, const char* format, ...);

#define io_JSON_as_text(strm,max_len,format,...) io_printf_to_string_stream(strm,max_len,format,__VA_ARGS__)

int io_JSON_key(io_string_stream *strm, const char* key);

int io_JSON_quoted(io_string_stream *strm, const char* value);

int io_JSON_quote(io_string_stream *strm);

int io_JSON_sep(io_string_stream *strm);

int io_JSON_begin_set(io_string_stream *strm);

int io_JSON_end_set(io_string_stream *strm);

int io_JSON_begin_array(io_string_stream *strm);

int io_JSON_end_array(io_string_stream *strm);

int io_JSON_complex_array(io_string_stream *json_re, io_string_stream *json_im, complex vec[], int length);

#endif
