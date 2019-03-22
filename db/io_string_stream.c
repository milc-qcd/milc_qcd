#include "io_string_stream.h"
#include <malloc.h>
#include <string.h>
#include <stdio.h>

char* io_string_stream_alloc(io_string_stream *strm, size_t buf_size)
{
  strm->length = buf_size + 1; // plus null termination
  strm->base = static_cast(char*,calloc(strm->length,sizeof(char)));
  strm->base[0] = '\0';
  return strm->base;
}

char* io_string_stream_realloc(io_string_stream *strm, size_t additional)
{
  size_t new_length = strm->length + additional + 1;
  strm->length = new_length;
  strm->base = static_cast(char*,realloc(strm->base, new_length * sizeof(char)));
  return strm->base;
}

int io_string_stream_free(io_string_stream *strm)
{
  if(strm->base)
    free(strm->base);
  strm->base = NULL;
  strm->length = 0;
  return(0);
}

int io_printf_to_string_stream(io_string_stream *strm, int max_len, const char* format, ...)
{
  int n;
  va_list ap;
  size_t eos = strlen(strm->base);
  size_t avail = strm->length - eos - 1;
  if( avail < max_len )
    {
      io_string_stream_realloc(strm, 2*max_len);
    }

  va_start(ap, format);
  n = vsnprintf(strm->base+eos, strm->length, format, ap);
  va_end(ap);

  return(n);
}

int io_JSON_key(io_string_stream *strm, const char* key)
{
  int n = io_printf_to_string_stream(strm, strlen(key)+3, "\"%s\":",key);
  return(n);
}

int io_JSON_sep(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, ",");
  return(n);
}

int io_JSON_quoted(io_string_stream *strm, const char* value)
{
  int n = io_printf_to_string_stream(strm, strlen(value)+2, "\"%s\"", value);
  return(n);
}

int io_JSON_quote(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, "\"");
  return(n);
}

int io_JSON_begin_set(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, "{");
  return(n);
}

int io_JSON_end_set(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, "}");
  return(n);
}

int io_JSON_begin_array(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, "[");
  return(n);
}

int io_JSON_end_array(io_string_stream *strm)
{
  int n = io_printf_to_string_stream(strm, 1, "]");
  return(n);
}

int io_JSON_complex_array(io_string_stream *json_re, io_string_stream *json_im, complex vec[], int length)
{
  const int precision = 6; // six digits real precision (format %e default)
  // each element of an array of reals: +f.pppppe+nn
  const int len_num = precision + 7;
  const int len_jvec = length * (len_num  + 1) + 2;
  if(json_re->length < len_jvec+1)
    io_string_stream_realloc(json_re,len_jvec-json_re->length);
  if(json_im->length < len_jvec+1)
    io_string_stream_realloc(json_im,len_jvec-json_im->length);
  size_t lr0 = strlen(json_re->base);
  size_t li0 = strlen(json_im->base);
  const int flen = 7;
  char fmt[flen+1];
  snprintf(fmt,flen,"%%.%de",precision);

  io_JSON_begin_array(json_re);
  io_JSON_begin_array(json_im);
  for(int j=0; j<length; ++j)
    {
      double re = vec[j].real;
      double im = vec[j].imag;
      io_JSON_as_text(json_re, len_num, fmt, re);
      io_JSON_as_text(json_im, len_num, fmt, im);
      if( j != length-1 )
	{
	  io_JSON_sep(json_re);
	  io_JSON_sep(json_im);
	}
    }
  io_JSON_end_array(json_re);
  io_JSON_end_array(json_im);
  int lr = static_cast(int,strlen(json_re->base)-lr0);
  int li = static_cast(int,strlen(json_im->base)-li0);
  if(lr>=li) { return(lr); } else { return(li); }
}

  
