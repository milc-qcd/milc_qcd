#include "io_string_stream.h"
#include <stdio.h>
#include <string.h>

int main()
{
  io_string_stream meta;
  int sz = 25;
  io_string_stream_alloc(&meta, sz);

  io_JSON_begin_set(&meta);
  io_JSON_key(&meta,"greet");
  io_JSON_quoted(&meta,"Hello World!");
  io_JSON_sep(&meta);
  io_JSON_key(&meta,"quark_mass");
  io_JSON_as_text(&meta,10,"%f",0.0031);
  io_JSON_sep(&meta);
  io_JSON_key(&meta,"antiquark_mass");
  io_JSON_as_text(&meta,10,"%f",0.631);
  io_JSON_end_set(&meta);
 
  printf("%s\n", meta.base);

  io_string_stream json_re; json_re.base=NULL;json_re.length=0;
  io_string_stream json_im; json_im.base=NULL;json_im.length=0;

  srand48(4141736);
  const int vlen = 8;
  complex vec[vlen];
  for(int j=0; j<vlen; ++j) { vec[j].real = (1.0 - 2.0*drand48()); vec[j].imag = (1.0-2.0*drand48()); }

  io_JSON_complex_array(&json_re, &json_im, vec, vlen);

  printf("re %d: %s\n",strlen(json_re.base),json_re.base);
  printf("im %d: %s\n",strlen(json_im.base),json_im.base);

  return(0);

}
