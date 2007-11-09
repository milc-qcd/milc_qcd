#ifndef CHECK_MALLOC_H
#define CHECK_MALLOC_H

#ifdef CHECK_MALLOC

static void * _malloc_ptr;

#define malloc(_size) \
  (( (_malloc_ptr = malloc(_size)), \
   (this_node == 0 ? \
   printf("%x = malloc(%d) %s:%d\n",_malloc_ptr,_size,__func__,__LINE__) \
   && fflush(stdout) : 0 )), _malloc_ptr)

#define realloc(_ptr,_size) \
  (( (_malloc_ptr = realloc(_ptr,_size)), \
   (this_node == 0 ? \
   printf("%x = realloc(%x,%d) %s:%d\n",_malloc_ptr,_ptr,_size,__func__,__LINE__) \
   && fflush(stdout) : 0 )), _malloc_ptr)

#define calloc(_nelem,_elsize) \
  (( (_malloc_ptr = calloc(_nelem,_elsize)), \
   (this_node == 0 ? \
   printf("%x = calloc(%d,%d) %s:%d\n",_malloc_ptr,_nelem,_elsize,__func__,__LINE__) \
   && fflush(stdout) : 0 )), _malloc_ptr)

#define free(_ptr) {\
    node0_printf("free(%x) %s:%d\n",_ptr,__func__,__LINE__); \
    fflush(stdout); \
  free(_ptr); }

#ifdef QCDOC

#define qcdoc_alloc(_size) \
  (( (_malloc_ptr = qcdoc_alloc(_size)), \
   (this_node == 0 ? \
   printf("%x = malloc(%d) %s:%d\n",_malloc_ptr,_size,__func__,__LINE__) \
   && fflush(stdout) : 0 )), _malloc_ptr)

#define qfree(_ptr) {\
    node0_printf("free(%x) %s:%d\n",_ptr,__func__,__LINE__); \
    fflush(stdout); \
  qfree(_ptr); }

#endif

#endif

#endif
