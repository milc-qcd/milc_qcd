#ifdef GB_BARYON /* To make sure other programs can compile properly */
#define _POSIX_C_SOURCE 200809
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "generic_includes.h"

// always map mmap cache to RAM
#undef GB_BARYON_MMAP

// Keep track of how big the mmap is
static long long tot_mmap_size; //accumulated total mmap cache size

int alloc_mmap_cache(mmap_cache* obj, unsigned nbuf, size_t buf_size, const char* cdir)
{
  obj->nbuffers = nbuf;
  obj->buffer_size = buf_size;
  obj->mem_page_size = sysconf(_SC_PAGESIZE); // system memory page size
  if(obj->mem_page_size == -1)
    {
      fprintf(stderr,"alloc_mmap_cache: invalid sysconf(PAGESIZE)\n");
      return -1;
    }
  unsigned nblocks = 1 + buf_size / obj->mem_page_size; // number of pages rounded up
  size_t block_size = nblocks * obj->mem_page_size;  // length of each buffer in bytes
  obj->cache_size = nbuf * block_size; // total size of backing cache
  char* fbase = "/back-store.XXXXXX";
  int slen = strlen(cdir) + strlen(fbase) + 1; // don't forget null terminator

  #ifdef GB_BARYON_MMAP
  obj->back_file = calloc(slen,sizeof(char));
  strcpy(obj->back_file,cdir); strcat(obj->back_file,fbase);
  mode_t oumask = umask(S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH); // set umask
  obj->fback = mkstemp(obj->back_file);  // O_RDWR|O_CREAT|O_TRUNC);
  umask(oumask); // reset mask
  if(obj->fback < 0)
    {
      fprintf(stderr,"alloc_mmap_cache: open failed, backing file %s\n",obj->back_file);
      return -1;
    }

  if(lseek(obj->fback,obj->cache_size-1,SEEK_SET) == -1) // seek to end of sparse file
    {
      fprintf(stderr,"alloc_mmap_cache: cannot seek backing file %s to offset %ld\n",obj->back_file,obj->cache_size-1);
      return -1;
    }
  {
  long word = 0;
  if (write(obj->fback, &word, sizeof(long)) != sizeof(long)) // write a null word
      {
        fprintf(stderr,"alloc_mmap_cache: cannot write backing file\n");
        return -1;
      }
  }
  fsync(obj->fback);

  // mmap file
  obj->base = mmap(static_cast(void*,NULL), obj->cache_size,
                   PROT_READ|PROT_WRITE,
                   MAP_SHARED,
                   obj->fback,static_cast(off_t,0));

  if(obj->base == MAP_FAILED)
    {
      fprintf(stderr,"alloc_mmap_cache: cannot mmap backing file\n");
      return -1;
    }
  if( unlink(obj->back_file) ){
    fprintf(stderr,"alloc_mmap_cache: cannot unlink %s\n",obj->back_file);
   }
  //pointers to buffers

  obj->buffer = calloc(obj->nbuffers,sizeof(void*));
  {
    unsigned j;
    for(j=0; j<obj->nbuffers; ++j)
      {
        obj->buffer[j] = static_cast(char*,obj->base) + j*block_size;
        int stat;
        if((stat = posix_madvise(obj->buffer[j],buf_size,POSIX_MADV_SEQUENTIAL)))
          {
            fprintf(stderr,"alloc_mmap_cache: advice MADV_SEQUENTIAL failed (%d)\n",stat);
          }

      }
  }

  #else
  node0_printf("No mmap created! Everything will be mapped to RAM!\n");
  obj->buffer = calloc(obj->nbuffers,sizeof(void*));
  {
    unsigned j;
    for(j=0; j<obj->nbuffers; ++j)
      {
        obj->buffer[j] = calloc(block_size,1);
  }
  }
  #endif
  node0_printf("(alloc_mmap_cache) previous mmap size on node0: %lld bytes\n", tot_mmap_size);
  tot_mmap_size += obj->cache_size;
  node0_printf("(alloc_mmap_cache) current mmap size on node0: %lld bytes\n", tot_mmap_size);
  return 0;
}


int free_mmap_cache(mmap_cache* obj)
{

#ifdef GB_BARYON_MMAP
  free(obj->buffer);
  obj->buffer = NULL;
  munmap(obj->base, obj->cache_size);
  obj->base = NULL;
  close(obj->fback);
  unlink(obj->back_file);
  free(obj->back_file);
#else
  for (int j = 0; j < obj->nbuffers; j ++)
    free(obj->buffer[j]);
  free(obj->buffer);
  obj->buffer = NULL;
#endif
  node0_printf("(free_mmap_cache) previous mmap size on node0: %lld bytes\n", tot_mmap_size);
  tot_mmap_size -= obj->cache_size;
  node0_printf("(free_mmap_cache) current mmap size on node0: %lld bytes\n", tot_mmap_size);
  return 0;
}

// NOT implemented! Do not use!
int msync_mmap_buffer(mmap_cache* obj,int nbuf)
{
exit(-1);
return -1;
}
#undef _POSIX_C_SOURCE
#endif /* GB baryon */
