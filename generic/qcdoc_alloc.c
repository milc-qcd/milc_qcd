/* This is a first stab at a wrapper function for the qalloc statement on 
   the QCDOC. Some include file needs to have a 
    #define malloc qcdoc_alloc
    and 
    #define free qfree

  this wrapper function should provide graceful behavior and a warning
  for when EDRAM is exceeded.  3/11/04 EBG */
#include <qmp.h>
#include <qalloc.h>
#include "generic_includes.h"

void *qcdoc_alloc(size_t nbytes){

  void* ptr;
  ptr = qalloc(QFAST|QCOMMS,nbytes);

  if(ptr==NULL){
    node0_printf("EDRAM exceeded. now using DDR!!!!!!!!!!!!!\n");
    fflush(stdout);
    ptr = qalloc( QCOMMS,nbytes);
  }
  return ptr;
}
