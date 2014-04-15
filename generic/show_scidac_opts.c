/****************** show_generic_opts.c ********************************/
/* MIMD Version 7 */

/* List options selected in the compilation */

extern int this_node;
#ifdef HAVE_QMP
#include <qmp.h>
#endif
#ifdef HAVE_QIO
#include <qio.h>
#endif
#ifdef HAVE_QLA
#include <qla.h>
#endif
#ifdef HAVE_QDP
#include <qdp.h>
#endif
#ifdef HAVE_QOP
#include <qop.h>
#endif

void show_scidac_opts(void){
#ifdef HAVE_QMP
  if(this_node==0)printf("QMP version %s\n",QMP_version_str());
#endif
#ifdef HAVE_QIO
  if(this_node==0)printf("QIO version %s\n",QIO_VERSION);
#endif
#ifdef HAVE_QLA
  if(this_node==0)printf("QLA version %s\n",QLA_version_str());
#endif
#ifdef HAVE_QDP
  if(this_node==0)printf("QDP version %s\n",QDP_version_str());
#endif
#ifdef HAVE_QOP
  if(this_node==0)printf("QOP version %s\n","Unknown");
#endif
}
