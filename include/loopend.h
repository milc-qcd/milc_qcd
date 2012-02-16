#ifndef _LOOPEND_H
#define _LOOPEND_H
/* Redefinition of macros in macros.h for slightly better looping */

/* macros.h must come first!! */
#include "../include/macros.h"

#define LOOPEND   /* We take this as the rule for now  -C.D. */

#ifdef LOOPEND

#undef FORSOMEPARITY
#define FORSOMEPARITY(i,s,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ), s= &(lattice[i]); \
i<loopend; i++,s++)

#undef FORSOMEFIELDPARITY
#define FORSOMEFIELDPARITY(i,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ); \
i<loopend; i++)

#define END_LOOP }

#else

#define END_LOOP 	/* define it to be nothing */

#endif

/* NOTE FORSOMEPARITYDOMAIN is still defined in terms of FORSOMEPARITY
   so this change gets applied to it as well */

#endif /* _LOOPEND_H */
