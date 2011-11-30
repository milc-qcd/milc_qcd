#ifndef _MACROS_H
#define _MACROS_H

#include <defines.h>
#include "../include/dirs.h"

/* Macros common to all applications */

/* ---------------------------------------------------------- */
/* Constants */

#define PI 3.14159265358979323846

/* ---------------------------------------------------------- */
/* Conventions for defining checkerboard parity of inversions */

#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03

/* ---------------------------------------------------------- */
/* Storage constants */

#define MAXFILENAME 256   /* ASCII string length for all file names */

/* ---------------------------------------------------------- */
/* Names of gauge fixing options */

#define NO_GAUGE_FIX 30
#define COULOMB_GAUGE_FIX 31
#define LANDAU_GAUGE_FIX 32
#define AXIAL_GAUGE_FIX 33

/* ---------------------------------------------------------- */
/* "field offset" and "field pointer" */

/* used when fields are arguments to subroutines */
/* Usage:  fo = F_OFFSET( field ), where "field" is the name of a field
  in lattice.
     address = F_PT( &site , fo ), where &site is the address of the
  site and fo is a field_offset.  Usually, the result will have to be
  cast to a pointer to the appropriate type. (It is naturally a char *).
*/
typedef int field_offset;
#define F_OFFSET(a) \
  ((field_offset)(((char *)&(lattice[0]. a ))-((char *)&(lattice[0])) ))
#define F_PT( site , fo )  ((char *)( site ) + (fo)) 

/* ---------------------------------------------------------- */
/* Macros for looping over directions */

#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)

#define FORALLUPDIRBUT(direction,dir) \
   for(dir=XUP; dir<= TUP; dir++)if(dir != direction)

#define OPP_PAR(parity) (0x03 ^ parity)	/* Switches EVEN and ODD. Nulls EVENANDOdd*/

/* ---------------------------------------------------------- */
/* printf on node zero only */
#define node0_printf if(this_node==0)printf

/* ---------------------------------------------------------- */
#define ON 1
#define OFF 0

/* ---------------------------------------------------------- */
/* Macros for looping over sites on a node */

/**********************************************************************/
/* WARNING: FORSOMEPARITY and FORSOMEPARITYDOMAIN is redefined in some
   routines if LOOPEND is specified.  See loopend.h */
/**********************************************************************/

#ifndef N_SUBL32   
/*--------------*/

/* Standard red-black checkerboard */

/* macros to loop over sites of a given parity.
   Usage:  
	int i;
	site *s;
	FOREVENSITES(i,s){
	    commands, where s is a pointer to the current site and i is
	    the index of the site on the node
	}
*/

#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<even_sites_on_node;i++,s++)
#define FOREVENFIELDSITES(i) \
    for(i=0;i<even_sites_on_node;i++)
#define FORODDSITES(i,s) \
    for(i=even_sites_on_node,s= &(lattice[i]);i<sites_on_node;i++,s++)
#define FORODDFIELDSITES(i) \
    for(i=even_sites_on_node;i<sites_on_node;i++)
#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]); \
    i< ( (choice)==EVEN ? even_sites_on_node : sites_on_node); \
    i++,s++)
#define FORSOMEFIELDPARITY(i,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 );	\
    i< ( (choice)==EVEN ? even_sites_on_node : sites_on_node); \
    i++)
#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)
#define FORALLFIELDSITES(i) \
    for(i=0;i<sites_on_node;i++)

/*--------------*/
#else	/* N_SUBL32 */
/*--------------*/

/*  32 sublattice checkerboard */

/* macros to loop over sites in a given sublattice.
   Usage:  
	int i, subl;
	site *s;
	FORSOMESUBLATTICE(i,s,subl){
	    commands, where s is a pointer to the current site and i is
	    the index of the site on the node
	}
*/
#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)
#define FORALLFIELDSITES(i) \
    for(i=0;i<sites_on_node;i++)
#define FORSOMESUBLATTICE(i,s,subl) \
    for( i=(subl*subl_sites_on_node), s= &(lattice[i]), \
	 last_in_loop=((subl+1)*subl_sites_on_node); \
	 i< last_in_loop; i++,s++)

/*--------------*/
#endif	/* N_SUBL32 */


#ifdef SCHROED_FUN
#define FOREVENSITESDOMAIN(i,s) \
    FOREVENSITES(i,s) if(s->t > 0)
#define FORODDSITESDOMAIN(i,s) \
    FORODDSITES(i,s) if(s->t > 0)
#define FORALLSITESDOMAIN(i,s) \
    FORALLSITES(i,s) if(s->t > 0)
#define FORSOMEPARITYDOMAIN(i,s,parity) \
    FORSOMEPARITY(i,s,parity) if(s->t > 0)
#else
#define FOREVENSITESDOMAIN FOREVENSITES
#define FORODDSITESDOMAIN FORODDSITES
#define FORALLSITESDOMAIN FORALLSITES
#define FORSOMEPARITYDOMAIN FORSOMEPARITY
#endif

#endif  /* _MACROS_H */
