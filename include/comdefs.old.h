#ifndef _COMDEFS_H
#define _COMDEFS_H
/************************ comdefs.h *************************************
*									*
*  Macros and declarations for communications routines com_*.c          *
*  MIMD version 6 							*
*									*
*/

/* Version that allows 32 sublattices for extended actions */

/* 
   
   Communications routines will assume that the lattice is stored as an
   array of structures of type "site".

   Modifications:
   1/27/00 Included option of 32 sublattices (for improved actions) -- UMH
   2/26/98 Changed int32type to signed
   8/5/97 ANSI prototyping for com_*.c routines C.D.
   8/10/96  Added defs for broadcast_bytes and wrappers for parallel file IO C.D.
*/

#include "../include/int32type.h"
#include "../include/dirs.h"
#include "../include/su3.h"
#include "../include/complex.h"
#include "../include/macros.h"
#include <lattice.h>     /* Needed ONLY for "site" in this header */

#ifdef MPI
/* For MPI, need to include mpi.h to define MPI_Request */
#include <mpi.h>
/* For MPI on paragon, mynode() and numnodes() can't replace nx routines */
#define mynode MILC_mynode
#define numnodes MILC_numnodes
#define dclock MILC_dclock
#endif

/* arguments to the make_gather() routine */
#define FORWARDS 1
#define BACKWARDS (-1)	/* BACKWARDS = -FORWARDS */
#define OWN_INVERSE 0
#define WANT_INVERSE 1
#define NO_INVERSE 2
#define ALLOW_EVEN_ODD 0
#define NO_EVEN_ODD 1
#define SAME_PARITY 0
#define SWITCH_PARITY 1
#define SCRAMBLE_PARITY 2

/* Structure to keep track of outstanding sends and receives */
typedef struct {
	/* node sending or receiving message */
    int msg_node;
	/* size of message in bytes */
    int msg_size;
	/* address of buffer malloc'd for message */
    char *msg_buf;
	/* message id returned by system call */
#ifdef MPI
    MPI_Request msg_id;
#else
    int msg_id;
#endif
} msg_tag;

/**********************************************************************/
/* Declarations for all routines called in the com_*.c files */

void start_handlers();
void initialize_machine(int argc, char **argv);
void make_nn_gathers();
void sort_eight_neighborlists(int index);
void sort_eight_special(void **pt);

void neighbor_coords_special(
 int x,int y,int z,int t, /* coordinates of site */
 int *dirpt,              /* direction (eg XUP) */
 int fb,                  /* "forwards/backwards"  */
 int *x2p,int *y2p,int *z2p,int *t2p);
                          /* pointers to coordinates of neighbor */
int make_gather(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
         		/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int inverse,		/* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int parity_conserve);	/* {SAME,SWITCH,SCRAMBLE}_PARITY */

void neighbor_coords(
 int x, int y, int z, int t,  /* coordinates of site */
 int dir,	              /* direction (eg XUP) */
 int *x2p, int *y2p, int *z2p, int *t2p);
                             /* pointers to coordinates of neighbor */
msg_tag * start_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest);		/* one of the vectors of pointers */

void restart_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf);        /* previously returned by start_gather */

msg_tag * start_gather_from_temp(
/* arguments */
 void * field,		/* which field? pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest);		/* one of the vectors of pointers */

void restart_gather_from_temp(
/* arguments */
 void * field,		/* which field? pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf);        /* previously returned by start_gather */

void wait_gather(msg_tag *mbuf);
void cleanup_gather(msg_tag *mbuf);

msg_tag * start_general_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest);		/* one of the vectors of pointers */

msg_tag * start_general_gather_from_temp(
/* arguments */
 void * field,	        /* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest);		/* one of the vectors of pointers */


void wait_general_gather(msg_tag *mbuf);
void cleanup_general_gather(msg_tag *mbuf);

char * field_pointer_at_coordinates(
/* arguments */
 int field,	/* offset of one of the fields in lattice[] */
 int size,	/* size of the field in bytes */
 int x,int y,int z,int t);	/* coordinates of point to get field from */

char * field_pointer_at_direction(
/* arguments */
 field_offset field,	/* offset of one of the fields in lattice[] */
 int size,	/* size of the field in bytes */
 site *s,	/* pointer to a site on this node */
 int direction);	/* direction of site's neighbor to get data from.
		   one of XUP,XDOWN,YUP... */

void cleanup_field_pointer(char * buf);
void send_field(char *buf, int size, int tonode);
void get_field(char *buf, int size, int fromnode);
char * machine_type();
int mynode();
int numnodes();
void g_sync();
void g_floatsum( Real *fpt);
void g_vecfloatsum( Real *fpt, int nReals);
void g_doublesum( double *dpt);
void g_vecdoublesum( double *dpt, int ndoubles);
void g_complexsum( complex *cpt);
void g_dcomplexsum( double_complex *cpt);
void g_veccomplexsum( complex *cpt, int ncomplex);
void g_vecdcomplexsum( double_complex *cpt, int ncomplex);
void g_wvectosumReal( wilson_vector *wvpt);
void g_xor32( u_int32type *pt );
void g_floatmax( Real *fpt);
void g_doublemax( double *dpt);
void broadcast_float(Real *fpt);
void broadcast_double(double *dpt);
void broadcast_complex(complex *cpt);
void broadcast_dcomplex(double_complex *cpt);
void broadcast_bytes(char *buf,int size);
void send_integer(int tonode, int *address);
void receive_integer(int fromnode, int *address);

/* On the Paragon dclock is a library routine with the
   same functionality as ours */
/* Either way, it needs to be declared double */
double dclock();
void time_stamp(char *msg);

void terminate(int status);

void normal_exit(int status);

#endif	/* _COMDEFS_H */
