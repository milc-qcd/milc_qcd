/******************  com_mpi.c *****************************************/
/* Communications routines for the SU3 program 
   MIMD version 6.
   This file is communications-scheme dependent.
   MPI version - allegedly machine independent
*/
#define NOWHERE -1	/* Not an index in array of fields */

/* Modifications   

   4/20/02  added start_general_gather_from_temp C.D.
   1/30/00  combined with Schroedinger functional and
            32 sublattice versions - UMH
   11/27/98 Corrected g_wvectosumReal and made independent of su3.h. C.D.
    8/05/97 ANSI prototyping for all routines C.D.
    9/2/97  Revised to allow gathers from temporary fields.  neighbor[]
	is now list of indices, add start/restart_gather_from_temp D.T.
   10/05/96 Moved parallel I/O wrappers to io_ansi.c C.D.
    9/23/96 Explicit void types for modules with empty returns C.D.
    9/20/96 Added restore/save_checkpoint C.D.  
    9/20/96 Improved sort_site_list C.D.  
    9/20/96 Added broadcast_bytes and wrappers for system-dependent
            parallel file system calls C.D.   
*/

/*
   initialize_machine() does any machine dependent setup at the
   very beginning.

   Some routines implement basic functions in a machine independent
   way.
   mynode() returns node number of this node.
   numnodes() returns number of nodes
   g_sync() provides a synchronization point for all nodes.
   g_floatsum() sums a Realing point number over all nodes.
   g_doublesum() sums a double over all nodes.
   g_vecdoublesum() sums a vector of doubles over all nodes.
   g_complexsum() sums a single precision complex number over all nodes.
   g_dcomplexsum() sums a double precision complex number over all nodes.
   g_veccomplexsum() sums a vector of single precision complex numbers over all
	nodes.
   g_wvectosumReal() sums a single precision wilson vector over all nodes.
   g_xor32() finds global exclusive or of 32-bit word
   g_floatmax() finds maximum of a Realing point number over all nodes.
   g_doublemax() finds maximum of a double over all nodes.
   broadcast_float()  broadcasts a single precision number from
	node 0 to all nodes.
   broadcast_double()  broadcasts a double precision number
   broadcast_complex()  broadcasts a single precision complex number
   broadcast_dcomplex()  broadcasts a double precision complex number
   broadcast_bytes()  broadcasts a number of bytes
   send_integer() sends an integer to one other node
   receive_integer() receives an integer
   dclock() returns a double precision time, with arbitrary zero
   terminate() kills the job on all processors

   make_nn_gathers()  makes all necessary lists for communications with
   nodes containing neighbor sites.

   field_pointer_at_coordinates() returns the address of a field in
   the lattice given its coordinates.
   field_pointer_at_direction() returns the address of a field in the
   lattice at a direction from a given site.
   cleanup_field_pointer() frees the buffers that field_pointer...
   allocated.

   start_gather() starts asynchronous sends and receives required
   to gather neighbors.
   start_gather_from_temp() starts asynchronous sends and receives required
   to gather neighbors from temporary array of fields.
   wait_gather()  waits for receives to finish, insuring that the
   data has actually arrived.
   cleanup_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   restart_gather() repeats the internode communications of a previous
   gather.
   restart_gather_from_temp() repeats the internode communications of a
   previous gather of temporary field.

   start_general_gather() starts asynchronous sends and receives required
   to gather fields at arbitrary displacement.
   start_general_gather_from_temp() starts asynchronous sends and receives required
   to gather neighbors from a temporary array of fields.
   wait_general_gather()  waits for receives to finish, insuring that the
   data has actually arrived, and sets pointers to received data.
   cleanup_general_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.

   send_field() sends a field to one other node.
   get_links() receives a field from some other node.
*/

/*  Note mpi.h is included in comdefs.h */

#include <assert.h>
#include "generic_includes.h"
#include "../include/config.h"

/* message types used here */
#define FIELD_TYPE 12	/* type of field sent from one node to another */
#define SEND_INTEGER_TYPE 17	/* send an integer to one other node */

#define GENERAL_GATHER_BASE 102	/* types from this to this+number_of_nodes
				are used by the general_gather routines */
#define GATHER_BASE 614	/* types greater than or equal to this are used
			   by the gather routines */

#define MAX_GATHERS 32	/* Maximum number of gather tables */

#include <sys/types.h>

/**********************************************************************/
/* "comlink" is the basic structure used in gathering neighboring sites.
   Each node will maintain one such structure for each direction for each
   (other) node that contains sites that are neighbors of the sites on
   this node.  For example, if the XUP neighbors of sites on this node
   are found on two other nodes, then this node will maintain a linked
   list of two comlink structures for gathering from the XUP direction.
*/
#ifdef N_SUBL32
struct comlink {
	/* pointer to next in list, NULL if this is last */
    struct comlink *nextcomlink;
	/* number of the node to which we connect */
    int othernode;
	/* number of sites in a sublattice on this node that have
	   neighbors on other node connected by this "comlink"
	   as well as the total number. */
    int n_subl_connected[N_SUBL32], n_connected;
	/* Address of list of indices of sites in one sublattice
	   (on this node) whose neighbors are found through this comlink.
	   The lists follow in order of the sublattice number,
	   so to get all sites you start at sitelist[0] and take
	   n_connected addresses. */
	/* When the sublattices are treated separately, the sublattice
	   index always corresponds to the the sublattice of the site
	   on this node. Different comlink structures may point to the
	   same list. For example, the receive list for one gather may
	   be a send list for the opposite gather. */
    int send_subl[N_SUBL32];
	/* The sublattice on the other node for ALLOW_EVEN_ODD.
	   Is not used otherwise. */
    int *sitelist[N_SUBL32];
};
#else
struct comlink {
	/* pointer to next in list, NULL if this is last */
    struct comlink *nextcomlink;
	/* number of the node to which we connect */
    int othernode;
	/* number of even sites on this node that have neighbors on
	   other node connected by this "comlink", and same for odd
	   sites on this node. */
    int n_even_connected, n_odd_connected;
	/* Address of list of indices of even sites (on this node)
	   whose neighbors are found through this comlink, same for odd.
	   The odd list follows the even list, so to get all sites you
	   start at esitelist and take n_even_connected+n_odd_connected
	   addresses. */
	/* When the comlink is describing sites to be sent, the "odd"
	   list lists sites whose neighbors are even.  This convention
	   is natural for the nearest neighbor gathers.  For gathers
	   which don't allow even and odd site gathers, the even list
	   is used for list of sites to be received and the odd
	   list for sites to be sent.  Different comlink structures
	   may point to the same list.  For example, the receive list
	   for one gather may be a send list for the opposite gather. */
    int *esitelist, *ositelist;
};
#endif
typedef struct comlink comlink;

/**********************************************************************/
/* Structure for requesting a field from another node */
typedef struct {
    int field;	/* offset of field in site */
    int size;	/* size of field */
    int index;	/* index of field on other node */
} msg_request;

/**********************************************************************/

/* Global variables for the communications stuff */
    /* message types for gather encode a sequence number for the gather
	so that if several gathers are going at once, you can read
	the message corresponding to the right one. */
    /* for computing message type in gather */
    int type_offset;	/* label gathers by round-robin */
    int gmt_offshift;	/* number of bits to shift field offset in
	general gather */
    /* macro to compute the message type */
#define GATHER_MSG_TYPE(type_offset) (GATHER_BASE + type_offset)
#define TYPE_OFFSET_MAX 1024
#define GENERAL_GATHER_MSG_TYPE (GENERAL_GATHER_BASE)

    /* Maintain lists of headers to lists of comlinks */
    /* One list for sites to be received, one for sites to be sent.
		These will often point to the same things */
	    comlink ** neighborlist, **neighborlist_send;
    /* addresses of neighboring sites, NULL if off-node */
    /* neighbor[X][i] is the index of the neighbor of site lattice[i] in
       gather number X */
    int ** neighbor;
    /* Number of gathers (mappings) that have been set up */
    int n_gathers;

/**********************************************************************/
/* Machine initialization */
void initialize_machine(int argc, char **argv){
  MPI_Init(&argc,&argv);
  type_offset = 0;
}


/**********************************************************************/
/* Set up "comlink" structures needed by gather routines.
   make_lattice() must be called first. */
void make_nn_gathers(){
int i;
void neighbor_coords_special(
 int x,int y,int z,int t, /* coordinates of site */
 int *dirpt,              /* direction (eg XUP) */
 int fb,                  /* "forwards/backwards"  */
 int *x2p,int *y2p,int *z2p,int *t2p);
                          /* pointers to coordinates of neighbor */
int max_type,tag_upper_bound,flag;

   /* Check for possible overflow of message type. */
   /* Note: with MPI-2 MPI_Comm_get_attr is preferred,
    but we keep MPI_Attr_get until MPI-2 is more widely 
    available */ 
   MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB,&tag_upper_bound,&flag);
   max_type = GATHER_BASE+TYPE_OFFSET_MAX;
    /* When this happens, fix the program */
    if( max_type >= tag_upper_bound ){
	if(this_node==0)printf(
            "make_nn_gathers: Possible overflow of gather message type\n %d > %d flag %d: Fix the program\n",max_type,tag_upper_bound,flag);
    }

    /* initialize neighborlist[] */
    neighborlist = (comlink **)malloc(NDIRS*sizeof(comlink *));
    neighborlist_send = (comlink **)malloc(NDIRS*sizeof(comlink *));
    /* Allocate space for lists of indexes of neighbor sites.
       (NOWHERE if neighbor not on this node) */
    neighbor = (int **)malloc(NDIRS*sizeof(int *));
    n_gathers=0;

    for(i=XUP;i<=TUP;i++)
	make_gather(neighbor_coords_special,&i,WANT_INVERSE,
	    ALLOW_EVEN_ODD,SWITCH_PARITY);

    /* Sort into the order we want for nearest neighbor gathers,
	so you can use XUP, XDOWN, etc. as argument in calling them. */
    /**    sort_eight_special((void **) neighbor );
    sort_eight_special((void **) neighborlist );
    sort_eight_special((void **) neighborlist_send ); **/
    sort_eight_neighborlists(XUP);
}

/**********************************************************************/
/* sort a list of eight pointers into the order we want for the
  nearest neighbor gathers:  XUP,YUP,ZUP,TUP,TDOWN,ZDOWN,YDOWN,XDOWN,
  starting from the index for the first pointer */

void sort_eight_neighborlists(int index){
  sort_eight_special((void **)&(neighbor[index]) );
  sort_eight_special((void **)&(neighborlist[index]) );
  sort_eight_special((void **)&(neighborlist_send[index]) );
}

/**********************************************************************/
/* sort a list of eight pointers into the order we want for the
  nearest neighbor gathers:  XUP,YUP,ZUP,TUP,TDOWN,ZDOWN,YDOWN,XDOWN */
void sort_eight_special(void **pt){
void *tt[8];
register int i;
    for(i=0;i<8;i++)tt[i]=pt[i];
    for(i=XUP;i<=TUP;i++){pt[i]=tt[2*i]; pt[OPP_DIR(i)]=tt[2*i+1];}
}

/**********************************************************************/
/* utility function for finding coordinates of neighbor */
/* This version for use by make_gather for nearest neighbor gathers */
void neighbor_coords_special(
 int x,int y,int z,int t, /* coordinates of site */
 int *dirpt,              /* direction (eg XUP) */
 int fb,                  /* "forwards/backwards"  */
 int *x2p,int *y2p,int *z2p,int *t2p)
                          /* pointers to coordinates of neighbor */
{
int dir;
    dir = (fb==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
    *x2p = x; *y2p = y; *z2p = z; *t2p = t;
    switch(dir){
	case XUP: *x2p = (x+1)%nx; break;
	case XDOWN: *x2p = (x+nx-1)%nx; break;
	case YUP: *y2p = (y+1)%ny; break;
	case YDOWN: *y2p = (y+ny-1)%ny; break;
	case ZUP: *z2p = (z+1)%nz; break;
	case ZDOWN: *z2p = (z+nz-1)%nz; break;
	case TUP: *t2p = (t+1)%nt; break;
	case TDOWN: *t2p = (t+nt-1)%nt; break;
	default: printf("BOTCH: bad direction\n"); exit(1);
    }
}


/**********************************************************************/
#define RECEIVE 0
#define SEND 1
/* add another gather to the list of tables */
int make_gather(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
        		/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int inverse,		/* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int parity_conserve)	/* {SAME,SWITCH,SCRAMBLE}_PARITY */
{
#ifdef N_SUBL32
/* copy linked list of comlinks, switching two sublattices */
comlink * copy_list_switch( comlink *old_compt, int *othersubl);
#else
/* copy linked list of comlinks, switching even and odd */
comlink * copy_list_switch( comlink *old_compt ); 
#endif

comlink *  make_send_receive_list(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
        		/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int forw_back,		/* FORWARDS or BACKWARDS */
 int send_recv);	/* SEND or RECEIVE list */
                        /* lists of sites to be sent or received */
register int i,j,k;	/* scratch */
register site *s;	/* scratch */
int dir;		/* direction */
int x,y,z,t;		/* coordinates */

#ifdef N_SUBL32
int subl, nsubl, *othersubl;

    if( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SCRAMBLE_PARITY ){
	printf("Unsupported combination for make_gather:");
	printf(" ALLOW_EVEN_ODD and SCRAMBLE_PARITY\n");
	terminate(1);
    }
#endif

    /* we will have one or two more gathers */
    if( inverse==WANT_INVERSE ) n_gathers += 2;
    else			n_gathers += 1;
    if(n_gathers > MAX_GATHERS){
	if(this_node==0)printf("Too many gathers! change MAX_GATHERS\n");
	exit(1);
    }
    /* lengthen neighborlist[] */
    neighborlist = (comlink **)realloc(neighborlist,
	n_gathers*sizeof(comlink *));
    neighborlist_send = (comlink **)realloc(neighborlist_send,
	n_gathers*sizeof(comlink *));
    /* Allocate more space for lists of pointers to neighbor sites.
       (NULL if neighbor not on this node) */
    neighbor = (int **)realloc(neighbor, n_gathers*sizeof(int *));
    if( inverse==WANT_INVERSE) {
	neighborlist[n_gathers-2] = neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-2] = neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-2] = (int *)malloc(sites_on_node*sizeof(int) );
        if(neighbor[n_gathers-2]==NULL){
            printf("make_gather: NODE %d no room for neighbor vector\n",this_node);
            terminate(1);
        }
        neighbor[n_gathers-1] = (int *)malloc(sites_on_node*sizeof(int) );
        if(neighbor[n_gathers-1]==NULL){
            printf("make_gather: NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-2;	/* index of gather we are working on */
    }
    else {
	neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-1] = (int *)malloc(sites_on_node*sizeof(int) );
        if(neighbor[n_gathers-1]==NULL){
            printf("make_gather: NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-1;
    }

#ifndef N_SUBL32

    /* Check to see if mapping has advertised parity and inverse properties */
    /* Also check to see if it returns legal values for coordinates */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
	k=(x+y+z+t)%2;

	if( x<0 || y<0 || z<0  || t<0 || x>=nx || y>=ny || z>=nz || t>=nt){
	    printf("DUMMY! Your gather mapping does not stay in lattice\n");
	    printf("It mapped %d %d %d %d to %d %d %d %d\n",
		s->x,s->y,s->z,s->t,x,y,z,t);
	    terminate(1);
	}
	if( ( parity_conserve==SAME_PARITY &&
		((k==0&&s->parity==ODD)||(k==1&&s->parity==EVEN))
	    )
	  ||( parity_conserve==SWITCH_PARITY &&
		((k==0&&s->parity==EVEN)||(k==1&&s->parity==ODD))
#ifdef SCHROED_FUN
	    ) ) if( *args!=TUP && *args!=TDOWN && t!=0 && t!=(nt-1) ) {
#else
	    ) ){
#endif
	    printf("DUMMY! Your gather mapping does not obey claimed parity\n");
	    printf("It mapped %d %d %d %d to %d %d %d %d\n",
		s->x,s->y,s->z,s->t,x,y,z,t);
	    terminate(1);
	}
	if( inverse==OWN_INVERSE ){
	    int x2,y2,z2,t2;
            func( x, y, z, t, args, FORWARDS, &x2,&y2,&z2,&t2);
	    if( s->x!=x2 || s->y!=y2 || s->z!=z2 || s->t!=t2 ){
	        printf(
		    "DUMMY! Your gather mapping is not its own inverse\n");
	        printf("It's square mapped %d %d %d %d to %d %d %d %d\n",
		    s->x,s->y,s->z,s->t,x2,y2,z2,t2);
	        terminate(1);
	    }
	}
    }

#else	/* N_SUBL32 */

    othersubl = (int *)malloc(N_SUBL32*sizeof(int));
    if(othersubl == NULL){
	printf("NODE %d: no room for othersubl\n",this_node);
	terminate(1);
    }
    for(subl=0;subl<N_SUBL32;subl++) othersubl[subl] = NOWHERE;

    /* Check to see if mapping has advertised parity and inverse properties */
    /* Also check to see if it returns legal values for coordinates */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
	subl = s->parity;
	nsubl = (x%2) + 2*(y%2) + 4*(z%2) + 8*(t%2);
	nsubl += 16*((x/2+y/2+z/2+t/2)%2);

	if( want_even_odd==ALLOW_EVEN_ODD ){
	    if( othersubl[subl] == NOWHERE ){
		othersubl[subl] = nsubl;
	    }
	    else if( othersubl[subl] != nsubl ){
		printf("DUMMY! Your gather mixes up sublattices: %d vs %d\n",
		    othersubl[subl], nsubl);
		terminate(1);
	    }
	}

	if( x<0 || y<0 || z<0  || t<0 || x>=nx || y>=ny || z>=nz || t>=nt){
	    printf("DUMMY! Your gather mapping does not stay in lattice\n");
	    printf("It mapped %d %d %d %d to %d %d %d %d\n",
		s->x,s->y,s->z,s->t,x,y,z,t);
	    terminate(1);
	}
	if( parity_conserve==SAME_PARITY && subl!=nsubl ){
	    printf("DUMMY! Your gather mapping does not obey claimed parity");
	    printf(", namely SAME_PARITY\n");
	    printf("It mapped %d %d %d %d with %d to %d %d %d %d with %d\n",
		s->x,s->y,s->z,s->t,subl,x,y,z,t,nsubl);
	    terminate(1);
	}
	if( parity_conserve==SWITCH_PARITY && subl==nsubl ){
	    printf("DUMMY! Your gather mapping does not obey claimed parity");
	    printf(", namely SWITCH_PARITY\n");
	    printf("It mapped %d %d %d %d with %d to %d %d %d %d with %d\n",
		s->x,s->y,s->z,s->t,subl,x,y,z,t,nsubl);
	    terminate(1);
	}
	if( inverse==OWN_INVERSE ){
	    int x2,y2,z2,t2;
            func( x, y, z, t, args, FORWARDS, &x2,&y2,&z2,&t2);
	    if( s->x!=x2 || s->y!=y2 || s->z!=z2 || s->t!=t2 ){
	        printf(
		    "DUMMY! Your gather mapping is not its own inverse\n");
	        printf("It's square mapped %d %d %d %d to %d %d %d %d\n",
		    s->x,s->y,s->z,s->t,x2,y2,z2,t2);
	        terminate(1);
	    }
	}
    }

#endif	/* N_SUBL32 */

    /* RECEIVE LISTS: */
    /* Fill in pointers to sites which are on this node, NOWHERE if
	they are off-node */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
        j = node_number(x,y,z,t);	/* node for neighbor site */
        /* if neighbor is on node, set up pointer */
        if( j == mynode() ) neighbor[dir][i]= node_index(x,y,z,t);
        else		    neighbor[dir][i]= NOWHERE;
    }

    /* make lists of sites which get data from other nodes.  */
    neighborlist[dir] = make_send_receive_list( func, args, want_even_odd,
        FORWARDS, RECEIVE );

#ifndef N_SUBL32

    /* SEND LISTS: */
    /* Now make lists of sites to which we send */
    /* Under some conditions, if mapping is its own inverse we can use
	the lists we have already made */
    if( inverse==OWN_INVERSE && want_even_odd==ALLOW_EVEN_ODD
        && parity_conserve==SWITCH_PARITY ){
        /* Under these conditions the send and receive comlinks are the same */
        /* Just plug head of send list into receive list. */
        neighborlist_send[dir]=neighborlist[dir];
    }  /* end code to use same comlinks for send as for receive */
    else if ( inverse==OWN_INVERSE && 
	( want_even_odd==NO_EVEN_ODD ||
	( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SAME_PARITY ) ) ){
        /* Under these conditions, the even list in a send comlink is the
	    odd list in the receive comlink, and vice versa. */
        neighborlist_send[dir] = copy_list_switch( neighborlist[dir] );
    }
    else{
        /* Make new linked list of comlinks for send lists */
        neighborlist_send[dir] = make_send_receive_list( func, args,
	    want_even_odd, FORWARDS, SEND );
    } /* End general case for send lists */

    if( inverse != WANT_INVERSE) return(dir);

#else	/* N_SUBL32 */

    /* SEND LISTS: */
    /* Now make lists of sites to which we send */
    /* If mapping is its own inverse we can use lists we have already made */
    if( inverse==OWN_INVERSE ){
	/* Under these conditions the send and receive comlinks are the same */
	/* Just plug head of send list into receive list. */
	neighborlist_send[dir]=neighborlist[dir];
    }
    else{
	/* Make new linked list of comlinks for send lists */
	neighborlist_send[dir] = make_send_receive_list( func, args,
	    want_even_odd, FORWARDS, SEND );
    } /* End general case for send lists */

    if( inverse != WANT_INVERSE){
	free(othersubl);
	return(dir);
    }

#endif	/* N_SUBL32 */

    /* INVERSE GATHER */
    /* Now, if necessary, make inverse gather */
    /* In most cases, we can use the same lists as the gather, in one
	form or another.  Of course, by the time you get to here
	you know that inverse = WANT_INVERSE */
    dir++;	/* inverse gather has direction one more than original */

    /* Always set up pointers to sites on this node */
    /* scan sites in lattice */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, BACKWARDS, &x,&y,&z,&t);
        j = node_number(x,y,z,t);	/* node for neighbor site */

        /* if neighbor is on node, set up pointer */
        if( j == mynode() ) neighbor[dir][i]= node_index(x,y,z,t);
        else 		    neighbor[dir][i]= NOWHERE;
    }

#ifndef N_SUBL32

    if( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SWITCH_PARITY ){
        /* Use same comlinks as inverse gather, switching send and receive.
           Nearest neighbor gathers are an example of this case. */
        neighborlist_send[dir]=neighborlist[dir-1];
        neighborlist[dir]=neighborlist_send[dir-1];
    }
    else if( (want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SAME_PARITY)
	|| want_even_odd==NO_EVEN_ODD  ){
        /* make new comlinks, but use same lists as inverse gather, switching
           send and receive, switching even and odd. */
        neighborlist_send[dir] = copy_list_switch( neighborlist[dir-1] );
        neighborlist[dir] = copy_list_switch( neighborlist_send[dir-1] );
    }
    else {  /* general case.  Really only get here if ALLOW_EVEN_ODD
		and SCRAMBLE_PARITY */

        /* RECEIVE LISTS */
        neighborlist[dir] = make_send_receive_list( func, args, want_even_odd,
	    BACKWARDS, RECEIVE );

        /* SEND LISTS: */
        /* Now make lists of sites to which we send */
        neighborlist_send[dir] = make_send_receive_list( func, args,
	    want_even_odd, BACKWARDS, SEND );
    } /* End making new lists for inverse gather */

#else	/* N_SUBL32 */

    if( (want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SAME_PARITY )
	|| want_even_odd==NO_EVEN_ODD  ){
        /* Use same comlinks as inverse gather, switching send and receive. */
        neighborlist_send[dir]=neighborlist[dir-1];
        neighborlist[dir]=neighborlist_send[dir-1];
    }
    else if( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SWITCH_PARITY){
        /* make new comlinks, but use same lists as inverse gather, switching
           send and receive, switching sublattices. */
        neighborlist_send[dir] = copy_list_switch( neighborlist[dir-1],
	    othersubl );
        neighborlist[dir] = copy_list_switch( neighborlist_send[dir-1],
	    othersubl );
    }
    else {  /* general case.  Really only get here if ALLOW_EVEN_ODD
		and SCRAMBLE_PARITY */
	printf("Should not get to HERE!\n");
	terminate(1);
/**	neighborlist[dir] = make_send_receive_list( func, args,
            want_even_odd, BACKWARDS, RECEIVE );
	neighborlist_send[dir] = make_send_receive_list( func, args,
            want_even_odd, BACKWARDS, SEND ); **/
    }  /* End making new lists for inverse gather */

    free(othersubl);
#endif	/* N_SUBL32 */

    return(dir-1);
}


/**********************************************************************/
#ifndef N_SUBL32

comlink *  make_send_receive_list(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
        		/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int forw_back,		/* FORWARDS or BACKWARDS */
 int send_recv)		/* SEND or RECEIVE list */
{
register int i,j;	/* scratch */
register site *s;	/* scratch */
int x,y,z,t;		/* coordinates */
int parity;		/*if send, parity of site on other node */
			/*if receive, parity of site on this node */
int *ebuf,*obuf;	/* to be malloc'd */
comlink **combuf;	/* to be malloc'd, remember where comlinks are */
register comlink *compt,**comptpt;
comlink *firstpt;

void sort_site_list(
 int n,		/* number of elements in list */
 int *list,	/* pointer to list */
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
                /* function which defines mapping */
 int *args,	/* arguments to pass to function */
 int forw_back);/* look forwards or backwards in map */

    /* make temporary buffers of numnodes() integers to count numbers of
       even and odd neighbors on each node */
    ebuf = (int *)malloc( numnodes()*sizeof(int) );
    obuf = (int *)malloc( numnodes()*sizeof(int) );
    combuf = (comlink **)malloc( numnodes()*sizeof(comlink *) );

    /* clear neighbor_numbers */
    for(i=0;i<numnodes();i++){
        ebuf[i] = obuf[i] = 0;
    }

    /* scan sites in lattice */
    FORALLSITES(i,s){
        /* find coordinates, node, and parity of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    parity = s->parity;
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    if( (x+y+z+t)%2==0 )parity=EVEN; else parity=ODD;
	}
	j = node_number(x,y,z,t);

        /* if site is off node, increment neighbor_counter */
        if( j != mynode() ){
	    if( send_recv==RECEIVE ){
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD) ebuf[j]++;
	        else   obuf[j]++;
	    }
	    else{
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD) obuf[j]++;
	        else   ebuf[j]++;
	    }
	}
    }

    firstpt=NULL;
    comptpt = &firstpt;
    /* for each neighbor_counter that is nonzero, create a comlink */
    for(j=0;j<numnodes();j++){
        if( j==mynode() )continue;	/* not for local node */
        if( ebuf[j]==0 && obuf[j]==0)continue;
	    /* no neighbors on this node */

        compt = (comlink *)malloc( sizeof(comlink) );
        *comptpt = compt;
        combuf[j] = compt;	/* to make it easy to find again */
        compt->nextcomlink=NULL;	/* currently terminates list */
        compt->othernode=j;
        compt->n_even_connected = ebuf[j];
        compt->n_odd_connected = obuf[j];
        compt->esitelist = (int *)malloc(
	    (ebuf[j]+obuf[j])*sizeof(int) );
        compt->ositelist = (compt->esitelist)+ebuf[j];
	    /*esitelist and ositelist must be filled in later */

        comptpt = &(compt->nextcomlink);	/* linked list, if we
	    extend it this will get address of next comlink. */
    }

    /* clear neighbor_numbers, to be used as counters now */
    for(j=0;j<numnodes();j++){
        ebuf[j]=obuf[j]=0;
    }

    /* scan sites in node again */
    FORALLSITES(i,s){
        /* find coordinates, node, and parity of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    parity = s->parity;
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    if( (x+y+z+t)%2==0 )parity=EVEN; else parity=ODD;
	}
	j = node_number(x,y,z,t);

        /* if neighbor is offnode, add to list in appropriate comlink */
        if( j != mynode() ){
	    if(send_recv==RECEIVE ){
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD ){
	            combuf[j]->esitelist[ebuf[j]] = i;
	            ebuf[j]++;
	        }
	        else{
	            combuf[j]->ositelist[obuf[j]] = i;
	            obuf[j]++;
	        }
	    }
	    else { /*SEND*/
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD ){
	            combuf[j]->ositelist[obuf[j]] = i;
	            obuf[j]++;
	        }
	        else{
	            combuf[j]->esitelist[ebuf[j]] = i;
	            ebuf[j]++;
	        }
	    }
        }
    }
    /* sort the lists of links according to the ordering of their
       even neighbors in the lower numbered node.  The list of sites
       on the lower numbered node is already in order. */
    for(compt=firstpt; compt != NULL; compt=compt->nextcomlink){
        if(compt->othernode > this_node)continue;
	    /* this is lower numbered node, so don't sort */
	if( send_recv==RECEIVE ) i = forw_back;
	else i = -forw_back;
	sort_site_list(compt->n_odd_connected, compt->ositelist,
	    func, args, i);
	sort_site_list(compt->n_even_connected, compt->esitelist,
	    func, args, i);
    }

    /* free temporary storage */
    free(ebuf); free(obuf); free(combuf);
    return(firstpt);
}
#else	/* N_SUBL32 */

comlink *  make_send_receive_list(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
			/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int forw_back,		/* FORWARDS or BACKWARDS */
 int send_recv)		/* SEND or RECEIVE list */
{
register int i,j;	/* scratch */
register site *s;	/* scratch */
int x,y,z,t;		/* coordinates */
int subl;		/* sublattice of site on this node */
int subl_s;		/* sublattice of site on sending node */
int subl_r;		/* sublattice of site on receiving node */
int *send_subl;		/* sublattice on sending node */
int *sbuf[N_SUBL32];	/* to be malloc'd */
int *tbuf;		/* to be malloc'd */
comlink **combuf;	/* to be malloc'd, remember where comlinks are */
register comlink *compt,**comptpt;
comlink *firstpt;

void sort_site_list(
 int n,		/* number of elements in list */
 int *list,	/* pointer to list */
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
		/* function which defines mapping */
 int *args,	/* arguments to pass to function */
 int forw_back);/* look forwards or backwards in map */

    /* make temporary buffers of numnodes() integers to count numbers of
       neighbors in each sublattice on each node */
    for(subl=0;subl<N_SUBL32;subl++){
	sbuf[subl] = (int *)malloc( numnodes()*sizeof(int) );
	/* clear neighbor_numbers */
	for(i=0;i<numnodes();i++) sbuf[subl][i] = 0;
    }
    tbuf = (int *)malloc( numnodes()*sizeof(int) );
    for(i=0;i<numnodes();i++) tbuf[i] = 0;
    combuf = (comlink **)malloc( numnodes()*sizeof(comlink *) );

    send_subl = (int *)malloc(N_SUBL32*sizeof(int));
    if(send_subl== NULL){
	printf("NODE %d: no room for send_subl\n",this_node);
	terminate(1);
    }
    for(subl=0;subl<N_SUBL32;subl++) send_subl[subl] = NOWHERE;

    /* scan sites in lattice */
    FORALLSITES(i,s){
        /* find coordinates, node, and sublattice of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    subl = subl_r = s->parity;
	    subl_s = (x%2) + 2*(y%2) + 4*(z%2) + 8*(t%2);
	    subl_s += 16*((x/2+y/2+z/2+t/2)%2);
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    subl = subl_s = s->parity;
	    subl_r = (x%2) + 2*(y%2) + 4*(z%2) + 8*(t%2);
	    subl_r += 16*((x/2+y/2+z/2+t/2)%2);
	}
	j = node_number(x,y,z,t);

	if( want_even_odd==ALLOW_EVEN_ODD ){
	    if( send_subl[subl_r] == NOWHERE ){
		send_subl[subl_r] = subl_s;
	    }
	    else if( send_subl[subl_r] != subl_s ){
		printf("Mixing up sublattices: %d vs %d\n",
		    send_subl[subl_r], subl_s);
		terminate(1);
	    }
	}

        /* if site is off node, increment neighbor_counter */
        if( j != mynode() ){
	    tbuf[j]++;
	    if( want_even_odd==ALLOW_EVEN_ODD) sbuf[subl][j]++;
	    else sbuf[0][j]++;
	}
    }

    firstpt=NULL;
    comptpt = &firstpt;
    /* for each neighbor_counter that is nonzero, create a comlink */
    for(j=0;j<numnodes();j++){
        if( j==mynode() )continue;	/* not for local node */
        if( tbuf[j]==0 )continue;
	    /* no neighbors on this node */

        compt = (comlink *)malloc( sizeof(comlink) );
        *comptpt = compt;
        combuf[j] = compt;	/* to make it easy to find again */
        compt->nextcomlink=NULL;	/* currently terminates list */
        compt->othernode=j;
	for(subl=0;subl<N_SUBL32;subl++){
	    compt->n_subl_connected[subl] = sbuf[subl][j];
	    compt->send_subl[subl] = send_subl[subl];
	}
        compt->n_connected = tbuf[j];
        compt->sitelist[0] = (int *)malloc(tbuf[j]*sizeof(int) );
	for(subl=1;subl<N_SUBL32;subl++)
	    compt->sitelist[subl] = (compt->sitelist[subl-1])+sbuf[subl-1][j];
	/*sitelist[...] must be filled in later */

        comptpt = &(compt->nextcomlink);	/* linked list, if we
	    extend it this will get address of next comlink. */
    }

    /* clear neighbor_numbers, to be used as counters now */
    for(subl=0;subl<N_SUBL32;subl++){
	for(i=0;i<numnodes();i++) sbuf[subl][i] = 0;
    }
    for(i=0;i<numnodes();i++) tbuf[i] = 0;

    /* scan sites in node again */
    FORALLSITES(i,s){
        /* find coordinates, node, and sublattice of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    subl = s->parity;
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    subl = s->parity;
	}
	j = node_number(x,y,z,t);

        /* if neighbor is offnode, add to list in appropriate comlink */
        if( j != mynode() ){
	    if( want_even_odd==ALLOW_EVEN_ODD ){
		combuf[j]->sitelist[subl][sbuf[subl][j]] = i;
		sbuf[subl][j]++;
	    }
	    else{
		combuf[j]->sitelist[0][sbuf[0][j]] = i;
		sbuf[0][j]++;
	    }
        }
    }
    /* sort the lists of links according to the ordering of their
       even neighbors in the lower numbered node.  The list of sites
       on the lower numbered node is already in order. */
    for(compt=firstpt; compt != NULL; compt=compt->nextcomlink){
        if(compt->othernode > this_node)continue;
	    /* this is lower numbered node, so don't sort */
	if( send_recv==RECEIVE ) i = forw_back;
	else i = -forw_back;
	for(subl=0;subl<N_SUBL32;subl++)
	    sort_site_list(compt->n_subl_connected[subl],
		compt->sitelist[subl], func, args, i);
    }

    /* free temporary storage */
    free(send_subl); free(tbuf); free(combuf);
    for(subl=0;subl<N_SUBL32;subl++) free(sbuf[subl]);

    return(firstpt);
}
#endif	/* N_SUBL32 */

/**********************************************************************/
#ifndef N_SUBL32

comlink * copy_list_switch( comlink *old_compt ){
comlink *firstpt,*compt;
    /* copy a linked list of comlinks, switching even and odd */
    if( old_compt==NULL )return(NULL);
    firstpt = compt = (comlink *)malloc( sizeof(comlink) );
    do{
        compt->othernode=old_compt->othernode;
        compt->n_even_connected = old_compt->n_odd_connected;
        compt->n_odd_connected = old_compt->n_even_connected;
        compt->esitelist = old_compt->ositelist;
        compt->ositelist = old_compt->esitelist;

        if( old_compt->nextcomlink != NULL)
	    compt->nextcomlink = (comlink *)malloc( sizeof(comlink) );
	else compt->nextcomlink = NULL;
	old_compt=old_compt->nextcomlink; 
	compt = compt->nextcomlink;
    } while( old_compt!=NULL );
    return(firstpt);
}
#else	/* N_SUBL32 */

comlink * copy_list_switch( comlink *old_compt, int *othersubl) {
comlink *firstpt,*compt;
int subl, subl_s, subl_r;
    /* copy a linked list of comlinks, switching even and odd */
    if( old_compt==NULL )return(NULL);
    firstpt = compt = (comlink *)malloc( sizeof(comlink) );
    do{
        compt->othernode=old_compt->othernode;
        compt->n_connected = old_compt->n_connected;
	for(subl=0;subl<N_SUBL32;subl++){
	    subl_s = old_compt->send_subl[subl];
	    if( subl_s == othersubl[subl] ){
		/* switch from receive to send */
		subl_r = subl_s;
		subl_s = subl;
	    }
	    else{
		/* switch from send to receive */
		subl_r = subl;
		subl_s = othersubl[subl];
	    }
	    compt->send_subl[subl_r] = subl_s;
	    compt->n_subl_connected[subl] = old_compt->n_subl_connected[subl];
	    compt->sitelist[subl] = old_compt->sitelist[subl];
	}

        if( old_compt->nextcomlink != NULL)
	    compt->nextcomlink = (comlink *)malloc( sizeof(comlink) );
	else compt->nextcomlink = NULL;
	old_compt=old_compt->nextcomlink; 
	compt = compt->nextcomlink;
    } while( old_compt!=NULL );
    return(firstpt);
}
#endif	/* N_SUBL32 */

/**********************************************************************/
/* sort a list of sites according to the order of the sites on the
   node with which they comunicate */
void sort_site_list(
 int n,		/* number of elements in list */
 int *list,	/* pointer to list */
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
                /* function which defines mapping */
 int *args,	/* arguments to pass to function */
 int forw_back)	/* look forwards or backwards in map */
{
register int j,k,in1,in2,flag;
register site *s;
int x,y,z,t;
int *key;

    if(n==0)return;
    key = (int *)malloc(n*sizeof(int));
    if(key == NULL)
      {
	printf("sort_site_list(%d): no room for key\n",mynode());
	terminate(1);
      }

    /* Construct sort key */
    for(j=0;j<n;j++)
      {
	s = &(lattice[list[j]]);
	func(s->x,s->y,s->z,s->t,args,forw_back,&x,&y,&z,&t);
	key[j] = node_index(x,y,z,t);
      }

    /* bubble sort, if this takes too long fix it later */
    for(j = n-1; j>0; j--){
        flag=0;
	for(k=0; k<j; k++){
	    in1 = key[k];
	    in2 = key[k+1];
	    if(in1>in2){
		flag=1;
		key[k]   = in2;
		key[k+1] = in1;
		in1 = list[k];
		list[k]=list[k+1];
		list[k+1]=in1;
	    }
	}
	if(flag==0)break;
    }
    free(key);
}

/**********************************************************************/
/* utility function for finding coordinates of neighbor */
void neighbor_coords(
 int x, int y, int z, int t,  /* coordinates of site */
 int dir,	              /* direction (eg XUP) */
 int *x2p, int *y2p, int *z2p, int *t2p)
                             /* pointers to coordinates of neighbor */
{
    *x2p = x; *y2p = y; *z2p = z; *t2p = t;
    switch(dir){
	case XUP: *x2p = (x+1)%nx; break;
	case XDOWN: *x2p = (x+nx-1)%nx; break;
	case YUP: *y2p = (y+1)%ny; break;
	case YDOWN: *y2p = (y+ny-1)%ny; break;
	case ZUP: *z2p = (z+1)%nz; break;
	case ZDOWN: *z2p = (z+nz-1)%nz; break;
	case TUP: *t2p = (t+1)%nt; break;
	case TDOWN: *t2p = (t+nt-1)%nt; break;
	default: printf("BOTCH: bad direction\n"); exit(1);
    }
}

/**********************************************************************/
/* Set up interrupt handlers to handle field_pointer routines
   make_lattice() must be called first. */
static msg_request mreqbuf;     /* global so handler can use it too */
void start_handlers(){
void fillfieldrequest(int type,int size,int node,int pid);
    /**
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    **/
    if(this_node==0)printf(
	"Can't start interrupt handler yet!!, no field_pointer\n");
}
/**********************************************************************/
/* The handler for field requests.  It fills the requests, then posts a
   receive for the next request. */
void fillfieldrequest(int type,int size,int node,int pid){
register int trapstate;
void fillfieldrequest();

    /**
    trapstate=masktrap(1);
    csend( FIELD_REPLY, F_PT( &(lattice[mreqbuf.index]), mreqbuf.field ),
	mreqbuf.size, node, pid);
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    masktrap(trapstate);
    **/
    printf("Oops: called fillfieldrequest\n"); terminate(0);
}


/**********************************************************************/
/* GATHER ROUTINES */
/* start_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_gather() and cleanup_gather() calls.
   This list contains msg_tags for all receive buffers, followed by
   a msg_tag with msg_id = MPI_REQUEST_NULL and msg_buf = NULL,
   followed by msg_tags for all send buffers, followed by a msg_tag
   with id=MPI_REQUEST_NULL, buf=NULL.
   If no messages at all are required, the routine will return NULL.
   msg_buf=NULL should be a reliable indicator of no message.

   usage:  tag = start_gather( source, size, direction, parity, dest )
   example:
	msg_tag *tag;
	tag = start_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0] );
	  ** do other stuff **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here.
	  **
	cleanup_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

    Under certain circumstances it is possible to efficiently gather
    a field that has been previously gathered.  This happens when the
    field being gathered has been modified, but the pointers (the
    destination of start_gather() ) have not been modified.  To use
    restart_gather, the original gather must have been waited for but
    not cleaned up.  The usage is:
	msg_tag *tag;
	tag = start_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0] );
	  ** do other stuff, but don't modify tag or gen_pt[0] **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here, but don't modify tag or
	   gen_pt[0].
	   Do modify the source field phi. **
	restart_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0], tag );
	  ** do other stuff **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the modified phi.
	   The restart-wait may be repeated as often as desired.  **
	cleanup_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

   Internal convention for message types is (schematically):
	( f_offset(src) | dir ) + GATHER_BASE
*/

#ifndef N_SUBL32

msg_tag * start_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
msg_tag *mbuf;		/* list of message tags, to be returned */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX; /* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - chase a linked list. */
    for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	n_recv_msgs++) compt = compt->nextcomlink;
    for(n_send_msgs=0, compt= neighborlist_send[index];
	compt != NULL; n_send_msgs++) compt = compt->nextcomlink;

    /* allocate a buffer for the msg_tags.  This is dynamically allocated
       because there may be an arbitrary number of gathers in progress
       in any direction. */
    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0, compt= neighborlist[index]; compt != NULL;
	i++,compt = compt->nextcomlink){

	/* allocate buffer to receive neighbors */
	switch (parity){
	    case EVEN:
		nsites = compt->n_even_connected;
		sitelist = compt->esitelist;
		break;
	    case ODD:
		nsites = compt->n_odd_connected;
		sitelist = compt->ositelist;
		break;
	    case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
		sitelist = compt->esitelist;
		break;
	}
	mbuf[i].msg_node = compt->othernode;
	mbuf[i].msg_size =  nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
	     dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = ((char *)(lattice+neighbor[index][j]))+field;
	    }}
	    break;
	case ODD:
	    FORODDSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = ((char *)(lattice+neighbor[index][j]))+field;
	    }}
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = ((char *)(lattice+neighbor[index][j]))+field;
	    }}
	    break;
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Allocate buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
        switch (parity){
            case EVEN:
		nsites = compt->n_odd_connected;
                break;
            case ODD:
		nsites = compt->n_even_connected;
                break;
            case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
                break;
        }
	tpt=(char *)malloc( nsites*size );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i].msg_node=compt->othernode;
	mbuf[i].msg_size=nsites*size;
	mbuf[i].msg_buf=tpt;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }
    /* return */
    return(mbuf);
}

/**********************************************************************/
/* Repeat a gather with the same source and destination as a
  previous gather.  The previous gather must have been waited for
  but not cleaned up.  Pointers to sites on the same node are not
  reset, and the same buffers are reused. */
void restart_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf)          /* previously returned by start_gather */
{
/* local variables */
register int i,j;	/* scratch */
register char *tpt;	/* scratch pointer in buffers */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX; /* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* use same buffer to receive neighbors */
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node, GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
}

/**********************************************************************/
msg_tag * start_gather_from_temp(
/* arguments */
 void * field,		/* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
msg_tag *mbuf;		/* list of message tags, to be returned */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX; /* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - chase a linked list. */
    for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	n_recv_msgs++) compt = compt->nextcomlink;
    for(n_send_msgs=0, compt= neighborlist_send[index];
	compt != NULL; n_send_msgs++) compt = compt->nextcomlink;

    /* allocate a buffer for the msg_tags.  This is dynamically allocated
       because there may be an arbitrary number of gathers in progress
       in any direction. */
    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0, compt= neighborlist[index]; compt != NULL;
	i++,compt = compt->nextcomlink){

	/* allocate buffer to receive neighbors */
	switch (parity){
	    case EVEN:
		nsites = compt->n_even_connected;
		sitelist = compt->esitelist;
		break;
	    case ODD:
		nsites = compt->n_odd_connected;
		sitelist = compt->ositelist;
		break;
	    case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
		sitelist = compt->esitelist;
		break;
	}
	mbuf[i].msg_node = compt->othernode;
	mbuf[i].msg_size =  nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
	     dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = (char *)field + neighbor[index][j]*size;
	    }}
	    break;
	case ODD:
	    FORODDSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = (char *)field + neighbor[index][j]*size;
	    }}
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ if(neighbor[index][j] != NOWHERE){
		dest[j] = (char *)field + neighbor[index][j]*size;
	    }}
	    break;
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Allocate buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
        switch (parity){
            case EVEN:
		nsites = compt->n_odd_connected;
                break;
            case ODD:
		nsites = compt->n_even_connected;
                break;
            case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
                break;
        }
	tpt=(char *)malloc( nsites*size );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i].msg_node=compt->othernode;
	mbuf[i].msg_size=nsites*size;
	mbuf[i].msg_buf=tpt;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, (char *)field + compt->ositelist[j]*size, size);
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, (char *)field + compt->esitelist[j]*size, size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }
    /* return */
    return(mbuf);
} /*start_gather_from_temp()*/

/**********************************************************************/
void restart_gather_from_temp(
/* arguments */
 void * field,		/* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf)          /* previously returned by start_gather */
{
/* local variables */
register int i,j;	/* scratch */
register char *tpt;	/* scratch pointer in buffers */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX; /* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* use same buffer to receive neighbors */
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, (char *)field + compt->ositelist[j]*size, size);
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, (char *)field + compt->esitelist[j]*size, size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node, GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
} /*restart_gather_from_temp() */

#else	/* N_SUBL32 */

msg_tag * start_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int subl_s;		/* sublattice on sending node */
int subl_r;		/* sublattice on receiving node */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
msg_tag *mbuf;		/* list of message tags, to be returned */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - chase a linked list. */
    if( subl == EVENANDODD ){
	for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	    n_recv_msgs++) compt = compt->nextcomlink;
	for(n_send_msgs=0, compt= neighborlist_send[index];
	    compt != NULL; n_send_msgs++) compt = compt->nextcomlink;
    }
    else{
	for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	    compt = compt->nextcomlink){
	    if( compt->n_subl_connected[subl] > 0) n_recv_msgs++;
	}
	for(n_send_msgs=0, compt= neighborlist_send[index];
	    compt != NULL; compt = compt->nextcomlink){
	    subl_s = compt->send_subl[subl];
	    if( compt->n_subl_connected[subl_s] > 0) n_send_msgs++;
	}
    }

    /* allocate a buffer for the msg_tags.  This is dynamically allocated
       because there may be an arbitrary number of gathers in progress
       in any direction. */
    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0, compt= neighborlist[index]; compt != NULL;
	compt = compt->nextcomlink){

	/* allocate buffer to receive neighbors */
	if( subl == EVENANDODD ){
	    nsites = compt->n_connected;
	    sitelist = compt->sitelist[0];
	}
	else{
	    nsites = compt->n_subl_connected[subl];
	    sitelist = compt->sitelist[subl];
	}
	if(nsites == 0)continue;	/* Nothing to receive for this sublattice */
	mbuf[i].msg_node = compt->othernode;
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
        MPI_Irecv( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
	i++;
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
	mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    if( subl == EVENANDODD ){
	FORALLSITES(j,s){ if(neighbor[index][j] != NOWHERE){
	    dest[j] = ((char *)(lattice+neighbor[index][j]))+field;
	}}
    }
    else{
	FORSOMESUBLATTICE(j,s,subl){ if(neighbor[index][j] != NOWHERE){
	    dest[j] = ((char *)(lattice+neighbor[index][j]))+field;
	}}
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL; 
	compt = compt->nextcomlink){
	/* Allocate buffer to gather data. */
	if( subl == EVENANDODD ){
	    nsites = compt->n_connected;
	}
	else{
	    subl_s = compt->send_subl[subl];
	    nsites = compt->n_subl_connected[subl_s];
	}
	if(nsites == 0)continue;	/* Nothing to send for this sublattice */
	tpt=(char *)malloc( nsites*size );
	if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i].msg_node=compt->othernode;
	mbuf[i].msg_size=nsites*size;
	mbuf[i].msg_buf=tpt;
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    s = &(lattice[compt->sitelist[0][j]]);
		    memcpy( tpt, ((char *)(s)+field), size);
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    s = &(lattice[compt->sitelist[subl_s][j]]);
		    memcpy( tpt, ((char *)(s)+field), size);
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		s = &(lattice[compt->sitelist[subl_s][j]]);
		memcpy( tpt, ((char *)(s)+field), size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
	i++;
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
}

/**********************************************************************/
/* Repeat a gather with the same source and destination as a
  previous gather.  The previous gather must have been waited for
  but not cleaned up.  Pointers to sites on the same node are not
  reset, and the same buffers are reused. */
void restart_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf)          /* previously returned by start_gather */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int subl_s;		/* sublattice on sending node */
int subl_r;		/* sublattice on receiving node */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* use same buffer to receive neighbors */
	/* post receive */
        MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	compt = compt->nextcomlink){
	if( subl != EVENANDODD ){
	    subl_s = compt->send_subl[subl];
	    if(compt->n_subl_connected[subl_s] == 0)continue;
	}
	/* Use same buffer to gather data. */
	tpt=mbuf[i].msg_buf;
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    s = &(lattice[compt->sitelist[0][j]]);
		    memcpy( tpt, ((char *)(s)+field), size);
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    s = &(lattice[compt->sitelist[subl_s][j]]);
		    memcpy( tpt, ((char *)(s)+field), size);
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		s = &(lattice[compt->sitelist[subl_s][j]]);
		memcpy( tpt, ((char *)(s)+field), size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node, GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
	i++;
    }
}

/**********************************************************************/
msg_tag * start_gather_from_temp(
/* arguments */
 void * field,		/* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int subl_s;		/* sublattice on sending node */
int subl_r;		/* sublattice on receiving node */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
msg_tag *mbuf;		/* list of message tags, to be returned */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - chase a linked list. */
    if( subl == EVENANDODD ){
	for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	    n_recv_msgs++) compt = compt->nextcomlink;
	for(n_send_msgs=0, compt= neighborlist_send[index];
	    compt != NULL; n_send_msgs++) compt = compt->nextcomlink;
    }
    else{
	for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	    compt = compt->nextcomlink){
	    if( compt->n_subl_connected[subl] > 0) n_recv_msgs++;
	}
	for(n_send_msgs=0, compt= neighborlist_send[index];
	    compt != NULL; compt = compt->nextcomlink){
	    subl_s = compt->send_subl[subl];
	    if( compt->n_subl_connected[subl_s] > 0) n_send_msgs++;
	}
    }

    /* allocate a buffer for the msg_tags.  This is dynamically allocated
       because there may be an arbitrary number of gathers in progress
       in any direction. */
    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0, compt= neighborlist[index]; compt != NULL;
	compt = compt->nextcomlink){

	/* allocate buffer to receive neighbors */
	if( subl == EVENANDODD ){
	    nsites = compt->n_connected;
	    sitelist = compt->sitelist[0];
	}
	else{
	    nsites = compt->n_subl_connected[subl];
	    sitelist = compt->sitelist[subl];
	}
	if(nsites == 0)continue;	/* Nothing to receive for this sublattice */
	mbuf[i].msg_node = compt->othernode;
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
        MPI_Irecv( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
	i++;
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
	mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    if( subl == EVENANDODD ){
	FORALLSITES(j,s){ if(neighbor[index][j] != NOWHERE){
	    dest[j] = (char *)field + neighbor[index][j]*size;
	}}
    }
    else{
	FORSOMESUBLATTICE(j,s,subl){ if(neighbor[index][j] != NOWHERE){
	    dest[j] = (char *)field + neighbor[index][j]*size;
	}}
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL; 
	compt = compt->nextcomlink){
	/* Allocate buffer to gather data. */
	if( subl == EVENANDODD ){
	    nsites = compt->n_connected;
	}
	else{
	    subl_s = compt->send_subl[subl];
	    nsites = compt->n_subl_connected[subl_s];
	}
	if(nsites == 0)continue;	/* Nothing to send for this sublattice */
	tpt=(char *)malloc( nsites*size );
	if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i].msg_node=compt->othernode;
	mbuf[i].msg_size=nsites*size;
	mbuf[i].msg_buf=tpt;
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    memcpy( tpt, (char *)field + compt->sitelist[0][j]*size,
			size);
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
			size);
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
		    size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, nsites*size, MPI_BYTE,
	    compt->othernode,GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
	i++;
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
} /*start_gather_from_temp()*/

/**********************************************************************/
void restart_gather_from_temp(
/* arguments */
 void * field,		/* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest,		/* one of the vectors of pointers */
 msg_tag *mbuf)          /* previously returned by start_gather */
{
/* local variables */
register int i,j;	/* scratch */
register char *tpt;	/* scratch pointer in buffers */
int subl_s;		/* sublattice on sending node */
int subl_r;		/* sublattice on receiving node */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* use same buffer to receive neighbors */
	/* post receive */
        MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node,  GATHER_MSG_TYPE(type_offset), MPI_COMM_WORLD,
	    &mbuf[i].msg_id );
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	compt = compt->nextcomlink){
	if( subl != EVENANDODD ){
	    subl_s = compt->send_subl[subl];
	    if(compt->n_subl_connected[subl_s] == 0)continue;
	}
	/* Use same buffer to gather data. */
	tpt=mbuf[i].msg_buf;
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    memcpy( tpt, (char *)field + compt->sitelist[0][j]*size,
			size);
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
			size);
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
		    size);
	    }
	}
	/* start the send */
	MPI_Issend( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE,
	    mbuf[i].msg_node, GATHER_MSG_TYPE(type_offset),
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
	i++;
    }
} /*restart_gather_from_temp() */

#endif	/* N_SUBL32 */

/**********************************************************************/
void wait_gather(msg_tag *mbuf) {
register int i,i0;
MPI_Status status;
    if(mbuf==NULL)return;
    /* wait for all receive messages */
    for(i=0; mbuf[i].msg_buf != NULL; i++)MPI_Wait( &mbuf[i].msg_id, &status );
    /* wait for all send messages */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	MPI_Wait( &mbuf[i].msg_id, &status );
    }
}

/**********************************************************************/
void cleanup_gather(msg_tag *mbuf) {
register int i,i0;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /*  free all send buffers */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	free( mbuf[i].msg_buf );
    }
    /* free the msg_tag buffer */
    free(mbuf);
}


/**********************************************************************/
/* GENERAL_GATHER ROUTINES */
/* start_general_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_general_gather() and
   cleanup_general_gather() calls.
   This list contains msg_tags for all receive buffers, followed by
   a msg_tag with msg_id = MPI_REQUEST_NULL and msg_buf = NULL,
   followed by msg_tags for all send buffers, followed by a msg_tag
   with id=MPI_REQUEST_NULL, buf=NULL.
   If no messages at all are required, the routine will return NULL.
   msg_buf=NULL should be a reliable indicator of no message.

   usage:  tag = start_general_gather( source, size, displacement, parity, dest)
   example:
	msg_tag *tag;
	int disp[4]; 
        disp[XUP]=1; disp[YUP]= -1; disp[ZUP] = disp[TUP] = 0;
	tag = start_general_gather( F_OFFSET(phi), sizeof(su3_vector), disp,
	    EVEN, gen_pt[0] );
	  ** do other stuff **
	wait_general_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here.
	  **
	cleanup_general_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

*/
struct msg_tmp { int node,count; }; /* temporary structure for keeping
	track of messages to be sent or received */
int g_gather_flag=0;	/* flag to tell if general gather in progress */
struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
int tsize;		/* size of entry in messages = 2*sizeof(int)+size */
char ** tdest;		/* tdest is copy of dest */
/* from_nodes, tsize and tdest are global because they are set in 
   start_general_gather() and used in wait_general_gather().  This
   works because we allow only one general_gather in progress at a
   time. */

#ifndef N_SUBL32

msg_tag * start_general_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int disp_parity;	/* parity of displacement vector */
int send_parity;	/* parity of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
msg_tag *mbuf;		/* list of message tags, to be returned */
int n_send_msgs, n_recv_msgs;


    /* check for gather already in progress */
    if(g_gather_flag!=0){
	fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
	exit(1);
    }
    n_recv_msgs = n_send_msgs = 0;
    tsize = 2*sizeof(int)+size;
	/* Use 2*sizeof int so pointer will be aligned to double word */
    tdest = dest;
    /* find parity of sites that may be sent */
    if( (displacement[XUP]+displacement[YUP]+displacement[ZUP]+
	displacement[TUP])%2 == 0)disp_parity=EVEN;
    else disp_parity=ODD;
    switch(parity){
	case EVENANDODD: send_parity=EVENANDODD; break;
	case EVEN:
	    if( disp_parity==EVEN )send_parity=EVEN;
	    else send_parity=ODD;
	    break;
	case ODD:
	    if( disp_parity==EVEN )send_parity=ODD;
	    else send_parity=EVEN;
	    break;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) Make
	list of nodes from whom we expect messages */
    FORSOMEPARITY(i,s,parity){
	if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode==this_node ){
	    dest[i] = (char *)( &lattice[node_index(tx,ty,tz,tt)]) + field;
	}
	else{
	    for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
	    if(j < n_recv_msgs){
		from_nodes[j].count++;
	    }
	    else {
	        if(n_recv_msgs==0){
		    from_nodes = (struct msg_tmp *)malloc(
			sizeof(struct msg_tmp) );
		    from_nodes[0].node = othernode;
		    from_nodes[0].count = 1;
		    n_recv_msgs++;
	        }
	        else{
		    from_nodes = (struct msg_tmp *)realloc( from_nodes,
			(n_recv_msgs+1)*sizeof(struct msg_tmp) );
		    from_nodes[j].node = othernode;
		    from_nodes[j].count = 1;
		    n_recv_msgs++;
	        }
	    }
	}
    }

    /* scan sites of parity we are sending, make list of nodes to which
	we must send messages and the number of messages to each. */
    FORSOMEPARITY(i,s,send_parity){
	if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    if(j < n_send_msgs){
		to_nodes[j].count++;
	    }
	    else {
	        if(n_send_msgs==0){
		    to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
		    to_nodes[0].node = othernode;
		    to_nodes[0].count = 1;
		    n_send_msgs++;
	        }
	        else{
		    to_nodes = (struct msg_tmp *)realloc( to_nodes,
			(n_send_msgs+1)*sizeof(struct msg_tmp) );
		    to_nodes[j].node = othernode;
		    to_nodes[j].count = 1;
		    n_send_msgs++;
	        }
	    }
	}
    }

    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = (char *)malloc( nsites*tsize );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*tsize, MPI_BYTE,
	    from_nodes[i].node, GENERAL_GATHER_MSG_TYPE, 
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=(char *)malloc( to_nodes[i].count*tsize );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
    }

    /* reset to_node counters */
    for(i=0;i<n_send_msgs;i++)to_nodes[i].count=0;
    /* gather data into the buffers. Each entry in the buffers consists
	of the index of the site to which the data is sent, followed by
	the actual data */
    FORSOMEPARITY(i,s,send_parity){
	tx = (s->x - displacement[XUP] + nx)%nx;
	ty = (s->y - displacement[YUP] + ny)%ny;
	tz = (s->z - displacement[ZUP] + nz)%nz;
	tt = (s->t - displacement[TUP] + nt)%nt;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		to_nodes[j].count*tsize;
	    *(int *)tpt = node_index(tx,ty,tz,tt);
	        /* index of site on other node */
	    memcpy( tpt+2*sizeof(int), F_PT(s,field), size);
	    to_nodes[j].count++;
	}
    }

    /* start the sends */
    for(i=0;i<n_send_msgs;i++){
	nsites = to_nodes[i].count;
	MPI_Issend( mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, MPI_BYTE,
	    to_nodes[i].node, GENERAL_GATHER_MSG_TYPE, 
	    MPI_COMM_WORLD, &mbuf[i+n_recv_msgs+1].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

#else	/* N_SUBL32 */

msg_tag * start_general_gather(
/* arguments */
 field_offset field,	/* which field? Some member of structure "site" */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int send_subl;		/* sublattice of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
msg_tag *mbuf;		/* list of message tags, to be returned */
int n_send_msgs, n_recv_msgs;

    /* check for gather already in progress */
    if(g_gather_flag!=0){
	fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
	exit(1);
    }
    n_recv_msgs = n_send_msgs = 0;
    tsize = 2*sizeof(int)+size;
        /* Use 2*sizeof int so pointer will be aligned to double word */
    tdest = dest;

    /* find sublattice of sites that may be sent */
    /* This is not needed for EVENANDODD */
    if( subl != EVENANDODD ){
	/* Displacments by multiples of 4 in any direction does
	   not change sublattice */
	tx = displacement[XUP]%4;
	ty = displacement[YUP]%4;
	tz = displacement[ZUP]%4;
	tt = displacement[TUP]%4;
	send_subl = subl;
	if( tx < 0 ){
	    for(i=0;i<(-tx);i++) send_subl = neighsubl[send_subl][XDOWN];
	}
	else
	    for(i=0;i<tx;i++) send_subl = neighsubl[send_subl][XUP];
	if( ty < 0 ){
	    for(i=0;i<(-ty);i++) send_subl = neighsubl[send_subl][YDOWN];
	}
	else
	    for(i=0;i<ty;i++) send_subl = neighsubl[send_subl][YUP];
	if( tz < 0 ){
	    for(i=0;i<(-tz);i++) send_subl = neighsubl[send_subl][ZDOWN];
	}
	else
	    for(i=0;i<tz;i++) send_subl = neighsubl[send_subl][ZUP];
	if( tt < 0 ){
	    for(i=0;i<(-tt);i++) send_subl = neighsubl[send_subl][TDOWN];
	}
	else
	    for(i=0;i<tt;i++) send_subl = neighsubl[send_subl][TUP];
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) Make
	list of nodes from whom we expect messages */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode==this_node ){
		dest[i] = (char *)( &lattice[node_index(tx,ty,tz,tt)]) + field;
	    }
	    else{
		for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
		if(j < n_recv_msgs){
		    from_nodes[j].count++;
		}
		else {
	            if(n_recv_msgs==0){
			from_nodes = (struct msg_tmp *)malloc(
			    sizeof(struct msg_tmp) );
			from_nodes[0].node = othernode;
			from_nodes[0].count = 1;
			n_recv_msgs++;
	            }
		    else{
			from_nodes = (struct msg_tmp *)realloc( from_nodes,
			    (n_recv_msgs+1)*sizeof(struct msg_tmp) );
			from_nodes[j].node = othernode;
			from_nodes[j].count = 1;
			n_recv_msgs++;
		    }
		}
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,subl){
	    if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode==this_node ){
		dest[i] = (char *)( &lattice[node_index(tx,ty,tz,tt)]) + field;
	    }
	    else{
		for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
		if(j < n_recv_msgs){
		    from_nodes[j].count++;
		}
		else {
	            if(n_recv_msgs==0){
			from_nodes = (struct msg_tmp *)malloc(
			    sizeof(struct msg_tmp) );
			from_nodes[0].node = othernode;
			from_nodes[0].count = 1;
			n_recv_msgs++;
	            }
		    else{
			from_nodes = (struct msg_tmp *)realloc( from_nodes,
			    (n_recv_msgs+1)*sizeof(struct msg_tmp) );
			from_nodes[j].node = othernode;
			from_nodes[j].count = 1;
			n_recv_msgs++;
		    }
		}
	    }
	}
    }

    /* scan sites of sublattice we are sending, make list of nodes to which
	we must send messages and the number of messages to each. */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		if(j < n_send_msgs){
		    to_nodes[j].count++;
		}
		else {
		    if(n_send_msgs==0){
			to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
			to_nodes[0].node = othernode;
			to_nodes[0].count = 1;
			n_send_msgs++;
		    }
		    else{
			to_nodes = (struct msg_tmp *)realloc( to_nodes,
			    (n_send_msgs+1)*sizeof(struct msg_tmp) );
			to_nodes[j].node = othernode;
			to_nodes[j].count = 1;
			n_send_msgs++;
		    }
		}
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,send_subl){
	    if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		if(j < n_send_msgs){
		    to_nodes[j].count++;
		}
		else {
		    if(n_send_msgs==0){
			to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
			to_nodes[0].node = othernode;
			to_nodes[0].count = 1;
			n_send_msgs++;
		    }
		    else{
			to_nodes = (struct msg_tmp *)realloc( to_nodes,
			    (n_send_msgs+1)*sizeof(struct msg_tmp) );
			to_nodes[j].node = othernode;
			to_nodes[j].count = 1;
			n_send_msgs++;
		    }
		}
	    }
	}
    }

    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = (char *)malloc( nsites*tsize );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*tsize, MPI_BYTE,
	    from_nodes[i].node, GENERAL_GATHER_MSG_TYPE,
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
	mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=(char *)malloc( to_nodes[i].count*tsize );
	if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
    }

    /* reset to_node counters */
    for(i=0;i<n_send_msgs;i++)to_nodes[i].count=0;
    /* gather data into the buffers. Each entry in the buffers consists
	of the index of the site to which the data is sent, followed by
	the actual data */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    tx = (s->x - displacement[XUP] + nx)%nx;
	    ty = (s->y - displacement[YUP] + ny)%ny;
	    tz = (s->z - displacement[ZUP] + nz)%nz;
	    tt = (s->t - displacement[TUP] + nt)%nt;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		    to_nodes[j].count*tsize;
		*(int *)tpt = node_index(tx,ty,tz,tt);
		/* index of site on other node */
		memcpy( tpt+2*sizeof(int), F_PT(s,field), size);
		to_nodes[j].count++;
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,send_subl){
	    tx = (s->x - displacement[XUP] + nx)%nx;
	    ty = (s->y - displacement[YUP] + ny)%ny;
	    tz = (s->z - displacement[ZUP] + nz)%nz;
	    tt = (s->t - displacement[TUP] + nt)%nt;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		    to_nodes[j].count*tsize;
		*(int *)tpt = node_index(tx,ty,tz,tt);
		/* index of site on other node */
		memcpy( tpt+2*sizeof(int), F_PT(s,field), size);
		to_nodes[j].count++;
	    }
	}
    }

    /* start the sends */
    for(i=0;i<n_send_msgs;i++){
	nsites = to_nodes[i].count;
	MPI_Issend( mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, MPI_BYTE,
	    to_nodes[i].node, GENERAL_GATHER_MSG_TYPE,
	    MPI_COMM_WORLD, &mbuf[i+n_recv_msgs+1].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

#endif	/* N_SUBL32 */

/**********************************************************************/

#ifndef N_SUBL32

msg_tag * start_general_gather_from_temp(
/* arguments */
 void * field,	        /* which field? Pointer returned by malloc() */
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int disp_parity;	/* parity of displacement vector */
int send_parity;	/* parity of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
msg_tag *mbuf;		/* list of message tags, to be returned */
int n_send_msgs, n_recv_msgs;


    /* check for gather already in progress */
    if(g_gather_flag!=0){
	fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
	exit(1);
    }
    n_recv_msgs = n_send_msgs = 0;
    tsize = 2*sizeof(int)+size;
	/* Use 2*sizeof int so pointer will be aligned to double word */
    tdest = dest;
    /* find parity of sites that may be sent */
    if( (displacement[XUP]+displacement[YUP]+displacement[ZUP]+
	displacement[TUP])%2 == 0)disp_parity=EVEN;
    else disp_parity=ODD;
    switch(parity){
	case EVENANDODD: send_parity=EVENANDODD; break;
	case EVEN:
	    if( disp_parity==EVEN )send_parity=EVEN;
	    else send_parity=ODD;
	    break;
	case ODD:
	    if( disp_parity==EVEN )send_parity=ODD;
	    else send_parity=EVEN;
	    break;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) Make
	list of nodes from whom we expect messages */
    FORSOMEPARITY(i,s,parity){
	if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode==this_node ){
	    dest[i] = (char *)field + node_index(tx,ty,tz,tt)*size;
}
	else{
	    for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
	    if(j < n_recv_msgs){
		from_nodes[j].count++;
	    }
	    else {
	        if(n_recv_msgs==0){
		    from_nodes = (struct msg_tmp *)malloc(
			sizeof(struct msg_tmp) );
		    from_nodes[0].node = othernode;
		    from_nodes[0].count = 1;
		    n_recv_msgs++;
	        }
	        else{
		    from_nodes = (struct msg_tmp *)realloc( from_nodes,
			(n_recv_msgs+1)*sizeof(struct msg_tmp) );
		    from_nodes[j].node = othernode;
		    from_nodes[j].count = 1;
		    n_recv_msgs++;
	        }
	    }
	}
    }

    /* scan sites of parity we are sending, make list of nodes to which
	we must send messages and the number of messages to each. */
    FORSOMEPARITY(i,s,send_parity){
	if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    if(j < n_send_msgs){
		to_nodes[j].count++;
	    }
	    else {
	        if(n_send_msgs==0){
		    to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
		    to_nodes[0].node = othernode;
		    to_nodes[0].count = 1;
		    n_send_msgs++;
	        }
	        else{
		    to_nodes = (struct msg_tmp *)realloc( to_nodes,
			(n_send_msgs+1)*sizeof(struct msg_tmp) );
		    to_nodes[j].node = othernode;
		    to_nodes[j].count = 1;
		    n_send_msgs++;
	        }
	    }
	}
    }

    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = (char *)malloc( nsites*tsize );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*tsize, MPI_BYTE,
	    from_nodes[i].node, GENERAL_GATHER_MSG_TYPE, 
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=(char *)malloc( to_nodes[i].count*tsize );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
    }

    /* reset to_node counters */
    for(i=0;i<n_send_msgs;i++)to_nodes[i].count=0;
    /* gather data into the buffers. Each entry in the buffers consists
	of the index of the site to which the data is sent, followed by
	the actual data */
    FORSOMEPARITY(i,s,send_parity){
	tx = (s->x - displacement[XUP] + nx)%nx;
	ty = (s->y - displacement[YUP] + ny)%ny;
	tz = (s->z - displacement[ZUP] + nz)%nz;
	tt = (s->t - displacement[TUP] + nt)%nt;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		to_nodes[j].count*tsize;
	    *(int *)tpt = node_index(tx,ty,tz,tt);
	        /* index of site on other node */
	    memcpy( tpt+2*sizeof(int), (char *)field+i*size, size);
	    to_nodes[j].count++;
	}
    }

    /* start the sends */
    for(i=0;i<n_send_msgs;i++){
	nsites = to_nodes[i].count;
	MPI_Issend( mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, MPI_BYTE,
	    to_nodes[i].node, GENERAL_GATHER_MSG_TYPE, 
	    MPI_COMM_WORLD, &mbuf[i+n_recv_msgs+1].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

#else	/* N_SUBL32 */

msg_tag * start_general_gather_from_temp(
/* arguments */
 void * field,	        /* which field? Pointer returned by malloc()*/
 int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement,	/* displacement to gather from. four components */
 int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
 char ** dest)		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int send_subl;		/* sublattice of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
msg_tag *mbuf;		/* list of message tags, to be returned */
int n_send_msgs, n_recv_msgs;

    /* check for gather already in progress */
    if(g_gather_flag!=0){
	fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
	exit(1);
    }
    n_recv_msgs = n_send_msgs = 0;
    tsize = 2*sizeof(int)+size;
        /* Use 2*sizeof int so pointer will be aligned to double word */
    tdest = dest;

    /* find sublattice of sites that may be sent */
    /* This is not needed for EVENANDODD */
    if( subl != EVENANDODD ){
	/* Displacments by multiples of 4 in any direction does
	   not change sublattice */
	tx = displacement[XUP]%4;
	ty = displacement[YUP]%4;
	tz = displacement[ZUP]%4;
	tt = displacement[TUP]%4;
	send_subl = subl;
	if( tx < 0 ){
	    for(i=0;i<(-tx);i++) send_subl = neighsubl[send_subl][XDOWN];
	}
	else
	    for(i=0;i<tx;i++) send_subl = neighsubl[send_subl][XUP];
	if( ty < 0 ){
	    for(i=0;i<(-ty);i++) send_subl = neighsubl[send_subl][YDOWN];
	}
	else
	    for(i=0;i<ty;i++) send_subl = neighsubl[send_subl][YUP];
	if( tz < 0 ){
	    for(i=0;i<(-tz);i++) send_subl = neighsubl[send_subl][ZDOWN];
	}
	else
	    for(i=0;i<tz;i++) send_subl = neighsubl[send_subl][ZUP];
	if( tt < 0 ){
	    for(i=0;i<(-tt);i++) send_subl = neighsubl[send_subl][TDOWN];
	}
	else
	    for(i=0;i<tt;i++) send_subl = neighsubl[send_subl][TUP];
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) Make
	list of nodes from whom we expect messages */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode==this_node ){
		dest[i] = (char *)field + node_index(tx,ty,tz,tt)*size;
	    }
	    else{
		for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
		if(j < n_recv_msgs){
		    from_nodes[j].count++;
		}
		else {
	            if(n_recv_msgs==0){
			from_nodes = (struct msg_tmp *)malloc(
			    sizeof(struct msg_tmp) );
			from_nodes[0].node = othernode;
			from_nodes[0].count = 1;
			n_recv_msgs++;
	            }
		    else{
			from_nodes = (struct msg_tmp *)realloc( from_nodes,
			    (n_recv_msgs+1)*sizeof(struct msg_tmp) );
			from_nodes[j].node = othernode;
			from_nodes[j].count = 1;
			n_recv_msgs++;
		    }
		}
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,subl){
	    if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode==this_node ){
		dest[i] = (char *)field + node_index(tx,ty,tz,tt)*size;
	    }
	    else{
		for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
		if(j < n_recv_msgs){
		    from_nodes[j].count++;
		}
		else {
	            if(n_recv_msgs==0){
			from_nodes = (struct msg_tmp *)malloc(
			    sizeof(struct msg_tmp) );
			from_nodes[0].node = othernode;
			from_nodes[0].count = 1;
			n_recv_msgs++;
	            }
		    else{
			from_nodes = (struct msg_tmp *)realloc( from_nodes,
			    (n_recv_msgs+1)*sizeof(struct msg_tmp) );
			from_nodes[j].node = othernode;
			from_nodes[j].count = 1;
			n_recv_msgs++;
		    }
		}
	    }
	}
    }

    /* scan sites of sublattice we are sending, make list of nodes to which
	we must send messages and the number of messages to each. */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		if(j < n_send_msgs){
		    to_nodes[j].count++;
		}
		else {
		    if(n_send_msgs==0){
			to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
			to_nodes[0].node = othernode;
			to_nodes[0].count = 1;
			n_send_msgs++;
		    }
		    else{
			to_nodes = (struct msg_tmp *)realloc( to_nodes,
			    (n_send_msgs+1)*sizeof(struct msg_tmp) );
			to_nodes[j].node = othernode;
			to_nodes[j].count = 1;
			n_send_msgs++;
		    }
		}
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,send_subl){
	    if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	    if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	    if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	    if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		if(j < n_send_msgs){
		    to_nodes[j].count++;
		}
		else {
		    if(n_send_msgs==0){
			to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
			to_nodes[0].node = othernode;
			to_nodes[0].count = 1;
			n_send_msgs++;
		    }
		    else{
			to_nodes = (struct msg_tmp *)realloc( to_nodes,
			    (n_send_msgs+1)*sizeof(struct msg_tmp) );
			to_nodes[j].node = othernode;
			to_nodes[j].count = 1;
			n_send_msgs++;
		    }
		}
	    }
	}
    }

    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = (char *)malloc( nsites*tsize );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	MPI_Irecv( mbuf[i].msg_buf, nsites*tsize, MPI_BYTE,
	    from_nodes[i].node, GENERAL_GATHER_MSG_TYPE,
	    MPI_COMM_WORLD, &mbuf[i].msg_id );
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=MPI_REQUEST_NULL;
	mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=(char *)malloc( to_nodes[i].count*tsize );
	if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
    }

    /* reset to_node counters */
    for(i=0;i<n_send_msgs;i++)to_nodes[i].count=0;
    /* gather data into the buffers. Each entry in the buffers consists
	of the index of the site to which the data is sent, followed by
	the actual data */
    if( subl == EVENANDODD ){
	FORALLSITES(i,s){
	    tx = (s->x - displacement[XUP] + nx)%nx;
	    ty = (s->y - displacement[YUP] + ny)%ny;
	    tz = (s->z - displacement[ZUP] + nz)%nz;
	    tt = (s->t - displacement[TUP] + nt)%nt;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		    to_nodes[j].count*tsize;
		*(int *)tpt = node_index(tx,ty,tz,tt);
		/* index of site on other node */
		memcpy( tpt+2*sizeof(int), field*i*size, size);
		to_nodes[j].count++;
	    }
	}
    }
    else{
	FORSOMESUBLATTICE(i,s,send_subl){
	    tx = (s->x - displacement[XUP] + nx)%nx;
	    ty = (s->y - displacement[YUP] + ny)%ny;
	    tz = (s->z - displacement[ZUP] + nz)%nz;
	    tt = (s->t - displacement[TUP] + nt)%nt;
	    othernode = node_number(tx,ty,tz,tt);
	    if( othernode != this_node ){
		for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
		tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		    to_nodes[j].count*tsize;
		*(int *)tpt = node_index(tx,ty,tz,tt);
		/* index of site on other node */
		memcpy( tpt+2*sizeof(int), field+i*size, size);
		to_nodes[j].count++;
	    }
	}
    }

    /* start the sends */
    for(i=0;i<n_send_msgs;i++){
	nsites = to_nodes[i].count;
	MPI_Issend( mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, MPI_BYTE,
	    to_nodes[i].node, GENERAL_GATHER_MSG_TYPE,
	    MPI_COMM_WORLD, &mbuf[i+n_recv_msgs+1].msg_id );
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=MPI_REQUEST_NULL;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

#endif	/* N_SUBL32 */

/**********************************************************************/
void wait_general_gather(msg_tag *mbuf) {
register int i,j,k;
MPI_Status status;
    g_gather_flag=0;
    if(mbuf==NULL)return;
    for(i=0; mbuf[i].msg_buf != NULL; i++){
	MPI_Wait( &mbuf[i].msg_id, &status );
	/* set pointers in sites to correct location */
	for(j=0;j<from_nodes[i].count;j++){
	    /* k = index of site on this node, sent in message */
	    k = *(int *)( mbuf[i].msg_buf + j*tsize );
	    tdest[k] = mbuf[i].msg_buf + j*tsize + 2*sizeof(int);
	}
    }
    if( i > 0)free(from_nodes);
    return;
}

/**********************************************************************/
void cleanup_general_gather(msg_tag *mbuf) {
register int i,i0;
MPI_Status status;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /* wait for all send messages, free all send buffers */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	MPI_Wait( &mbuf[i].msg_id, &status );
	free( mbuf[i].msg_buf );
    }
    /* free the msg_tag buffer */
    free(mbuf);
}


/**********************************************************************/
/* FIELD POINTER ROUTINES */

/* Return a pointer to a field in the lattice at some coordinates.
   If the site is on this node, just return its address.  If it is
   on another node, make a buffer, get the data, and return address of
   buffer.
   Example:  
	su3_matrix *pt;
	pt = (su3_matrix *)field_pointer_at_coordinates(
	    F_OFFSET(xlink), sizeof(su3_matrix), x,y,z,t );
	... do stuff here ...
	cleanup_field_pointer( (char *)pt );
*/
char * field_pointer_at_coordinates(
/* arguments */
 int field,	/* offset of one of the fields in lattice[] */
 int size,	/* size of the field in bytes */
 int x,int y,int z,int t)	/* coordinates of point to get field from */
{
    /* if on my node, return address */
printf("Oops: tried to use field_pointer()\n"); terminate(1);
return NULL;  /* To satisfy compilers looking for a return value */
}

/**********************************************************************/
/* Return a pointer to a field in the lattice at some direction from
   a site on this node. (Usually the "current" site.)
   Works like field_pointer_at_coordinates. */
char * field_pointer_at_direction(
/* arguments */
 field_offset field,	/* offset of one of the fields in lattice[] */
 int size,	/* size of the field in bytes */
 site *s,	/* pointer to a site on this node */
 int direction)	/* direction of site's neighbor to get data from.
		   one of XUP,XDOWN,YUP... */
{
printf("Oops: tried to use field_pointer()\n"); terminate(1);
return NULL;  /* To satisfy compilers looking for a return value */
}

/**********************************************************************/
/* free any buffers allocated by above routines.  Argument is the
address returned by field_pointer...() */
void cleanup_field_pointer(char * buf) {
    /* if buf is an address in the lattice, leave it alone.  Otherwise
	it was created by malloc and should be freed */
    if( buf < (char *)lattice || buf >= (char *)(lattice+sites_on_node))
	free(buf);
}

/* Copy "size" bytes from "src" to "dest".  This should be assembler
   coded. */
/**
void memcpy( char *src, char *dest, int size) {
register int j;
    for(j=0;j<size;j++)dest[j]=src[j];
}
**/
/**
void memcpy( char *src, char *dest, int size) {
register int j;
    if( j%4==0){ ** transfer words **
	size >>=2;
        for(j=0;j<size;j++)((int *)dest)[j]=((int *)src)[j];
    }
    else{
        for(j=0;j<size;j++)dest[j]=src[j];
    }
}
**/


/**********************************************************************/
/* SEND AND RECEIVE FIELD */
/* send_field is to be called only by the node doing the sending */
/* get_field is to be called only by the node to which the field was sent */
void send_field(char *buf, int size, int tonode) {
    MPI_Send(buf,size,MPI_BYTE,tonode,FIELD_TYPE,MPI_COMM_WORLD);
}
/**********************************************************************/
void get_field(char *buf, int size, int fromnode){
    MPI_Status status;
    MPI_Recv(buf,size,MPI_BYTE,fromnode,FIELD_TYPE,
	MPI_COMM_WORLD,&status);
}

/**********************************************************************/
/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="MPI (portable)";
char * machine_type(){
    return(name);
}

/* Return my node number */
int mynode(){
    int node;
    MPI_Comm_rank( MPI_COMM_WORLD, &node );
    return(node);
}

/* Return number of nodes */
int numnodes(){
    int nodes;
    MPI_Comm_size( MPI_COMM_WORLD, &nodes );
    return(nodes);
}

/* Synchronize all nodes */
void g_sync(){
   MPI_Barrier( MPI_COMM_WORLD );
}

/* Sum Real over all nodes */
void g_floatsum( Real *fpt) {
Real work;
    MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
    *fpt = work;
}

/* Sum a vector of Reals over all nodes */
void g_vecfloatsum( Real *fpt, int nReals) {
register Real *work;
register int i;
    work = (Real *)malloc(nReals*sizeof(Real));
    MPI_Allreduce( fpt, work, nReals, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
    for(i=0;i<nReals;i++)fpt[i]=work[i];
    free(work);
}

/* Sum double over all nodes */
void g_doublesum( double *dpt) {
double work;
    MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    *dpt = work;
}
 
/* Sum a vector of doubles over all nodes */
void g_vecdoublesum( double *dpt, int ndoubles) {
register double *work;
register int i;
    work = (double *)malloc(ndoubles*sizeof(double));
    MPI_Allreduce( dpt, work, ndoubles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    for(i=0;i<ndoubles;i++)dpt[i]=work[i];
    free(work);
}

/* Sum complex over all nodes */
void g_complexsum( complex *cpt) {
complex work;
    MPI_Allreduce( cpt, &work, 2, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
    *cpt = work;
}

/* Sum double_complex over all nodes */
void g_dcomplexsum( double_complex *cpt) {
double_complex work;
    MPI_Allreduce( cpt, &work, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    *cpt = work;
}

/* Sum a vector of complex over all nodes */
void g_veccomplexsum( complex *cpt, int ncomplex){
register complex *work;
register int i;
    work = (complex *)malloc(ncomplex*sizeof(complex));
    MPI_Allreduce( cpt, work, 2*ncomplex, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
    for(i=0;i<ncomplex;i++)cpt[i]=work[i];
    free(work);
}

/* Sum a vector of double_complex over all nodes */
void g_vecdcomplexsum( double_complex *cpt, int ncomplex){
register double_complex *work;
register int i;
    work = (double_complex *)malloc(ncomplex*sizeof(double_complex));
    MPI_Allreduce( cpt, work, 2*ncomplex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    for(i=0;i<ncomplex;i++)cpt[i]=work[i];
    free(work);
}

/* Sum wilson_vector over all nodes */
void g_wvectosumReal( wilson_vector *wvpt) {
    g_veccomplexsum((complex *)wvpt, 12);
}

/* Global exclusive or acting on u_int32type */
void g_xor32( u_int32type *pt) {
  u_int32type work;
#ifdef SHORT_IS_32BIT
  assert(sizeof(unsigned short) == 4);
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED_SHORT, 
		 MPI_BXOR, MPI_COMM_WORLD );
#else
  assert(sizeof(unsigned int) == 4);
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED, 
		 MPI_BXOR, MPI_COMM_WORLD );
#endif
  *pt = work;
}

/* Find maximum of Real over all nodes */
void g_floatmax( Real *fpt) {
Real work;
    MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_MAX, MPI_COMM_WORLD );
    *fpt = work;
}

/* Find maximum of double over all nodes */
void g_doublemax( double *dpt) {
double work;
    MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    *dpt = work;
}

/* Broadcast Realing point number from node zero */
void broadcast_float(Real *fpt) {
    MPI_Bcast( fpt, 1, MILC_MPI_REAL, 0, MPI_COMM_WORLD );
}

/* Broadcast double precision Realing point number from node zero */
void broadcast_double(double *dpt) {
    MPI_Bcast( dpt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(complex *cpt) {
    MPI_Bcast( cpt, 2, MILC_MPI_REAL, 0, MPI_COMM_WORLD );
}

/* Broadcast double precision complex number from node zero */
void broadcast_dcomplex(double_complex *cpt) {
    MPI_Bcast( cpt, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}

/* Broadcast bytes from node 0 to all others */
void broadcast_bytes(char *buf, int size) {
    MPI_Bcast( buf, size, MPI_BYTE, 0, MPI_COMM_WORLD );
}

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(int tonode, int *address) {
    MPI_Send( address, 1, MPI_INT, tonode, SEND_INTEGER_TYPE,
	MPI_COMM_WORLD);
}

/* Receive an integer from another node */
/* Note we do not check if this was really meant for us */
void receive_integer(int fromnode, int *address) {
    MPI_Status status;
    MPI_Recv( address, 1, MPI_INT, fromnode, SEND_INTEGER_TYPE,
	MPI_COMM_WORLD, &status);
}

/* Double precision time */
double dclock(){ 
   /*return( MPI_Wtime() ); */
    long fine;
    fine = clock();
    return( ((double)fine)/1000000.0 );
}

/* Print time stamp */
void time_stamp(char *msg){
  time_t time_stamp;
  
  if(mynode()==0){
    time(&time_stamp);
    printf("%s: %s\n",msg,ctime(&time_stamp));
  }
}

/* version of exit for multinode processes -- kill all nodes */
void terminate(int status) {
  time_stamp("termination");
  printf("Termination: node %d, status = %d\n",this_node,status);
  MPI_Abort( MPI_COMM_WORLD, 0);
  exit(status);
}

/* version of normal exit for multinode processes */
void normal_exit(int status) {
  time_stamp("exit");
  MPI_Finalize();
  exit(status);
}
