/******************  com_t3e.c **************************************/
/* Communications routines for the SU3 program
   MIMD version 6.
   This file is communications-scheme dependent.
   T3E PVM shared-memory implementation */

/* Modifications

   4/11/00 Added g_xor32
   8/08/97 ANSI prototyping for all routines C.D.
   8/08/97 Corrected error in g_doublemax (broadcast_float-> double)
   SG 9/30/96 pragmas commented out for t3e
   9/23/96 Explicit void types for modules with empty returns C.D.
   8/30/96 Improved sort_site_list C.D.
   8/05/96 Added broadcast_bytes for system-dependent
      parallel file system calls C.D.
   DT 11/22/93 version for T3D
   DT 7/29/93: First T3D implementation.  Use pvm_get_PE as node
      number, use this number in all PVM communication calls?
   DT 4/02/93: PVM version 3.0
   DT 3/16/93: PVM version 2.4, node 0 starts other node programs 
   2/18/93 Port to T3D, PVM version 2.3 C.D.

 */

/*
   This version implements asynchronous sends and receives.
   No interrupts that I know of, so field_pointer is a problem.

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
   wait_gather()  waits for receives to finish, insuring that the
   data has actually arrived.
   cleanup_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   restart_gather() repeats the internode communications of a previous
   gather.

   start_general_gather() starts asynchronous sends and receives required
   to gather fields at arbitrary displacement.
   wait_general_gather()  waits for receives to finish, insuring that the
   data has actually arrived, and sets pointers to received data.
   cleanup_general_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.

   send_field() sends a field to one other node.
   get_links() receives a field from some other node.
*/

#include "generic_includes.h"
#include <pvm3.h>

/* For shared memory calls.  this is currently /mpp/local/include/shmem.h*/
/*#include <mpp/shmem.h>*/
#include "shmem_t3e.h"

#define TRUE 1
#define FALSE 0

/* Global variables for the communications stuff */
    /* message types for gather encode the direction, the node sending
	the message, and the field offset of the source.  The base is
	added to this so that gather message types are higher than any
	others we use. The direction is in the lowest order bits, then
	the node number, then the field offset. */
    /* for computing message type in gather */
/*#pragma _CRI taskcommon mt_nodeshift,mt_offshift */
    int mt_nodeshift;	/* number of bits to shift node number */
    int mt_offshift;	/* number of bits to shift field offset */
    /* macro to compute the message type */
#define GATHER_MSG_TYPE(offset,node,dir) \
   (GATHER_BASE + ((offset)<<mt_offshift) + ((node)<<mt_nodeshift) + (dir))
#define GENERAL_GATHER_MSG_TYPE(node) (GENERAL_GATHER_BASE + node)

    /* Maintain lists of headers to lists of comlinks */
    /* One list for sites to be received, one for sites to be sent.
	These will often point to the same things */
/* #pragma _CRI taskcommon neighborlist,neighborlist_send */
    comlink ** neighborlist, **neighborlist_send;
    /* addresses of neighboring sites, NULL if off-node */
    /* neighbor[X][i] is a pointer to the neighbor of site lattice[i] in
       gather number X */
/* #pragma _CRI taskcommon neighbor */
    site *** neighbor;
    /* Number of gathers (mappings) that have been set up */
/* #pragma _CRI taskcommon n_gathers */
    int n_gathers;

/* #pragma _CRI taskcommon pvm_numnodes */
    int pvm_numnodes;

    /* For shared memory, need workspace and synchronization buffers */
    /* workspace buffers must be maximum of _SHMEM_REDUCE_MIN_WRKDATA_SIZE
       and nreduce/2+1 times the size of the object (Right now, our
	largest object is double.  So make it large enough to do a
	double wilson_vector. Syncronization buffers must be
	maximum of _SHMEM_BARRIER_SYNC_SIZE and _SHMEM_REDUCE_SYNC_SIZE
	 */
#define SHMEM_SYNCSIZE ( _SHMEM_BARRIER_SYNC_SIZE > \
_SHMEM_REDUCE_SYNC_SIZE ? _SHMEM_BARRIER_SYNC_SIZE : _SHMEM_REDUCE_SYNC_SIZE )
#define SHMEM_WORKSIZE ( _SHMEM_REDUCE_MIN_WRKDATA_SIZE > \
(48/2+4) ? _SHMEM_REDUCE_MIN_WRKDATA_SIZE : (48/2+4) )
/* #pragma _CRI cache_align  s_tmp_buf0,s_tmp_buf1,d_tmp_buf0,d_tmp_buf1 */
    static long s_tmp_buf0[SHMEM_SYNCSIZE];
    static long s_tmp_buf1[SHMEM_SYNCSIZE];
    static long * s_sh_buf[2];
    static double d_tmp_buf0[ SHMEM_WORKSIZE ];
    static double d_tmp_buf1[ SHMEM_WORKSIZE ];
    static double * d_sh_buf[2];

struct msg_tmp { int node,count; }; /* temporary structure for keeping
	track of messages to be sent or received */
/* #pragma _CRI taskcommon g_gather_flag,to_nodes,from_nodes,tsize,tdest */
int g_gather_flag;	/* flag to tell if general gather in progress */
struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
int tsize;		/* size of entry in messages = sizeof(int)+size */
char ** tdest;		/* tdest is copy of dest */
/* from_nodes, tsize and tdest are global because they are set in 
   start_general_gather() and used in wait_general_gather().  This
   works because we allow only one general_gather in progress at a
   time. */

#include <signal.h>
#ifdef PSC
/* PSC feature: might get signal 27 before a rollout */
void rollout_handler( int sig );
#endif


/**********************************************************************/
/* Machine initialization */
void initialize_machine(int argc, char **argv){
    int i,code;
    char * envpt;

    g_gather_flag=0; /* For emulator, put back outside for T3D */

    /* Figure out how many nodes we have */
	pvm_numnodes = pvm_gsize(0);

	/* set up shared memory stuff */
	s_sh_buf[0] = s_tmp_buf0;
	s_sh_buf[1] = s_tmp_buf1;
	for(i=0;i<SHMEM_SYNCSIZE;i++){
	    s_sh_buf[0][i] = _SHMEM_SYNC_VALUE;
	    s_sh_buf[1][i] = _SHMEM_SYNC_VALUE;
	}
	shmem_set_cache_inv();
	d_sh_buf[0] = d_tmp_buf0;
	d_sh_buf[1] = d_tmp_buf1;

#ifdef PSC
	/* PSC feature:  you might get signal 27 before a rollout */
	signal(27,rollout_handler);
#endif
}

#ifdef PSC
/* PSC feature:  you might get signal 27 before a rollout */
void rollout_handler( int sig ){
    if(this_node==0)printf("ROLLOUT COMING!! quitting\n");
    fflush(stdout);
    exit(0);
}
#endif


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

   /* Set up variables for constructing message types */
   for( mt_nodeshift=0,i=MAX_GATHERS-1; i>0; i>>=1,mt_nodeshift++);
   for( mt_offshift=mt_nodeshift,i=numnodes()-1; i>0; i>>=1,mt_offshift++);
   /* Assume that fields will be ints or bigger, offset divisible by 4 */
    mt_offshift -= 2;
   /* Check for possible overflow of message type.  I think 2^30-1 is
	the largest allowed on ncube */
    /* When this happens, fix the program */
    if( GATHER_MSG_TYPE(sizeof(site),numnodes(),MAX_GATHERS-1) >= 0x3fffffff ){
	if(this_node==0)printf(
	    "Possible overflow of gather message type - Fix the program\n");
	exit(1);
    }

    /* initialize neighborlist[] */
    neighborlist = (comlink **)malloc(NDIRS*sizeof(comlink *));
    neighborlist_send = (comlink **)malloc(NDIRS*sizeof(comlink *));
    /* Allocate space for lists of pointers to neighbor sites.
       (NULL if neighbor not on this node) */
    neighbor = (site ***)malloc(NDIRS*sizeof(site **));
    n_gathers=0;

    for(i=XUP;i<=TUP;i++)
	make_gather(neighbor_coords_special,&i,WANT_INVERSE,
	    ALLOW_EVEN_ODD,SWITCH_PARITY);

    /* Sort into the order we want for nearest neighbor gathers,
	so you can use XUP, XDOWN, etc. as argument in calling them. */
    sort_eight_special((void **) neighbor );
    sort_eight_special((void **) neighborlist );
    sort_eight_special((void **) neighborlist_send );
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
/* copy linked list of comlinks, switching even and odd */
comlink * copy_list_switch( comlink *old_compt ); 

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
    neighbor = (site ***)realloc(neighbor, n_gathers*sizeof(site **));
    if( inverse==WANT_INVERSE) {
	neighborlist[n_gathers-2] = neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-2] = neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-2] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-2]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-2;	/* index of gather we are working on */
    }
    else {
	neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-1;
    }

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
	    ) ){
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

    /* RECEIVE LISTS: */
    /* Fill in pointers to sites which are on this node, NULL if
	they are off-node */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
        j = node_number(x,y,z,t);	/* node for neighbor site */
        /* if neighbor is on node, set up pointer */
        if( j == mynode() ) neighbor[dir][i]= &(lattice[node_index(x,y,z,t)]);
        else		    neighbor[dir][i]= NULL;
    }

    /* make lists of sites which get data from other nodes.  */
    neighborlist[dir] = make_send_receive_list( func, args, want_even_odd,
        FORWARDS, RECEIVE );

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
        if( j == mynode() ) neighbor[dir][i]= &(lattice[node_index(x,y,z,t)]);
        else 		    neighbor[dir][i]= NULL;
    }

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

    return(dir-1);
}


/**********************************************************************/
comlink *  make_send_receive_list(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
        		/* function which defines sites to gather from */
 int *args,		/* list of arguments, to be passed to function */
 int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
 int forw_back,		/* FORWARDS or BACKWARDS */
 int send_recv)		/* SEND or RECEIVE list */
{
register int i,j,k;	/* scratch */
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

/**********************************************************************/
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
/* #pragma _CRI taskcommon mreqbuf */
static msg_request mreqbuf;	/* global so handler can use it too */
void start_handlers(){
void fillfieldrequest(int type,int size,int node,int pid);
    /**
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    **/
if(this_node==0)printf("Can't start interrupt handler yet!!\n");
}
/**********************************************************************/
/* The handler for field requests.  It fills the requests, then posts a
   receive for the next request. */
void fillfieldrequest(int type,int size,int node,int pid){
register char *buf;
register int trapstate;
void fillfieldrequest();

if(type != FIELD_REQUEST){
printf("BOTCH node %d from node %d type %d size %d pid %d\n",
mynode(),node,type,size,pid);
printf("BOTCH: request: size = %d, index = %d, offset = %d\n",
mreqbuf.size,mreqbuf.index,mreqbuf.field);
}
    /* Use "masktraps() for simulator, "masktrap()" for real machine */
    /**
    trapstate=masktrap(1);
    csend( FIELD_REPLY, F_PT( &(lattice[mreqbuf.index]), mreqbuf.field ),
	mreqbuf.size, node, pid);
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    masktrap(trapstate);
    **/
    printf("Oops: called fillfieldrequest\n");
}


/**********************************************************************/
/* GATHER ROUTINES */
/* start_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_gather() and cleanup_gather() calls.
   This list contains msg_tags for all receive buffers, followed by
   a msg_tag with msg_id = 0 and msg_buf = NULL, followed by msg_tags
   for all send buffers, followed by a msg_tag with id=0, buf=NULL.
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
	   vector (or a copy therof) on the neighbor of site i in the
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
	( f_offset(src) | sending_node | dir ) + GATHER_BASE
*/

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
char *tbuf;		/* temporary buffer pointer */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

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

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
	    }}
	    break;
	case ODD:
	    FORODDSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
	    }}
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
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
	mbuf[i].msg_OK=FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(field,mynode(),index);
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);
                /**pvm_pkbyte( ((char *)(&lattice[compt->ositelist[j]])+field),
                    size, 1 );**/
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);
                /**pvm_pkbyte( ((char *)(&lattice[compt->esitelist[j]])+field),
                    size, 1 );**/
	    }
	}
	/* start the send */
	pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );
	pvm_send( mbuf[i].msg_node, GATHER_MSG_TYPE(field,mynode(),index) );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
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
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(field, mbuf[i].msg_node, index);
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=0;
        mbuf[n_recv_msgs].msg_buf=NULL;
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
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
char *tbuf;		/* temporary buffer pointer */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	mbuf[i].msg_id = GATHER_MSG_TYPE(field,mynode(),index);
	mbuf[i].msg_OK = FALSE;
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);
                /**pvm_pkbyte( ((char *)(&lattice[compt->ositelist[j]])+field),
                    size, 1 );**/
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);
                /**pvm_pkbyte( ((char *)(&lattice[compt->esitelist[j]])+field),
                    size, 1 );**/
	    }
	}
	/* start the send */
	pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );
	pvm_send( compt->othernode, GATHER_MSG_TYPE(field,mynode(),index));
    }
    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
      /* use same buffer to receive neighbors */
      mbuf[i].msg_OK = FALSE;
    }
}

/**********************************************************************/
void wait_gather(msg_tag *mbuf) {
  /* For asynchronous receives with random order of arrival */
  register int i,bufid;
  int all_done;

  if(mbuf==NULL)return;

  /* for each node which has neighbors of my sites */
  /* wait for all receive messages */
  all_done = FALSE;

  while(!all_done) {
      all_done = TRUE;
      for(i=0; mbuf[i].msg_buf !=NULL; i++) {
	  /* If we have processed this message already, go on. */
	  if(!mbuf[i].msg_OK) {
	      /* If it has arrived, grab it.  Otherwise, go on. */
	      if( pvm_nrecv( ANY_NODE, mbuf[i].msg_id ) != 0){
	      /**if( pvm_recv( ANY_NODE, mbuf[i].msg_id ) != 0){**/
		  pvm_upkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );
		  /* Check off this receive */
		  mbuf[i].msg_OK = TRUE;
	      }
	      else {
		/* This message not received yet, so we are not done */
		  all_done = FALSE;
	      }
	  }
      }
  }
  return;
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
   a msg_tag with msg_id = 0 and msg_buf = NULL, followed by msg_tags
   for all send buffers, followed by a msg_tag with id=0, buf=NULL.
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
register int i,j,k;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int disp_parity;	/* parity of displacement vector */
int send_parity;	/* parity of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
msg_tag *mbuf;		/* list of message tags, to be returned */
char *tbuf;		/* temporary buffer pointer */
int n_send_msgs, n_recv_msgs;
int type;		/* message type - nread needs its address */

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
	    dest[i] = F_PT( &lattice[node_index(tx,ty,tz,tt)], field );
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

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=(char *)malloc( to_nodes[i].count*tsize );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
	mbuf[i+n_recv_msgs+1].msg_OK=FALSE;
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
	mbuf[i+n_recv_msgs+1].msg_id = GENERAL_GATHER_MSG_TYPE(this_node);
	pvm_initsend( PvmDataRaw );
	pvm_pkbyte( (char *)mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, 1 );
	pvm_send( to_nodes[i].node, GENERAL_GATHER_MSG_TYPE(this_node));
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = (char *)malloc( nsites*tsize );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GENERAL_GATHER_MSG_TYPE(from_nodes[i].node);
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=0;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }
    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

/**********************************************************************/
void wait_general_gather(msg_tag *mbuf) {
  register int i,j,k;
  int all_done;

  g_gather_flag=0;
  if(mbuf==NULL)return;

  /* for each node which has neighbors of my sites */
  /* wait for all receive messages */
  all_done = FALSE;

  while(!all_done) {
      all_done = TRUE;
      for(i=0; mbuf[i].msg_buf !=NULL; i++) {
	  /* If we have processed this message already, go on. */
	  if(!mbuf[i].msg_OK) {
	      /* If it has arrived, grab it.  Otherwise, go on. */
	      if( pvm_nrecv( ANY_NODE, mbuf[i].msg_id ) != 0){
	      /**if( pvm_recv( ANY_NODE, mbuf[i].msg_id ) != 0){**/
		  pvm_upkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );
		  /* Check off this receive */
		  mbuf[i].msg_OK = TRUE;
		  
		  /* set pointers in sites to correct location */
		  for(j=0;j<from_nodes[i].count;j++){
		    /* k = index of site on this node, sent in message */
		    k = *(int *)( mbuf[i].msg_buf + j*tsize );
		    tdest[k] = mbuf[i].msg_buf + j*tsize + 2*sizeof(int);
		  }
		}
	      else {
		/* This message not received yet, so we are not done */
		  all_done = FALSE;
	      }
	  }
      }
  }
  if( i > 0)free(from_nodes);
  return;
}

/**********************************************************************/
void cleanup_general_gather(msg_tag *mbuf) {
register int i,i0;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /* wait for all send messages, free all send buffers */
    /* In simulator, messages were synchronous so don't wait */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
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
register int node,index;
register char *buf;
msg_request mreq;
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
register int node,index,trapstate,id;
int x,y,z,t;
register char *buf;
msg_request mreq;
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

/**********************************************************************/
/* SEND AND RECEIVE FIELD */
void send_field(char *buf, int size, int tonode) {

  pvm_initsend( PvmDataRaw );
  pvm_pkbyte( (char *)buf,size,1);
  pvm_send(tonode,FIELD_TYPE);
}
/**********************************************************************/
void get_field(char *buf, int size){

  pvm_recv( ANY_NODE, FIELD_TYPE);
  pvm_upkbyte( (char *)buf, size, 1 );
}

/**********************************************************************/
/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="T3E (PVM-shmem)";
char * machine_type(){
    return(name);
}

/* Return my node number */
int mynode()
{
  return( pvm_get_PE( pvm_mytid() ) );
}

/* Return number of nodes */
int numnodes(){
  return( pvm_numnodes );
}

/* Synchronize all nodes */
void g_sync(){
  barrier();
}

/* #pragma _CRI taskcommon bufswitch */
int bufswitch=0;

/* Sum Real over all nodes */
/* Use shared memory T3d library */
void g_floatsum( Real *fpt) {
Real t;
    shmem_Real_sum_to_all( fpt,fpt,1,0,0,pvm_numnodes,
	(Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    bufswitch=1-bufswitch;
}

void old_floatsum( Real *fpt) {
  Real work,sum;
  int othernode,type;
  
/*THIS IS INDETERMINATE AT THE MOMENT BECAUSE THINGS ARE ADDED
IN UNKNOWN ORDER */
    type=SUM_MILC_REAL_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkReal( fpt,1,1);
      pvm_send( 0, type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *fpt; /* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkReal( &work, 1, 1 );
	  sum += work;
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_float(&sum);

  /* Answer goes in fpt */
  *fpt = sum;
}

/* Sum double over all nodes */
/* Use shared memory T3d library */
void g_doublesum( double *dpt) {
    shmem_double_sum_to_all( dpt,dpt,1,0,0,pvm_numnodes,
	d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    bufswitch=1-bufswitch;
}

/* Sum double over all nodes */
void old_doublesum( double *dpt) {
  double work,sum;
  int othernode,type;
double t1,t2;
t1 = *dpt;

  type=SUM_DOUBLE_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkdouble( dpt,1,1 );
      pvm_send( 0, type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *dpt; /* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkdouble( &work, 1, 1);
	  sum += work;
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_double(&sum);

  /* Answer goes in dpt */
  *dpt = sum;
}

/* Sum a vector of doubles over all nodes */
/* Use shared memory T3d library */
void g_vecdoublesum( double *dpt, int ndoubles) {
    register int i,j; /* 24 at a time, from size of SHMEM_WORKSIZE */
	/* Actually, batches of more than 8 don't seem to work 11/13/95 */
    for(i=0;i<ndoubles;i+=8){
	j= ndoubles-i < 8 ? ndoubles-i : 8 ;
        shmem_double_sum_to_all( ((double *)dpt)+i,((double *)dpt)+i,
	    j,0,0,pvm_numnodes,
	    (double *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
        bufswitch=1-bufswitch;
    }
    /*register int i;
    for(i=0;i<ndoubles;i++)g_doublesum( dpt+i );*/
}

/* Sum complex over all nodes */
/* Use shared memory T3d library */
void g_complexsum( complex *cpt) {
    shmem_Real_sum_to_all( (Real *)cpt,(Real *)cpt,2,0,0,pvm_numnodes,
	(Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    bufswitch=1-bufswitch;
}

/* Sum complex over all nodes */
void old_complexsum( complex *cpt) {
  complex work,sum;
  int othernode,type;

    type=SUM_COMPLEX_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkcplxHELP( (Real *)cpt, 1, 1 );
      pvm_send( 0, type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *cpt; /* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkcplxHELP( (Real *)&work, 1, 1 );
	  CSUM(sum,work);
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_complex(&sum);

  /* Answer goes in cpt */
  *cpt = sum;
}

/* Use shared memory T3d library */
void g_dcomplexsum( double_complex *cpt) {
    shmem_double_sum_to_all( (double *)cpt,(double *)cpt,
	2,0,0,pvm_numnodes, d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    bufswitch=1-bufswitch;
}

/* Sum double complex over all nodes */
void old_dcomplexsum( double_complex *cpt) {
  double_complex work,sum;
  int othernode,type;

    type=SUM_DCOMPLEX_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkdcplx( (double *)cpt, 1, 1 );
      pvm_send( 0, type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *cpt; /* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkdcplx( (double *)&work, 1, 1 );
	  CSUM(sum,work);
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_dcomplex(&sum);

  /* Answer goes in cpt */
  *cpt = sum;
}

/* Sum a vector of complex over all nodes */
void g_veccomplexsum( complex *cpt, int ncomplex) {
    register int i,j; /* 24 at a time, from size of SHMEM_WORKSIZE */
	/* Actually, batches of more than 8 don't seem to work 11/13/95 */
    for(i=0;i<2*ncomplex;i+=8){
	j= 2*ncomplex-i < 8 ? 2*ncomplex-i : 8 ;
        shmem_Real_sum_to_all( ((Real *)cpt)+i,((Real *)cpt)+i,
	    j,0,0,pvm_numnodes,
	    (Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
        bufswitch=1-bufswitch;
    }
}

/* Sum wilson_vector over all nodes */
void g_wvectosumReal( wilson_vector *wvpt) {
    /* vector lengths > 8 don't seem to work 11/13/95 */
    shmem_Real_sum_to_all( ((Real *)wvpt)+0, ((Real *)wvpt)+0
	,8,0,0,pvm_numnodes, (Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    shmem_Real_sum_to_all( ((Real *)wvpt)+8, ((Real *)wvpt)+8
	,8,0,0,pvm_numnodes, (Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    shmem_Real_sum_to_all( ((Real *)wvpt)+16, ((Real *)wvpt)+16
	,8,0,0,pvm_numnodes, (Real *)d_sh_buf[bufswitch],s_sh_buf[bufswitch]);
    bufswitch=1-bufswitch;
}

/* Global exclusive or acting on int32type */
void g_xor32( u_int32type *pt) {
#ifdef SHORT32
  shmem_short_xor_to_all((short *)pt, (short *)pt, 1, 0, 0, pvm_numnodes,
		 (short *)d_sh_buf[bufswitch],(long *)s_sh_buf[bufswitch]);
#else
  shmem_int_xor_to_all((int *)pt, (int *)pt, 1, 0, 0, pvm_numnodes,
		 (int *)d_sh_buf[bufswitch],(long *)s_sh_buf[bufswitch]);
#endif
  bufswitch=1-bufswitch;
}

/* Find maximum of Real over all nodes */
void g_floatmax( Real *fpt) {
  Real work,high;
  int othernode,type;
  
    type=MAX_MILC_REAL_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkReal( fpt,1,1 );
      pvm_send( 0, type);
    }
    else{
	/* node zero finds the max value */
	high = *fpt; /*  from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkReal( &work, 1, 1);
	  if(work>high)high = work;
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_float(&high);

  /* Answer goes in fpt */
  *fpt = high;
}


/* Find maximum of double over all nodes */
void g_doublemax( double *dpt) {
  double work,high;
  int othernode,type;
  
    type=MAX_DOUBLE_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkdouble( dpt,1,1 );
      pvm_send( 0, type);
    }
    else{
	/* node zero finds the max value */
	high = *dpt; /*  from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkdouble( &work, 1, 1 );
	  if(work>high)high = work;
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_double(&high);

  /* Answer goes in dpt */
  *dpt = high;
}


/* Broadcast Realing point number from node zero */
void broadcast_float(Real *fpt) {
    int i;
  
    if(mynode()==0) {
        pvm_initsend( PvmDataRaw );
        pvm_pkReal( fpt, 1,1 );
        /**for(i = 0; i < numnodes(); i++) if(i != mynode())
	  pvm_send( i, BROADCAST_MILC_REAL_TYPE);**/
        pvm_bcast( NULL, BROADCAST_MILC_REAL_TYPE);
    }
    else {
        pvm_recv( ANY_NODE, BROADCAST_MILC_REAL_TYPE);
        pvm_upkReal( fpt,1,1 );
    }
}

/* Broadcast double precision Realing point number from node zero */
void broadcast_double(double *dpt) {
    int i;

    if(mynode()==0) {
        pvm_initsend( PvmDataRaw );
        pvm_pkdouble( dpt, 1, 1 );
        /**for(i = 0; i < numnodes(); i++) if(i != mynode())
	  pvm_send( i, BROADCAST_DOUBLE_TYPE);**/
	pvm_bcast( NULL, BROADCAST_DOUBLE_TYPE);
    }
    else {
        pvm_recv( ANY_NODE, BROADCAST_DOUBLE_TYPE);
        pvm_upkdouble( dpt,1,1 );
    }
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(complex *cpt) {
    int i;

    if(mynode()==0) {
        pvm_initsend( PvmDataRaw );
        pvm_pkcplxHELP( (Real *)cpt, 1, 1);
        /**for(i = 0; i < numnodes(); i++) if(i != mynode())
	  pvm_send( i, BROADCAST_COMPLEX_TYPE);**/
	pvm_bcast( NULL, BROADCAST_COMPLEX_TYPE);
    }
    else {
        pvm_recv( ANY_NODE, BROADCAST_COMPLEX_TYPE);
        pvm_upkcplxHELP( (Real *)cpt,1,1 );
    }
}

/* Broadcast double precision complex number from node zero */
void broadcast_dcomplex(double_complex *cpt) {
    int i;

    if(mynode()==0) {
        pvm_initsend( PvmDataRaw );
        pvm_pkdcplx( (double *)cpt, 1, 1);
        /**for(i = 0; i < numnodes(); i++) if(i != mynode())
	  pvm_send( i, BROADCAST_DCOMPLEX_TYPE);**/
	pvm_bcast( NULL, BROADCAST_DCOMPLEX_TYPE);
    }
    else {
        pvm_recv( ANY_NODE, BROADCAST_DCOMPLEX_TYPE);
        pvm_upkdcplx( (double *)cpt,1,1 );
    }
}

/* Broadcast bytes */
void broadcast_bytes(char *buf, int size) {
     int i;

    if(mynode()==0) {
        pvm_initsend( PvmDataRaw );
        pvm_pkbyte( buf, size, 1);
        /**for(i = 0; i < numnodes(); i++) if(i != mynode())
	  pvm_send( i, BROADCAST_BYTES_TYPE);**/
	pvm_bcast( NULL, BROADCAST_BYTES_TYPE);
    }
    else {
        pvm_recv( ANY_NODE, BROADCAST_BYTES_TYPE);
        pvm_upkbyte( buf, size ,1 );
    }
}
 

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(int tonode, int *address) {

  pvm_initsend( PvmDataRaw );
  pvm_pkint( (int *)address, 1, 1 );
  pvm_send( tonode, SEND_INTEGER_TYPE);

}

/* Receive an integer from another node */
/* Note we do not check if this was really meant for us */
void receive_integer(int *address) {
  
  pvm_recv( ANY_NODE, SEND_INTEGER_TYPE);
  pvm_upkint( (int *)address,1,1 );
}

/* Double precision time */
/* This one wraps around after 36 minutes!! It gives the cpu time,
   not the wall clock time */
double dclock(){
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
  int i;

  time_stamp("termination");
  printf("Termination: node %d, status = %d\n",this_node,status);
  for( i = 0; i < numnodes(); i++){
    if(i != mynode())pvm_kill(i);
  }
  pvm_exit();
  globalexit(status);  /* Suggested by C.McNeile */
  exit(status);
}
