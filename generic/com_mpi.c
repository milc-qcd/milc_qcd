/******************  com_mpi.c *****************************************/
/* Communications routines for the SU3 program
   MIMD version 6.
   This file is communications-scheme dependent.
   MPI version - allegedly machine independent
*/
/* Modifications

    4/20/02 added start_general_gather_from_temp C.D.
   10/15/01 condensed and modified to use multiple gathers - JCO
    1/30/00 combined with Schroedinger functional and
            32 sublattice versions - UMH
   11/27/98 Corrected g_wvectorsumfloat and made independent of su3.h. C.D.
    9/02/97 Revised to allow gathers from temporary fields.  neighbor[]
	    is now list of indices, add start/restart_gather_field D.T.
    8/05/97 ANSI prototyping for all routines C.D.
   10/05/96 Moved parallel I/O wrappers to io_ansi.c C.D.
    9/23/96 Explicit void types for modules with empty returns C.D.
    9/20/96 Added restore/save_checkpoint C.D.  
    9/20/96 Improved sort_site_list C.D.  
    9/20/96 Added broadcast_bytes and wrappers for system-dependent
            parallel file system calls C.D.   
*/
/*
  Exported Functions:

   initialize_machine()  does any machine dependent setup at the
                           very beginning.
   normal_exit()         closes communications and exits
   terminate()           halts program abruptly and exits
   machine_type()        returns string describing communications architecture
   mynode()              returns node number of this node.
   numnodes()            returns number of nodes
   g_sync()              provides a synchronization point for all nodes.
   g_floatsum()          sums a floating point number over all nodes.
   g_vecfloatsum()       sums a vector of generic floats over all nodes 
   g_doublesum()         sums a double over all nodes.
   g_vecdoublesum()      sums a vector of doubles over all nodes.
   g_complexsum()        sums a generic precision complex number over all nodes.
   g_veccomplexsum()     sums a vector of generic precision complex numbers
                           over all nodes.
   g_dcomplexsum()       sums a double precision complex number over all nodes.
   g_vecdcomplexsum()    sums a vector of double_complex over all nodes 
   g_wvectorsumfloat()   sums a generic precision wilson vector over all nodes.
   g_xor32()             finds global exclusive or of 32-bit word
   g_floatmax()          finds maximum floating point number over all nodes.
   g_doublemax()         finds maximum double over all nodes.
   broadcast_float()     broadcasts a generic precision number from
	                   node 0 to all nodes.
   broadcast_double()    broadcasts a double precision number
   broadcast_complex()   broadcasts a generic precision complex number
   broadcast_dcomplex()  broadcasts a double precision complex number
   broadcast_bytes()     broadcasts a number of bytes
   send_integer()        sends an integer to one other node
   receive_integer()     receives an integer
   send_field()          sends a field to one other node.
   get_field()           receives a field from some other node.
   dclock()              returns a double precision time, with arbitrary zero
   time_stamp()          print wall clock time with message
   sort_eight_gathers()  sorts eight contiguous gathers from order
                           XUP,XDOWN,YUP,YDOWN... to XUP,YUP,...XDOWN,YDOWN...
   make_nn_gathers()     makes all necessary lists for communications with
                           nodes containing neighbor sites.
   make_gather()         calculates and stores necessary communications lists
                           for a given gather mapping
   declare_gather()      creates a message tag that defines specific details
                           of a gather to be used later
   declare_gather_from_temp()  creates a message tag that defines specific
                                 details of a gather from temp to be used later
   prepare_gather()      optional call that allocates buffers for a previously
                           declared gather.  will automatically be called from
                           do_gather() if not done before.
   do_gather()           executes a previously declared gather
   wait_gather()         waits for receives to finish, insuring that the
                           data has actually arrived.
   cleanup_gather()      frees all the buffers that were allocated, WHICH
                           MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   accumulate_gather()   combines gathers into single message tag
   declare_accumulate_gather()  does declare_gather() and accumulate_gather()
                                  in single step.
   declare_accumulate_gather_from_temp()  does declare_gather_from_temp() and
                                            accumulate_gather() in single step.
   start_gather_site()        older function which does declare/prepare/do_gather
                           in a single step
   start_gather_field()  older function which does
                               declare/prepare/do_gather_from_temp
   restart_gather_site()      older function which is obsoleted by do_gather()
   restart_gather_field() older function which is obsoleted by do_gather() 
   start_general_gather()  starts asynchronous sends and receives required
                             to gather fields at arbitrary displacement.
   start_general_gather_from_temp() starts asynchronous sends and receives 
                             required to gather neighbors from a temporary 
			     array of fields.
   wait_general_gather()   waits for receives to finish, insuring that the
                             data has actually arrived, and sets pointers to
			     received data.
   cleanup_general_gather()  frees all the buffers that were allocated, WHICH
                               MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.

*/

#include <time.h>
#include "generic_includes.h"
#include <mpi.h>
#if PRECISION == 1
#define MILC_MPI_REAL MPI_FLOAT
#else
#define MILC_MPI_REAL MPI_DOUBLE
#endif

#define NOWHERE -1	/* Not an index in array of fields */

/* message types used here */
#define SEND_INTEGER_ID    1  /* send an integer to one other node */
#define SEND_FIELD_ID      2  /* id of field sent from one node to another */
#define GENERAL_GATHER_ID  3  /* id used by general_gather routines */
#define GATHER_BASE_ID     4  /* ids greater than or equal to this are used
                                 by the gather routines */

/* macro to compute the message id */
#define GATHER_ID(x) (GATHER_BASE_ID+(x))

/* hacks needed to unify even/odd and 32 sublattice cases */
#ifdef N_SUBL32
#define NUM_SUBL 32
#define FORSOMEPARITY FORSOMESUBLATTICE
#else
#define NUM_SUBL 2
#endif


/**********************************************************************
 *                      INTERNAL DATA TYPES                           *
 **********************************************************************/

/* "comlink" is the basic structure used in gathering neighboring sites.
   Each node will maintain one such structure for each direction for each
   (other) node that contains sites that are neighbors of the sites on
   this node.  For example, if the XUP neighbors of sites on this node
   are found on two other nodes, then this node will maintain a linked
   list of two comlink structures for gathering from the XUP direction.
*/
typedef struct comlink {
  struct comlink *nextcomlink;  /* pointer to next in list, NULL if last */
  int othernode;                /* number of the node to which we connect */
  int n_subl_connected[NUM_SUBL+1];
  /* Number of sites on this node that have neighbors on other node connected
     by this "comlink" of certain parity of the receiver.
     The indicies 0..NUM_SUBL-1 refer to a specific parity and the
     index NUM_SUBL refers to all parities */
  int *sitelist[NUM_SUBL+1];
  /* Address of list of indices of a certain receiver parity whose
     neighbors are found through this comlink.  The index is the same as for
     n_subl_connected above. */
  /* Different comlink structures may point to the same list.
     For example, the receive list for one gather may be a send list for
     the opposite gather. */
} comlink;

/* Linked list type to store id offsets for the sender.
   Needed to match the id that receiver is expecting */
typedef struct id_list_t {
  int id_offset;           /* id offset */
  struct id_list_t *next;  /* linked list */
} id_list_t;

/* structure to hold all necessary info for a gather */
typedef struct gather_t {
  int *neighbor;    /* keeps track if gather neighbor is on our node or not */
  comlink *neighborlist;         /* comlink for receiving messages */
  comlink *neighborlist_send;    /* comlink for sending messages */
  id_list_t *id_list;            /* list of id offsets for sending */
  int n_recv_msgs, n_send_msgs;  /* number of messages to receive and send */ 
  int offset_increment; /* total number of message ids used for this gather */
} gather_t;

/* structure to keep track of details of a declared gather */
typedef struct gmem_t {
  char *mem;            /* source (destination) address for send (receive) */
  int size;             /* size of sent field */
  int stride;           /* stride of source/destination field */
  int num;              /* number of sites in sitelist */
  int *sitelist;        /* sites gathered to/from */
  struct gmem_t *next;  /* linked list */
} gmem_t;

/* Structure to keep track of outstanding sends and receives */
typedef struct {
  int msg_node;         /* node sending or receiving message */
  int id_offset;        /* id offset for this message */
  int msg_size;         /* size of message in bytes */
  char *msg_buf;        /* address of buffer malloc'd for message */
  gmem_t *gmem;         /* linked list explaining detailed usage for buffer */
  MPI_Request msg_req;  /* message handle returned by system call */
} msg_sr_t;

/* structure to store declared gathers
   this is the actual structure used internally
   it has the same name as the typedef which contains this structure which
   the user sees */
struct msg_tag {
  int *ids;          /* array of message ids used in gather */
  int nids;          /* number of message ids used in gather */
  int nrecvs;        /* number of messages to receive in gather */
  int nsends;        /* number of messages to send in gather */
  msg_sr_t *recv_msgs;  /* array of messages to receive */
  msg_sr_t *send_msgs;  /* array of messages to send */
};


/***************************************************
 *  Global variables for the communications stuff  *
 ***************************************************/

/* message ids for gather encode a sequence number for the gather
   so that if several gathers are going at once, you can read
   the message corresponding to the right one. */
/* for computing message id in gather */
/* not needed anymore, but may be used for a check later */
static int id_offset;	    /* label gathers by round-robin */
static int num_gather_ids;  /* number of id offsets allowed */

/* keep track of used ids */
static int *id_array;

/* array storing gather setup info */
static gather_t *gather_array;

/* Number of gathers (mappings) that have been set up */
static int n_gathers, gather_array_len;


/**********************************************************************
 *                BASIC COMMUNICATIONS FUNCTIONS                      *
 **********************************************************************/

void
err_func(MPI_Comm *comm, int *stat, ...)
{
  int len;
  char err_string[MPI_MAX_ERROR_STRING];

  printf("MPI error number: %i\n", *stat);
  MPI_Error_string(*stat, err_string, &len);
  printf("%s\n", err_string);
  terminate(*stat);
}

/*
**  Machine initialization
*/
void
initialize_machine(int argc, char **argv)
{
  int i, flag, *tag_ub;
  MPI_Comm comm;
  MPI_Errhandler errhandler;

  flag = MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  if(flag) err_func(&comm, &flag);
  flag = MPI_Errhandler_create(err_func, &errhandler);
  if(flag) err_func(&comm, &flag);
  flag = MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);
  if(flag) err_func(&comm, &flag);

  /* check if 32 bit int is set correctly */
#ifdef SHORT_IS_32BIT
  if(sizeof(unsigned short)!=4) {
    printf("node %i: SHORT_IS_32BIT is set but sizeof(unsigned short)=%i\n",
	   mynode(), sizeof(unsigned short));
    terminate(1);
  }
#else
  if(sizeof(unsigned int)!=4) {
    printf("node %i: SHORT_IS_32BIT is not set but sizeof(unsigned int)=%i\n",
	   mynode(), sizeof(unsigned int));
    terminate(1);
  }
#endif

  /* get the number of message types */
  /* Note: with MPI-2 MPI_Comm_get_attr is preferred,
     but we keep MPI_Attr_get until MPI-2 is more widely available */ 
  MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
  num_gather_ids = *tag_ub + 1 - GATHER_BASE_ID;
  if(num_gather_ids>1024) num_gather_ids = 1024;

  id_offset = 0;
  id_array = (int *)malloc(num_gather_ids*sizeof(int));
  for(i=0; i<num_gather_ids; ++i) id_array[i] = 0;

  n_gathers = 0;
  gather_array_len = 0;
  gather_array = NULL;
}

/*
**  version of normal exit for multinode processes
*/
void
normal_exit(int status)
{
  time_stamp("exit");
  g_sync();
  MPI_Finalize();
  exit(status);
}

/*
**  version of exit for multinode processes -- kill all nodes
*/
void
terminate(int status)
{
  time_stamp("termination");
  printf("Termination: node %d, status = %d\n", this_node, status);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(status);
}

/*
**  Tell what kind of machine we are on
*/
static char name[]="MPI (portable)";
char *
machine_type(void)
{
  return(name);
}

/*
**  Return my node number
*/
int
mynode(void)
{
  int node;
  MPI_Comm_rank( MPI_COMM_WORLD, &node );
  return(node);
}

/*
**  Return number of nodes
*/
int
numnodes(void)
{
  int nodes;
  MPI_Comm_size( MPI_COMM_WORLD, &nodes );
  return(nodes);
}

/*
**  Synchronize all nodes
*/
void
g_sync(void)
{
  MPI_Barrier( MPI_COMM_WORLD );
}

/*
**  Sum generic floating type over all nodes
*/
void
g_floatsum(Real *fpt)
{
  Real work;
  MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
  *fpt = work;
}

/*
**  Sum a vector of generic floating point types over all nodes
*/
void
g_vecfloatsum(Real *fpt, int length)
{
  Real *work;
  int i;
  work = (Real *)malloc(length*sizeof(Real));
  MPI_Allreduce( fpt, work, length, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
  for(i=0; i<length; i++) fpt[i] = work[i];
  free(work);
}

/*
**  Sum double over all nodes
*/
void
g_doublesum(double *dpt)
{
  double work;
  MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  *dpt = work;
}

/*
**  Sum a vector of doubles over all nodes
*/
void
g_vecdoublesum(double *dpt, int ndoubles)
{
  double *work;
  int i;
  work = (double *)malloc(ndoubles*sizeof(double));
  MPI_Allreduce( dpt, work, ndoubles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  for(i=0; i<ndoubles; i++) dpt[i] = work[i];
  free(work);
}

/*
**  Sum the generic precision complex type over all nodes
*/
void
g_complexsum(complex *cpt)
{
  complex work;
  MPI_Allreduce( cpt, &work, 2, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
  *cpt = work;
}

/*
**  Sum a vector of the generic precision complex type over all nodes
*/
void
g_veccomplexsum(complex *cpt, int ncomplex)
{
  complex *work;
  int i;
  work = (complex *)malloc(ncomplex*sizeof(complex));
  MPI_Allreduce( cpt, work, 2*ncomplex, MILC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD );
  for(i=0; i<ncomplex; i++) cpt[i] = work[i];
  free(work);
}

/*
**  Sum double_complex over all nodes
*/
void
g_dcomplexsum(double_complex *cpt)
{
  double_complex work;
  MPI_Allreduce( cpt, &work, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  *cpt = work;
}

/*
**  Sum a vector of double_complex over all nodes
*/
void
g_vecdcomplexsum(double_complex *cpt, int ncomplex)
{
  double_complex *work;
  int i;
  work = (double_complex *)malloc(ncomplex*sizeof(double_complex));
  MPI_Allreduce( cpt, work, 2*ncomplex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  for(i=0; i<ncomplex; i++) cpt[i] = work[i];
  free(work);
}

/*
**  Sum wilson_vector over all nodes
*/
void
g_wvectorsumfloat(wilson_vector *wvpt)
{
  g_veccomplexsum((complex *)wvpt, 12);
}

/*
**  Global exclusive or acting on u_int32type
*/
void
g_xor32(u_int32type *pt)
{
  u_int32type work;
#ifdef SHORT_IS_32BIT
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED_SHORT, 
		 MPI_BXOR, MPI_COMM_WORLD );
#else
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED, 
		 MPI_BXOR, MPI_COMM_WORLD );
#endif
  *pt = work;
}

/*
**  Find maximum of the generic precision floating point type over all nodes
*/
void
g_floatmax(Real *fpt)
{
  Real work;
  MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_MAX, MPI_COMM_WORLD );
  *fpt = work;
}

/*
**  Find maximum of double over all nodes
*/
void
g_doublemax(double *dpt)
{
  double work;
  MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  *dpt = work;
}

/*
**  Broadcast generic precision floating point number from node zero
*/
void
broadcast_float(Real *fpt)
{
  MPI_Bcast( fpt, 1, MILC_MPI_REAL, 0, MPI_COMM_WORLD );
}

/*
**  Broadcast double precision floating point number from node zero
*/
void
broadcast_double(double *dpt)
{
  MPI_Bcast( dpt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}

/*
**  Broadcast generic precision complex number from node zero
*/
void
broadcast_complex(complex *cpt)
{
  MPI_Bcast( cpt, 2, MILC_MPI_REAL, 0, MPI_COMM_WORLD );
}

/*
**  Broadcast double precision complex number from node zero
*/
void
broadcast_dcomplex(double_complex *cpt)
{
  MPI_Bcast( cpt, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}

/*
**  Broadcast bytes from node 0 to all others
*/
void
broadcast_bytes(char *buf, int size)
{
  MPI_Bcast( buf, size, MPI_BYTE, 0, MPI_COMM_WORLD );
}


/******************************
 *  SEND AND RECEIVE INTEGER  *
 ******************************/

/*
**  Send an integer to one other node
**  This is to be called only by the node doing the sending
*/
void
send_integer(int tonode, int *address)
{
  MPI_Send( address, 1, MPI_INT, tonode, SEND_INTEGER_ID, MPI_COMM_WORLD );
}

/*
**  Receive an integer from another node
*/
void
receive_integer(int fromnode, int *address)
{
  MPI_Status status;
  MPI_Recv( address, 1, MPI_INT, fromnode, SEND_INTEGER_ID,
	    MPI_COMM_WORLD, &status );
}


/****************************
 *  SEND AND RECEIVE FIELD  *
 ****************************/

/*
**  send_field is to be called only by the node doing the sending
*/
void
send_field(char *buf, int size, int tonode)
{
  MPI_Send( buf, size, MPI_BYTE, tonode, SEND_FIELD_ID, MPI_COMM_WORLD );
}

/*
**  get_field is to be called only by the node to which the field was sent
*/
void
get_field(char *buf, int size, int fromnode)
{
  MPI_Status status;
  MPI_Recv( buf, size, MPI_BYTE, fromnode, SEND_FIELD_ID, MPI_COMM_WORLD,
	    &status );
}


/*********************
 *  TIMING ROUTINES  *
 *********************/

/*
**  Double precision CPU time in seconds
*/
double
dclock_cpu(void)
{ 
  long fine;
  fine = clock();
  return( ((double)fine)/CLOCKS_PER_SEC);
}

/*
**  Double precision wall clock time in seconds
*/
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
double dclock(void){
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
#else
double dclock(void){
  return dclock_cpu();
}
#endif

/*
**  Print time stamp
*/
void
time_stamp(char *msg)
{
  time_t time_stamp;

  if(mynode()==0){
    time(&time_stamp);
    printf("%s: %s\n", msg, ctime(&time_stamp));
    fflush(stdout);
  }
}


/**********************************************************************
 *                  FUNCTIONS USED FOR GATHERS                        *
 **********************************************************************/

/*
**  sort a list of eight gather_t structures into the order we want for the
**  nearest neighbor gathers:  XUP,YUP,ZUP,TUP,TDOWN,ZDOWN,YDOWN,XDOWN,
**  starting from the index for the first pointer
*/
void
sort_eight_gathers(int index)
{
  gather_t tt[8];
  int i;

  for(i=0; i<8; i++) memcpy(&tt[i], &gather_array[index+i], sizeof(gather_t));
  for(i=XUP; i<=TUP; i++) {
    memcpy(&gather_array[index+i], &tt[2*i], sizeof(gather_t));
    memcpy(&gather_array[index+OPP_DIR(i)], &tt[2*i+1], sizeof(gather_t));
  }
}

/*
**  utility function for finding coordinates of neighbor
**  This version for use by make_gather for nearest neighbor gathers
*/
static void
neighbor_coords_special(
  int x, int y, int z, int t,       /* coordinates of site */
  int *dirpt,                       /* direction (eg XUP) */
  int fb,                           /* "forwards/backwards"  */
  int *x2p, int *y2p, int *z2p, int *t2p)
                                    /* pointers to coordinates of neighbor */
{
  int dir;

  dir = (fb==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *x2p = x; *y2p = y; *z2p = z; *t2p = t;
  switch(dir) {
    case XUP   : *x2p = (x+1)%nx; break;
    case XDOWN : *x2p = (x+nx-1)%nx; break;
    case YUP   : *y2p = (y+1)%ny; break;
    case YDOWN : *y2p = (y+ny-1)%ny; break;
    case ZUP   : *z2p = (z+1)%nz; break;
    case ZDOWN : *z2p = (z+nz-1)%nz; break;
    case TUP   : *t2p = (t+1)%nt; break;
    case TDOWN : *t2p = (t+nt-1)%nt; break;
    default: printf("BOTCH: bad direction\n"); terminate(1);
  }
}

/*
**  Set up "comlink" structures needed by nearest neighbor gather routines.
**  make_lattice() must be called first.
*/
void
make_nn_gathers(void)
{
  int i, gather_parity;

  if(n_gathers!=0) {
    printf("error: make_nn_gathers must come before any make_gather\n");
    terminate(1);
  }

  gather_array_len = 8;
  gather_array = (gather_t *)malloc(gather_array_len*sizeof(gather_t));
  if(gather_array==NULL) {
    printf("error: not enough room for gather_array in make_nn_gathers\n");
    terminate(1);
  }

  if((nx&1)||(ny&1)||(nz&1)||(nt&1)) gather_parity = SCRAMBLE_PARITY;
  else gather_parity = SWITCH_PARITY;

  for(i=XUP; i<=TUP; i++)
    make_gather( neighbor_coords_special, &i, WANT_INVERSE,
                 ALLOW_EVEN_ODD, gather_parity );

  /* Sort into the order we want for nearest neighbor gathers,
     so you can use XUP, XDOWN, etc. as argument in calling them. */
  sort_eight_gathers( 0 );
}


/**********************************************************************
 *                  FUNCTIONS USED TO MAKE GATHERS                    *
 **********************************************************************/

#define RECEIVE 0
#define SEND    1

static int
parity_function(int x, int y, int z, int t)
{
#ifndef N_SUBL32
  return (x+y+z+t)&1;
#else
  return (x%2) + 2*(y%2) + 4*(z%2) + 8*(t%2) + 16*((x/2+y/2+z/2+t/2)%2);
#endif
}

/*
** copy a linked list of comlinks, switching send and receive parity
*/
static comlink *
copy_list_switch(comlink *old_compt, int *send_subl)
{
  comlink *firstpt, *compt;
  int r_subl, s_subl;

  if(old_compt==NULL) return(NULL);

  firstpt = compt = (comlink *)malloc( sizeof(comlink) );
  do{
    compt->othernode = old_compt->othernode;
    for(r_subl=0; r_subl<NUM_SUBL; r_subl++) {
      s_subl = send_subl[r_subl];
      compt->n_subl_connected[s_subl] = old_compt->n_subl_connected[r_subl];
      compt->sitelist[s_subl] = old_compt->sitelist[r_subl];
    }
    compt->n_subl_connected[NUM_SUBL] = old_compt->n_subl_connected[NUM_SUBL];
    compt->sitelist[NUM_SUBL] = old_compt->sitelist[NUM_SUBL];
    if( old_compt->nextcomlink != NULL)
      compt->nextcomlink = (comlink *)malloc( sizeof(comlink) );
    else compt->nextcomlink = NULL;
    old_compt = old_compt->nextcomlink; 
    compt = compt->nextcomlink;
  } while( old_compt!=NULL );
  return(firstpt);
}

/*
**  sort a list of sites according to the order of the sites on the
**  node with which they comunicate
*/
static void
sort_site_list(
  int n,		/* number of elements in list */
  int *list,	        /* pointer to list */
  void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
                        /* function which defines mapping */
  int *args,	        /* arguments to pass to function */
  int forw_back)	/* look forwards or backwards in map */
{
  register int j,k,in1,in2,flag;
  register site *s;
  int x,y,z,t;
  int *key;

  if(n==0) return;
  key = (int *)malloc(n*sizeof(int));
  if(key == NULL) {
    printf("sort_site_list(%d): no room for key\n",mynode());
    terminate(1);
  }

  /* Construct sort key */
  for(j=0; j<n; j++) {
    s = &(lattice[list[j]]);
    func(s->x,s->y,s->z,s->t,args,forw_back,&x,&y,&z,&t);
    key[j] = node_index(x,y,z,t);
  }

  /* bubble sort, if this takes too long fix it later */
  for(j = n-1; j>0; j--) {
    flag=0;
    for(k=0; k<j; k++){
      in1 = key[k];
      in2 = key[k+1];
      if(in1>in2){
	flag=1;
	key[k]   = in2;
	key[k+1] = in1;
	in1 = list[k];
	list[k] = list[k+1];
	list[k+1] = in1;
      }
    }
    if(flag==0)break;
  }
  free(key);
}

/*
**  make comlink for send or receive
*/
static comlink *
make_send_receive_list(
  void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
        		/* function which defines sites to gather from */
  int *args,		/* list of arguments, to be passed to function */
  int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int forw_back,	/* FORWARDS or BACKWARDS */
  int send_recv,        /* SEND or RECEIVE list */
  int *n_msgs)          /* returns number of messages in list */
{
  int i,j,subl;	        /* scratch */
  site *s;	        /* scratch */
  int x,y,z,t;		/* coordinates */
  int *sbuf[NUM_SUBL];	/* to be malloc'd */
  int *tbuf;	        /* to be malloc'd */
  comlink **combuf;	/* to be malloc'd, remember where comlinks are */
  comlink *compt,**comptpt;
  comlink *firstpt;

  /* make temporary buffers of numnodes() integers to count numbers of
     neighbors in each sublattice on each node */
  for(subl=0; subl<NUM_SUBL; subl++) {
    sbuf[subl] = (int *)malloc( numnodes()*sizeof(int) );
    /* clear neighbor_numbers */
    for(i=0; i<numnodes(); i++) sbuf[subl][i] = 0;
  }
  tbuf = (int *)malloc( numnodes()*sizeof(int) );
  for(i=0; i<numnodes(); i++) tbuf[i] = 0;
  combuf = (comlink **)malloc( numnodes()*sizeof(comlink *) );

  /* scan sites in lattice */
  FORALLSITES(i,s) {
    /* find coordinates, node, and sublattice of receiving site */
    if( send_recv==RECEIVE ) {
      func( s->x, s->y, s->z, s->t, args, forw_back, &x, &y, &z, &t );
      subl = parity_function(s->x,s->y,s->z,s->t);
    }
    else {  /* SEND */
      func( s->x, s->y, s->z, s->t, args, -forw_back, &x, &y, &z, &t );
      subl = parity_function(x,y,z,t);
    }
    j = node_number(x,y,z,t);

    /* if site is off node, increment neighbor_counter */
    if( j != mynode() ) {
      ++tbuf[j];
      if(want_even_odd==NO_EVEN_ODD) subl = 0;
      ++sbuf[subl][j];
    }
  }

  *n_msgs = 0;
  firstpt = NULL;
  comptpt = &firstpt;
  /* for each neighbor_counter that is nonzero, create a comlink */
  for(j=0; j<numnodes(); j++) {
    if( j==mynode() ) continue;  /* not for local node */
    if( tbuf[j]==0 ) continue;   /* no neighbors on this node */

    compt = (comlink *)malloc( sizeof(comlink) );
    *comptpt = compt;
    combuf[j] = compt;	/* to make it easy to find again */
    compt->nextcomlink = NULL;	/* currently terminates list */
    compt->othernode = j;
    compt->n_subl_connected[NUM_SUBL] = tbuf[j];
    for(subl=0; subl<NUM_SUBL; subl++) {
      compt->n_subl_connected[subl] = sbuf[subl][j];
    }
    compt->sitelist[0] = compt->sitelist[NUM_SUBL] =
      (int *)malloc( tbuf[j]*sizeof(int) );
    for(subl=1; subl<NUM_SUBL; subl++)
      compt->sitelist[subl] = (compt->sitelist[subl-1]) + sbuf[subl-1][j];
    /* sitelist[...] must be filled in later */
    comptpt = &(compt->nextcomlink);	/* linked list, if we
		      extend it this will get address of next comlink. */
    ++(*n_msgs);
  }

  /* clear neighbor_numbers, to be used as counters now */
  for(subl=0; subl<NUM_SUBL; subl++) {
    for(i=0; i<numnodes(); i++) sbuf[subl][i] = 0;
  }

  /* scan sites in node again */
  FORALLSITES(i,s){
    /* find coordinates, node, and sublattice of receiving site */
    if( send_recv==RECEIVE ){
      func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
      subl = parity_function(s->x,s->y,s->z,s->t);
    }
    else {  /* SEND */
      func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
      subl = parity_function(x,y,z,t);
    }
    j = node_number(x,y,z,t);

    /* if neighbor is offnode, add to list in appropriate comlink */
    if( j != mynode() ){
      if(want_even_odd==NO_EVEN_ODD) subl = 0;
      combuf[j]->sitelist[subl][sbuf[subl][j]] = i;
      ++sbuf[subl][j];
    }
  }
  /* sort the lists of links according to the ordering of their
     even neighbors in the lower numbered node.  The list of sites
     on the lower numbered node is already in order. */
  for(compt=firstpt; compt != NULL; compt=compt->nextcomlink) {
    if(compt->othernode > this_node)continue;
    /* this is lower numbered node, so don't sort */
    if( send_recv==RECEIVE ) i = forw_back;
    else i = -forw_back;
    for(subl=0; subl<NUM_SUBL; subl++)
      sort_site_list( compt->n_subl_connected[subl],
		      compt->sitelist[subl], func, args, i );
  }

  /* free temporary storage */
  free(combuf);
  free(tbuf);
  for(subl=0; subl<NUM_SUBL; subl++) free(sbuf[subl]);

  return(firstpt);
}

/*
**  determine tag offsets needed by sender
*/
static id_list_t *
make_id_list(
  comlink *recv,       /* neighborlist */
  int n_recv,          /* number of receives */
  comlink *send)       /* neighborlist_send */
{
  int i, *buf;
  id_list_t *tol_top, *tol, **tol_next;
  MPI_Request *req, sreq;
  MPI_Status stat;

  buf = (int *)malloc(n_recv*sizeof(int));
  req = (MPI_Request *)malloc(n_recv*sizeof(MPI_Request));

  for(i=0; recv!=NULL; ++i, recv=recv->nextcomlink) {
    buf[i] = i;
    MPI_Isend( &buf[i], 1, MPI_INT, recv->othernode, 0, MPI_COMM_WORLD,
	       &req[i] );
  }
  if(i!=n_recv) {printf("error i!=n_recv\n"); terminate(1);}

  tol_next = &tol_top;
  while(send!=NULL) {
    tol = *tol_next = (id_list_t *)malloc(sizeof(id_list_t));
    MPI_Irecv( &i, 1, MPI_INT, send->othernode, 0, MPI_COMM_WORLD, &sreq );
    MPI_Wait( &sreq, &stat );
    tol->id_offset = i;
    tol_next = &(tol->next);
    send = send->nextcomlink;
  }
  *tol_next = NULL;

  for(i=0; i<n_recv; ++i) {
    MPI_Wait( &req[i], &stat );
  }

  free(req);
  free(buf);

  return tol_top;
}

/*
**  determine max number of ids needed for gather
*/
static int
get_max_receives(int n_recv)
{
  int work;
  MPI_Allreduce( &n_recv, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  return work;
}

/*
**  add another gather to the list of tables
*/
int
make_gather(
  void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
                        /* function which defines sites to gather from */
  int *args,		/* list of arguments, to be passed to function */
  int inverse,		/* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
  int want_even_odd,	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int parity_conserve)	/* {SAME,SWITCH,SCRAMBLE}_PARITY */
{
  int i,j,subl;	        /* scratch */
  site *s;	        /* scratch */
  int dir;		/* direction */
  int x,y,z,t;		/* coordinates */
  int *send_subl;       /* sublist of sender for a given receiver */

  /* we will have one or two more gathers */
  if( inverse==WANT_INVERSE ) n_gathers += 2;
  else			      n_gathers += 1;

  if(n_gathers>gather_array_len) {
    gather_array_len = n_gathers;
    /* lengthen gather array to add more gathers */
    gather_array =
      (gather_t *)realloc(gather_array, gather_array_len*sizeof(gather_t));
  }

  dir = n_gathers - 1;	/* index of gather we are working on */
  gather_array[dir].neighbor = (int *)malloc( sites_on_node*sizeof(int) );
  if( gather_array[dir].neighbor==NULL ) {
    printf("make_gather: NODE %d: no room for neighbor vector\n",this_node);
    terminate(1);
  }
  if( inverse==WANT_INVERSE ) {
    dir = n_gathers - 2;	/* index of gather we are working on */
    gather_array[dir].neighbor = (int *)malloc( sites_on_node*sizeof(int) );
    if( gather_array[dir].neighbor==NULL ) {
      printf("make_gather: NODE %d no room for neighbor vector\n",this_node);
      terminate(1);
    }
  }

  if( want_even_odd==ALLOW_EVEN_ODD && parity_conserve!=SCRAMBLE_PARITY ) {
    send_subl = (int *)malloc(NUM_SUBL*sizeof(int));
    if(send_subl==NULL){
      printf("NODE %d: no room for send_subl\n",this_node);
      terminate(1);
    }
    for(subl=0; subl<NUM_SUBL; subl++) send_subl[subl] = NOWHERE;
  } else {
    send_subl = NULL;
  }

  /* Check to see if mapping has advertised parity and inverse properties */
  /* Also check to see if it returns legal values for coordinates */
  FORALLSITES(i,s) {
    /* find coordinates of neighbor who sends us data */
    func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);

    if( x<0 || y<0 || z<0  || t<0 || x>=nx || y>=ny || z>=nz || t>=nt){
      printf("DUMMY! Your gather mapping does not stay in lattice\n");
      printf("It mapped %d %d %d %d to %d %d %d %d\n",
	     s->x,s->y,s->z,s->t,x,y,z,t);
      terminate(1);
    }

    if(parity_conserve!=SCRAMBLE_PARITY) {
      int r_subl, s_subl;

      r_subl = parity_function(s->x,s->y,s->z,s->t);
      s_subl = parity_function(x,y,z,t);

      if( want_even_odd==ALLOW_EVEN_ODD ) {
	if( send_subl[r_subl] == NOWHERE ) {
	  send_subl[r_subl] = s_subl;
	}
	else if( send_subl[r_subl] != s_subl ){
	  printf("DUMMY! Your gather mixes up sublattices: %d vs %d\n",
		 send_subl[r_subl], s_subl);
	  printf("on mapping %i %i %i %i -> %i %i %i %i\n",
		 s->x,s->y,s->z,s->t, x,y,z,t);
	  terminate(1);
	}
      }

      if( parity_conserve==SAME_PARITY && s_subl!=r_subl ){
	printf("DUMMY! Your gather mapping does not obey claimed parity");
	printf(", namely SAME_PARITY\n");
	printf("It mapped %d %d %d %d with %d to %d %d %d %d with %d\n",
	       s->x,s->y,s->z,s->t,r_subl,x,y,z,t,s_subl);
	terminate(1);
      }
      if( parity_conserve==SWITCH_PARITY && s_subl==r_subl ){
	printf("DUMMY! Your gather mapping does not obey claimed parity");
	printf(", namely SWITCH_PARITY\n");
	printf("It mapped %d %d %d %d with %d to %d %d %d %d with %d\n",
	       s->x,s->y,s->z,s->t,r_subl,x,y,z,t,s_subl);
	terminate(1);
      }

      if( inverse==OWN_INVERSE ) {
	int x2,y2,z2,t2;
	func( x, y, z, t, args, FORWARDS, &x2,&y2,&z2,&t2);
	if( s->x!=x2 || s->y!=y2 || s->z!=z2 || s->t!=t2 ) {
	  printf("DUMMY! Your gather mapping is not its own inverse\n");
	  printf("It's square mapped %d %d %d %d to %d %d %d %d\n",
		 s->x,s->y,s->z,s->t,x2,y2,z2,t2);
	  terminate(1);
	}
      }
    }
  }

  /* RECEIVE LISTS: */
  /* Fill in pointers to sites which are on this node, NOWHERE if
     they are off-node */
  FORALLSITES(i,s){
    /* find coordinates of neighbor who sends us data */
    func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
    j = node_number(x,y,z,t);	/* node for neighbor site */
    /* if neighbor is on node, set up pointer */
    if( j == mynode() ) gather_array[dir].neighbor[i] = node_index(x,y,z,t);
    else		gather_array[dir].neighbor[i] = NOWHERE;
  }

  /* make lists of sites which get data from other nodes.  */
  gather_array[dir].neighborlist =
    make_send_receive_list( func, args, want_even_odd, FORWARDS, RECEIVE,
			    &gather_array[dir].n_recv_msgs );

  /* SEND LISTS: */
  /* Now make lists of sites to which we send */
  /* Under some conditions, if mapping is its own inverse we can use
     the lists we have already made */
  if( inverse==OWN_INVERSE && 
      ( want_even_odd!=ALLOW_EVEN_ODD || parity_conserve!=SCRAMBLE_PARITY ) ) {
    if( want_even_odd==NO_EVEN_ODD || parity_conserve==SAME_PARITY ) {
      gather_array[dir].neighborlist_send = gather_array[dir].neighborlist;
      gather_array[dir].n_send_msgs = gather_array[dir].n_recv_msgs;
    } else {
      gather_array[dir].neighborlist_send =
	copy_list_switch( gather_array[dir].neighborlist, send_subl );
      gather_array[dir].n_send_msgs = gather_array[dir].n_recv_msgs;
    }
  } else {
    /* Make new linked list of comlinks for send lists */
    gather_array[dir].neighborlist_send =
      make_send_receive_list( func, args, want_even_odd, FORWARDS, SEND,
			      &gather_array[dir].n_send_msgs );
  } /* End general case for send lists */

  gather_array[dir].id_list =
    make_id_list( gather_array[dir].neighborlist,
		  gather_array[dir].n_recv_msgs,
		  gather_array[dir].neighborlist_send );

  gather_array[dir].offset_increment =
    get_max_receives( gather_array[dir].n_recv_msgs );

  if( inverse != WANT_INVERSE ) {
    free(send_subl);
    return(dir);
  }

  /******************
   * INVERSE GATHER *
   ******************/

  /* Now, if necessary, make inverse gather */
  /* In most cases, we can use the same lists as the gather, in one
     form or another.  Of course, by the time you get to here
     you know that inverse = WANT_INVERSE */
  dir++;	/* inverse gather has direction one more than original */

  /* Always set up pointers to sites on this node */
  /* scan sites in lattice */
  FORALLSITES(i,s) {
    /* find coordinates of neighbor who sends us data */
    func( s->x, s->y, s->z, s->t, args, BACKWARDS, &x,&y,&z,&t);
    j = node_number(x,y,z,t);	/* node for neighbor site */

    /* if neighbor is on node, set up pointer */
    if( j == mynode() ) gather_array[dir].neighbor[i] = node_index(x,y,z,t);
    else 		gather_array[dir].neighbor[i] = NOWHERE;
  }

  if( parity_conserve==SAME_PARITY || want_even_odd==NO_EVEN_ODD ) {
    /* Use same comlinks as inverse gather, switching send and receive.
       Nearest neighbor gathers are an example of this case. */
    gather_array[dir].neighborlist = gather_array[dir-1].neighborlist_send;
    gather_array[dir].neighborlist_send = gather_array[dir-1].neighborlist;
    gather_array[dir].n_recv_msgs = gather_array[dir-1].n_send_msgs;
    gather_array[dir].n_send_msgs = gather_array[dir-1].n_recv_msgs;
  } else if( parity_conserve==SWITCH_PARITY ) {
    /* make new comlinks, but use same lists as inverse gather, switching
       send and receive, switching even and odd. */
    gather_array[dir].neighborlist = 
      copy_list_switch( gather_array[dir-1].neighborlist_send, send_subl );
    gather_array[dir].neighborlist_send = 
      copy_list_switch( gather_array[dir-1].neighborlist, send_subl );
    gather_array[dir].n_recv_msgs = gather_array[dir-1].n_send_msgs;
    gather_array[dir].n_send_msgs = gather_array[dir-1].n_recv_msgs;
  } else {  /* general case.  Really only get here if ALLOW_EVEN_ODD
	       and SCRAMBLE_PARITY */
    /* RECEIVE LISTS */
    gather_array[dir].neighborlist =
      make_send_receive_list( func, args, want_even_odd, BACKWARDS, RECEIVE,
			      &gather_array[dir].n_recv_msgs );
    /* SEND LISTS */
    gather_array[dir].neighborlist_send =
      make_send_receive_list( func, args, want_even_odd, BACKWARDS, SEND,
			      &gather_array[dir].n_send_msgs );
  } /* End making new lists for inverse gather */

  gather_array[dir].id_list =
    make_id_list( gather_array[dir].neighborlist,
		  gather_array[dir].n_recv_msgs,
		  gather_array[dir].neighborlist_send );

  gather_array[dir].offset_increment =
    get_max_receives(gather_array[dir].n_recv_msgs);

  free(send_subl);
  return(dir-1);
}


/**********************************************************************
 *                         GATHER ROUTINES                            *
 **********************************************************************

 declare_strided_gather() returns a pointer to msg_tag which will
   be used as input to subsequent prepare_gather() (optional), do_gather(),
   wait_gather() and cleanup_gather() calls.

   This handles gathers from both field offset and temp and is not called
   directly by the user.  Instead they should call declare_gather() or
   declare_gather_from_temp().

 prepare_gather() allocates buffers needed for the gather.  This call is
   optional since it will automatically be called from do_gather() if
   not explicitly called before.

 do_gather() starts the actual gather.  This may be repeated after a
    wait_gather() to repeat the exact same gather.

 wait_gather() waits for the gather to finish.

 cleanup_gather() frees memory allocated for the gather including the msg_tag.

   example:
	msg_tag *tag;
	tag = declare_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	                      EVEN, gen_pt[0] );
        prepare_gather(tag);  ** this step is optional **
        do_gather(tag);
	  ** do other stuff, but don't modify tag or gen_pt[0] **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here, but don't modify tag or
	   gen_pt[0].
	   Do modify the source field phi. **
	do_gather(tag);
	  ** do other stuff **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the modified phi.
	   The restart-wait may be repeated as often as desired.  **
	cleanup_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

*/

/*
**  returns msg_tag containing details for specific gather
**  handles gathers from both field offset and temp
*/
static msg_tag *
declare_strided_gather(
  void *field,	        /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
  char ** dest)		/* one of the vectors of pointers */
{
  int i;	        /* scratch */
  site *s;	        /* scratch pointer to site */
  msg_tag *mtag;	/* message tag structure we will return a pointer to */
  msg_sr_t *mrecv, *msend; /* arrays for send and receive lists */
  gmem_t *gmem;
  comlink *compt;	/* pointer to current comlink */
  gather_t *gt;         /* pointer to current gather */
  id_list_t *idl;

  gt = &gather_array[index];

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) */
  if(subl==EVENANDODD) {
    FORALLSITES(i,s){ if(gt->neighbor[i] != NOWHERE){
      dest[i] = (char *)field + gt->neighbor[i]*stride;
    }}
  } else {
    FORSOMEPARITY(i,s,subl){ if(gt->neighbor[i] != NOWHERE){
      dest[i] = (char *)field + gt->neighbor[i]*stride;
    }}
  }

#ifndef N_SUBL32
  switch(subl) {
    case EVEN:        subl = 0; break;
    case ODD:         subl = 1; break;
    case EVENANDODD:  subl = 2; break;
    default:  printf("ERROR: bad sublattice\n"); terminate(subl);
  }
#else
  if(subl==EVENANDODD) subl = NUM_SUBL;
#endif

  /*  allocate the message tag */
  mtag = (msg_tag *)malloc(sizeof(msg_tag));

  mtag->nids = gt->offset_increment;
  mtag->ids = NULL;

  /* allocate a buffer for the msg_sr_t's.  This is dynamically allocated
     because there may be an arbitrary number of gathers in progress
     in any direction. */

  for( i=0, compt = gt->neighborlist; compt != NULL;
       compt = compt->nextcomlink ) {
    if(compt->n_subl_connected[subl]!=0) ++i;
  }
  mtag->nrecvs = i;
  if( gt->n_recv_msgs==0 ) mrecv = NULL;
  else {
    mrecv = (msg_sr_t *)malloc(gt->n_recv_msgs*sizeof(msg_sr_t));
    if(mrecv==NULL) {
      printf("NO ROOM for mrecv, node %d\n", mynode());
      terminate(1);
    }
  }
  mtag->recv_msgs = mrecv;

  for( i=0, compt = gt->neighborlist_send; compt != NULL;
       compt = compt->nextcomlink ) {
    if(compt->n_subl_connected[subl]!=0) ++i;
  }
  mtag->nsends = i;
  if( gt->n_send_msgs==0 ) msend = NULL;
  else {
    msend = (msg_sr_t *)malloc(gt->n_send_msgs*sizeof(msg_sr_t));
    if(msend==NULL) {
      printf("NO ROOM for msend, node %d\n", mynode());
      terminate(1);
    }
  }
  mtag->send_msgs = msend;

  /* for each node which has neighbors of my sites */
  for( i=0, compt = gt->neighborlist; compt != NULL;
       i++, compt = compt->nextcomlink ) {
    if(compt->n_subl_connected[subl]==0) continue;
    mrecv[i].msg_node = compt->othernode;
    mrecv[i].id_offset = i;
    mrecv[i].msg_size = size*compt->n_subl_connected[subl];
    mrecv[i].msg_buf = NULL;
    gmem = (gmem_t *)malloc(sizeof(gmem_t));
    mrecv[i].gmem = gmem;
    gmem->num = compt->n_subl_connected[subl];
    gmem->sitelist = compt->sitelist[subl];
    gmem->mem = (char *)dest;
    gmem->stride = sizeof(char *);
    gmem->size = size;
    gmem->next = NULL;
  }

  /* for each node whose neighbors I have */
  idl = gt->id_list;
  for( i=0, compt = gt->neighborlist_send; compt != NULL;
       i++, compt = compt->nextcomlink, idl = idl->next ) {
    if(compt->n_subl_connected[subl]==0) continue;
    msend[i].msg_node = compt->othernode;
    msend[i].id_offset = idl->id_offset;
    msend[i].msg_size = size*compt->n_subl_connected[subl];
    msend[i].msg_buf = NULL;
    gmem = (gmem_t *)malloc(sizeof(gmem_t));
    msend[i].gmem = gmem;
    gmem->num = compt->n_subl_connected[subl];
    gmem->sitelist = compt->sitelist[subl];
    gmem->mem = field;
    gmem->stride = stride;
    gmem->size = size;
    gmem->next = NULL;
  }

  return mtag;
}

/*
**  allocate buffers for gather
*/
void
prepare_gather(msg_tag *mtag)
{
  int i, j, nids;
  int *ids;
  msg_sr_t *mrecv,*msend;
  gmem_t *gmem;
  char *tpt;

  if(mtag->ids!=NULL) {
    printf("error: already prepared\n");
    terminate(1);
  }

  nids = mtag->nids;
  if(nids!=0) {
    mtag->ids = ids = (int *)malloc(nids*sizeof(int));
    for(i=0, j=id_offset; i<nids; i++, j=(j+1)%num_gather_ids) {
      /* find next available type */
      while(id_array[j]!=0) {
	j = (j+1)%num_gather_ids;
	if(j==id_offset) {
	  printf("error: not enough message ids\n");
	  terminate(1);
	}
      }
      ids[i] = j;
      id_array[j] = 1;
    }
    id_offset = j;
  }

  mrecv = mtag->recv_msgs;
  /* for each node which has neighbors of my sites */
  for(i=0; i<mtag->nrecvs; ++i) {
    if(mrecv[i].msg_size==0) {
      node0_printf("error: unexpected zero msg_size\n");
      terminate(1);
    }
    mrecv[i].msg_buf = tpt = (char *)malloc( mrecv[i].msg_size );
    if(tpt==NULL) {
      printf("NO ROOM for msg_buf, node %d\n", mynode());
      terminate(1);
    }
    /* set pointers in sites to correct location */
    gmem = mrecv[i].gmem;
    do {
      for(j=0; j<gmem->num; ++j,tpt+=gmem->size) {
	((char **)gmem->mem)[gmem->sitelist[j]] = tpt;
      }
    } while((gmem=gmem->next)!=NULL);
  }

  msend = mtag->send_msgs;
  /* for each node whose neighbors I have */
  for(i=0; i<mtag->nsends; ++i) {
    msend[i].msg_buf = (char *)malloc( msend[i].msg_size );
    if(msend[i].msg_buf==NULL) {
      printf("NO ROOM for msg_buf, node %d\n",mynode());
      terminate(1);
    }
  }
}

/*
**  actually execute the gather
*/
void
do_gather(msg_tag *mtag)  /* previously returned by start_gather_site */
{
  register int i,j;	/* scratch */
  register char *tpt;	/* scratch pointer in buffers */
  msg_sr_t *mbuf;
  gmem_t *gmem;

  if((mtag->ids==NULL)&&(mtag->nids!=0)) prepare_gather(mtag);

  mbuf = mtag->recv_msgs;
  /* for each node which has neighbors of my sites */
  for(i=0; i<mtag->nrecvs; i++) {
    /* post receive */
    MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE, MPI_ANY_SOURCE,
	       GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_WORLD,
	       &mbuf[i].msg_req );
  }

  mbuf = mtag->send_msgs;
  /* for each node whose neighbors I have */
  for(i=0; i<mtag->nsends; ++i) {
    /* gather data into the buffer */
    tpt = mbuf[i].msg_buf;
    gmem = mbuf[i].gmem;
    do {
      for(j=0; j<gmem->num; ++j,tpt+=gmem->size) {
	memcpy( tpt, gmem->mem + gmem->sitelist[j]*gmem->stride, gmem->size );
      }
    } while((gmem=gmem->next)!=NULL);
    /* start the send */
    MPI_Isend( mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE, mbuf[i].msg_node,
	       GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_WORLD,
	       &mbuf[i].msg_req );
  }
}

/*
**  wait for gather to finish
*/
void
wait_gather(msg_tag *mtag)
{
  MPI_Status status;
  int i;

  /* wait for all receive messages */
  for(i=0; i<mtag->nrecvs; i++) {
    MPI_Wait( &mtag->recv_msgs[i].msg_req, &status );
  }

  /* wait for all send messages */
  for(i=0; i<mtag->nsends; i++) {
    MPI_Wait( &mtag->send_msgs[i].msg_req, &status );
  }
}

/*
**  free buffers associated with message tag
*/
void
cleanup_gather(msg_tag *mtag)
{
  int i;
  gmem_t *gmem, *next;

  if(mtag->ids!=NULL)
    for(i=0; i<mtag->nids; ++i) id_array[mtag->ids[i]] = 0;

  /* free all receive buffers */
  for(i=0; i<mtag->nrecvs; i++) {
    free( mtag->recv_msgs[i].msg_buf );
    gmem = mtag->recv_msgs[i].gmem;
    do {
      next = gmem->next;
      free(gmem);
      gmem = next;
    } while(gmem!=NULL);
  }
  /*  free all send buffers */
  for(i=0; i<mtag->nsends; i++) {
    free( mtag->send_msgs[i].msg_buf );
    gmem = mtag->send_msgs[i].gmem;
    do {
      next = gmem->next;
      free(gmem);
      gmem = next;
    } while(gmem!=NULL);
  }
  /* free the msg_tag buffer */
  free(mtag->recv_msgs);
  free(mtag->send_msgs);
  free(mtag->ids);
  free(mtag);
}


/***********************************************************************
 *                 Convenience Routines for Gathers                    *
 ***********************************************************************/

/*
**  declare gather with a field offset
*/
msg_tag *
declare_gather(
  field_offset field,	/* which field? Some member of structure "site" */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  return declare_strided_gather( (char *)lattice + field, sizeof(site), size,
				 index, parity, dest );
}

/*
**  old style gather routine which declares and starts in one call
*/
msg_tag *
start_gather_site(
  field_offset field,	/* which field? Some member of structure "site" */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  msg_tag *mt;

  mt = declare_strided_gather( (char *)lattice + field, sizeof(site), size,
			       index, parity, dest );
  prepare_gather(mt);
  do_gather(mt);

  return mt;
}

/*
**  old style routine used to restart a previously waited gather
**  this finction is now depreciated and users should call do_gather()
**  instead
*/
void
restart_gather_site(
  field_offset field,	/* which field? Some member of structure "site" */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest,		/* one of the vectors of pointers */
  msg_tag *mtag)        /* previously returned by start_gather_site */
{
  msg_sr_t *mbuf;

  if(mtag->nsends!=0) mbuf = mtag->send_msgs;
  else mbuf = NULL;

  /* sanity checks for improper usage */
  if(mbuf!=NULL) {
    if(((char *)lattice+field)!=mbuf->gmem->mem) {
      printf("error: wrong field in restart gather\n");
      terminate(1);
    }
    if(sizeof(site)!=mbuf->gmem->stride) {
      printf("error: wrong stride in restart gather\n");
      terminate(1);
    }
    if(size!=mbuf->gmem->size) {
      printf("error: wrong size in restart gather\n");
      terminate(1);
    }
    if(((char *)lattice+field)!=mbuf->gmem->mem) {
      printf("error: wrong field in restart gather\n");
      terminate(1);
    }
  }

  do_gather(mtag);
}

/*****************************
 * gather routines from temp *
 *****************************/

/*
**  declares a gather from a temp array
*/
msg_tag *
declare_gather_from_temp(
  void * field,		/* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  return declare_strided_gather( field, size, size, index, parity, dest );
}

/*
**  old style gather routine which declares and starts in one call
*/
msg_tag *
start_gather_field(
  void * field,		/* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  msg_tag *mt;

  mt = declare_strided_gather( field, size, size, index, parity, dest );
  prepare_gather(mt);
  do_gather(mt);

  return mt;
}

/*
**  old style routine used to restart a previously waited gather
**  this finction is now depreciated and users should call do_gather()
**  instead
*/
void
restart_gather_field(
  void *field,		/* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest,		/* one of the vectors of pointers */
  msg_tag *mtag)          /* previously returned by start_gather_field */
{
  msg_sr_t *mbuf;

  if(mtag->nsends!=0) mbuf = mtag->send_msgs;
  else mbuf = NULL;

  /* sanity checks for improper usage */
  if(mbuf!=NULL) {
    if(field!=mbuf->gmem->mem) {
      printf("error: wrong field in restart gather\n");
      terminate(1);
    }
    if(size!=mbuf->gmem->stride) {
      printf("error: wrong stride in restart gather\n");
      terminate(1);
    }
    if(size!=mbuf->gmem->size) {
      printf("error: wrong size in restart gather\n");
      terminate(1);
    }
    if(field!=mbuf->gmem->mem) {
      printf("error: wrong field in restart gather\n");
      terminate(1);
    }
  }

  do_gather(mtag);
}


/**********************************************************************
 *                      MULTI-GATHER ROUTINES                         *
 **********************************************************************

 accumulate_gather(msg_tag **mtag, msg_tag *tag)
   Joins declared gathers together under a single msg_tag.
   The second argument (tag) would be merged with the first (mtag).
   If mtag is NULL then this just copies tag into mtag.

 declare_accumulate_gather() declares and joins gathers.

 example:

   msg_tag *tag1, *tag2, *mtag;

   tag1 = declare_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	                  EVEN, gen_pt1 );
   tag2 = declare_gather( F_OFFSET(phi), sizeof(su3_vector), XDOWN,
	                  EVEN, gen_pt2 );
   mtag = NULL;
   accumulate_gather( &mtag, tag1 );
   accumulate_gather( &mtag, tag2 );
   prepare_gather( mtag );  ** optional **
   do_gather( mtag );
   wait_gather( mtag );
   ** stuff **
   do_gather( tag1 );     ** this is valid as long as the combined gather
   wait_gather( tag1 );      (mtag) has been waited on **
   ** stuff **
   do_gather( mtag );
   wait_gather( mtag );
   cleanup_gather( mtag );
   cleanup_gather( tag1 );
   cleanup_gather( tag2 );

 Note that mtag must be set to NULL first in this case.
 If there is no need to use the single gathers alone one could do:

   msg_tag *mtag;

   mtag = NULL;
   declare_accumulate_gather( &mtag, F_OFFSET(phi), sizeof(su3_vector), XUP,
	                      EVEN, gen_pt1 );
   declare_accumulate_gather( &mtag, F_OFFSET(phi), sizeof(su3_vector), XDOWN,
	                      EVEN, gen_pt2 );
   prepare_gather( mtag );  ** optional **
   do_gather( mtag );
   wait_gather( mtag );
   ** stuff **
   do_gather( mtag );
   wait_gather( mtag );
   cleanup_gather( mtag );

 one coule also replace
   mtag = NULL;
   declare_accumulate_gather( &mtag, F_OFFSET(phi), sizeof(su3_vector), XUP,
	                      EVEN, gen_pt1 );
 with
   mtag = declare_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	                  EVEN, gen_pt1 );
 since they do the same thing, however the first form is a bit more uniform
 in the given example.
*/

/*
**  helper function to copy the gmem_t structure
*/
static void
copy_gmem(gmem_t **dest, gmem_t *src)
{
  while(*dest!=NULL) dest = &((*dest)->next);
  do {
    *dest = (gmem_t *)malloc(sizeof(gmem_t));
    if(*dest==NULL) {
      printf("error copy_gmem malloc node:%i\n",mynode());
      terminate(1);
    }
    memcpy(*dest, src, sizeof(gmem_t));
    dest = &((*dest)->next);
    src = src->next;
  } while(src!=NULL);
  *dest = NULL;
}

/*
**  helper function that merges a source msg_sr_t structure into the dest
*/
static void
add_msgt(msg_sr_t **dest, int *ndest, msg_sr_t *src, int nsrc, int nids)
{
  int i, j, n;

  n = 0;
  for(i=0; i<nsrc; ++i) {
    for(j=0; j<*ndest; ++j) {
      if((*dest)[j].msg_node==src[i].msg_node) {
	++n;
	break;
      }
    }
  }
  n = *ndest + nsrc - n;

  if(n!=0) {
    *dest = (msg_sr_t *)realloc(*dest, n*sizeof(msg_sr_t));
    if(*dest==NULL) {
      printf("error add_msgt malloc node:%i\n",mynode());
      terminate(1);
    }
    for(i=0; i<nsrc; ++i) {
      for(j=0; j<*ndest; ++j) {
	if((*dest)[j].msg_node==src[i].msg_node) break;
      }
      if(j<*ndest) {
	(*dest)[j].msg_size += src[i].msg_size;
	copy_gmem(&((*dest)[j].gmem), src[i].gmem);
      } else {
	(*dest)[*ndest+i].msg_node = src[i].msg_node;
	(*dest)[*ndest+i].id_offset = nids + src[i].id_offset;
	(*dest)[*ndest+i].msg_size = src[i].msg_size;
	(*dest)[*ndest+i].msg_buf = NULL;
	(*dest)[*ndest+i].gmem = NULL;
	copy_gmem(&((*dest)[*ndest+i].gmem), src[i].gmem);
      }
    }
  }
  *ndest = n;
}

/*
**  merges already declared gather
*/
void
accumulate_gather(msg_tag **mmtag, msg_tag *mtag)
{
  msg_tag *amtag;

  if(*mmtag==NULL) {
    amtag = (msg_tag *)malloc(sizeof(msg_tag));
    if(amtag==NULL) {
      printf("error accumulate_gather malloc node:%i\n",mynode());
      terminate(1);
    }
    amtag->nids = 0;
    amtag->ids = NULL;
    amtag->nrecvs = 0;
    amtag->recv_msgs = NULL;
    amtag->nsends = 0;
    amtag->send_msgs = NULL;
    *mmtag = amtag;
  } else {
    amtag = *mmtag;
  }

  add_msgt( &(amtag->recv_msgs), &(amtag->nrecvs),
	    mtag->recv_msgs, mtag->nrecvs, amtag->nids );
  add_msgt( &(amtag->send_msgs), &(amtag->nsends),
	    mtag->send_msgs, mtag->nsends, amtag->nids );
  amtag->nids += mtag->nids;
}

/*
**  declares and merges gather
**  handles both field offset and temp
*/
static void
declare_accumulate_strided_gather(
  msg_tag **mmtag,      /* tag to accumulate gather into */
  void *field,	        /* which field? Some member of structure "site" */
  int stride,           /* bytes between fields in source buffer */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  msg_tag *mtag;

  mtag = declare_strided_gather( field, stride, size, index, parity, dest );
  if(*mmtag==NULL) {
    *mmtag = mtag;
  } else {
    accumulate_gather( mmtag, mtag );
    cleanup_gather( mtag );
  }
}

/*
**  declares and merges gather from field offset
*/
void
declare_accumulate_gather(
  msg_tag **mmtag,
  field_offset field,	/* which field? Some member of structure "site" */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  declare_accumulate_strided_gather( mmtag, (char *)lattice + field,
				     sizeof(site), size, index, parity, dest );
}

/*
**  declares and merges gather from temp
*/
void
declare_accumulate_gather_from_temp(
  msg_tag **mmtag,
  void * field,		/* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  declare_accumulate_strided_gather( mmtag, field, size, size, index, parity,
				     dest );
}


/**********************************************************************
 *                     GENERAL GATHER ROUTINES                        *
 **********************************************************************

 start_general_gather() returns a msg_tag which will
   be used as input to subsequent wait_general_gather() and
   cleanup_general_gather() calls.

   usage: tag = start_general_gather( source, size, displacement, parity, dest)
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

struct msg_tmp { int node, count; }; /* temporary structure for keeping track
					of messages to be sent or received */
static struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
static int g_gather_flag=0; /* flag to tell if general gather in progress */
static int tsize;	    /* size of entry in messages =2*sizeof(int)+size */
static char ** tdest;	    /* tdest is copy of dest */
/* from_nodes, tsize and tdest are global because they are set in 
   start_general_gather() and used in wait_general_gather().  This
   works because we allow only one general_gather in progress at a
   time. */

#ifndef N_SUBL32

msg_tag *
start_general_strided_gather(
  char *field,	        /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int *displacement,	/* displacement to gather from. four components */
  int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  register int i,j;	/* scratch */
  register site *s;	/* scratch pointer to site */
  register char *tpt;	/* scratch pointer in buffers */
  int nsites;		/* number of sites in this receive or send */
  int disp_parity;	/* parity of displacement vector */
  int send_parity;	/* parity of sites that may be sent */
  int tx,ty,tz,tt;	/* temporary coordinates */
  int othernode;		/* node sent to or received from */
  msg_sr_t *mrecv,*msend;
  msg_tag *mtag; 		/* message tag, to be returned */
  int n_send_msgs, n_recv_msgs;

  /* check for gather already in progress */
  if(g_gather_flag!=0){
    printf("ERROR: node %d, two general_gathers() at once!\n", mynode());
    terminate(1);
  }
  n_recv_msgs = n_send_msgs = 0;
  tsize = 2*sizeof(int)+size;
  /* Use 2*sizeof int so pointer will be aligned to double word */
  tdest = dest;
  /* find parity of sites that may be sent */
  if( (displacement[XUP]+displacement[YUP]+displacement[ZUP]+
       displacement[TUP])%2 == 0 ) disp_parity = EVEN;
  else disp_parity = ODD;
  switch(parity) {
    case EVEN:
      if( disp_parity==EVEN ) send_parity = EVEN;
      else send_parity = ODD;
      break;
    case ODD:
      if( disp_parity==EVEN ) send_parity = ODD;
      else send_parity = EVEN;
      break;
    default: /* EVENANDODD */
      if(parity!=EVENANDODD) {
	printf("ERROR: bad parity\n");
	terminate(parity);
      }
      send_parity = EVENANDODD;
      break;
  }

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) Make
     list of nodes from whom we expect messages */
  FORSOMEPARITY(i,s,parity){
    if(displacement[XUP]!=0) tx = (s->x + displacement[XUP] + nx)%nx;
    else                     tx = s->x;
    if(displacement[YUP]!=0) ty = (s->y + displacement[YUP] + ny)%ny;
    else                     ty = s->y;
    if(displacement[ZUP]!=0) tz = (s->z + displacement[ZUP] + nz)%nz;
    else                     tz = s->z;
    if(displacement[TUP]!=0) tt = (s->t + displacement[TUP] + nt)%nt;
    else                     tt = s->t;
    othernode = node_number(tx,ty,tz,tt);
    if( othernode==this_node ) {
      dest[i] = field + node_index(tx,ty,tz,tt) * stride;
    }
    else{
      for(j=0;j<n_recv_msgs;j++) if(from_nodes[j].node==othernode) break;
      if(j < n_recv_msgs) {
	from_nodes[j].count++;
      }
      else {
	if(n_recv_msgs==0) {
	  from_nodes = (struct msg_tmp *)malloc( sizeof(struct msg_tmp) );
	  from_nodes[0].node = othernode;
	  from_nodes[0].count = 1;
	  n_recv_msgs++;
	}
	else{
	  from_nodes = (struct msg_tmp *)
	    realloc( from_nodes, (n_recv_msgs+1)*sizeof(struct msg_tmp) );
	  from_nodes[j].node = othernode;
	  from_nodes[j].count = 1;
	  n_recv_msgs++;
	}
      }
    }
  }

  /* scan sites of parity we are sending, make list of nodes to which
     we must send messages and the number of messages to each. */
  FORSOMEPARITY(i,s,send_parity) {
    if(displacement[XUP]!=0) tx = (s->x - displacement[XUP] + nx)%nx;
    else                     tx = s->x;
    if(displacement[YUP]!=0) ty = (s->y - displacement[YUP] + ny)%ny;
    else                     ty = s->y;
    if(displacement[ZUP]!=0) tz = (s->z - displacement[ZUP] + nz)%nz;
    else                     tz = s->z;
    if(displacement[TUP]!=0) tt = (s->t - displacement[TUP] + nt)%nt;
    else                     tt = s->t;
    othernode = node_number(tx,ty,tz,tt);
    if( othernode != this_node ) {
      for(j=0;j<n_send_msgs;j++) if(to_nodes[j].node==othernode) break;
      if(j < n_send_msgs) {
	to_nodes[j].count++;
      }
      else {
	if(n_send_msgs==0) {
	  to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
	  to_nodes[0].node = othernode;
	  to_nodes[0].count = 1;
	  n_send_msgs++;
	}
	else{
	  to_nodes = (struct msg_tmp *)
	    realloc( to_nodes, (n_send_msgs+1)*sizeof(struct msg_tmp) );
	  to_nodes[j].node = othernode;
	  to_nodes[j].count = 1;
	  n_send_msgs++;
	}
      }
    }
  }

  mtag = (msg_tag *)malloc(sizeof(msg_tag));

  if( n_recv_msgs==0 ) mrecv = NULL;
  else {
    mrecv = (msg_sr_t *)malloc(n_recv_msgs*sizeof(msg_sr_t));
    if(mrecv==NULL) {
      printf("NO ROOM for mrecv, node %d\n",mynode());
      terminate(1);
    }
  }
  if( n_send_msgs==0 ) msend = NULL;
  else {
    msend = (msg_sr_t *)malloc(n_send_msgs*sizeof(msg_sr_t));
    if(msend==NULL) {
      printf("NO ROOM for msend, node %d\n",mynode());
      terminate(1);
    }
  }
  mtag->recv_msgs = mrecv;
  mtag->send_msgs = msend;

  mtag->nrecvs = n_recv_msgs;
  mtag->nsends = n_send_msgs;

  /* for each node which has neighbors of my sites */
  for(i=0; i<n_recv_msgs; i++) {
    /* allocate buffer to receive neighbors */
    nsites = from_nodes[i].count;
    mrecv[i].msg_node = from_nodes[i].node;
    mrecv[i].msg_size = nsites*tsize;
    mrecv[i].msg_buf = (char *)malloc( nsites*tsize );
    if(mrecv[i].msg_buf==NULL){
      printf("NO ROOM for msg_buf, node %d\n",mynode());
      terminate(1);
    }
    /* post receive */
    MPI_Irecv( mrecv[i].msg_buf, nsites*tsize, MPI_BYTE,
	       from_nodes[i].node, GENERAL_GATHER_ID, 
	       MPI_COMM_WORLD, &mrecv[i].msg_req );
  }

  /* for each node whose neighbors I have */
  for(i=0; i<n_send_msgs; i++) {
    /* Allocate buffer to gather data. */
    tpt=(char *)malloc( to_nodes[i].count*tsize );
    if(tpt==NULL) {
      printf("NO ROOM for tpt, node %d\n",mynode());
      terminate(1);
    }
    msend[i].msg_node = to_nodes[i].node;
    msend[i].msg_size = to_nodes[i].count*tsize;
    msend[i].msg_buf = tpt;
  }

  /* reset to_node counters */
  for(i=0; i<n_send_msgs; i++) to_nodes[i].count = 0;
  /* gather data into the buffers. Each entry in the buffers consists
     of the index of the site to which the data is sent, followed by
     the actual data */
  FORSOMEPARITY(i, s, send_parity) {
    tx = (s->x - displacement[XUP] + nx)%nx;
    ty = (s->y - displacement[YUP] + ny)%ny;
    tz = (s->z - displacement[ZUP] + nz)%nz;
    tt = (s->t - displacement[TUP] + nt)%nt;
    othernode = node_number(tx,ty,tz,tt);
    if( othernode != this_node ) {
      for(j=0; j<n_send_msgs; j++) if(to_nodes[j].node==othernode) break;
      tpt = msend[j].msg_buf + to_nodes[j].count*tsize;
      *(int *)tpt = node_index(tx,ty,tz,tt);
      /* index of site on other node */
      memcpy( tpt+2*sizeof(int), field+i*stride, size);
      to_nodes[j].count++;
    }
  }

  /* start the sends */
  for(i=0; i<n_send_msgs; i++) {
    nsites = to_nodes[i].count;
    MPI_Isend( msend[i].msg_buf, nsites*tsize, MPI_BYTE,
	       to_nodes[i].node, GENERAL_GATHER_ID, 
	       MPI_COMM_WORLD, &msend[i].msg_req );
  }

  /* free temporary arrays */
  if(n_send_msgs>0) free(to_nodes);
  /* mark gather in progress and return */
  g_gather_flag = 1;

  return mtag;
}

#else	/* N_SUBL32 */

msg_tag *
start_general_strided_gather(
  char *field,	        /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int *displacement,	/* displacement to gather from. four components */
  int subl,		/* subl of sites whose neighbors we gather.
			   It is EVENANDODD, if all sublattices are done. */
  char ** dest)		/* one of the vectors of pointers */
{
  register int i,j;	/* scratch */
  register site *s;	/* scratch pointer to site */
  register char *tpt;	/* scratch pointer in buffers */
  int nsites;		/* number of sites in this receive or send */
  int send_subl;		/* sublattice of sites that may be sent */
  int tx,ty,tz,tt;	/* temporary coordinates */
  int othernode;		/* node sent to or received from */
  msg_sr_t *mrecv,*msend;
  msg_tag *mtag;		/* message tag, to be returned */
  int n_send_msgs, n_recv_msgs;

  /* check for gather already in progress */
  if(g_gather_flag!=0) {
    printf("ERROR: node %d, two general_gathers() at once!\n", mynode());
    terminate(1);
  }
  n_recv_msgs = n_send_msgs = 0;
  tsize = 2*sizeof(int)+size;
  /* Use 2*sizeof int so pointer will be aligned to double word */
  tdest = dest;

  /* find sublattice of sites that may be sent */
  /* This is not needed for EVENANDODD */
  send_subl = subl;
  if( subl != EVENANDODD ) {
    /* Displacments by multiples of 4 in any direction does
       not change sublattice */
    tx = displacement[XUP]%4;
    ty = displacement[YUP]%4;
    tz = displacement[ZUP]%4;
    tt = displacement[TUP]%4;
    if( tx < 0 ) {
      for(i=0;i<(-tx);i++) send_subl = neighsubl[send_subl][XDOWN];
    }
    else
      for(i=0;i<tx;i++) send_subl = neighsubl[send_subl][XUP];
    if( ty < 0 ) {
      for(i=0;i<(-ty);i++) send_subl = neighsubl[send_subl][YDOWN];
    }
    else
      for(i=0;i<ty;i++) send_subl = neighsubl[send_subl][YUP];
    if( tz < 0 ) {
      for(i=0;i<(-tz);i++) send_subl = neighsubl[send_subl][ZDOWN];
    }
    else
      for(i=0;i<tz;i++) send_subl = neighsubl[send_subl][ZUP];
    if( tt < 0 ) {
      for(i=0;i<(-tt);i++) send_subl = neighsubl[send_subl][TDOWN];
    }
    else
      for(i=0;i<tt;i++) send_subl = neighsubl[send_subl][TUP];
  }

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) Make
     list of nodes from whom we expect messages */
  if( subl == EVENANDODD ) {
    FORALLSITES(i,s) {
      if(displacement[XUP]!=0) tx = (s->x + displacement[XUP] + nx)%nx;
      else                     tx = s->x;
      if(displacement[YUP]!=0) ty = (s->y + displacement[YUP] + ny)%ny;
      else                     ty = s->y;
      if(displacement[ZUP]!=0) tz = (s->z + displacement[ZUP] + nz)%nz;
      else                     tz = s->z;
      if(displacement[TUP]!=0) tt = (s->t + displacement[TUP] + nt)%nt;
      else                     tt = s->t;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode==this_node ) {
	dest[i] = field + node_index(tx,ty,tz,tt) * stride;
      }
      else{
	for(j=0;j<n_recv_msgs;j++) if(from_nodes[j].node==othernode) break;
	if(j < n_recv_msgs) {
	  from_nodes[j].count++;
	}
	else {
	  if(n_recv_msgs==0) {
	    from_nodes = (struct msg_tmp *)malloc( sizeof(struct msg_tmp) );
	    from_nodes[0].node = othernode;
	    from_nodes[0].count = 1;
	    n_recv_msgs++;
	  }
	  else{
	    from_nodes = (struct msg_tmp *)
	      realloc( from_nodes, (n_recv_msgs+1)*sizeof(struct msg_tmp) );
	    from_nodes[j].node = othernode;
	    from_nodes[j].count = 1;
	    n_recv_msgs++;
	  }
	}
      }
    }
  }
  else {
    FORSOMESUBLATTICE(i,s,subl) {
      if(displacement[XUP]!=0) tx = (s->x + displacement[XUP] + nx)%nx;
      else                     tx = s->x;
      if(displacement[YUP]!=0) ty = (s->y + displacement[YUP] + ny)%ny;
      else                     ty = s->y;
      if(displacement[ZUP]!=0) tz = (s->z + displacement[ZUP] + nz)%nz;
      else                     tz = s->z;
      if(displacement[TUP]!=0) tt = (s->t + displacement[TUP] + nt)%nt;
      else                     tt = s->t;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode==this_node ) {
	dest[i] = field + node_index(tx,ty,tz,tt) * stride;
      }
      else {
	for(j=0;j<n_recv_msgs;j++) if(from_nodes[j].node==othernode) break;
	if(j < n_recv_msgs) {
	  from_nodes[j].count++;
	}
	else {
	  if(n_recv_msgs==0) {
	    from_nodes = (struct msg_tmp *)malloc( sizeof(struct msg_tmp) );
	    from_nodes[0].node = othernode;
	    from_nodes[0].count = 1;
	    n_recv_msgs++;
	  }
	  else{
	    from_nodes = (struct msg_tmp *)
	      realloc( from_nodes, (n_recv_msgs+1)*sizeof(struct msg_tmp) );
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
  if( subl == EVENANDODD ) {
    FORALLSITES(i,s) {
      if(displacement[XUP]!=0) tx = (s->x - displacement[XUP] + nx)%nx;
      else                     tx = s->x;
      if(displacement[YUP]!=0) ty = (s->y - displacement[YUP] + ny)%ny;
      else                     ty = s->y;
      if(displacement[ZUP]!=0) tz = (s->z - displacement[ZUP] + nz)%nz;
      else                     tz = s->z;
      if(displacement[TUP]!=0) tt = (s->t - displacement[TUP] + nt)%nt;
      else                     tt = s->t;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode != this_node ) {
	for(j=0;j<n_send_msgs;j++) if(to_nodes[j].node==othernode) break;
	if(j < n_send_msgs) {
	  to_nodes[j].count++;
	}
	else {
	  if(n_send_msgs==0) {
	    to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
	    to_nodes[0].node = othernode;
	    to_nodes[0].count = 1;
	    n_send_msgs++;
	  }
	  else {
	    to_nodes = (struct msg_tmp *)
	      realloc( to_nodes, (n_send_msgs+1)*sizeof(struct msg_tmp) );
	    to_nodes[j].node = othernode;
	    to_nodes[j].count = 1;
	    n_send_msgs++;
	  }
	}
      }
    }
  }
  else {
    FORSOMESUBLATTICE(i,s,send_subl) {
      if(displacement[XUP]!=0) tx = (s->x - displacement[XUP] + nx)%nx;
      else                     tx = s->x;
      if(displacement[YUP]!=0) ty = (s->y - displacement[YUP] + ny)%ny;
      else                     ty = s->y;
      if(displacement[ZUP]!=0) tz = (s->z - displacement[ZUP] + nz)%nz;
      else                     tz = s->z;
      if(displacement[TUP]!=0) tt = (s->t - displacement[TUP] + nt)%nt;
      else                     tt = s->t;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode != this_node ) {
	for(j=0;j<n_send_msgs;j++) if(to_nodes[j].node==othernode) break;
	if(j < n_send_msgs) {
	  to_nodes[j].count++;
	}
	else {
	  if(n_send_msgs==0) {
	    to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
	    to_nodes[0].node = othernode;
	    to_nodes[0].count = 1;
	    n_send_msgs++;
	  }
	  else {
	    to_nodes = (struct msg_tmp *)
	      realloc( to_nodes, (n_send_msgs+1)*sizeof(struct msg_tmp) );
	    to_nodes[j].node = othernode;
	    to_nodes[j].count = 1;
	    n_send_msgs++;
	  }
	}
      }
    }
  }

  mtag = (msg_tag *)malloc(sizeof(msg_tag));

  if( n_recv_msgs==0 ) mrecv = NULL;
  else {
    mrecv = (msg_sr_t *)malloc( n_recv_msgs*sizeof(msg_sr_t) );
    if(mrecv==NULL) {
      printf("NO ROOM for mrecv, node %d\n",mynode());
      terminate(1);
    }
  }
  if( n_send_msgs==0 ) msend=NULL;
  else {
    msend = (msg_sr_t *)malloc( n_send_msgs*sizeof(msg_sr_t) );
    if(msend==NULL) {
      printf("NO ROOM for msend, node %d\n",mynode());
      terminate(1);
    }
  }

  mtag->recv_msgs = mrecv;
  mtag->send_msgs = msend;

  mtag->nrecvs = n_recv_msgs;
  mtag->nsends = n_send_msgs;

  /* for each node which has neighbors of my sites */
  for(i=0; i<n_recv_msgs; i++) {
    /* allocate buffer to receive neighbors */
    nsites = from_nodes[i].count;
    mrecv[i].msg_node = from_nodes[i].node;
    mrecv[i].msg_size = nsites*tsize;
    mrecv[i].msg_buf = (char *)malloc( nsites*tsize );
    if(mrecv[i].msg_buf==NULL){
      printf("NO ROOM for msg_buf, node %d\n",mynode());
      terminate(1);
    }
    /* post receive */
    MPI_Irecv( mrecv[i].msg_buf, nsites*tsize, MPI_BYTE,
	       from_nodes[i].node, GENERAL_GATHER_ID,
	       MPI_COMM_WORLD, &mrecv[i].msg_req );
  }

  /* for each node whose neighbors I have */
  for(i=0; i<n_send_msgs; i++) {
    /* Allocate buffer to gather data. */
    tpt = (char *)malloc( to_nodes[i].count*tsize );
    if(tpt==NULL) {
      printf("NO ROOM for tpt, node %d\n",mynode());
      terminate(1);
    }
    msend[i].msg_node = to_nodes[i].node;
    msend[i].msg_size = to_nodes[i].count*tsize;
    msend[i].msg_buf = tpt;
  }

  /* reset to_node counters */
  for(i=0; i<n_send_msgs; i++) to_nodes[i].count = 0;
  /* gather data into the buffers. Each entry in the buffers consists
     of the index of the site to which the data is sent, followed by
     the actual data */
  if( subl == EVENANDODD ) {
    FORALLSITES(i, s) {
      tx = (s->x - displacement[XUP] + nx)%nx;
      ty = (s->y - displacement[YUP] + ny)%ny;
      tz = (s->z - displacement[ZUP] + nz)%nz;
      tt = (s->t - displacement[TUP] + nt)%nt;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode != this_node ) {
	for(j=0; j<n_send_msgs; j++) if(to_nodes[j].node==othernode) break;
	tpt = msend[j].msg_buf + to_nodes[j].count*tsize;
	*(int *)tpt = node_index(tx,ty,tz,tt);
	/* index of site on other node */
	memcpy( tpt+2*sizeof(int), field+i*stride, size);
	to_nodes[j].count++;
      }
    }
  }
  else {
    FORSOMESUBLATTICE(i, s, send_subl) {
      tx = (s->x - displacement[XUP] + nx)%nx;
      ty = (s->y - displacement[YUP] + ny)%ny;
      tz = (s->z - displacement[ZUP] + nz)%nz;
      tt = (s->t - displacement[TUP] + nt)%nt;
      othernode = node_number(tx,ty,tz,tt);
      if( othernode != this_node ) {
	for(j=0; j<n_send_msgs; j++) if(to_nodes[j].node==othernode) break;
	tpt = msend[j].msg_buf + to_nodes[j].count*tsize;
	*(int *)tpt = node_index(tx,ty,tz,tt);
	/* index of site on other node */
	memcpy( tpt+2*sizeof(int), field+i*stride, size);
	to_nodes[j].count++;
      }
    }
  }

  /* start the sends */
  for(i=0; i<n_send_msgs; i++) {
    nsites = to_nodes[i].count;
    MPI_Isend( msend[i].msg_buf, nsites*tsize, MPI_BYTE,
	       to_nodes[i].node, GENERAL_GATHER_ID,
	       MPI_COMM_WORLD, &msend[i].msg_req );
  }

  /* free temporary arrays */
  if( n_send_msgs > 0) free(to_nodes);
  /* mark gather in progress and return */
  g_gather_flag = 1;

  return mtag;
}

#endif	/* N_SUBL32 */

msg_tag *
start_general_gather(
  field_offset field,	/* which field? Some member of structure "site" */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int *displacement,	/* displacement to gather from. four components */
  int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  return start_general_strided_gather( (char *)lattice + field, sizeof(site),
				       size, displacement, parity, dest );
}

msg_tag *
start_general_gather_from_temp(
  void * field,	        /* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int *displacement,	/* displacement to gather from. four components */
  int parity,		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest)		/* one of the vectors of pointers */
{
  return start_general_strided_gather( field, size, size,
				       displacement, parity, dest );
}

/*
**  wait for a general gather to complete
*/
void
wait_general_gather(msg_tag *mtag)
{
  int i,j,k;
  MPI_Status status;

  g_gather_flag=0;

  for(i=0; i<mtag->nrecvs; i++) {
    MPI_Wait( &mtag->recv_msgs[i].msg_req, &status );
    /* set pointers in sites to correct location */
    for(j=0; j<from_nodes[i].count; j++) {
      /* k = index of site on this node, sent in message */
      k = *(int *)( mtag->recv_msgs[i].msg_buf + j*tsize );
      tdest[k] = mtag->recv_msgs[i].msg_buf + j*tsize + 2*sizeof(int);
    }
  }
  if(i>0) free(from_nodes);
}

/*
**  free memory associated with general gather
*/
void
cleanup_general_gather(msg_tag *mtag)
{
  int i;
  MPI_Status status;

  /* free all receive buffers */
  for(i=0; i<mtag->nrecvs; i++) {
    free( mtag->recv_msgs[i].msg_buf );
  }
  /* wait for all send messages, free all send buffers */
  for(i=0; i<mtag->nsends; i++) {
    MPI_Wait( &mtag->send_msgs[i].msg_req, &status );
    free( mtag->send_msgs[i].msg_buf );
  }
  /* free the msg_tag buffer */
  free(mtag->recv_msgs);
  free(mtag->send_msgs);
  free(mtag);
}
