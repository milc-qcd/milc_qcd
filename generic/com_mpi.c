/******************  com_mpi.c *****************************************/
/* Communications routines for the SU3 program
   MIMD version 7.
   This file is communications-scheme dependent.
   MPI version - allegedly machine independent
   This version breaks the MPI machine into a number of separate lattices
*/
/* Modifications

    4/20/02 added start_general_gather_field C.D.
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
   myjobid()             returns the jobid of this node.
   numjobs()             returns number of jobs
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
   get_utc_datetime()    get GM time as ASCII string
   sort_eight_gathers()  sorts eight contiguous gathers from order
                           XUP,XDOWN,YUP,YDOWN... to XUP,YUP,...XDOWN,YDOWN...
   make_nn_gathers()     makes all necessary lists for communications with
                           nodes containing neighbor sites.
   make_gather()         calculates and stores necessary communications lists
                           for a given gather mapping
   declare_gather_site()      creates a message tag that defines specific details
                           of a gather to be used later
   declare_gather_field()  creates a message tag that defines specific
                               details of a gather from field to be used later
   prepare_gather()      optional call that allocates buffers for a previously
                           declared gather.  will automatically be called from
                           do_gather() if not done before.
   do_gather()           executes a previously declared gather
   wait_gather()         waits for receives to finish, insuring that the
                           data has actually arrived.
   cleanup_gather()      frees all the buffers that were allocated, WHICH
                           MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   accumulate_gather()   combines gathers into single message tag
   declare_accumulate_gather_site()  does declare_gather_site() and 
                                  accumulate_gather() in single step.
   declare_accumulate_gather_field()  does declare_gather_field() and
                                            accumulate_gather() in single step.
   start_gather_site()        older function which does declare/prepare/do_gather
                           in a single step
   start_gather_field()  older function which does
                               declare/prepare/do_gather_field
   restart_gather_site()      older function which is obsoleted by do_gather()
   restart_gather_field() older function which is obsoleted by do_gather() 
   start_general_gather_site()  starts asynchronous sends and receives required
                             to gather fields at arbitrary displacement.
   start_general_gather_field() starts asynchronous sends and receives 
                             required to gather neighbors from an
			     array of fields.
   wait_general_gather()   waits for receives to finish, insuring that the
                             data has actually arrived, and sets pointers to
			     received data.
   cleanup_general_gather()  frees all the buffers that were allocated, WHICH
                               MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   myjobid()                 The index number of this job
   numjobs()                 Number of jobs in multijob execution
   jobgeom()                 Dimensions of the multijob layout.  Product = numjobs
   ionodegeom()              Dimensions of the I/O partition layout.  Product =
                              number of files.
   nodegeom()                Allocated dimensions of the nodes.

*/

#include <time.h>
#include "generic_includes.h"
#include <mpi.h>
#include <ctype.h>
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

/* If we want to do our own checksums */
#ifdef COM_CRC
u_int32type crc32(u_int32type crc, const unsigned char *buf, size_t len);
#define CRCBYTES 8
#else
#define CRCBYTES 0
#endif

/* hacks needed to unify even/odd and 32 sublattice cases */
#ifdef N_SUBL32
#define NUM_SUBL 32
#define FORSOMEPARITY FORSOMESUBLATTICE
#else
#define NUM_SUBL 2
#endif

static int jobid = 0;
static int num_jobs = 1;
static int *geom = NULL;
static int *jobgeomvals = NULL;
static int *worldcoord = NULL;
static MPI_Comm  MPI_COMM_THISJOB;


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
#ifdef CRC_DEBUG
  int index;
#endif
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

static void
get_arg(int argc, char **argv, char *tag, int *first, int *last,
	char **c, int **a)
{
  int i;
  *first = -1;
  *last = -1;
  *c = NULL;
  *a = NULL;
  for(i=1; i<argc; i++) {
    if(strcmp(argv[i], tag)==0) {
      *first = i;
      //printf("%i %i\n", i, argc);
      if( ((i+1)<argc) && !(isdigit(argv[i+1][0])) ) {
	//printf("c %i %s\n", i+1, argv[i+1]);
	*c = argv[i+1];
	*last = i+1;
      } else {
	//printf("a %i %s\n", i+1, argv[i+1]);
	while( (++i<argc) && isdigit(argv[i][0]) );
	*last = i-1;
	int n = *last - *first;
	if(n) {
	  int j;
	  *a = (int *) malloc(n*sizeof(int));
	  //printf("%i %p\n", n, *a);
	  for(j=0; j<n; j++) {
	    (*a)[j] = atoi(argv[*first+1+j]);
	    //printf(" %i", (*a)[j]);
	  }
	  //printf("\n");
	}
      }
    }
  }
}

static void
remove_from_args(int *argc, char ***argv, int first, int last)
{
  int n = last - first;
  if(first>=0) {
    int i;
    for(i=last+1; i<*argc; i++) (*argv)[i-n-1] = (*argv)[i];
    *argc -= n + 1;
  }
}

static int 
lex_rank(const int coords[], int dim, int size[])
{
  int d;
  int rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/* Create partitions of equal size from the allocated machine, based
   on num_jobs */
static void
repartition_switch_machine(void){
  int localnodeid;
  int num_nodes = numnodes();
  int nodeid = mynode();
  MPI_Comm jobcomm;
  int localgeom;
  int flag;

  /* localgeom gives the number of nodes in the job partition */
  if(num_nodes % num_jobs != 0){
    printf("num_jobs %i must divide number of nodes %i\n",
	   num_jobs, num_nodes);
    terminate(1);
  }
  localgeom = num_nodes/num_jobs;
  jobid = nodeid/localgeom;

  /* Split the communicator */

  flag = MPI_Comm_split(MPI_COMM_THISJOB, jobid, 0, &jobcomm);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);

  flag = MPI_Comm_rank(jobcomm, &localnodeid);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);
  
  /* Make MPI on this node think I live in just this one job partition */
  
  flag = MPI_Comm_free(&MPI_COMM_THISJOB);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);

  MPI_COMM_THISJOB = jobcomm;
}

/* Create partitions of equal size from the allocated machine, based
   on jobgeom */
static void
repartition_mesh_machine(void){
  int i;
  int localnodeid;
  int nd = 4;
  int flag;
  MPI_Comm jobcomm;
  int *jobcoord, *localgeom, *localcoord;

  if(jobgeomvals == NULL)return;

  /* localgeom gives the node dimensions of the job partition */
  localgeom = (int *)malloc(sizeof(int)*nd);
  for(i=0; i<nd; i++){
    if(geom[i] % jobgeomvals[i] != 0){
      printf( "job partition[%i] = %i must divide machine geometry %i\n",
	      i, jobgeomvals[i], geom[i]);fflush(stdout);
      terminate(1);
    }
    localgeom[i] = geom[i]/jobgeomvals[i];
  }

  /* jobcoord locates my job partition in the world of job partitions */
  /* localcoord locates my node within the job partition */
  jobcoord = (int *)malloc(sizeof(int)*nd);
  localcoord = (int *)malloc(sizeof(int)*nd);

  for(i=0; i<nd; i++){
    localcoord[i] = worldcoord[i]%localgeom[i];
    jobcoord[i]   = worldcoord[i]/localgeom[i];
  }

  jobid = lex_rank(jobcoord, nd, jobgeomvals);

  /* Split the communicator */

  flag = MPI_Comm_split(MPI_COMM_THISJOB, jobid, 0, &jobcomm);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);
  flag = MPI_Comm_rank(jobcomm, &localnodeid);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);

  //printf("node %d jobid %d\n", localnodeid, jobid); fflush(stdout);

  /* Make MPI on this node think I live in just this one job partition */

  flag = MPI_Comm_free(&MPI_COMM_THISJOB);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);
  MPI_COMM_THISJOB = jobcomm;
  for(i=0; i<nd; i++)
    geom[i]  = localgeom[i];

  free(localcoord);
  free(jobcoord);
  free(localgeom);
}

/*
**  Machine initialization
**  This version breaks the MPI machine into a number of
**  separate lattices.
*/
void
initialize_machine(int *argc, char ***argv)
{
  int i, n, flag, found, *tag_ub;
  int nj, nd;
  int first, last, *a = NULL;
  char *c = NULL;
  char myname[] = "initialize_machine";
  
  MPI_Comm comm;
  MPI_Errhandler errhandler;

  flag = MPI_Init(argc, argv);
  flag = MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_THISJOB);
  comm = MPI_COMM_THISJOB;
  if(flag != MPI_SUCCESS) err_func(&comm, &flag);

  /* check if 32 bit int is set correctly */
#ifdef SHORT_IS_32BIT
  if(sizeof(unsigned short)!=4) {
    printf("node %d: SHORT_IS_32BIT is set but sizeof(unsigned short)=%d\n",
	   mynode(), sizeof(unsigned short));
    terminate(1);
  }
#else
  if(sizeof(unsigned int)!=4) {
    printf("node %d: SHORT_IS_32BIT is not set but sizeof(unsigned int)=%d\n",
	   mynode(), (int)sizeof(unsigned int));
    terminate(1);
  }
#endif

  /* Process and remove our command-line arguments */

  /* process -geom */
  get_arg(*argc, *argv, "-geom", &first, &last, &c, &a);
  if( c != 0){
    node0_printf("%s: unknown argument to -geom: %s\n",myname,c);
    terminate(1);
  }
  nd = last - first;
  if(nd <= 0){
    geom = NULL;
  } else {
    if (nd != 4){
      node0_printf("%s: found %d -geom values, but wanted 4\n",myname, nd);
      terminate(1);
    }
    
    geom = (int *)malloc(4*sizeof(int));
    for(i = 0; i < 4; i++)
      geom[i] = a[i];

    worldcoord = (int *)malloc(4*sizeof(int));
    n = mynode();
    for(i=0; i<4; i++) {
      worldcoord[i] = n % geom[i];
      n /= geom[i];
    }
  }
    
  remove_from_args(argc, argv, first, last);

  /* process -jobs */
  /* This option causes the allocated machine to be subdivided into independent
     partitions in which separate jobs run with the same executable.
     -geom must first be specified. */

  /* The integer a[i] specifies the number of divisions of the ith geom dimension */
  /* The default a[i] = 1 for all i implies no subdivision */

  /* default settings */
  jobid = 0;
  num_jobs = 1;

  get_arg(*argc, *argv, "-jobs", &first, &last, &c, &a);
  if( c ) {
    printf("%s: unknown argument to -jobs: %s\n", myname, c);
    terminate(1);
  }
  nj = last - first;
  if(nj) {
    int i;
    jobgeomvals = a;
    /* Check sanity of job partition divisions */
    if(nj != 1 && geom == NULL){
      fprintf(stderr, "-jobs requires -geom\n");
      terminate(1);
    }
    if(geom != NULL && nj!=4) {
      printf("%s: allocated number dimensions %d != job partition dimensions %d\n", myname, 4, nj);
      terminate(1);
    }
    for(i=0; i<nj; i++){
      if(jobgeomvals[i]<=0){
	printf("%s: job partition division[%i] = %d <= 0\n", myname,
	       i, jobgeomvals[i]);
      }
      num_jobs *= jobgeomvals[i];
    }
    
    if(nj==1)
      repartition_switch_machine();
    else
      repartition_mesh_machine();
  }

  remove_from_args(argc, argv, first, last);

  /* process -ionodes a[0] a[1] a[2] a[3] flag */
  get_arg(*argc, *argv, "-ionodes", &first, &last, &c, &a);
  if(last - first > 0){
    if(mynode()==0)printf("-ionodes option requires QIO\n");
  }
  remove_from_args(argc, argv, first, last);


  /* Set error handler for this job */

  /* Note: with MPI-2 MPI_Comm_create_errhandler and
     MPI_Comm_set_errorhandler are preferred, but we keep MPI_Attr_get
     until MPI-2 is more widely available */ 
  flag = MPI_Errhandler_create(err_func, &errhandler);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);
  flag = MPI_Errhandler_set(MPI_COMM_THISJOB, errhandler);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);

  /* get the number of message types */
  /* Note: with MPI-2 MPI_Comm_get_attr is preferred,
     but we keep MPI_Attr_get until MPI-2 is more widely available */ 
  flag = MPI_Attr_get(MPI_COMM_THISJOB, MPI_TAG_UB, &tag_ub, &found);
  if(flag != MPI_SUCCESS) err_func(&MPI_COMM_THISJOB, &flag);
  if(found == 0){
    num_gather_ids = 1024;
    if(mynode() == 0){
      printf("%s: MPI won't give me an upper limit on the number of message types\n", 
	     myname);
      printf("%s: setting the limit to %d\n",myname, num_gather_ids);
    
    } 
  } else {
    num_gather_ids = *tag_ub + 1 - GATHER_BASE_ID;
  }
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
  // g_sync();
  MPI_Barrier( MPI_COMM_WORLD );  // wait for all lattices to finish?
  MPI_Finalize();
  fflush(stdout);
  exit(status);
}

/*
** UTC time as ASCII string
*/

void 
get_utc_datetime(char *time_string)
{
  time_t time_stamp;
  struct tm *gmtime_stamp;

  time(&time_stamp);
  gmtime_stamp = gmtime(&time_stamp);
  strncpy(time_string,asctime(gmtime_stamp),64);
  
  /* Remove trailing end-of-line character */
  if(time_string[strlen(time_string) - 1] == '\n')
    time_string[strlen(time_string) - 1] = '\0';
}


/*
**  version of exit for multinode processes -- kill all nodes
*/
// MPI_Abort has implementation dependent effects - sometimes kills
// all processes in MPI_COMM_WORLD, sometimes not.   Also, we can't
// decide if this is desirable.  So, for the moment, terminate()
// really works the same as normal_exit()
void
terminate(int status)
{
  time_stamp("termination");
  printf("Termination: node %d, status = %d\n", this_node, status);
  fflush(stdout);
  MPI_Barrier( MPI_COMM_WORLD );  // wait for all jobs to finish?
  //MPI_Abort(MPI_COMM_THISLATTICE, 0);
  MPI_Finalize();
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
  MPI_Comm_rank( MPI_COMM_THISJOB, &node );
  return(node);
}

/*
**  Return number of nodes
*/
int
numnodes(void)
{
  int nodes;
  MPI_Comm_size( MPI_COMM_THISJOB, &nodes );
  return(nodes);
}

/*
** Return the allocated dimensions (node geometry) if a grid is being used
*/
int const *
nodegeom(void)
{
  return geom;
}

/*
**  Return my jobid
*/
int
myjobid(void)
{
  return jobid;
}

/*
**  Return number of jobs
*/
int
numjobs(void)
{
  return num_jobs;
}

/*
** Return the job geometry
*/
int const *
jobgeom(void)
{
  return jobgeomvals;
}


/*
** Return the ionode geometry (supported only for QIO/QMP)
*/
int *
ionodegeom(void)
{
  return NULL;
}

/*
**  Synchronize all nodes
*/
void
g_sync(void)
{
  MPI_Barrier( MPI_COMM_THISJOB );
}

/*
**  Sum signed integer over all nodes
*/
void
g_intsum(int *ipt)
{
  int work;
  MPI_Allreduce( ipt, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_THISJOB );
  *ipt = work;
}

/*
**  Sum unsigned 32-bit integer type
*/
void
g_uint32sum(u_int32type *pt)
{
  u_int32type work;
#ifdef SHORT_IS_32BIT
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED_SHORT, 
		 MPI_SUM, MPI_COMM_THISJOB );
#else
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED, 
		 MPI_SUM, MPI_COMM_THISJOB );
#endif
  *pt = work;
}


/*
**  Sum generic floating type over all nodes
*/
void
g_floatsum(Real *fpt)
{
  Real work;
  MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( fpt, work, length, MILC_MPI_REAL, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( dpt, work, ndoubles, MPI_DOUBLE, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( cpt, &work, 2, MILC_MPI_REAL, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( cpt, work, 2*ncomplex, MILC_MPI_REAL, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( cpt, &work, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_THISJOB );
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
  MPI_Allreduce( cpt, work, 2*ncomplex, MPI_DOUBLE, MPI_SUM, MPI_COMM_THISJOB );
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
		 MPI_BXOR, MPI_COMM_THISJOB );
#else
  MPI_Allreduce( pt, &work, 1, MPI_UNSIGNED, 
		 MPI_BXOR, MPI_COMM_THISJOB );
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
  MPI_Allreduce( fpt, &work, 1, MILC_MPI_REAL, MPI_MAX, MPI_COMM_THISJOB );
  *fpt = work;
}

/*
**  Find maximum of double over all nodes
*/
void
g_doublemax(double *dpt)
{
  double work;
  MPI_Allreduce( dpt, &work, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_THISJOB );
  *dpt = work;
}

/*
**  Broadcast generic precision floating point number from node zero
*/
void
broadcast_float(Real *fpt)
{
  MPI_Bcast( fpt, 1, MILC_MPI_REAL, 0, MPI_COMM_THISJOB );
}

/*
**  Broadcast double precision floating point number from node zero
*/
void
broadcast_double(double *dpt)
{
  MPI_Bcast( dpt, 1, MPI_DOUBLE, 0, MPI_COMM_THISJOB );
}

/*
**  Broadcast generic precision complex number from node zero
*/
void
broadcast_complex(complex *cpt)
{
  MPI_Bcast( cpt, 2, MILC_MPI_REAL, 0, MPI_COMM_THISJOB );
}

/*
**  Broadcast double precision complex number from node zero
*/
void
broadcast_dcomplex(double_complex *cpt)
{
  MPI_Bcast( cpt, 2, MPI_DOUBLE, 0, MPI_COMM_THISJOB );
}

/*
**  Broadcast bytes from node 0 to all others
*/
void
broadcast_bytes(char *buf, int size)
{
  MPI_Bcast( buf, size, MPI_BYTE, 0, MPI_COMM_THISJOB );
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
  MPI_Send( address, 1, MPI_INT, tonode, SEND_INTEGER_ID, MPI_COMM_THISJOB );
}

/*
**  Receive an integer from another node
*/
void
receive_integer(int fromnode, int *address)
{
  MPI_Status status;
  MPI_Recv( address, 1, MPI_INT, fromnode, SEND_INTEGER_ID,
	    MPI_COMM_THISJOB, &status );
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
  MPI_Send( buf, size, MPI_BYTE, tonode, SEND_FIELD_ID, MPI_COMM_THISJOB );
}

/*
**  get_field is to be called only by the node to which the field was sent
*/
void
get_field(char *buf, int size, int fromnode)
{
  MPI_Status status;
  MPI_Recv( buf, size, MPI_BYTE, fromnode, SEND_FIELD_ID, MPI_COMM_THISJOB,
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
    MPI_Isend( &buf[i], 1, MPI_INT, recv->othernode, 0, MPI_COMM_THISJOB,
	       &req[i] );
  }
  if(i!=n_recv) {printf("error i!=n_recv\n"); terminate(1);}

  tol_next = &tol_top;
  while(send!=NULL) {
    tol = *tol_next = (id_list_t *)malloc(sizeof(id_list_t));
    MPI_Irecv( &i, 1, MPI_INT, send->othernode, 0, MPI_COMM_THISJOB, &sreq );
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
  MPI_Allreduce( &n_recv, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_THISJOB );
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

   This handles gathers from both the site structure and an array of
   fields and is not called directly by the user.  Instead they should
   call declare_gather_site() or declare_gather_field().

 prepare_gather() allocates buffers needed for the gather.  This call is
   optional since it will automatically be called from do_gather() if
   not explicitly called before.

 do_gather() starts the actual gather.  This may be repeated after a
    wait_gather() to repeat the exact same gather.

 wait_gather() waits for the gather to finish.

 cleanup_gather() frees memory allocated for the gather including the msg_tag.

   example:
	msg_tag *tag;
	tag = declare_gather_site( F_OFFSET(phi), sizeof(su3_vector), XUP,
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
**  handles gathers from both site structure and arrays (field)
*/
msg_tag *
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
#ifdef CRC_DEBUG
  mtag->index = index;
#endif

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
    mrecv[i].msg_buf = tpt = (char *)malloc( mrecv[i].msg_size+CRCBYTES );
    if(tpt==NULL) {
      printf("NO ROOM for msg_buf, node %d\n", mynode());
      terminate(1);
    }
#ifdef CRC_DEBUG
    memset(tpt, '\0', mrecv[i].msg_size+CRCBYTES);
#endif
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
    msend[i].msg_buf = (char *)malloc( msend[i].msg_size+CRCBYTES );
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
    MPI_Irecv( mbuf[i].msg_buf, mbuf[i].msg_size+CRCBYTES, MPI_BYTE, MPI_ANY_SOURCE,
	       GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_THISJOB,
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
#ifdef COM_CRC
    {
      int msg_size;
      char *crc_pt;
      u_int32type *crc;

      tpt = mbuf[i].msg_buf;
      msg_size = mbuf[i].msg_size;
      crc_pt = tpt + msg_size;
      crc = (u_int32type *)crc_pt;

      *crc = crc32(0, tpt, msg_size );
#ifdef CRC_DEBUG
      {
	char filename[128];
	FILE *dump;
	sprintf(filename,"/tmp/send.%d.to.%d.dir%d.msg%d",
		mynode(),mbuf[i].msg_node,mtag->index,i);
	dump = fopen(filename,"w");
	fwrite(tpt, 1, msg_size + CRCBYTES, dump);
	fclose(dump);
      }
#endif      
    }
#endif
    MPI_Isend( mbuf[i].msg_buf, mbuf[i].msg_size+CRCBYTES, MPI_BYTE, mbuf[i].msg_node,
	       GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_THISJOB,
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
#ifdef COM_CRC
  int fail = 0, work = 0;
#endif

  /* wait for all receive messages */
  for(i=0; i<mtag->nrecvs; i++) {
    MPI_Wait( &mtag->recv_msgs[i].msg_req, &status );
  }

  /* wait for all send messages */
  for(i=0; i<mtag->nsends; i++) {
    MPI_Wait( &mtag->send_msgs[i].msg_req, &status );
  }
#if COM_CRC
  /* Verify the checksums received */
  for(i=0; i<mtag->nrecvs; i++) {
    {
      u_int32type crcgot;
      msg_sr_t *mbuf;
      char *tpt;
      int msg_size;
      char *crc_pt;
      u_int32type *crc;

      mbuf = mtag->recv_msgs;
      tpt = mbuf[i].msg_buf;
      msg_size = mbuf[i].msg_size;
      crc_pt = tpt + msg_size;
      crc = (u_int32type *)crc_pt;
      crcgot = crc32(0, tpt, msg_size );

      if(*crc != crcgot){
	fprintf(stderr,
		"Node %d received checksum %x != node %d sent checksum %x\n",
		mynode(),*crc, mbuf[i].msg_node, crcgot);
	fflush(stdout);
	fail = 1;
#ifdef CRC_DEBUG
	{
	  char filename[128];
	  FILE *dump;
	  sprintf(filename,"/tmp/receive.%d.from.%d.dir%d.msg%d",mynode(),
		  mbuf[i].msg_node,mtag->index,i);
	  dump = fopen(filename,"w");
	  fwrite(tpt, 1, msg_size + CRCBYTES, dump);
	  fclose(dump);
	}
#endif      
      }
    }
  }
  MPI_Allreduce( &fail, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_THISJOB );
  fail = work;
  if(fail > 0)terminate(1);
#endif
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
declare_gather_site(
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
 * gather routines from arrays of fields *
 *****************************/

/*
**  declares a gather from a field
*/
msg_tag *
declare_gather_field(
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

 declare_accumulate_gather_site() declares and joins gathers.

 example:

   msg_tag *tag1, *tag2, *mtag;

   tag1 = declare_gather_site( F_OFFSET(phi), sizeof(su3_vector), XUP,
	                  EVEN, gen_pt1 );
   tag2 = declare_gather_site( F_OFFSET(phi), sizeof(su3_vector), XDOWN,
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
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), XUP,
	                      EVEN, gen_pt1 );
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), XDOWN,
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
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), XUP,
	                      EVEN, gen_pt1 );
 with
   mtag = declare_gather_site( F_OFFSET(phi), sizeof(su3_vector), XUP,
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
**  handles both the site structure and an array of fields
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
declare_accumulate_gather_site(
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
**  declares and merges gather from field
*/
void
declare_accumulate_gather_field(
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

 start_general_gather_site() returns a msg_tag which will
   be used as input to subsequent wait_general_gather() and
   cleanup_general_gather() calls.

   usage: tag = start_general_gather_site( source, size, displacement, parity, dest)
   example:
	msg_tag *tag;
	int disp[4]; 
        disp[XUP]=1; disp[YUP]= -1; disp[ZUP] = disp[TUP] = 0;
	tag = start_general_gather_site( F_OFFSET(phi), sizeof(su3_vector), disp,
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
   start_general_gather_site() and used in wait_general_gather().  This
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
	       MPI_COMM_THISJOB, &mrecv[i].msg_req );
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
	       MPI_COMM_THISJOB, &msend[i].msg_req );
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
	       MPI_COMM_THISJOB, &mrecv[i].msg_req );
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
	       MPI_COMM_THISJOB, &msend[i].msg_req );
  }

  /* free temporary arrays */
  if( n_send_msgs > 0) free(to_nodes);
  /* mark gather in progress and return */
  g_gather_flag = 1;

  return mtag;
}

#endif	/* N_SUBL32 */

msg_tag *
start_general_gather_site(
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
start_general_gather_field(
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
#ifdef COM_CRC
/*
** compute crc32 checksum
*/
/* Taken from the GNU CVS distribution and
   modified for SciDAC use  C. DeTar 10/11/2003 
   and MILC use 5/3/2005 */

/* crc32.c -- compute the CRC-32 of a data stream
 * Copyright (C) 1995-1996 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h 
 */

/* Copyright notice reproduced from zlib.h -- (C. DeTar)

  version 1.0.4, Jul 24th, 1996.

  Copyright (C) 1995-1996 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  gzip@prep.ai.mit.edu    madler@alumni.caltech.edu


  The data format used by the zlib library is described by RFCs (Request for
  Comments) 1950 to 1952 in the files ftp://ds.internic.net/rfc/rfc1950.txt
  (zlib format), rfc1951.txt (deflate format) and rfc1952.txt (gzip format).
*/

typedef u_int32type uLong;            /* At least 32 bits */
typedef unsigned char Byte;
typedef Byte Bytef;
typedef uLong uLongf;
#define Z_NULL  0  /* for initializing zalloc, zfree, opaque */

#define local static

#ifdef DYNAMIC_CRC_TABLE

local int crc_table_empty = 1;
local uLongf crc_table[256];
local void make_crc_table OF((void));

/*
  Generate a table for a byte-wise 32-bit CRC calculation on the polynomial:
  x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1.

  Polynomials over GF(2) are represented in binary, one bit per coefficient,
  with the lowest powers in the most significant bit.  Then adding polynomials
  is just exclusive-or, and multiplying a polynomial by x is a right shift by
  one.  If we call the above polynomial p, and represent a byte as the
  polynomial q, also with the lowest power in the most significant bit (so the
  byte 0xb1 is the polynomial x^7+x^3+x+1), then the CRC is (q*x^32) mod p,
  where a mod b means the remainder after dividing a by b.

  This calculation is done using the shift-register method of multiplying and
  taking the remainder.  The register is initialized to zero, and for each
  incoming bit, x^32 is added mod p to the register if the bit is a one (where
  x^32 mod p is p+x^32 = x^26+...+1), and the register is multiplied mod p by
  x (which is shifting right by one and adding x^32 mod p if the bit shifted
  out is a one).  We start with the highest power (least significant bit) of
  q and repeat for all eight bits of q.

  The table is simply the CRC of all possible eight bit values.  This is all
  the information needed to generate CRC's on data a byte at a time for all
  combinations of CRC register values and incoming bytes.
*/
local void 
make_crc_table()
{
  uLong c;
  int n, k;
  uLong poly;            /* polynomial exclusive-or pattern */
  /* terms of polynomial defining this crc (except x^32): */
  static Byte p[] = {0,1,2,4,5,7,8,10,11,12,16,22,23,26};

  /* make exclusive-or pattern from polynomial (0xedb88320L) */
  poly = 0L;
  for (n = 0; n < sizeof(p)/sizeof(Byte); n++)
    poly |= 1L << (31 - p[n]);
 
  for (n = 0; n < 256; n++)
  {
    c = (uLong)n;
    for (k = 0; k < 8; k++)
      c = c & 1 ? poly ^ (c >> 1) : c >> 1;
    crc_table[n] = c;
  }
  crc_table_empty = 0;
}
#else
/* ========================================================================
 * Table of CRC-32's of all single-byte values (made by make_crc_table)
 */
local uLongf crc_table[256] = {
  0x00000000L, 0x77073096L, 0xee0e612cL, 0x990951baL, 0x076dc419L,
  0x706af48fL, 0xe963a535L, 0x9e6495a3L, 0x0edb8832L, 0x79dcb8a4L,
  0xe0d5e91eL, 0x97d2d988L, 0x09b64c2bL, 0x7eb17cbdL, 0xe7b82d07L,
  0x90bf1d91L, 0x1db71064L, 0x6ab020f2L, 0xf3b97148L, 0x84be41deL,
  0x1adad47dL, 0x6ddde4ebL, 0xf4d4b551L, 0x83d385c7L, 0x136c9856L,
  0x646ba8c0L, 0xfd62f97aL, 0x8a65c9ecL, 0x14015c4fL, 0x63066cd9L,
  0xfa0f3d63L, 0x8d080df5L, 0x3b6e20c8L, 0x4c69105eL, 0xd56041e4L,
  0xa2677172L, 0x3c03e4d1L, 0x4b04d447L, 0xd20d85fdL, 0xa50ab56bL,
  0x35b5a8faL, 0x42b2986cL, 0xdbbbc9d6L, 0xacbcf940L, 0x32d86ce3L,
  0x45df5c75L, 0xdcd60dcfL, 0xabd13d59L, 0x26d930acL, 0x51de003aL,
  0xc8d75180L, 0xbfd06116L, 0x21b4f4b5L, 0x56b3c423L, 0xcfba9599L,
  0xb8bda50fL, 0x2802b89eL, 0x5f058808L, 0xc60cd9b2L, 0xb10be924L,
  0x2f6f7c87L, 0x58684c11L, 0xc1611dabL, 0xb6662d3dL, 0x76dc4190L,
  0x01db7106L, 0x98d220bcL, 0xefd5102aL, 0x71b18589L, 0x06b6b51fL,
  0x9fbfe4a5L, 0xe8b8d433L, 0x7807c9a2L, 0x0f00f934L, 0x9609a88eL,
  0xe10e9818L, 0x7f6a0dbbL, 0x086d3d2dL, 0x91646c97L, 0xe6635c01L,
  0x6b6b51f4L, 0x1c6c6162L, 0x856530d8L, 0xf262004eL, 0x6c0695edL,
  0x1b01a57bL, 0x8208f4c1L, 0xf50fc457L, 0x65b0d9c6L, 0x12b7e950L,
  0x8bbeb8eaL, 0xfcb9887cL, 0x62dd1ddfL, 0x15da2d49L, 0x8cd37cf3L,
  0xfbd44c65L, 0x4db26158L, 0x3ab551ceL, 0xa3bc0074L, 0xd4bb30e2L,
  0x4adfa541L, 0x3dd895d7L, 0xa4d1c46dL, 0xd3d6f4fbL, 0x4369e96aL,
  0x346ed9fcL, 0xad678846L, 0xda60b8d0L, 0x44042d73L, 0x33031de5L,
  0xaa0a4c5fL, 0xdd0d7cc9L, 0x5005713cL, 0x270241aaL, 0xbe0b1010L,
  0xc90c2086L, 0x5768b525L, 0x206f85b3L, 0xb966d409L, 0xce61e49fL,
  0x5edef90eL, 0x29d9c998L, 0xb0d09822L, 0xc7d7a8b4L, 0x59b33d17L,
  0x2eb40d81L, 0xb7bd5c3bL, 0xc0ba6cadL, 0xedb88320L, 0x9abfb3b6L,
  0x03b6e20cL, 0x74b1d29aL, 0xead54739L, 0x9dd277afL, 0x04db2615L,
  0x73dc1683L, 0xe3630b12L, 0x94643b84L, 0x0d6d6a3eL, 0x7a6a5aa8L,
  0xe40ecf0bL, 0x9309ff9dL, 0x0a00ae27L, 0x7d079eb1L, 0xf00f9344L,
  0x8708a3d2L, 0x1e01f268L, 0x6906c2feL, 0xf762575dL, 0x806567cbL,
  0x196c3671L, 0x6e6b06e7L, 0xfed41b76L, 0x89d32be0L, 0x10da7a5aL,
  0x67dd4accL, 0xf9b9df6fL, 0x8ebeeff9L, 0x17b7be43L, 0x60b08ed5L,
  0xd6d6a3e8L, 0xa1d1937eL, 0x38d8c2c4L, 0x4fdff252L, 0xd1bb67f1L,
  0xa6bc5767L, 0x3fb506ddL, 0x48b2364bL, 0xd80d2bdaL, 0xaf0a1b4cL,
  0x36034af6L, 0x41047a60L, 0xdf60efc3L, 0xa867df55L, 0x316e8eefL,
  0x4669be79L, 0xcb61b38cL, 0xbc66831aL, 0x256fd2a0L, 0x5268e236L,
  0xcc0c7795L, 0xbb0b4703L, 0x220216b9L, 0x5505262fL, 0xc5ba3bbeL,
  0xb2bd0b28L, 0x2bb45a92L, 0x5cb36a04L, 0xc2d7ffa7L, 0xb5d0cf31L,
  0x2cd99e8bL, 0x5bdeae1dL, 0x9b64c2b0L, 0xec63f226L, 0x756aa39cL,
  0x026d930aL, 0x9c0906a9L, 0xeb0e363fL, 0x72076785L, 0x05005713L,
  0x95bf4a82L, 0xe2b87a14L, 0x7bb12baeL, 0x0cb61b38L, 0x92d28e9bL,
  0xe5d5be0dL, 0x7cdcefb7L, 0x0bdbdf21L, 0x86d3d2d4L, 0xf1d4e242L,
  0x68ddb3f8l, 0x1fda836eL, 0x81be16cdL, 0xf6b9265bL, 0x6fb077e1L,
  0x18b74777L, 0x88085ae6L, 0xff0f6a70L, 0x66063bcaL, 0x11010b5cL,
  0x8f659effL, 0xf862ae69L, 0x616bffd3L, 0x166ccf45L, 0xa00ae278L,
  0xd70dd2eeL, 0x4e048354L, 0x3903b3c2L, 0xa7672661L, 0xd06016f7L,
  0x4969474dL, 0x3e6e77dbL, 0xaed16a4aL, 0xd9d65adcL, 0x40df0b66L,
  0x37d83bf0L, 0xa9bcae53L, 0xdebb9ec5L, 0x47b2cf7fL, 0x30b5ffe9L,
  0xbdbdf21cL, 0xcabac28aL, 0x53b39330L, 0x24b4a3a6L, 0xbad03605L,
  0xcdd70693L, 0x54de5729L, 0x23d967bfL, 0xb3667a2eL, 0xc4614ab8L,
  0x5d681b02L, 0x2a6f2b94L, 0xb40bbe37L, 0xc30c8ea1L, 0x5a05df1bL,
  0x2d02ef8dL
};
#endif

/* =========================================================================
 * This function can be used by asm versions of crc32()
 */
static uLongf *get_crc_table()
{
#ifdef DYNAMIC_CRC_TABLE
  if (crc_table_empty) make_crc_table();
#endif
  return (uLongf *)crc_table;
}

/* ========================================================================= */
#define DO1(buf) crc = crc_table[((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
#define DO2(buf)  DO1(buf); DO1(buf);
#define DO4(buf)  DO2(buf); DO2(buf);
#define DO8(buf)  DO4(buf); DO4(buf);

/* ========================================================================= */
u_int32type 
crc32(u_int32type crc, const unsigned char *buf, size_t len)
{
    if (buf == Z_NULL) return 0L;
#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
      make_crc_table();
#endif
    crc = crc ^ 0xffffffffL;
    while (len >= 8)
    {
      DO8(buf);
      len -= 8;
    }
    if (len) do {
      DO1(buf);
    } while (--len);
    return crc ^ 0xffffffffL;
}

#endif
