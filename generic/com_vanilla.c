/******************  com_vanilla.c **************************************/
/* Communications routines for the SU3 program
   MIMD version 7.
   This file is communications-scheme dependent.
   Version for single processor machines.
*/
/* Modifications

   4/20/02 added start_general_gather_field C.D.
   1/28/00 combined with Schroedinger functional and
           32 sublattice versions - UMH
   9/2/97  Revised to allow gathers from temporary fields.  neighbor[]
           is now list of indices, add start/restart_gather_field D.T.
   8/05/97 ANSI prototyping for all routines C.D.
   10/5/96 Removed parallel I/O wrappers. Use io_ansi.c now.
   8/30/96 Added restore/save_checkpoint for compatibility
   8/05/96 Added broadcast_bytes and wrappers for system-dependent
           parallel file system calls C.D.
*/
/*
   This is a trivial version, where there is only one node, every site
   is on node ...
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
   g_floatsum()          sums a Realing point number over all nodes.
   g_vecfloatsum()       sums a vector of Reals over all nodes 
   g_doublesum()         sums a double over all nodes.
   g_vecdoublesum()      sums a vector of doubles over all nodes.
   g_complexsum()        sums a single precision complex number over all nodes.
   g_veccomplexsum()     sums a vector of single precision complex numbers
                           over all nodes.
   g_dcomplexsum()       sums a double precision complex number over all nodes.
   g_vecdcomplexsum()    sums a vector of double_complex over all nodes 
   g_wvectorsumfloat()   sums a single precision wilson vector over all nodes.
   g_xor32()             finds global exclusive or of 32-bit word
   g_floatmax()          finds maximum Realing point number over all nodes.
   g_doublemax()         finds maximum double over all nodes.
   broadcast_float()     broadcasts a single precision number from
	                   node 0 to all nodes.
   broadcast_double()    broadcasts a double precision number
   broadcast_complex()   broadcasts a single precision complex number
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
                               declare_gather_field/prepare_gather/do_gather
   restart_gather_site()      older function which is obsoleted by do_gather()
   restart_gather_field() older function which is obsoleted by do_gather() 
   start_general_gather_site()  starts asynchronous sends and receives required
                             to gather fields at arbitrary displacement.
   start_general_gather_field() starts asynchronous sends and receives 
                             required to gather neighbors from a temporary 
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

#define NOWHERE -1      /* Not an index in array of fields */

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

/* structure to hold all necessary info for a gather */
typedef struct gather_t {
  int *neighbor;    /* keeps track of where the neighbor is */
} gather_t;

/* Structure to keep track of outstanding sends and receives */
struct msg_tag {
 int dummy;
};  /* don't need anything */

/***************************************************
 *  Global variables for the communications stuff  *
 ***************************************************/

/* array storing gather setup info */
static gather_t *gather_array;

/* Number of gathers (mappings) that have been set up */
static int n_gathers, gather_array_len;


/**********************************************************************
 *                BASIC COMMUNICATIONS FUNCTIONS                      *
 **********************************************************************/

/*
**  Machine initialization
*/
void
initialize_machine(int *argc, char ***argv)
{
  /* check if 32 bit int is set correctly */
#ifdef SHORT_IS_32BIT
  if(sizeof(unsigned short)!=4) {
    printf("node %i: SHORT_IS_32BIT is set but sizeof(unsigned short)=%i\n",
	   mynode(), sizeof(unsigned short));
    terminate(1);
  }
#else
  if(sizeof(unsigned int)!=4) {
    printf("node %i: SHORT_IS_32BIT is not set but sizeof(unsigned int)=%d\n",
	   mynode(), (int)sizeof(unsigned int));
    terminate(1);
  }
#endif

  n_gathers = 0;
  gather_array_len = 0;
  gather_array = NULL;
}

/*
**  version of normal exit for scalar processes
*/
void
normal_exit(int status)
{
  time_stamp("exit");
  fflush(stdout);
  exit(status);
}

/*
**  version of exit for scalar processes
*/
void
terminate(int status)
{
  time_stamp("termination");
  printf("Termination: node %d, status = %d\n", this_node, status);
  fflush(stdout);
  exit(status);
}

/*
**  Tell what kind of machine we are on
*/
static char name[]="Scalar processor";
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
  return(0);
}

/*
**  Return number of nodes
*/
int
numnodes(void)
{
  return(1);
}

/*
**  Return my jobid
*/
int
myjobid(void)
{
  return 0;
}

/*
**  Return number of jobs
*/
int
numjobs(void)
{
  return 1;
}

/*
** Return the job geometry
*/
int const *
jobgeom(void)
{
  static int ones[4] = {1,1,1,1};
  return ones;
}

/*
** Return the ionode geometry
*/
int *
ionodegeom(void)
{
  static int ones[4] = {1,1,1,1};
  return ones;
}

/*
** Return the allocated dimensions (node geometry) if a grid is being used
*/
const int *
nodegeom(void)
{
  static int ones[4] = {1,1,1,1};
  return ones;
}

/*
**  Synchronize all nodes
*/
void
g_sync(void)
{
}

/*
**  Sum signed integer over all nodes
*/
void
g_intsum(int *ipt)
{
}

/*
**  Sum unsigned 32-bit integer type
*/
void
g_uint32sum(u_int32type *pt)
{
}

/*
**  Sum Real over all nodes
*/
void
g_floatsum(Real *fpt)
{
}

/*
**  Sum a vector of Reals over all nodes
*/
void
g_vecfloatsum(Real *fpt, int nReals)
{
}

/*
**  Sum double over all nodes
*/
void
g_doublesum(double *dpt)
{
}

/*
**  Sum a vector of doubles over all nodes
*/
void
g_vecdoublesum(double *dpt, int ndoubles)
{
}

/*
**  Sum complex over all nodes
*/
void
g_complexsum(complex *cpt)
{
}

/*
**  Sum a vector of complex over all nodes
*/
void
g_veccomplexsum(complex *cpt, int ncomplex)
{
}

/*
**  Sum double_complex over all nodes
*/
void
g_dcomplexsum(double_complex *cpt)
{
}

/*
**  Sum a vector of double_complex over all nodes
*/
void
g_vecdcomplexsum(double_complex *cpt, int ncomplex)
{
}

/*
**  Sum wilson_vector over all nodes
*/
void
g_wvectorsumfloat(wilson_vector *wvpt)
{
}

/*
**  Global exclusive or acting on u_int32type
*/
void
g_xor32(u_int32type *pt)
{
}

/*
**  Find maximum of Real over all nodes
*/
void
g_floatmax(Real *fpt)
{
}

/*
**  Find maximum of double over all nodes
*/
void
g_doublemax(double *dpt)
{
}

/*
**  Broadcast Realing point number from node zero
*/
void
broadcast_float(Real *fpt)
{
}

/*
**  Broadcast double precision Realing point number from node zero
*/
void
broadcast_double(double *dpt)
{
}

/*
**  Broadcast single precision complex number from node zero
*/
void
broadcast_complex(complex *cpt)
{
}

/*
**  Broadcast double precision complex number from node zero
*/
void
broadcast_dcomplex(double_complex *cpt)
{
}

/*
**  Broadcast bytes from node 0 to all others
*/
void
broadcast_bytes(char *buf, int size)
{
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
  printf("ERROR: called send_integer() in com_vanilla.c\n");
  terminate(1);
}

/*
**  Receive an integer from another node
*/
void
receive_integer(int fromnode, int *address)
{
  printf("ERROR: called receive_integer() in com_vanilla.c\n");
  terminate(1);
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
  printf("ERROR: called send_field() in com_vanilla.c\n");
  terminate(1);
}

/*
**  get_field is to be called only by the node to which the field was sent
*/
void
get_field(char *buf, int size, int fromnode)
{
  printf("ERROR: called get_field() in com_vanilla.c\n");
  terminate(1);
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

  time(&time_stamp);
  printf("%s: %s\n", msg, ctime(&time_stamp));
  fflush(stdout);
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
    printf("error: not enought room for gather_array in make_nn_gathers\n");
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
  int i,subl;	        /* scratch */
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
  /* Fill in pointers to sites */
  FORALLSITES(i,s){
    /* find coordinates of neighbor who sends us data */
    func( s->x, s->y, s->z, s->t, args, FORWARDS, &x, &y, &z, &t );
    gather_array[dir].neighbor[i] = node_index(x,y,z,t);
  }

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
    func( s->x, s->y, s->z, s->t, args, BACKWARDS, &x, &y, &z, &t );
    /* set up pointer */
    gather_array[dir].neighbor[i] = node_index(x,y,z,t);
  }

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
**  handles gathers from both the site structure and an array of fields
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
  gather_t *gt;         /* pointer to current gather */

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

  return(NULL);
}

/*
**  allocate buffers for gather
*/
void
prepare_gather(msg_tag *mtag)
{
}

/*
**  actually execute the gather
*/
void
do_gather(msg_tag *mtag)  /* previously returned by start_gather_site */
{
}

/*
**  wait for gather to finish
*/
void
wait_gather(msg_tag *mtag)
{
}

/*
**  free buffers associated with message tag
*/
void
cleanup_gather(msg_tag *mtag)
{
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
  do_gather(mtag);
}

/*****************************
 * gather routines from an array of fields *
 *****************************/

/*
**  declares a gather from an array of fields
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
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), 
                              XUP, EVEN, gen_pt1 );
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), 
                              XDOWN, EVEN, gen_pt2 );
   prepare_gather( mtag );  ** optional **
   do_gather( mtag );
   wait_gather( mtag );
   ** stuff **
   do_gather( mtag );
   wait_gather( mtag );
   cleanup_gather( mtag );

 one coule also replace
   mtag = NULL;
   declare_accumulate_gather_site( &mtag, F_OFFSET(phi), sizeof(su3_vector), 
                              XUP,EVEN, gen_pt1 );
 with
   mtag = declare_gather_site( F_OFFSET(phi), sizeof(su3_vector), XUP,
	                  EVEN, gen_pt1 );
 since they do the same thing, however the first form is a bit more uniform
 in the given example.
*/

/*
**  merges already declared gather
*/
void
accumulate_gather(msg_tag **mmtag, msg_tag *mtag)
{
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
**  declares and merges gather from an array of fields
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
  declare_accumulate_strided_gather( mmtag, (char *)field, size, size, index, parity,
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

static int g_gather_flag=0; /* flag to tell if general gather in progress */

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
  site *s;	    /* scratch pointer to site */
  int i;     	    /* scratch */
  int tx,ty,tz,tt;  /* temporary coordinates */

  /* check for gather already in progress */
  if(g_gather_flag!=0) {
    fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
    exit(1);
  }

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) Make
     list of nodes from whom we expect messages */
  if( subl == EVENANDODD ) {
    FORALLSITES(i,s) {
      if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
      else    tx = s->x;
      if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
      else    ty = s->y;
      if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
      else    tz = s->z;
      if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
      else    tt = s->t;
      dest[i] = field + stride*node_index(tx,ty,tz,tt);
    }
  } else {
    FORSOMEPARITY(i,s,subl) {
      if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
      else    tx = s->x;
      if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
      else    ty = s->y;
      if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
      else    tz = s->z;
      if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
      else    tt = s->t;
      dest[i] = field + stride*node_index(tx,ty,tz,tt);
    }
  }

  /* mark gather in progress and return */
  g_gather_flag = 1;

  return(NULL);
}

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
  return start_general_strided_gather( (char *)field, size, size,
				       displacement, parity, dest );
}

/*
**  wait for a general gather to complete
*/
void
wait_general_gather(msg_tag *mtag)
{
  g_gather_flag=0;
}

/*
**  free memory associated with general gather
*/
void
cleanup_general_gather(msg_tag *mtag)
{
}
