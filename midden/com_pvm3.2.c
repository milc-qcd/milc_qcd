/******************  com_pvm3.2.c **************************************/
/* Communications routines for the SU3 program */
/*   MIMD version 6. */

/* C. DeTar port 18 February 1993 */
/* DT 3/16/93: PVM version 2.4, node 0 starts other node programs 
   DT 4/02/93: PVM version 3.0 
   DT 9/16/93: PVM version 3.2
   9/5/96 CED sort_site_list uses sort key.
*/
#define NOWHERE -1      /* Not an index in array of fields */

/* Modifications

   12/18/97 ANSI prototyping for all routines;
	    added broadcast_bytes and start/restart_gather_from_temp:
	      neighbor[] is now list of indices;
	    added routine g_xor32.
   1/8/98:  added routines g_vecdoublesum(), g_veccomplexsum()
            and g_wvectosumReal().
   1/28/00  combined with Schroedinger functional and
            32 sublattice versions - UMH
*/

/*
   This file is machine dependent.
        FOR PVM VERSION 3.2
        This version implements asynchronous sends and receives.
	Field pointer not done yet.

   The lattice is stored as an array of sites, where "site" is a
   structure defined in "<lattice.h>".

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

   send_field() sends a field to one other node.
   get_links() receives a field from some other node.
*/

#include "generic_includes.h"

/*#include <pvm3.2.h>*/
/* #include "/usr/local/include/pvm3.2.h" UMH */
#include "/usr/local/include/pvm3.h"	/* UMH */
/**#include "/usr/local/include/xab3.h"**/ /* TEMPORARY */

#define TRUE 1
#define FALSE 0

/* message types  (not all are used on all machines) */
#define PARAM_TYPE 11	/* type of parameter message to nodes */
#define FIELD_TYPE 12	/* type of field sent from one node to another */
#define BROADCAST_MILC_REAL_TYPE 13	/* broadcast of Realing point number */
#define BROADCAST_DOUBLE_TYPE 14 /* broadcast of double */
#define BROADCAST_COMPLEX_TYPE 15 /* broadcast of single precision complex */
#define BROADCAST_DCOMPLEX_TYPE 16 /* broadcast of double precision complex */
#define SEND_INTEGER_TYPE 17	/* send an integer to one other node */
#define CM_GATHER_TYPE 18	/* type for CM5 (cooperative) messages */
#define BROADCAST_BYTES_TYPE 19 /* send n raw bytes from node 0 to all nodes */
#define SYNC_TYPE 50		/* Synchronize all nodes */
#define SUM_MILC_REAL_TYPE 51	/* Sum Real over all nodes */
#define SUM_DOUBLE_TYPE 52	/* Sum double over all nodes */
#define SUM_COMPLEX_TYPE 53	/* Sum complex over all nodes */
#define SUM_DCOMPLEX_TYPE 54	/* Sum double_complex over all nodes */
#define MAX_MILC_REAL_TYPE 55	/* Maximum Real over all nodes */
#define MAX_DOUBLE_TYPE 56	/* Maximum double over all nodes */
#define SUM_VECDOUBLE_TYPE 61   /* Sum a vector of doubles over all nodes */
#define SUM_VECDCOMPLEX_TYPE 62 /* Sum a vector of double complex over nodes */
#define FIELD_REQUEST 100	/* used by field_pointer...() */
#define FIELD_REPLY 101		/* used by field_pointer...() */

#define GENERAL_GATHER_BASE 102	/* types from this to this+number_of_nodes
				are used by the general_gather routines */
#define GATHER_BASE 614	/* types greater than or equal to this are used
			   by the gather routines */

#define MAX_GATHERS 32	/* Maximum number of gather tables */

char pvm_node_name[32] = "Dummy";
int pvm_my_instance=99;	/* instance number supplied by PVM */
int pvm_tid;		/* process identifier supplied by PVM  */
int pvm_number_nodes;
int pvm_nformats;
int pvm_node_tid[MAX_NUMBER_NODES]; /*  task id number for 
					     logical node number */
int my_logical_node;	/* used as mynode value */

/* Global variables for the communications stuff */
    /* message types for gather encode a sequence number for the gather
        so that if several gathers are going at once, you can read
        the message corresponding to the right one. */
    /* for computing message type in gather */
    int mt_nodeshift;	/* number of bits to shift node number */
    int type_offset;	/* label gathers by round-robin */
    /* macro to compute the message type */
#define TYPE_OFFSET_MAX 1024
#define GATHER_MSG_TYPE(node,type_offset) \
   (GATHER_BASE + ((node)<<mt_nodeshift) + (type_offset))
#define GENERAL_GATHER_MSG_TYPE(node) (GENERAL_GATHER_BASE + node)

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

#include <signal.h>
/* Machine initialization */
void initialize_machine(int argc, char **argv){
  int i,j,n,nodes_started;
  char hostlist_name[128]; /* filename, list of hosts and number of procs */
  char other_host[128],this_host[16];
  FILE *fp;
  void gethostname(char *this_host, int length );

  /* Get process identification number.  If this is node 0, enroll in PVM */
  pvm_tid = pvm_mytid();

  /* Use direct task-to-task channels */
  pvm_advise( PvmRouteDirect );

  /* Find out if this is node 0.  If so, start other nodes. */
  if( pvm_parent() == PvmNoParent ){ /* No parent, so I must be node zero */

    /* Find program name and number of nodes, check arguments */
    sscanf(argv[0],"%s",pvm_node_name);
    if( argc > 3 || argc < 2 || (sscanf(argv[1],"%d",&pvm_number_nodes)!=1 )
	|| (argc==3 && (sscanf(argv[2],"%s",hostlist_name)!=1 ) ) ){
        fprintf(stderr,"PVM 3 Usage: %s <processes> [<hostlist>]\n", argv[0]);
	pvm_exit(); exit(1);
    }
    if(pvm_number_nodes > MAX_NUMBER_NODES) {
        printf("%s: Number of processes must not exceed MAX_NUMBER_NODES %d\n",
             pvm_node_name,MAX_NUMBER_NODES);
        pvm_exit(); exit(1);
    }
    if(pvm_number_nodes <= 0) {
        printf("%s: Number of processes must be positive\n",pvm_node_name);
        pvm_exit(); exit(1);
    }

    /* Set up my own bookkeeping variables */
    pvm_node_tid[0] = pvm_tid;
    my_logical_node = 0;

    /* Start up all the other nodes */
    /* Note this is node 0, so start N-1 other nodes.  Put TID's into
	array starting at element 1, since element 0 is already set. */
printf("INIT_M: nodes = %d\n",pvm_number_nodes);
    if( argc==2 ){	/* let PVM decide where to put them */
      if(pvm_number_nodes > 1){
        j = pvm_spawn( pvm_node_name, NULL, PvmTaskDefault, NULL,
	    pvm_number_nodes-1, &(pvm_node_tid[1]) );
        if( j != pvm_number_nodes-1 ){
	    /* something went wrong.  Try to clean up */
	    printf("PVM spawn failed, returned %d\n",j);
	    for(i=1;i<=j;i++)pvm_kill( pvm_node_tid[i] );
	    pvm_exit(); exit(1);
        }
      }
    }
    else {	/* read a file to find out where they go */
	fp = fopen(hostlist_name, "r" );
	if( fp==NULL ){
	    fprintf(stderr,"Can't open file %s\n", hostlist_name );
	    pvm_exit(); exit(1);
	}

	nodes_started = 0;
	while( nodes_started < pvm_number_nodes ){
	    if( (fscanf(fp,"%s%d\n",other_host,&n))!=2 ){
		fprintf(stderr,"Can't understand host file %s\n",
		    hostlist_name); 
		for(i=1;i<nodes_started;i++)pvm_kill( pvm_node_tid[i] );
		pvm_exit(); exit(1);
	    }

	    if( nodes_started==0 ){
		/* First line should be this host */
		gethostname( this_host, 16 );
		if( strcmp( this_host, other_host ) != 0 ){
		    fprintf(stderr,"First host must be %s\n", this_host); 
		    /**for(i=1;i<nodes_started;i++)pvm_kill( pvm_node_tid[i] );
		    pvm_exit(); exit(1);**/
		}
		nodes_started++;
		n--;	/* This node is already started */
	    }
	    if( n > 0 ){
/*printf("Starting %d nodes on %s\n",n,other_host);*/
		if( nodes_started + n > pvm_number_nodes )
		    n = pvm_number_nodes - nodes_started;
        	j = pvm_spawn( pvm_node_name, NULL, PvmTaskHost, other_host,
		    n, &(pvm_node_tid[nodes_started]) );
        	if( j != n ){
		    /* something went wrong.  Try to clean up */
		    printf("PVM spawn failed, returned %d\n",j);
		    for(i=1;i<=nodes_started+j;i++)pvm_kill( pvm_node_tid[i] );
		    pvm_exit(); exit(1);
		}
		else { nodes_started += j; }
	    }
	}
	fclose(fp);
    }

    /* Send environmental information to each of the nodes */

    hcs_ident.msg_type = NODE_IDENT_CALL;
	/* Copy of message type inside message tests readability */
    hcs_ident.number_nodes = pvm_number_nodes;
    for( i = 0; i < pvm_number_nodes; i++)
      hcs_ident.node_tid[i] = pvm_node_tid[i];

    for( i = 1; i < pvm_number_nodes; i++) {
        hcs_ident.your_node = i;
        pvm_initsend(PvmDataRaw);
        pvm_pkbyte( (char *)&hcs_ident, HCS_IDENT_SIZE, 1 );
        if(pvm_send(pvm_node_tid[i],NODE_IDENT_CALL) < 0) {
	    printf("%s: Error sending ident to node %d\n",pvm_node_name,i);
	    terminate(1);
	}
    }

    /* Once started, control-C should kill all processes */
    signal(SIGHUP,terminate);
    signal(SIGINT,terminate);
    signal(SIGQUIT,terminate);

  } /* End if node 0 */

  else {	/* I am not node zero */
    /* At the very start, we expect a message from the host */
    pvm_recv(ANY_NODE,NODE_IDENT_CALL);
    pvm_upkbyte( (char *)&hcs_ident, HCS_IDENT_SIZE, 1 );
  
    /* Verify that this message is of correct type and readable */
    if(hcs_ident.msg_type != NODE_IDENT_CALL) {
        fprintf(stderr,
  	"com_pvm: Missed ident message at node %s instance %d, got %d = %x\n",
  	   pvm_node_name,pvm_my_instance,hcs_ident.msg_type,hcs_ident.msg_type);
        pvm_exit(); exit(1);
    }
    else{ /* just for debugging */
      fprintf(stderr, "com_pvm: Got ident message at node %d tid %d\n",
          hcs_ident.your_node,hcs_ident.node_tid[hcs_ident.your_node]);
    }
  
    /* Unpack message */
    my_logical_node = hcs_ident.your_node;
    pvm_number_nodes = hcs_ident.number_nodes;
    for(i = 0; i < pvm_number_nodes; i++)
      pvm_node_tid[i] = hcs_ident.node_tid[i];
  }	/* End if not node 0 */

  type_offset = 0;
}


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

   /* Set up variables for constructing message types */
   for( mt_nodeshift=0,i=TYPE_OFFSET_MAX-1; i>0; i>>=1,mt_nodeshift++);

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
    sort_eight_special((void **) neighbor );
    sort_eight_special((void **) neighborlist );
    sort_eight_special((void **) neighborlist_send );
}

/* sort a list of eight pointers into the order we want for the
  nearest neighbor gathers:  XUP,YUP,ZUP,TUP,TDOWN,ZDOWN,YDOWN,XDOWN */
void sort_eight_special(void **pt){
void *tt[8];
register int i;
    for(i=0;i<8;i++)tt[i]=pt[i];
    for(i=XUP;i<=TUP;i++){pt[i]=tt[2*i]; pt[OPP_DIR(i)]=tt[2*i+1];}
}

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
			/* SWITCH_PARITY for the 32 sublattice case means,
			   that for all sites in a sublattice, subl, the
			   target sublattice is the same, and different
			   from sublattice subl. */
{
#ifdef N_SUBL32
/* copy linked list of comlinks, switching two sublattices */
comlink * copy_list_switch( comlink *old_compt, int *othersubl);
#else
/* copy linked list of comlinks, switching even and odd */
comlink * copy_list_switch( comlink *old_compt );
#endif

comlink * make_send_receive_list(
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
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
        neighbor[n_gathers-1] = (int *)malloc(sites_on_node*sizeof(int) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-2;	/* index of gather we are working on */
    }
    else {
	neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-1] = (int *)malloc(sites_on_node*sizeof(int) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
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

#ifndef N_SUBL32

comlink * copy_list_switch( comlink *old_compt ) {
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
    if(key == NULL) {
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

/* utility function for finding coordinates of neighbor */
void neighbor_coords(
 int x, int y, int z, int t,  /* coordinates of site */
 int dir,                     /* direction (eg XUP) */
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


/* Set up interrupt handlers to handle field_pointer routines
   make_lattice() must be called first. */
static msg_request mreqbuf;	/* global so handler can use it too */
void start_handlers(){
void fillfieldrequest(int type,int size,int node,int pid);
    /**
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    **/
if(this_node==0)printf("Can't start interrupt handler yet!!\n");
}
/* The handler for field requests.  It fills the requests, then posts a
   receive for the next request. */
void fillfieldrequest(int type,int size,int node,int pid){
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

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

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
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(compt->othernode,type_offset);
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
	mbuf[i].msg_OK=FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
        pvm_initsend(PvmDataRaw);
   /**pvm_initsend(PvmDataInPlace);**/
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        /**memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);**/
	        pvm_pkbyte( ((char *)(&lattice[compt->ositelist[j]])+field),
		    size, 1 );
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        /**memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);**/
	        pvm_pkbyte( ((char *)(&lattice[compt->esitelist[j]])+field),
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[mbuf[i].msg_node],
	    GATHER_MSG_TYPE(mynode(),type_offset) );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
}

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
      mbuf[i].msg_OK = FALSE;
      mbuf[i].msg_id = GATHER_MSG_TYPE(mbuf[i].msg_node,type_offset);
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	mbuf[i].msg_OK = FALSE;
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        /**memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);**/
	        pvm_pkbyte( ((char *)(&lattice[compt->ositelist[j]])+field),
		    size, 1 );
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        /**memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);**/
	        pvm_pkbyte( ((char *)(&lattice[compt->esitelist[j]])+field),
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[compt->othernode], 
	    GATHER_MSG_TYPE(mynode(),type_offset));
    }
}

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

    type_offset = (type_offset+1)%TYPE_OFFSET_MAX;	/* label this gather */

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
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = (char *)malloc( nsites*size );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(compt->othernode,type_offset);
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
	mbuf[i].msg_OK=FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
        pvm_initsend(PvmDataRaw);
   /**pvm_initsend(PvmDataInPlace);**/
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        /**memcpy( tpt, (char *)field + compt->ositelist[j]*size,
		    size);**/
	        pvm_pkbyte( (char *)field + compt->ositelist[j]*size,
		    size, 1 );
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        /**memcpy( tpt, (char *)field + compt->esitelist[j]*size,
		    size);**/
	        pvm_pkbyte( (char *)field + compt->esitelist[j]*size,
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[mbuf[i].msg_node],
	    GATHER_MSG_TYPE(mynode(),type_offset) );
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
} /*start_gather_from_temp()*/

/* Repeat a gather with the same source and destination as a
  previous gather.  The previous gather must have been waited for
  but not cleaned up.  Pointers to sites on the same node are not
  reset, and the same buffers are reused. */
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
      mbuf[i].msg_OK = FALSE;
      mbuf[i].msg_id = GATHER_MSG_TYPE(mbuf[i].msg_node,type_offset);
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	mbuf[i].msg_OK = FALSE;
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        /**memcpy( tpt, (char *)field + compt->ositelist[j]*size,
		    size);**/
	        pvm_pkbyte( (char *)field + compt->ositelist[j]*size,
		    size, 1 );
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        /**memcpy( tpt, (char *)field + compt->esitelist[j]*size,
		    size);**/
	        pvm_pkbyte( (char *)field + compt->esitelist[j]*size,
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[compt->othernode], 
	    GATHER_MSG_TYPE(mynode(),type_offset));
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
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(compt->othernode,type_offset);
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
	i++;
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=0;
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
	mbuf[i].msg_OK=FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	pvm_initsend(PvmDataRaw);
	/**pvm_initsend(PvmDataInPlace);**/
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    s = &(lattice[compt->sitelist[0][j]]);
		    /**memcpy( tpt, ((char *)(s)+field), size);**/
		    pvm_pkbyte( ((char *)(s)+field), size, 1 );
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    s = &(lattice[compt->sitelist[subl_s][j]]);
		    /**memcpy( tpt, ((char *)(s)+field), size);**/
		    pvm_pkbyte( ((char *)(s)+field), size, 1 );
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		s = &(lattice[compt->sitelist[subl_s][j]]);
		/**memcpy( tpt, ((char *)(s)+field), size);**/
		pvm_pkbyte( ((char *)(s)+field), size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[mbuf[i].msg_node],
	    GATHER_MSG_TYPE(mynode(),type_offset) );
	i++;
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
}

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
	mbuf[i].msg_OK = FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mbuf[i].msg_node,type_offset);
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
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	mbuf[i].msg_OK = FALSE;
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    s = &(lattice[compt->sitelist[0][j]]);
		    /**memcpy( tpt, ((char *)(s)+field), size);**/
		    pvm_pkbyte( ((char *)(s)+field), size, 1 );
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    s = &(lattice[compt->sitelist[subl_s][j]]);
		    /**memcpy( tpt, ((char *)(s)+field), size);**/
		    pvm_pkbyte( ((char *)(s)+field), size, 1 );
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		s = &(lattice[compt->sitelist[subl_s][j]]);
		/**memcpy( tpt, ((char *)(s)+field), size);**/
		pvm_pkbyte( ((char *)(s)+field), size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[compt->othernode], 
	    GATHER_MSG_TYPE(mynode(),type_offset));
	i++;
    }
}

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
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for pvm the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(compt->othernode,type_offset);
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
	i++;
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
	mbuf[n_recv_msgs].msg_id=0;
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
	mbuf[i].msg_OK=FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	pvm_initsend(PvmDataRaw);
	/**pvm_initsend(PvmDataInPlace);**/
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    /**memcpy( tpt, (char *)field + compt->sitelist[0][j]*size,
			size);**/
		    pvm_pkbyte( (char *)field + compt->sitelist[0][j]*size,
			size, 1 );
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    /**memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
			size);**/
		    pvm_pkbyte( (char *)field + compt->sitelist[subl_s][j]*size,
			size, 1 );
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		/**memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
		    size);**/
		pvm_pkbyte( (char *)field + compt->sitelist[subl_s][j]*size,
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[mbuf[i].msg_node],
	    GATHER_MSG_TYPE(mynode(),type_offset) );
	i++;
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* return */
    return(mbuf);
} /*start_gather_from_temp()*/

/* Repeat a gather with the same source and destination as a
  previous gather.  The previous gather must have been waited for
  but not cleaned up.  Pointers to sites on the same node are not
  reset, and the same buffers are reused. */
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
	mbuf[i].msg_OK = FALSE;
	mbuf[i].msg_id = GATHER_MSG_TYPE(mbuf[i].msg_node,type_offset);
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
	mbuf[i].msg_id = GATHER_MSG_TYPE(mynode(),type_offset);
	mbuf[i].msg_OK = FALSE;
	pvm_initsend(PvmDataRaw);
	/* gather data into the buffer */
	if( subl == EVENANDODD ){
	    if( compt->send_subl[0] == NOWHERE ){
		/* All is done by "subl 0" */
		for(j=0;j<compt->n_connected;j++,tpt += size){
		    /**memcpy( tpt, (char *)field + compt->sitelist[0][j]*size,
			size);**/
		    pvm_pkbyte( (char *)field + compt->sitelist[0][j]*size,
			size, 1 );
		}
	    }
	    else for(subl_r=0;subl_r<N_SUBL32;subl_r++){
		subl_s = compt->send_subl[subl_r];
		for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		    /**memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
			size);**/
		    pvm_pkbyte( (char *)field + compt->sitelist[subl_s][j]*size,
			size, 1 );
		}
	    }
	}
	else{
	    subl_s = compt->send_subl[subl];
	    for(j=0;j<compt->n_subl_connected[subl_s];j++,tpt += size){
		/**memcpy( tpt, (char *)field + compt->sitelist[subl_s][j]*size,
		    size);**/
		pvm_pkbyte( (char *)field + compt->sitelist[subl_s][j]*size,
		    size, 1 );
	    }
	}
	/* start the send */
	/**pvm_pkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );**/
	pvm_send( pvm_node_tid[compt->othernode], 
	    GATHER_MSG_TYPE(mynode(),type_offset));
	i++;
    }
} /*restart_gather_from_temp() */

#endif	/* N_SUBL32 */


void wait_gather(msg_tag *mbuf) {
  /* For asynchronous receives with random order of arrival */
  register int i,flag;
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
	      flag = pvm_nrecv( ANY_NODE, mbuf[i].msg_id );
	      if( flag > 0 ){
		  pvm_upkbyte( (char *)mbuf[i].msg_buf, mbuf[i].msg_size, 1 );
		  /* Check off this receive */
		  mbuf[i].msg_OK = TRUE;
	      }
	      else if ( flag < 0 ){	/* some PVM error */
		  printf("PVM error %d, node %d\n",flag,mynode() );
		  fflush(stdout);
		  terminate(-1);
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
struct msg_tmp { int node,count; }; /* temporary structure for keeping
	track of messages to be sent or received */
int g_gather_flag=0;	/* flag to tell if general gather in progress */
struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
int tsize;		/* size of entry in messages = sizeof(int)+size */
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
	pvm_send( pvm_node_tid[to_nodes[i].node],
	    GENERAL_GATHER_MSG_TYPE(this_node));
    }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
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
	mbuf[i+n_recv_msgs+1].msg_id = GENERAL_GATHER_MSG_TYPE(this_node);
	pvm_initsend( PvmDataRaw );
	pvm_pkbyte( (char *)mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, 1 );
	pvm_send( pvm_node_tid[to_nodes[i].node],
	    GENERAL_GATHER_MSG_TYPE(this_node));
    }
    /* terminate list */
    if(mbuf != NULL){
	mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
	mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

#endif	/* N_SUBL32 */


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
printf("Oops: tried to use field_pointer()\n"); terminate(1);
return(NULL);
}

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
return(NULL);
}

/* free any buffers allocated by above routines.  Argument is the
address returned by field_pointer...() */
void cleanup_field_pointer(char * buf) {
    /* if buf is an address in the lattice, leave it alone.  Otherwise
	it was created by malloc and should be freed */
    if( buf < (char *)lattice || buf >= (char *)(lattice+sites_on_node))
	free(buf);
}

/* SEND AND RECEIVE FIELD */
void send_field(char *buf, int size, int tonode) {

  pvm_initsend( PvmDataRaw );
  pvm_pkbyte( (char *)buf,size,1);
  pvm_send(pvm_node_tid[tonode],FIELD_TYPE);
}
void get_field(char *buf, int size){

  pvm_recv( ANY_NODE, FIELD_TYPE);
  pvm_upkbyte( (char *)buf, size, 1 );
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="pvm";
char * machine_type(){
    return(name);
}

/* Return my node number */
int mynode()
{
  return(my_logical_node );
}

/* Return number of nodes */
int numnodes(){
  return(pvm_number_nodes );
}

/* Synchronize all nodes */
void g_sync(){
  int i;
  /* Each node sends SYNC_TYPE message to node 0, waits to get
     reply */

  if(my_logical_node==0){
    for(i=1;i<pvm_number_nodes;i++){
      pvm_recv( pvm_node_tid[i], SYNC_TYPE );
    }
    pvm_initsend( PvmDataRaw );
    pvm_mcast( &(pvm_node_tid[1]), pvm_number_nodes-1, SYNC_TYPE );
  }
  else{  /* generic node != 0 */
      pvm_initsend( PvmDataRaw );
      pvm_send( pvm_node_tid[0], SYNC_TYPE );
      pvm_recv( pvm_node_tid[0], SYNC_TYPE );
  }
}

/* Sum Real over all nodes */
void g_floatsum( Real *fpt) {
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
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *fpt;	/* contribution from node 0 */
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
void g_doublesum( double *dpt) {
  double work,sum;
  int othernode,type;

  type=SUM_DOUBLE_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkdouble( dpt,1,1 );
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *dpt;	/* contribution from node 0 */
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
void g_vecdoublesum( double *dpt, int ndoubles) {
  double *work,*sum;
  int othernode,type,size,i;

  size = ndoubles * sizeof(double);
  work = (double *)malloc(size);
  sum = (double *)malloc(size);
  type = SUM_VECDOUBLE_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkbyte( (char *)dpt, size, 1);
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	for(i=0;i<ndoubles;i++)sum[i] = dpt[i];	/* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkbyte( (char *)work, size, 1);
	  for(i=0;i<ndoubles;i++)sum[i] += work[i];
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_bytes((char *)sum, size);

  /* Answer goes in dpt */
  for(i=0;i<ndoubles;i++)dpt[i]=sum[i];
  free(work);
  free(sum);
}

/* Sum complex over all nodes */
void g_complexsum( complex *cpt) {
  complex work,sum;
  int othernode,type;

    type=SUM_COMPLEX_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkcplxHELP( (Real *)cpt, 1, 1 );
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *cpt;	/* contribution from node 0 */
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

/* Sum double_complex over all nodes */
void g_dcomplexsum( double_complex *cpt) {
  double_complex work,sum;
  int othernode,type;

    type=SUM_DCOMPLEX_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkdcplx( (double *)cpt, 1, 1 );
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	sum = *cpt;	/* contribution from node 0 */
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
void g_veccomplexsum( complex *cpt, int ncomplex){
  complex *work,*sum;
  int othernode,type,size,i;

  size = ncomplex * sizeof(complex);
  work = (complex *)malloc(size);
  sum = (complex *)malloc(size);
  type = SUM_VECDCOMPLEX_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkbyte( (char *)cpt, size, 1);
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero sums up all the numbers */
	for(i=0;i<ncomplex;i++)sum[i] = cpt[i];	/* contribution from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkbyte( (char *)work, size, 1);
	  for(i=0;i<ncomplex;i++)CSUM(sum[i], work[i]);
	}
    }

  /* Broadcast sum from node zero to all others */
  broadcast_bytes((char *)sum, size);

  /* Answer goes in cpt */
  for(i=0;i<ncomplex;i++)cpt[i]=sum[i];
  free(work);
  free(sum);
}

/* Sum wilson_vector over all nodes */
void g_wvectosumReal( wilson_vector *wvpt) {
  g_veccomplexsum((complex *)wvpt, 12);
}

/* 32 bit exclusive or over all nodes */
void g_xor32( u_int32type *pt ) {
  u_int32type work,result;
  int othernode,type;
  
    type=SEND_INTEGER_TYPE;
    /* For the moment, do this the dumb way.  Node 0 does all the work. */
    if(mynode()!=0){
	/* send my value to node 0 */
      pvm_initsend( PvmDataRaw );
      pvm_pkbyte( (char *)pt, sizeof(u_int32type), 1);
      pvm_send( pvm_node_tid[0], type);
    }
    else{
	/* node zero finds the result */
	result = *pt; /*  from node 0 */
	for(othernode=1;othernode<numnodes();othernode++){
	  /* This receive does not specify the node of the sender */
	  /* Thus it waits until exactly numnodes()-1 values are received */
	  pvm_recv( ANY_NODE, type);
	  pvm_upkbyte( (char *)&work, sizeof(u_int32type), 1);
	  result ^= work;
	}
    }

  /* Broadcast result from node zero to all others */
  broadcast_bytes((char *)&result, sizeof(u_int32type));

  /* Answer goes in pt */
  *pt = result;
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
      pvm_send( pvm_node_tid[0], type);
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
      pvm_send( pvm_node_tid[0], type);
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
  
  if(mynode()==0)
    {
      pvm_initsend( PvmDataRaw );
      pvm_pkReal( fpt, 1,1 );
      for(i = 0; i < pvm_number_nodes; i++) if(i != mynode())
	pvm_send( pvm_node_tid[i], BROADCAST_MILC_REAL_TYPE);
    }

  else
    {
      pvm_recv( ANY_NODE, BROADCAST_MILC_REAL_TYPE);
      pvm_upkReal( fpt,1,1 );
    }
}

/* Broadcast double precision Realing point number from node zero */
void broadcast_double(double *dpt) {
  int i;

  if(mynode()==0)
    {
      pvm_initsend( PvmDataRaw );
      pvm_pkdouble( dpt, 1, 1 );
      for(i = 0; i < pvm_number_nodes; i++) if(i != mynode())
	pvm_send( pvm_node_tid[i], BROADCAST_DOUBLE_TYPE);
    }

  else
    {
      pvm_recv( ANY_NODE, BROADCAST_DOUBLE_TYPE);
      pvm_upkdouble( dpt,1,1 );
    }
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(complex *cpt) {
  int i;

  if(mynode()==0)
    {
      pvm_initsend( PvmDataRaw );
      pvm_pkcplxHELP( (Real *)cpt, 1, 1);
      for(i = 0; i < pvm_number_nodes; i++) if(i != mynode())
	pvm_send( pvm_node_tid[i], BROADCAST_COMPLEX_TYPE);
    }

  else
    {
      pvm_recv( ANY_NODE, BROADCAST_COMPLEX_TYPE);
      pvm_upkcplxHELP( (Real *)cpt,1,1 );
    }
}

/* Broadcast double precision complex number from node zero */
void broadcast_dcomplex(double_complex *cpt) {
  int i;

  if(mynode()==0)
    {
      pvm_initsend( PvmDataRaw );
      pvm_pkdcplx( (double *)cpt, 1, 1);
      for(i = 0; i < pvm_number_nodes; i++) if(i != mynode())
	pvm_send( pvm_node_tid[i], BROADCAST_DCOMPLEX_TYPE);
    }

  else
    {
      pvm_recv( ANY_NODE, BROADCAST_DCOMPLEX_TYPE);
      pvm_upkdcplx( (double *)cpt,1,1 );
    }
}

/* Broadcast bytes */
void broadcast_bytes(char *buf, int size) {
  int i;

  if(mynode()==0)
    {
      pvm_initsend( PvmDataRaw );
      pvm_pkbyte( (char *)buf, size, 1);
      for(i = 0; i < pvm_number_nodes; i++) if(i != mynode())
	pvm_send( pvm_node_tid[i], BROADCAST_BYTES_TYPE);
    }
  else
    {
      pvm_recv( ANY_NODE, BROADCAST_BYTES_TYPE);
      pvm_upkbyte( (char *)buf, size ,1 );
    }
}

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(int tonode, int *address) {

  pvm_initsend( PvmDataRaw );
  pvm_pkint( (int *)address, 1, 1 );
  pvm_send( pvm_node_tid[tonode], SEND_INTEGER_TYPE);

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
  for( i = 0; i < pvm_number_nodes; i++){
    if(i != my_logical_node)pvm_kill(pvm_node_tid[i]);
  }
  pvm_exit();
  exit(status);
}

