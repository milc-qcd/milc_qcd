/* Added for pvm */
#ifdef PVM
#define ANY_MSG -1              /* Any message */

/* Message structures for communication between host and nodes */
/* The first two fields must always be the same for each type */
/* The basic structure must be the shortest */

/* For most messages */
struct hcs_basic {
  int msg_type;
  int node;	/* integer identifies caller's instance */
  int arg1,arg2,arg3;  /* Use depends on which routine */
} ;
#define HCS_BASIC_SIZE (sizeof(struct hcs_basic))
int put_hcs_basic( struct hcs_basic * hcs );
int get_hcs_basic( struct hcs_basic * hcs );

/* For printf and scanf calls */
#define STRING_LENGTH 256
struct hcs_stdio {
  int msg_type;
  int node;	/* integer identifies caller's instance */
  int length;
  char s[STRING_LENGTH];
} ;
#define HCS_STDIO_SIZE (sizeof(struct hcs_stdio))
int put_hcs_stdio( struct hcs_stdio * hcs );
int get_hcs_stdio( struct hcs_stdio * hcs );

/* For initialization call */
#define MAX_NUMBER_NODES 512
#define HOST_NAME_LENGTH 128
struct hcs_ident {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_instance[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident))
int put_hcs_ident( struct hcs_ident * hcs );
int get_hcs_ident( struct hcs_ident * hcs );

union {
  struct hcs_basic basic;
  struct hcs_stdio stdio;
  struct hcs_ident ident;
} hcs;

#define HOST_CALL 77	/* pvm message type for call to host for service */
#define HOST_REPLY 87

/* Message subtypes internal to host-node service calls */
/* (Not used by pvm to identify messages) */
#define PRINTF_HOST_CALL 11
#define SCANF_HOST_CALL 12
#define FPRINTF_HOST_CALL 13 /* Not used */ 
#define FSCANF_HOST_CALL 14 /* Not used */
#define FFLUSH_HOST_CALL 20
#define FOPEN_HOST_CALL 30 /* Not used */
#define FCLOSE_HOST_CALL 31 /* Not used */
#define OPEN_HOST_CALL 32 /* Not used */
#define CLOSE_HOST_CALL 33 /* Not used */
#define CREAT_HOST_CALL 34 /* Not used */
#define READ_HOST_CALL 40 /* Not used */
#define WRITE_HOST_CALL 41 /* Not used */
#define NODE_IDENT_CALL 78	
#define NODES_DONE_HOST_CALL 99

/* end of pvm additions */
#endif	/* end ifdef PVM */

/* Added for pvm */
#ifdef PVM24
#define ANY_MSG -1              /* Any message */

/* Message structures for communication between host and nodes */
/* The first two fields must always be the same for each type */
/* The basic structure must be the shortest */

/* For initialization call */
#define MAX_NUMBER_NODES 512
#define HOST_NAME_LENGTH 128
struct hcs_ident_struct {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_instance[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} hcs_ident ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident_struct))
int put_hcs_ident( struct hcs_ident_struct * hcs );
int get_hcs_ident( struct hcs_ident_struct * hcs );

#define NODE_IDENT_CALL 77	
#define terminate g_terminate   /* Because of name conflict */

/* end of pvm version 2.4 additions */
#endif	/* end ifdef PVM24 */


/* Added for pvm */
#ifdef PVM3
#define ANY_MSG -1              /* Any message */
#define ANY_NODE -1              /* Any node */

/* Message structures for communication between host and nodes */

/* For initialization call */
#define MAX_NUMBER_NODES 512
#define HOST_NAME_LENGTH 128
struct hcs_ident_struct {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_tid[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} hcs_ident ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident_struct))
int put_hcs_ident( struct hcs_ident_struct * hcs );
int get_hcs_ident( struct hcs_ident_struct * hcs );

#define NODE_IDENT_CALL 77	

/* end of pvm version 3 additions */
#endif	/* end ifdef PVM3 */

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
#if defined(PVM) || defined(PVM24) || defined(PVM3) || defined(MPL)
    int msg_OK;
        /* flag to track the asynchronous arrival of messages */
#endif
#ifdef MPL
    int mpl_msgid;
        /* MPL assigned message id for checking status with mpc_status */
#endif
} msg_tag;

#ifdef PVM
#ifndef HOST_CODE
#define terminate g_terminate   /* Because of name conflict */
#define scanf host_scanf
#define printf host_printf
#define fflush host_fflush
#endif
#endif

