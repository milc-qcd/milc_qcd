#ifndef QPHIXJ_INT_H
#define QPHIXJ_INT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  QPHIXJ_SUCCESS = 0,
  QPHIXJ_FAIL = 1,
  QPHIXJ_MEM_ERROR, /* memory errors such as failure to allocate, etc */
} QPHIXJ_status_t;

typedef enum {
  QPHIXJ_EVEN = 0x02,
  QPHIXJ_ODD = 0x01,
  QPHIXJ_EVENODD = 0x03
} QPHIXJ_evenodd_t;

typedef float QPHIXJ_F_Real;
typedef double QPHIXJ_D_Real;

typedef struct {
  int (*node_number)(const int coords[]); /* node no for given latt coord */
  int (*node_index)(const int coords[]);  /* site rank for given latt coord */
  int latdim;                             /* number of lattice dimensions */
  int *latsize;                           /* physical lattice lengths */
  int machdim;                            /* number of logical machine dims */
  int *machsize;                          /* logical grid lengths */
  int this_node;                          /* lexicographic node number */
  int sites_on_node;
  int even_sites_on_node;                 /* If needed */
} QPHIXJ_layout_t;
#define QPHIXJ_LAYOUT_ZERO ((QPHIXJ_layout_t){NULL,NULL,0,NULL,0,NULL,0,0})

typedef struct {
  int By;
  int Bz;
  int PadXY;
  int PadXYZ;
  int NCores;
  int Sy;
  int Sz;
  int MinCt;
  int compress12;
} QPHIXJ_vec_layout_t;
  //#define QPHIXJ_VEC_LAYOUT_DEFAULT ((QPHIXJ_vec_layout_t){4,4,0,0,1,1,2,1,0})
#define QPHIXJ_VEC_LAYOUT_DEFAULT ((QPHIXJ_vec_layout_t){4,4,0,0,1,1,1,1,0})

typedef struct {
  double final_sec;        /* (out) total number of seconds used */
  double final_flop;       /* (out) total number of flops performed */
  QPHIXJ_status_t status;     /* (out) error status */
  int count1, count2;      /* (out) generic counters */
} QPHIXJ_info_t;
#define QPHIXJ_INFO_ZERO ((QPHIXJ_info_t){0,0,QPHIXJ_SUCCESS,0,0})

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  int max;                  /* (in) max number of iterations */
  int restart;              /* (in) number of iterations before restart */
  int nrestart;             /* (in) number of restarts allowed */
  QPHIXJ_evenodd_t evenodd;  /* (in) subset of source vector */
} QPHIXJ_invert_arg_t;
#define QPHIXJ_INVERT_ARG_DEFAULT ((QPHIXJ_invert_arg_t){2000,1000,5,QPHIXJ_EVENODD})

  /* these are quantities that vary for each mass in the multi inverter */
typedef struct {
  double resid;          /* (in) desired squared residual. Ignored if 0 */
  double final_rsq;      /* (out) actual squared residual */
  double relresid;       /* (in) desired squared relative norm. Ignored if 0 */
  double final_rel;      /* (out) actual squared relative norm */
  double size_r;         /* resulting cumulative residual. Same
			    normalization as final_rsq. */
  double size_relr;      /* resulting cumulative relative
			    residual. Same normalization as
			    final_rel. */
  int final_iter;        /* (out) number of iterations done */
  int final_restart;     /* (out) number of restarts done */
} QPHIXJ_resid_arg_t;
#define QPHIXJ_RESID_ARG_DEFAULT ((QPHIXJ_resid_arg_t){1e-6,0,0,0,0,0})

  /**********************/
  /*  General routines  */
  /**********************/

void QPHIXJ_init(QPHIXJ_layout_t *layout, QPHIXJ_vec_layout_t *params);
void QPHIXJ_finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* QPHIXJ_INT_H */
