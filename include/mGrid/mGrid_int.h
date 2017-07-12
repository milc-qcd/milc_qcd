#ifndef _MGRID_INT_H
#define _MGRID_INT_H

#ifdef __cplusplus
extern "C" {
#endif

/* This file defines miscellaneous types that may be used in all compilations */
/* It is included in mGrid.h so it doesn't need to be included separately */

typedef enum {
  GRID_SUCCESS = 0,
  GRID_FAIL = 1,
  GRID_MEM_ERROR, /* memory errors such as failure to allocate, etc */
} GRID_status_t;

typedef enum {
  GRID_EVEN = 0,
  GRID_ODD = 1,
  GRID_EVENODD = 2
} GRID_evenodd_t;

typedef float GRID_F_Real;
typedef double GRID_D_Real;
//typedef Real  GRID_Real;

#define GRID_LAYOUT_ZERO ((GRID_layout_t){NULL,NULL,0,NULL,0,NULL,0,0,0,0})

typedef struct {
  double final_sec;        /* (out) total number of seconds used */
  double final_flop;       /* (out) total number of flops performed */
  GRID_status_t status;     /* (out) error status */
  int count1, count2;      /* (out) generic counters */
} GRID_info_t;
#define GRID_INFO_ZERO ((GRID_info_t){0,0,GRID_SUCCESS,0,0})

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  int max;			/* (in) max number of iterations */
  int restart;			/* (in) number of iterations before restart */
  int nrestart;			/* (in) number of restarts allowed */
  GRID_evenodd_t parity;	/* (in) subset of source vector */
} GRID_invert_arg_t;
#define GRID_INVERT_ARG_DEFAULT ((GRID_invert_arg_t){2000,1000,5,GRID_EVENODD})

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
} GRID_resid_arg_t;
#define GRID_RESID_ARG_DEFAULT ((GRID_resid_arg_t){1e-6,0,0,0,0,0})

  /**********************/
  /*  General routines  */
  /**********************/

int GRID_verbose(int level);
int GRID_profcontrol(int level);

/*
GRID_status_t GRID_asqtad_invert_set_opts(GRID_opt_t opts[], int nopts);
GRID_status_t GRID_asqtad_force_set_opts(GRID_opt_t opts[], int nopts);

GRID_status_t GRID_hisq_links_set_opts(GRID_opt_t opts[], int nopts);
GRID_status_t GRID_hisq_force_set_opts(GRID_opt_t opts[], int nopts);

GRID_status_t GRID_wilson_invert_set_opts(GRID_opt_t opts[], int nopts);
GRID_status_t GRID_wilson_force_set_opts(GRID_opt_t opts[], int nopts);

GRID_status_t GRID_dw_invert_set_opts(GRID_opt_t opts[], int nopts);
GRID_status_t GRID_dw_force_set_opts(GRID_opt_t opts[], int nopts);
*/
#ifdef __cplusplus
}
#endif

#endif
