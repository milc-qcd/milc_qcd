#ifndef _MGRID_INT_H
#define _MGRID_INT_H

#ifdef __cplusplus
extern "C" {
#endif

// The following strategem is intended to provide a strong chain of typedefs for both C and C++ compilation
// Above the API the data types are opaque (incomplete).
// If mGrid_internal.h is used, it must precede this header file

#if defined( _MGRID_INTERNAL_H) && defined( __cplusplus)
#define STRUCT
#else
#define STRUCT struct
#endif

  /**********************************************/
  /*  Opaque container for a 5D Grid definition */
  /**********************************************/

typedef STRUCT GRID_4Dgrid_struct   GRID_4Dgrid;
typedef STRUCT GRID_4DRBgrid_struct GRID_4DRBgrid;
typedef STRUCT GRID_5Dgrid_struct   GRID_5Dgrid;
typedef STRUCT GRID_5DRBgrid_struct GRID_5DRBgrid;

GRID_4Dgrid *GRID_create_grid(void);
GRID_4DRBgrid *GRID_create_RBgrid(GRID_4Dgrid *full_grid);
GRID_5Dgrid *GRID_create_5Dgrid(int n, GRID_4Dgrid *full_grid);
GRID_5DRBgrid *GRID_create_5DRBgrid(int n, GRID_4Dgrid *full_grid);
void GRID_destroy_4Dgrid(GRID_4Dgrid *g);
void GRID_destroy_4DRBgrid(GRID_4DRBgrid *g);
void GRID_destroy_5Dgrid(GRID_5Dgrid *g);
void GRID_destroy_5DRBgrid(GRID_5DRBgrid *g);

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

#define GRID_LAYOUT_ZERO ((GRID_layout_t){NULL,NULL,0,NULL,0,NULL,0,0,0,0})

typedef struct {
  double final_sec;        /* (out) number of seconds in inverter */
  double misc_sec;         /* (out) number of seconds instantiating operators */
  double final_flop;       /* (out) total number of flops performed */
  GRID_status_t status;     /* (out) error status */
  int count1, count2;      /* (out) generic counters */
} GRID_info_t;
#define GRID_INFO_ZERO ((GRID_info_t){0,0,0,GRID_SUCCESS,0,0})

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  int max;			/* (in) max number of iterations */
  int maxInner;			/* (in) max number of iterations for inner solve (mixed prec) */
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

/* Parameters for Chebyshev polynomial */
typedef struct
{
  double alpha;
  double beta;
  int Npoly;
} GRID_ChebyParams;

/* Grid diagonalization algorithms */
typedef enum
{
  GRID_IRLdiagonaliseWithDSTEGR,
  GRID_IRLdiagonaliseWithQR,
  GRID_IRLdiagonaliseWithEigen
} GRID_IRLdiagonalisation;

/* Parameters for Grid implicitly restarted Lanczos */
typedef struct
{
  GRID_evenodd_t parity;
  int Nstop;
  int Nk;
  int Nm;
  double tol;
  int maxIter;
  int restartMin;
  int reorth_period;
  GRID_ChebyParams chebyParams;
  GRID_IRLdiagonalisation diag;
} GRID_eig_arg_t;

/* Parameters for the HISQ links and force */

/* Only GRID_HISQ_REUNIT_U3 is currently supported */
enum GRID_hisq_unitarize_group_t {
  GRID_HISQ_REUNIT_U3,
  GRID_HISQ_REUNIT_SU3,
  GRID_HISQ_REUNIT_SVD_REL_ERROR,
  GRID_HISQ_REUNIT_SVD_ABS_ERROR
};

/* Only GRID_HISQ_UNITARIZE_ANALYTIC is currently supported */
enum GRID_hisq_unitarize_method_t {
  GRID_HISQ_UNITARIZE_NONE,
  GRID_HISQ_UNITARIZE_APE,
  GRID_HISQ_UNITARIZE_ROOT,
  GRID_HISQ_UNITARIZE_RATIONAL,
  GRID_HISQ_UNITARIZE_HISQ,
  GRID_HISQ_UNITARIZE_ANALYTIC,
  GRID_HISQ_UNITARIZE_STOUT,
  GRID_HISQ_UNITARIZE_ITERATIVE,
  GRID_HISQ_REUNIT_ALLOW_SVD
  };

typedef  struct {
    enum GRID_hisq_unitarize_group_t ugroup;
    enum GRID_hisq_unitarize_method_t umethod;
  } GRID_hisq_coeffs_t;

#define GRID_HISQ_COEFFS_DEFAULT {GRID_HISQ_REUNIT_U3, GRID_HISQ_UNITARIZE_ANALYTIC}
  
  
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
