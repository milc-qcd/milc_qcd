#include "../../include/int32type.h"

typedef struct {
  int spect;
  int other;
  int mom;
  int oper;
  int copy;
  Real wt;
} twopt_select;

typedef struct {
  int forwback;
  int nfile;
  char filename[MAX_NO_FILE][128];
  int nselect;
  twopt_select select[MAX_SELECT_PER_FILE];
} two_list;

typedef struct {
  char filename[128];
  int nselect;
  twopt_select select[MAX_SELECT_PER_FILE];
} twopt_oneselect;

typedef struct {
  int spect;
  int zonked;
  int seq;
  int q;
  int p;
  int oper;
  int copy;
  Real wt;
} threept_select;

typedef struct {
  char filename[128];                  /* Output dump file name */
  int nselect;
  threept_select select[MAX_SELECT_PER_FILE];
} threept_oneselect;

typedef struct {
  int forwback;
  int nfile;
  char filename[MAX_NO_FILE][128];
  int nselect;
  threept_select select[MAX_SELECT_PER_FILE];
} three_list;

enum forwback_type { FORWARD, BACKWARD, FOLD };

void read_input_param(
  char filename[80],           /* Name of parameter file */
  three_list *threept,         /* List of three-pt functions */
  two_list *twopt_recoil,      /* List of recoil meson propagators */
  two_list *twopt_sequential   /* List of incoming B meson propagators */
);

void read_propagating_form_corr(
  complex **corr,           /* the three point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies             /* rotations */
);

void read_3pt_onemom(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies,            /* rotations */
  int p_mom_select          /* Momentum index selected */
);

void read_merge1_form_corr(
  complex **corr,           /* the three point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies             /* rotations */
);

void read_merge2_form_corr(
  complex **corr,           /* the three point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies             /* rotations */
);


void read_prop_twopt(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  char filename[80],        /* the name of the disk file */
  int32type **q_momstore,      /* list of momenta */
  int *hl_flag,             /* flag identifying type of correlator */
  int *no_k_one ,           /* the number of  kappa values */
  int *no_k_two,            /* the number of  kappa values */
  int *numer_of_operators,  /* number of operators */
  int *nt ,                 /* time values */
  int *no_q_values,         /* momentum values */
  int *nocopies             /* rotations */
  );

void read_select_twopt(
  complex **corr,      /* pointer to complex array: the two point correlator */
  char filename[80],   /* the name of the disk file */
  int nselect,         /* Number of selection parameter sets */
  twopt_select *twoselect, /* Selection parameter sets */
  int32type **q_momselect,   /* momentum list */		       
  int *hl_flag,        /* identifies meson source/sink smearing */
  int *nt              /* number of time steps */
  );

void read_multiselect_twopt_param(
  int *nfile,             /* Number of files to write */
  twopt_oneselect *fp  /* two point selection parameter sets */
);

void write_prop_twopt(
  complex *corr,            /* the two point correlators */
  char filename[80],        /* the name of the disk file */
  int32type *q_momstore,       /* list of momenta */
  int hl_flag,              /* flag identifying type of correlator */
  int no_spectator,         /* the number of  kappa values */
  int no_zonked,            /* the number of  kappa values */
  int number_of_operators,  /* number of operators */
  int nt ,                  /* time values */
  int no_q_values,          /* momentum values */
  int no_copies             /* rotations */
  );

void write_prop_twopt_onemom(
  complex *corr,            /* the two point correlator */
  char filename[80],        /* the name of the output file */
  int q_mom_select,         /* Index of selected momentum */
  int32type *q_momstore,       /* Array of momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked,            /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_q_corr_values,     /* Number of momentum values in corr */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  );

void write_prop_twopt_onespect(
  complex *corr,            /* the two point correlator */
  char filename[80],        /* the name of the output file */
  int spect_select,         /* Index of selected spectator kappa */
  int32type *q_momstore,       /* Array of momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,           /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_q_corr_values,     /* Number of momentum values in corr */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  );

void write_3pt(
  complex *corr,            /* the two point correlator */
  char filename[80],        /* the name of the output file */
  int32type *q_momstore,       /* Array of q momentum values */
  int32type *p_momstore,       /* Array of p momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,           /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_sequential,        /* Sequential quarks */
  int no_p_values,          /* Number of p momentum values in corr */
  int no_q_values,          /* Number of q momentum values */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  );

void write_3pt_onemom(
  complex *corr,            /* the three point correlator */
  char filename[80],        /* the name of the output file */
  int p_mom_select,         /* Index of selected momentum */
  int32type *q_momstore,       /* Array of q momentum values */
  int32type *p_momstore,       /* Array of p momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked,            /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_sequential,        /* Sequential quarks */
  int no_p_corr_values,     /* Number of p momentum values in corr */
  int no_q_values,          /* Number of q momentum values */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  );


void write_3pt_onespect(
  complex *corr,            /* the three point correlator */
  char filename[80],        /* the name of the output file */
  int spect_select,         /* Index of selected spectator */
  int32type *q_momstore,       /* Array of q momentum values */
  int32type *p_momstore,       /* Array of p momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,           /* Zonked quarks */
  int no_spectator_corr,    /* Spectator quarks */
  int no_sequential,        /* Sequential quarks */
  int no_p_values,          /* Number of p momentum values in corr */
  int no_q_values,          /* Number of q momentum values */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  );

void write_data(Real *data,Real *data_err, int dim ,char filename[80] );

void read_select_twopt_param(
  int *nselect,            /* Number of selection sets */
  twopt_select *twoselect  /* two point selection parameter sets */
);

void read_select_form_param(
  int *nselect,                 /* Number of selection sets */
  threept_select *threeselect   /* Selection parameter sets */
);

void read_multiselect_form_param(
  int *nfile,             /* Number of files to write */
  threept_oneselect *fp   /* Selection parameter set for each file */
);

void read_select_form_corr(
  complex **corr,      /* pointer to complex array: the two point correlator */
  char filename[80],   /* the name of the disk file */
  int nselect,         /* Number of sets of selection parameters */
  threept_select *threeselect, /* Selection parameter sets */
  int32type **p_momselect,   /* momentum list */	
  int32type **q_momselect,   /* momentum list */	
  int *hl_flag,        /* identifies meson source/sink smearing */
  int *ntime           /* number of time steps */
  );

void read_multiselect_form_corr(
  char filename[80],   /* the name of the input correlator file */
  int nfile,           /* Number of output files */
  threept_oneselect *fp /* Selection parameter sets */
  );

void read_multiselect_twopt(
  char filename[80],   /* the name of the input correlator file */
  int nfile,           /* Number of output files */
  twopt_oneselect *fp /* Selection parameter sets */
  );

void dump_twopt(complex *corr,  int corr_oper_list[], int corr_copy_list[],
  int no_zonked, int no_spectator, 
  int no_q_values, int no_oper, int nt, int nocopies) ;

void dump_prop_form(complex *corr,  int corr_oper_list[],
      int corr_copy_list[],
      int no_zonked, int no_spectator, int no_sequential,
      int no_q_values, int no_p_values, int no_oper, int nt, int nocopies);

int bsd_sum (char *data,int total_bytes) ;
void byte_rev_array(int32type buf[], int words) ;

void titles_hl(int nt, int no_p_mom, int no_q_mom, int no_oper,
int32type  *p_momentum, int32type  *q_momentum, int hl_flag);

void titles_2pt(int nt, int no_q_mom, int no_oper,int32type *q_momentum, int hl_flag, int no_copies);

char *three_oper_name(int n, int copy ) ;
char *two_oper_name(int n, int copy ) ;
int oper_chirality(int oper) ;
int oper_phase(int oper);

complex jack_sample_twopt(
  int jomit, 
  int nosample, 
  complex *two_corr[MAX_NO_FILE ] ,
  int t);

complex jack_sample_threept(
  int jomit, 
  int nosample, 
  complex *three_corr[MAX_NO_FILE ] ,
  int t);

Real jackmean(Real *data,int nodata,int jomit);
Real jackerror(Real mean,Real *data,int nodata);
