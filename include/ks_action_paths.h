#ifndef _KS_ACTION_PATHS_H
#define _KS_ACTION_PATHS_H

#include "../include/complex.h"

/* Structure specifying each rotation and reflection of each kind of
	path.  */
#define MAX_PATH_LENGTH 16
typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  Real coeff;	        /* coefficient, including minus sign if backwards */
  Real forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

/* Structure defining the fermion action using paths or optimized
   coefficients */

#define MAX_NAIK 10 // max number of quarks that require Naik epsilon correction
                    // without c quark this is normally 1, however this
                    // constant should be set to at least 2

typedef struct {
  Real one_link;
  Real naik;
  Real three_staple;
  Real five_staple;
  Real seven_staple;
  Real lepage;
} asqtad_coeffs_t;

typedef struct {
  asqtad_coeffs_t act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;                   /* For all actions */
  Q_path *q_paths;                   /* For all actions */
} ks_component_paths;

typedef struct {
  int n_naiks;
  double eps_naik[MAX_NAIK];
  ks_component_paths p1, p2, p3;
  int umethod;
  int ugroup;
  int constructed;       /* Boolean */
} ks_action_paths_hisq;

typedef struct {
  ks_component_paths p;
  int constructed;         /* Boolean */
} ks_action_paths;


/* ks_action_paths.c */
ks_action_paths *create_path_table(void);
void destroy_path_table(ks_action_paths *ap);
int make_path_table(ks_action_paths *ap, ks_action_paths *ap_dmdu0);
char *get_ap_string(ks_action_paths *ap);

/* ks_action_paths_hisq.c */
ks_action_paths_hisq *create_path_table_hisq(void);
void destroy_path_table_hisq(ks_action_paths_hisq *ap);
void load_act_path_coeff_hisq(ks_action_paths_hisq *ap, int n_naiks,
			      double *eps_naik);
int make_path_table_hisq(ks_action_paths_hisq *ap, int n_naiks, double *eps_naik);
int get_n_naiks(ks_action_paths_hisq *ap);
double *get_eps_naik(ks_action_paths_hisq *ap);
int get_umethod(ks_action_paths_hisq *ap);
int get_ugroup(ks_action_paths_hisq *ap);
char *get_ap_string_hisq(ks_action_paths_hisq *ap);

#endif /* _KS_ACTION_PATHS_H */
