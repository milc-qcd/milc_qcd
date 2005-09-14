#ifndef RG_include_h
#define RG_include_h

#define NRG 3
#define RG_Nd 4
#define RG_Ncn 16 
#define RG_Nf 4 
#define RG_Ns 4 

#define node0_printf if(this_node==0)printf
#define FIX_GAUGE  /*Allow to fix the gauge */
/* Please note that the checks are better to be done without fixing the gauge */
//#define PRINT
//#define READ

//#define CHECK_TRACE_UNIT
//#define CHECK_TRACE
//#define CHECK_INV
//#define CHECK  


#ifdef CHECK
//#define CHECK_LINK
//#define CHECK_TRANS
//#define CHECK_PATH

#ifdef CHECK_PATH
 #define CHECK_DEGRAND_WO_SMEAR
#endif

#ifdef CHECK_LINK
//#define CHECK_PLAQ   
//#define CHECK_SMEAR_GAUGE  /* smearing with qdp: s(l)=gs(l) */
//#define CHECK_SMEAR_GAUGE_2 /* smearing with qdp: s(l) = sg(l) */
//#define CHECK_SMEAR_MILC_2  /* smearing with milc: s(l) = sg(l) */
/* Please note that both with qdp and mild sg(l) is equal to s(l) only
   up to 10e-6 */
 
//#define CHECK_SMEAR_QDP_MILC  /* Check smearing qdp equal to milc */


#define CHECK_DEGRAND  /* Check trick of DeGrand */
/* please note that with smearing the accuracy will be up to 10e-5 for t
he same reason why s(l) = sg(l) are equal up to 10e-5 */
#ifdef CHECK_DEGRAND
//#define CHECK_DEGRAND_WO_SMEAR 
#define CHECK_DEGRAND_W_SMEAR 
#endif
#endif
#endif


typedef struct
{
int *x;
int s[RG_Nd];
int rv;
}
shift_v;

typedef int RG_gamma;

typedef struct {
  RG_gamma column;
  int phase;
} st_row;

typedef struct {
  st_row row[16];
} st_matrix;

typedef struct {
  RG_gamma g;
  int sign;
} gamma_phase;

typedef struct {
  QDP_Subset sub;
  int fact;
} QDP_Sub_Block;

int intpow(int base, int exp);

int hyp_func(int x[], void *arg);

void check_cv(QLA_ColorVector *df, int coords[]);

void print_cv(QLA_ColorVector *df, int coords[]);

void print_cv_node(QLA_ColorVector *df, int coords[]);

void print_cv_node_tot(QLA_ColorVector *df, int coords[]);

void print_df(QLA_DiracFermion *df, int coords[]);

void print_gl(QLA_ColorMatrix *gl, int coords[]);

void check_gl(QLA_ColorMatrix *gl, int coords[]);

void print_dp(QLA_DiracPropagator *df, int coords[]);

void check_df(QLA_DiracFermion *df, int coords[]);

void point_V(QLA_ColorVector *s, int coords[]);


/******************** Operation for RGT blocks *******************/
void RG_create_block(QDP_Sub_Block *block, int n);

void SQDP_M_eq_M(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Sub_Block s);

void SQDP_M_peq_M(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Sub_Block s);

void SQDP_M_eq_c(QDP_ColorMatrix *a ,QLA_Complex *c, QDP_Sub_Block s);

void SQDP_M_eq_sM(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s);

void SQDP_M_eq_func(QDP_ColorMatrix *a, void (*func)(QLA_ColorMatrix *gl, int coord[]), QDP_Sub_Block s);

void SQDP_M_eq_zero(QDP_ColorMatrix *a, QDP_Sub_Block s);

void SQDP_M_eq_M_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_eq_M_times_Ma(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_peq_M_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_eq_Ma_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_eq_M_minus_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_eq_M_times_sM(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Shift shift, QDP_ShiftDir dir,QDP_Sub_Block s);

void SQDP_M_peq_M_times_sM(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Shift shift, QDP_ShiftDir dir,QDP_Sub_Block s);

void SQDP_M_eq_r_times_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_peq_r_times_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c, QDP_Sub_Block s);

void SQDP_M_eq_r_times_M_plus_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c,QDP_ColorMatrix *d,QDP_Sub_Block s);

void SQDP_R_eq_re_M_dot_M(QLA_Real *r, QDP_ColorMatrix *a, QDP_ColorMatrix *b,QDP_Sub_Block s);

void SQDP_V_eq_V(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Sub_Block s);

void SQDP_V_peq_V(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Sub_Block s);

void SQDP_V_eq_sV(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s);

void SQDP_V_eq_func(QDP_ColorVector *a, void (*func)(QLA_ColorVector *gl, int coord[]), QDP_Sub_Block s);

void SQDP_V_eq_zero(QDP_ColorVector *a, QDP_Sub_Block s);

void SQDP_V_eq_M_times_sV(QDP_ColorVector *a, QDP_ColorMatrix *b, QDP_ColorVector *c, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s);

void SQDP_V_peq_r_times_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_Sub_Block s);

void SQDP_V_eq_r_times_V_minus_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_ColorVector *c, QDP_Sub_Block s);

void SQDP_V_eq_Ma_times_V(QDP_ColorVector *a, QDP_ColorMatrix *b, QDP_ColorVector *c, QDP_Sub_Block s);

void SQDP_V_eq_r_times_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_Sub_Block s);

void SQDP_r_eq_re_M_dot_M(QLA_Real *r, QDP_ColorMatrix *a, QDP_ColorMatrix *b, QDP_Sub_Block s);
/******************** Operation for RGT blocks *******************/



void RG_check_subset(QDP_Sub_Block block[NRG+1]);

void RG_smear_dir (QDP_ColorMatrix *sm_link, QDP_ColorMatrix *link[], QLA_Real
w_l, QLA_Real w_s, QLA_Int dir, QDP_Sub_Block s, int len);

void RG_smearing_qdp (QDP_ColorMatrix *sm_link[], QDP_ColorMatrix *link[], QLA_Real *s_w, QLA_Real *l0, QDP_Sub_Block s, int len);

void project_qdp(QDP_ColorMatrix *link_src[],QDP_ColorMatrix *link_dest[],int *space_only);

void RG_check_hermicity(QDP_ColorMatrix *w[RG_Ncn],QDP_Sub_Block s);

void RG_smearing(QDP_ColorMatrix *res[RG_Nd], QDP_ColorMatrix *src[RG_Nd],QDP_Sub_Block s,int len);

void RG_setup(QDP_Sub_Block QDP_block[NRG+1],QDP_ColorMatrix *wlink[NRG][RG_Ncn]);

int RG_check_for_setup();

void RG_gauge(QDP_ColorMatrix *rg_link[NRG][RG_Nd], QDP_ColorMatrix *link[RG_Nd], QDP_Sub_Block s[NRG+1]);

void RG_create_gauge(QDP_ColorMatrix *res[RG_Nd], QDP_ColorMatrix *src[RG_Nd], QDP_Sub_Block s,int len);

void RG_transf_field(QDP_ColorVector *dest, QDP_ColorVector *src[RG_Ncn], int len);

void RG_inv_transf_field(QDP_ColorVector *dest, QDP_ColorVector *src, QDP_ColorMatrix *wlink[RG_Ncn],QDP_Sub_Block s,int len);

void RG_create_field(QDP_ColorVector *phi_c[RG_Ncn], QDP_ColorVector *phi, QDP_ColorMatrix *wlink[RG_Ncn],QDP_Sub_Block s);

void RG_coarse_to_fine(QDP_ColorVector *dest,QDP_Sub_Block QDP_block[NRG+1],QDP_ColorVector *src,QDP_ColorMatrix *wlink[NRG][RG_Ncn]);

void RG_fine_to_coarse(QDP_ColorVector *dest,QDP_Sub_Block QDP_block[NRG+1],QDP_ColorVector *src,QDP_ColorMatrix *wlink[NRG][RG_Ncn]);

QLA_Real plaquette_qdp(QDP_ColorMatrix *link_qdp[RG_Nd],QDP_Sub_Block s,int n);

void RG_create_path(QDP_ColorMatrix *pr_wlink[RG_Ncn], QDP_ColorMatrix *link_qdp[RG_Nd], QDP_Sub_Block s,int len);

void RG_transformation(QDP_Sub_Block s[NRG+1]);

int RG_check(QDP_Sub_Block s[NRG+1]);

int RG_check_smear(QDP_Sub_Block s[NRG+1]);

void RG_M_inv(QDP_ColorVector *dest,QDP_ColorVector *src);

void RG_decompose(QDP_ColorVector *dest, QDP_ColorVector *src[RG_Ncn],int cmp[3],QDP_ColorMatrix *wlink[RG_Ncn],QDP_Sub_Block s);

void RG_bulk(QDP_ColorVector *dest[RG_Ncn], QDP_ColorVector *src);

void RG_check_inversion(QDP_ColorVector *src,QDP_ColorVector *chi[RG_Ncn]);

void RG_trace(int *sign,int r[4], int cmp[4]);

void point_2V(QLA_ColorVector *s, int coords[]);

void gamma_to_int(RG_gamma g, int r[4]);

RG_gamma int_to_gamma(RG_gamma g);



#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int rvs[RG_Ncn][RG_Nd];
EXTERN int S;
EXTERN int T;


#ifdef RG_BLOCK_SUB
#define EXTERN_OP
#else
#define EXTERN_OP extern
#endif

EXTERN_OP int fact;

#endif /* RG_include_h */
