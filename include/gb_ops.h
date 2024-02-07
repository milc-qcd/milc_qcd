#ifndef _GOLTERMAN_BAILEY_OPS_H
#define _GOLTERMAN_BAILEY_OPS_H

#include <string.h>
#include "../generic_ks/generic_ks_includes.h"

#ifdef GB_BARYON

enum gb_sym_type
{
  GBSYM_UNDEF = -1,
  GBSYM_SYM,
  GBSYM_MIXED,
  GBSYM_ANTISYM
};

enum gb_baryon_op
{
  /* Organized by the operator construction
     operators which may have different strangeness/isospin but that use
     the same construction are grouped together */

  GB_UNDEFINED = -1,
  /* Symmetric 8 */
  GB_8_S_C1, /* Only local operator */
  GB_8_S_C2,
  GB_8_S_C3,
  GB_8_S_C5,
  GB_8_S_C6,

  /* Symmetric 8' */
  GB_8P_S_C4, // #5
  GB_8P_S_C7,

  /* Symmetric 16 up/down (or +,-) */
  GB_16U_S_C2, // #7
  GB_16U_S_C3,
  GB_16U_S_C4,
  GB_16U_S_C6,
  GB_16D_S_C2,
  GB_16D_S_C3,
  GB_16D_S_C4,
  GB_16D_S_C6,

  /* Antisymmetric 8 */
  GB_8_A_C4, // #15
  GB_8_A_C6,
  GB_8_A_C7,

  /* Antisymmetric 16 up/down */
  GB_16U_A_C4, // #18
  GB_16U_A_C6,
  GB_16D_A_C4,
  GB_16D_A_C6,

  /* Mixed 8 with s=0 I_3=1/2 or s=-1 I_3=1
     "I_3=1 -like" */
  GB_8_M_I1_C2, // #22
  GB_8_M_I1_C3,
  GB_8_M_I1_C4,
  GB_8_M_I1_C5,
  GB_8_M_I1_C6_1,
  GB_8_M_I1_C6_2,

  /* Mixed 8', both 0,1-like isospin */
  GB_8P_M_I1_C4, // #28

  /* Mixed 16 up/down
     "I_3=1 -like"*/
  GB_16U_M_I1_C2, // #29
  GB_16U_M_I1_C3,
  GB_16U_M_I1_C4_1,
  GB_16U_M_I1_C4_2,
  GB_16U_M_I1_C6_1,
  GB_16U_M_I1_C6_2,
  GB_16U_M_I1_C7,
  GB_16D_M_I1_C2,
  GB_16D_M_I1_C3,
  GB_16D_M_I1_C4_1,
  GB_16D_M_I1_C4_2,
  GB_16D_M_I1_C6_1,
  GB_16D_M_I1_C6_2,
  GB_16D_M_I1_C7,
  MAX_GB_BARYON     // #43
};

/* Extra spin-taste index for doing 2-point tie-ups with 3-point infrastructure */
#define GB_2POINT_BACKPROP 127

/* generic_ks/gb_baryon_snk.c declaration */
void gb_baryon(ks_prop_field *qko0[], ks_prop_field *qko1[], ks_prop_field *qko2[],
               su3_matrix *links, enum gb_baryon_op src_op[],
               enum gb_baryon_op snk_op[],
               int stIdx, short dowall[], short docube[], int num_d, int num_s, int r0[],
               int mom[], char par[], complex *momfld, int flip_snk[],
               int num_corr_gb, int phase[], Real fact[], complex *prop[]);

/* generic_ks/gb_ops.c function declarations */
int gb_baryon_op(char *label);
char *gb_baryon_label(enum gb_baryon_op gbop);
enum  gb_baryon_op decode_gb_op (char* sym_label, char* gts_irrep, int ns, int cls);

enum gb_sym_type gb_get_sym_type(enum gb_baryon_op gbop);
void triplet_to_singlet_index(int t_idx, int s_idx[]);
int singlet_to_triplet_index(int s_idx[]);
int offset_singlet_index(int s_idx, int offset);
int offset_triplet_index(int t_idx, int offset);
int singlet_index_to_disp(int s_idx);
void singlet_index_to_dir(int s_idx, int dir[]);
int gb_get_class(enum gb_baryon_op gbop);
int gb_get_num_terms(enum gb_baryon_op gbop);
int gb_get_num_permutations(enum gb_baryon_op gbop, int num_d, int num_s);
int gb_get_isospin(enum gb_baryon_op gbop);
int gb_get_class_variant(enum gb_baryon_op gbop);
void gb_get_triplet_index(enum gb_baryon_op gbop, int t_idx[]);
void gb_get_prefactors(enum gb_baryon_op gbop, int pf[]);
void gb_get_permutation_order(enum gb_baryon_op gbop, int idx, int num_d, int num_s,
  int in[], int perm[]);
Real gb_get_permutation_sign(enum gb_baryon_op gbop, int idx, int num_d, int num_s);
#endif

#endif /* GB_BARYON */
