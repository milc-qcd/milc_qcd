
//#include "ks_spectrum_includes.h"
#include "../include/generic.h"
#include "../include/gb_ops.h"
#include <assert.h>
#include <string.h>

#ifdef GB_BARYON

/*-------------------------------------------------------------*/
/**
  Array of strings which correspond to the Golterman-Bailey baryon operators
 */
static char *gb_baryon_string[MAX_GB_BARYON] = {
  "GB_8_S_C1", // only local operator
  "GB_8_S_C2",
  "GB_8_S_C3",
  "GB_8_S_C5",
  "GB_8_S_C6",
  // symmetric 8'
  "GB_8P_S_C4",
  "GB_8P_S_C7",
  // symmetric 16 up/down (or +,-)
  "GB_16U_S_C2",
  "GB_16U_S_C3",
  "GB_16U_S_C4",
  "GB_16U_S_C6",
  //
  "GB_16D_S_C2",
  "GB_16D_S_C3",
  "GB_16D_S_C4",
  "GB_16D_S_C6",
  // antisymmetric" 8
  "GB_8_A_C4",
  "GB_8_A_C6",
  "GB_8_A_C7",
  // antisymmetric 16 up/down
  "GB_16U_A_C4",
  "GB_16U_A_C6",
  //
  "GB_16D_A_C4",
  "GB_16D_A_C6",
  // mixed 8 with s=0 I_3=1/2 or s=-1 I_3=1
  //  "I_3=1 -like"
  "GB_8_M_I1_C2",
  "GB_8_M_I1_C3",
  "GB_8_M_I1_C4",
  "GB_8_M_I1_C5",
  "GB_8_M_I1_C6_1",
  "GB_8_M_I1_C6_2",
  // mixed 8', both 0,1-like isospin
  "GB_8P_M_I1_C4",
  // mixed 16 up/down
  //  "I_3=1 -like"
  "GB_16U_M_I1_C2",
  "GB_16U_M_I1_C3",
  "GB_16U_M_I1_C4_1",
  "GB_16U_M_I1_C4_2",
  "GB_16U_M_I1_C6_1",
  "GB_16U_M_I1_C6_2",
  "GB_16U_M_I1_C7",
  //
  "GB_16D_M_I1_C2",
  "GB_16D_M_I1_C3",
  "GB_16D_M_I1_C4_1",
  "GB_16D_M_I1_C4_2",
  "GB_16D_M_I1_C6_1",
  "GB_16D_M_I1_C6_2",
  "GB_16D_M_I1_C7"
};

/*-------------------------------------------------------------*/
/**
  Convert a baryon string into the corresponding enum index
 */
int gb_baryon_op(char *label){
  int i;
  for(i = 0; i < MAX_GB_BARYON; i++){
    if(strcmp(label,gb_baryon_string[i]) == 0) return i;
  }
  return -1;
}

/*-------------------------------------------------------------*/
/**
  Convert a Golterman-Bailey baryon operator enum into a string
 */
char *gb_baryon_label(enum gb_baryon_op gbop){
  return gb_baryon_string[(int)gbop];
}

/*------------------------------------------------------------------*/
/**
  Return the symmetry type of the Golterman-Bailey baryon operator.
  Choices are Symmetric, Mixed, Antisymmetric
 */
enum gb_sym_type gb_get_sym_type(enum gb_baryon_op gbop){
  int gbi = (int) gbop;
  if(0<=gbi && gbi<(int)GB_8_A_C4 )                          return GBSYM_SYM;
  else if((int)GB_8_A_C4    <=gbi && gbi<(int)GB_8_M_I1_C2)  return GBSYM_ANTISYM;
  else if((int)GB_8_M_I1_C2 <=gbi && gbi<(int)MAX_GB_BARYON) return GBSYM_MIXED;
  else return GBSYM_UNDEF;
}

/*------------------------------------------------------------------*/
/**
  Convert a three-quark index into an array of one-quark indices.
  The three-quark index corresponds to the point-splitting of the
  three quarks which make up the baryon operator.
 */
void triplet_to_singlet_index(int t_idx, int s_idx[]){
  s_idx[0] = (int)t_idx   %8;
  s_idx[1] =((int)t_idx/8)%8;
  s_idx[2] = (int)t_idx  /64;
}

/*------------------------------------------------------------------*/
/**
  Convert an array of three one-quark indices into a three-quark index.
  The three indices correspond to the point-splitting of the
  three quarks which make up the baryon operator.
 */
int singlet_to_triplet_index(int s_idx[]){
  return s_idx[2]*64 + s_idx[1]*8 + s_idx[0];
}

/*------------------------------------------------------------------*/
/**
  Offset a singlet index by another singlet index. Effectively
  adding the two positions mod 2.
 */
int offset_singlet_index(int s_idx, int offset){
  int disp[3];
  disp[2] = ((s_idx/4)   + (offset/4))   % 2;
  disp[1] =(((s_idx/2)%2)+((offset/2)%2))% 2;
  disp[0] = ((s_idx%2)   + (offset%2))   % 2;
  return disp[2]*4 + disp[1]*2 + disp[0];
}

/*------------------------------------------------------------------*/
/**
  Offset a singlet index by another singlet index. Effectively
  adding the two positions mod 2.
 */
int offset_triplet_index(int t_idx, int offset){
  int s_idx[3];
  triplet_to_singlet_index(t_idx,s_idx);
  s_idx[0] = offset_singlet_index(s_idx[0],offset);
  s_idx[1] = offset_singlet_index(s_idx[1],offset);
  s_idx[2] = offset_singlet_index(s_idx[2],offset);
  return singlet_to_triplet_index(s_idx);
}

/*------------------------------------------------------------------*/
/**
   Convert a singlet (one quark) position index into a displacement.
   Displacement means N for an N-link point-splitting. See gb_ops.c
 */
int
singlet_index_to_disp(int s_idx){
  switch(s_idx){
    case 0: return 0; break;
    case 1:
    case 2:
    case 4: return 1; break;
    case 3:
    case 5:
    case 6: return 2; break;
    case 7: return 3; break;
    default: return -1; break;
  }
}

/*------------------------------------------------------------------*/
/**
   Convert a singlet (one quark) position index into an array of
   directions. Directions correspond to the point-splitting necessary
   to construct the baryon operator. See gb_ops.c
 */
void
singlet_index_to_dir(int s_idx, int dir[]){
  switch(s_idx){
    case 0: break;
    case 1: dir[0] = XUP; break;
    case 2: dir[0] = YUP; break;
    case 3: dir[0] = XUP; dir[1] = YUP; break;
    case 4: dir[0] = ZUP; break;
    case 5: dir[0] = XUP; dir[1] = ZUP; break;
    case 6: dir[0] = YUP; dir[1] = ZUP; break;
    case 7: dir[0] = XUP; dir[1] = YUP; dir[2] = ZUP; break;
    default: break;
  }
}

/*------------------------------------------------------------------*/
/**
  Get the operator class defined in the Bailey thesis. Classes are
  defined by their structure on the unit cube, so all operators
  within a class will have similar operator structure.
 */
int gb_get_class(enum gb_baryon_op gbop){
  switch(gbop){
    case GB_8_S_C1       :
      return 1; break;
    case GB_8_S_C2       :
    case GB_16U_S_C2     :
    case GB_16D_S_C2     :
    case GB_8_M_I1_C2    :
    case GB_16U_M_I1_C2  :
    case GB_16D_M_I1_C2  :
      return 2; break;
    case GB_8_S_C3       :
    case GB_16U_S_C3     :
    case GB_16D_S_C3     :
    case GB_8_M_I1_C3    :
    case GB_16U_M_I1_C3  :
    case GB_16D_M_I1_C3  :
      return 3; break;
    case GB_8P_S_C4      :
    case GB_16U_S_C4     :
    case GB_16D_S_C4     :
    case GB_8_A_C4       :
    case GB_16U_A_C4     :
    case GB_16D_A_C4     :
    case GB_8_M_I1_C4    :
    case GB_8P_M_I1_C4   :
    case GB_16U_M_I1_C4_1:
    case GB_16U_M_I1_C4_2:
    case GB_16D_M_I1_C4_1:
    case GB_16D_M_I1_C4_2:
      return 4;
    case GB_8_S_C5       :
    case GB_8_M_I1_C5    :
      return 5;
    case GB_8_S_C6       :
    case GB_16U_S_C6     :
    case GB_16D_S_C6     :
    case GB_8_A_C6       :
    case GB_16U_A_C6     :
    case GB_16D_A_C6     :
    case GB_8_M_I1_C6_1  :
    case GB_8_M_I1_C6_2  :
    case GB_16U_M_I1_C6_1:
    case GB_16U_M_I1_C6_2:
    case GB_16D_M_I1_C6_1:
    case GB_16D_M_I1_C6_2:
      return 6;
    case GB_8P_S_C7      :
    case GB_8_A_C7       :
    case GB_16U_M_I1_C7  :
    case GB_16D_M_I1_C7  :
      return 7;
    default:
      return -1;
  }
}

/*-------------------------------------------------------------*/
/**
  Get the number of terms (spin-taste terms specifically) required
  to construct the Golterman-Bailey baryon operator.
 */
int gb_get_num_terms(enum gb_baryon_op gbop){
  switch(gbop){
    /* n = 1 */
    case GB_8_S_C1        :
    case GB_8_S_C5        :
    case GB_8P_S_C7       :
    case GB_8_A_C7        :
    case GB_8_M_I1_C5    :
      return 1; break;
    /* n = 2 */
    case GB_16D_S_C2      :
    case GB_16D_S_C3      :
    case GB_16U_S_C4      :
    case GB_16D_S_C6      :
    case GB_16D_A_C4      :
    case GB_16D_A_C6      :
    case GB_16D_M_I1_C2   :
    case GB_16D_M_I1_C3   :
    case GB_16D_M_I1_C4_1 :
    case GB_16D_M_I1_C6_1 :
    case GB_16D_M_I1_C6_2 :
    case GB_16D_M_I1_C7  :
    case GB_16U_M_I1_C4_2 :
    case GB_16U_M_I1_C7  :
      return 2; break;
    /* n = 3 */
    case GB_8_S_C2        :
    case GB_8_S_C3        :
    case GB_8_S_C6        :
    case GB_8P_S_C4       :
    case GB_16U_S_C2      :
    case GB_16U_S_C3      :
    case GB_16U_S_C6      :
    case GB_16D_S_C4      :
    case GB_8_A_C4        :
    case GB_8_A_C6        :
    case GB_16U_A_C4      :
    case GB_16U_A_C6      :
    case GB_16D_M_I1_C4_2 :
    case GB_8_M_I1_C3     :
    case GB_8_M_I1_C2     :
    case GB_8_M_I1_C4     :
    case GB_8_M_I1_C6_1   :
    case GB_8_M_I1_C6_2   :
    case GB_8P_M_I1_C4    :
    case GB_16U_M_I1_C2   :
    case GB_16U_M_I1_C3   :
    case GB_16U_M_I1_C4_1 :
    case GB_16U_M_I1_C6_1 :
    case GB_16U_M_I1_C6_2 :
      return 3; break;
    default:
      return -1;
  }
}

/*------------------------------------------------------------------*/
/**
  Return the number of permutations of quarks needed for constructing
  Golterman-Bailey baryon operators. Only relevant for mixed-symmetry
  operators in the isospin-1-like operators. Number returned is based
  on the number of d and s quarks in the baryon for mixed symmetry operators.
 */
int gb_get_num_permutations(enum gb_baryon_op gbop, int num_d, int num_s){
  //int fact = (num_d == 1 && num_s == 1 ? 2 : 1);
  int fact = 1;
  if(gb_get_sym_type(gbop) == GBSYM_MIXED){
    /* only class 4 has double the number */
    if (gb_get_class(gbop) == 4) { return 2*fact; }
    else                         { return   fact; }

  } else return 1; // symmetric, antisymmetric
}

/*------------------------------------------------------------------*/
/**
   Return the isospin symmetry behavior type
   I=1/2, I=1 lumped together
  */
int
gb_get_isospin(enum gb_baryon_op gbop){
  /* only care about I quantum number for mixed symmetry */
  /* Isospin1-like, not necessarily isospin 1 */
  switch(gbop){
    case GB_8_M_I1_C2    :
    case GB_8_M_I1_C3    :
    case GB_8_M_I1_C4    :
    case GB_8_M_I1_C5    :
    case GB_8_M_I1_C6_1  :
    case GB_8_M_I1_C6_2  :
    case GB_8P_M_I1_C4   :
    case GB_16U_M_I1_C2  :
    case GB_16U_M_I1_C3  :
    case GB_16U_M_I1_C4_1:
    case GB_16U_M_I1_C4_2:
    case GB_16U_M_I1_C6_1:
    case GB_16U_M_I1_C6_2:
    case GB_16U_M_I1_C7  :
    case GB_16D_M_I1_C2  :
    case GB_16D_M_I1_C3  :
    case GB_16D_M_I1_C4_1:
    case GB_16D_M_I1_C4_2:
    case GB_16D_M_I1_C6_1:
    case GB_16D_M_I1_C6_2:
    case GB_16D_M_I1_C7  :
      return 1; break;
    default:
      return 0; break;
  }
}

/*------------------------------------------------------------------*/
/**
   No distinct names for 1,2 variants. Need to have some handle on these
   for mixed symmetry. Only care about classes 4 and 6.

   1 variant has -sign on second flavor term for class 4 and
   iij flavor content for class 6
  */
int
gb_get_class_variant(enum gb_baryon_op gbop){
  /* class 4 - gives sign for second term */
  /* class 6 - tells which flavor combo to use (iij or iji) */
  switch(gbop){
    case GB_8_M_I1_C4    :
    case GB_8_M_I1_C6_1  :
    case GB_16U_M_I1_C4_1:
    case GB_16U_M_I1_C6_1:
    case GB_16D_M_I1_C4_1:
    case GB_16D_M_I1_C6_1:
      return -1; break; // variant -
    case GB_8P_M_I1_C4   :
    case GB_8_M_I1_C6_2  :
    case GB_16U_M_I1_C4_2:
    case GB_16U_M_I1_C6_2:
    case GB_16D_M_I1_C4_2:
    case GB_16D_M_I1_C6_2:
      return 1; break; // variant +
    default:
      return 0; break;
  }
}

/*------------------------------------------------------------------*/
/**
   Permutes the quarks from the argument in and saves the objects to a
   3-dimensional array perm. If more than one permutation is needed,
   idx can specify which permutation is requested. Excluding class 4,
   idx==1 is the other permutation needed to generate the I_3=0 isospin
   singlet or triplet variant.
 */
void
gb_get_permutation_order(enum gb_baryon_op gbop, int idx, int num_d, int num_s,
                         int in[], int perm[]){
  // always assumed to be uuu or uud quark order
  int tperm[3];
  ///* prepermute the order in select cases, others should already be ordered */
  //if ((num_s == 0 && num_d == 2) || // udd -> ddu , uss -> ssu , dss -> ssd
  //     num_s == 2)
  //      { tperm[0] = in[2]; tperm[1] = in[1]; tperm[2] = in[0]; }
  //else  { tperm[0] = in[0]; tperm[1] = in[1]; tperm[2] = in[2]; }

  if (gb_get_sym_type(gbop) == GBSYM_MIXED){
    tperm[0] = in[0]; tperm[1] = in[1]; tperm[2] = in[2]; // no prepermuting
    if (idx == 0){ // index 0 of permutations
      switch(gb_get_class(gbop)){
        case 6:
          if (gb_get_class_variant(gbop) == 1) // 6_2 has different order
               { perm[0] = tperm[0]; perm[1] = tperm[2]; perm[2] = tperm[1]; return; } // iji
          else { perm[0] = tperm[0]; perm[1] = tperm[1]; perm[2] = tperm[2]; return; } // iij
        case 5:
        case 7:  // 1->2->3 permuted based on mixed operator direction
           perm[0] = tperm[2]; perm[1] = tperm[0]; perm[2] = tperm[1]; return; // jii
        case 2:
        case 3:
        case 4:
        default: // normal order
           perm[0] = tperm[0]; perm[1] = tperm[1]; perm[2] = tperm[2]; return; // iij
      }
    } else { // idx != 0 , must be class 4 (has 2 terms)
      if (gb_get_class(gbop) == 4 && idx == 1) {
        perm[0] = tperm[0]; perm[1] = tperm[2]; perm[2] = tperm[1]; return; // iji
      } else {
        printf("invalid permutation order!\n");
        assert(0);
      }
    } // idx == 0
  } else { // sym_type
    perm[0] = in[0]; perm[1] = in[1]; perm[2] = in[2]; return; // ijk
  }
}

/*------------------------------------------------------------------*/
/**
  Get the sign of the taste permutation for mixed symmetry operators.
  Number of terms can be different based on operator.
 */
Real
gb_get_permutation_sign(enum gb_baryon_op gbop, int idx, int num_d, int num_s){
  Real fact;
  if (gb_get_sym_type(gbop) == GBSYM_MIXED){
    if (idx == 0){
      return 1.;
    } else { // idx != 0 , must be class 4
      if(gb_get_class(gbop) == 4){
        return (Real) gb_get_class_variant(gbop);
      } else {
        printf("invalid permutation sign request!\n");
        assert(0);
      } // gb_get_class == 4
    }
  } else {
    printf("gb_get_permutation_sign called with non-mixed symmetry operator!\n");
    assert(0);
    //return 1.;
  }
}

/*------------------------------------------------------------------*/
/**
  Return an array of up to four sets of triplet indices which correspond
  to the Golterman-Bailey baryon operator. Triplet indices are organized
  by their class, which dictates their position on the unit cube.
 */
void gb_get_triplet_index(enum gb_baryon_op gbop, int t_idx[]){
  int cls = gb_get_class(gbop);
  switch(cls){
    /* 0,0,0 */
    case 1: t_idx[0]=   0; break;
    /* 0,1,1; 0,2,2; 0,3,3; */
    case 2: t_idx[0]=  72; t_idx[1]= 144; t_idx[2]= 288; break;
    /* 0,23,23; 0,13,13; 0,12,12; */
    case 3: t_idx[0]= 432; t_idx[1]= 360; t_idx[2]= 216; break;
    /* NOT EXACTLY THE SAME AS BAILEY PAPER */
    /* 23,3,2; 13,3,1; 12,2,1; */
    case 4: t_idx[0]= 166; t_idx[1]= 101; t_idx[2]=  83; break;
    /* 123,123,0; 0,123,123; */
    case 5: t_idx[0]=  63; t_idx[1]= 504; break;
    /* 123,23,1; 123,13,2; 123,12,3; */
    case 6: t_idx[0]= 119; t_idx[1]= 175; t_idx[2]= 287; break;
    /* 23,13,12; 13,23,12; OUT OF ORDER FROM BAILEY THESIS */
    case 7: t_idx[0]= 238; t_idx[1]= 245; break;
    /* WRONG TASTE ! */
    ///* 1,12,13; 2,12,23; 3,13,23; */
    //case 4: t_idx[0]= 345; t_idx[1]= 410; t_idx[2]= 428; break;
    ///* 0,0,123; 123,0,0 */
    //case 5: t_idx[0]= 448; t_idx[1]=   7; break;
    ///* 0,1,23; 0,2,13; 0,3,12; */
    //case 6: t_idx[0]= 392; t_idx[1]= 336; t_idx[2]= 224; break;
    ///* 1,2,3; 3,1,2; 2,1,3; 1,3,2; OUT OF ORDER FROM BAILEY THESIS */
    //case 7: t_idx[0]= 273; t_idx[1]= 140; t_idx[2]= 266; t_idx[3]= 161; break;
    default:t_idx[0]=  -1; break;
  }
}

/*------------------------------------------------------------------*/
/**
  Get the array of prefactors used to contruct a Golterman-Bailey baryon operator.
  Operators with similar prefactors are grouped together.
 */
void gb_get_prefactors(enum gb_baryon_op gbop, int pf[]){
  switch(gbop){
    case GB_8_S_C1       :
    case GB_8_S_C2       :
    case GB_8_S_C3       :
    case GB_8_S_C5       :
    case GB_8P_S_C7      :
    case GB_16D_S_C6     :
    case GB_8_A_C4       :
    case GB_8_A_C7       :
    case GB_16D_A_C6     :
    case GB_8_M_I1_C2    :
    case GB_8_M_I1_C3    :
    case GB_8_M_I1_C4    :
      pf[0]= 1; pf[1]= 1; pf[2]= 1; break;
    case GB_16U_S_C4     :
    case GB_16U_M_I1_C4_2:
      pf[0]=-1; pf[1]=-1; pf[2]= 0; break;
    case GB_8_M_I1_C5    :
      pf[0]= 1; pf[1]= 0; pf[2]= 0; break;
    case GB_16D_M_I1_C6_1:
    case GB_16D_M_I1_C6_2:
    case GB_16D_M_I1_C7  :
      pf[0]= 1; pf[1]= 1; pf[2]= 0; break;
    case GB_8_S_C6       :
    case GB_8P_S_C4      :
    case GB_16D_S_C2     :
    case GB_16D_S_C3     :
    case GB_8_A_C6       :
    case GB_16D_A_C4     :
    case GB_8P_M_I1_C4   :
    case GB_8_M_I1_C6_1  :
    case GB_8_M_I1_C6_2  :
      pf[0]= 1; pf[1]=-1; pf[2]= 1; break;
    case GB_16D_M_I1_C2  :
    case GB_16D_M_I1_C3  :
    case GB_16D_M_I1_C4_1:
    case GB_16U_M_I1_C7  :
      pf[0]= 1; pf[1]=-1; pf[2]= 0; break;
    case GB_16U_S_C2     :
    case GB_16U_S_C3     :
    case GB_16U_A_C4     :
    case GB_16U_M_I1_C2  :
    case GB_16U_M_I1_C3  :
    case GB_16U_M_I1_C4_1:
      pf[0]= 1; pf[1]= 1; pf[2]=-2; break;
    case GB_16U_S_C6     :
    case GB_16D_S_C4     :
    case GB_16U_A_C6     :
    case GB_16U_M_I1_C6_1:
    case GB_16U_M_I1_C6_2:
    case GB_16D_M_I1_C4_2:
      pf[0]= 1; pf[1]=-1; pf[2]=-2; break;
    default:
      assert(0);
      break;
  }
}

/*-------------------------------------------------------------*/
enum gb_baryon_op decode_gb_op (char* sym_label, char* gts_irrep, int ns, int cls)
{
  /*
   gts_irrep in {"8","8'","16+","16-"}
   sym_label in {"S","A","M0","M1/2","M1"}
   ns = number of strange quarks in {0,1,2,3}
   cls = class in {1,2,3,41,42,5,61,62,7}
   */
  if(strcmp("S",sym_label) == 0 || strcmp("S*",sym_label) == 0){
    /* all strangeness accepted */
    if (strcmp("8",gts_irrep) == 0) {
      switch (cls)
      {
        case 1:  return GB_8_S_C1; break;
        case 2:  return GB_8_S_C2; break;
        case 3:  return GB_8_S_C3; break;
        case 5:  return GB_8_S_C5; break;
        case 61: return GB_8_S_C6; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("8'",gts_irrep) == 0){
      switch (cls)
      {
        case 41: return GB_8P_S_C4; break;
        case 7:  return GB_8P_S_C7; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("16+",gts_irrep) == 0){
      switch (cls)
      {
        case 2:  return GB_16U_S_C2; break;
        case 3:  return GB_16U_S_C3; break;
        case 41: return GB_16U_S_C4; break;
        case 61: return GB_16U_S_C6; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("16-",gts_irrep) == 0){
      switch (cls)
      {
        case 2:  return GB_16D_S_C2; break;
        case 3:  return GB_16D_S_C3; break;
        case 41: return GB_16D_S_C4; break;
        case 61: return GB_16D_S_C6; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }

    else{
      printf("cannot parse input: unknown GTS irrep %s\n",gts_irrep);
      return GB_UNDEFINED;
    } /* GTS irrep */
  } /* Symmetric */

  else if(strcmp("A",sym_label) == 0 || strcmp("A*",sym_label) == 0){
    /* only ns=1 accepted */
    if (ns == 1){
      if (strcmp("8",gts_irrep) == 0){
        switch (cls)
        {
          case 41: return GB_8_A_C4; break;
          case 61: return GB_8_A_C6; break;
          case 7:  return GB_8_A_C7; break;
          default:
            printf("invalid class %i for symmetry %s and GTS irrep %s\n",
              cls,sym_label,gts_irrep);
            return GB_UNDEFINED; break;
        }
      }
      else if (strcmp("8'",gts_irrep) == 0){
        printf("there are no antisymmetric 8' GTS irreps for baryons\n");
        return GB_UNDEFINED;
      }
      else if (strcmp("16+",gts_irrep) == 0){
        switch (cls)
        {
          case 41: return GB_16U_A_C4; break;
          case 61: return GB_16U_A_C6; break;
          default:
            printf("invalid class %i for symmetry %s and GTS irrep %s\n",
              cls,sym_label,gts_irrep);
            return GB_UNDEFINED; break;
        }
      }
      else if (strcmp("16-",gts_irrep) == 0){
        switch (cls)
        {
          case 41: return GB_16D_A_C4; break;
          case 61: return GB_16D_A_C6; break;
          default:
            printf("invalid class %i for symmetry %s and GTS irrep %s\n",
              cls,sym_label,gts_irrep);
            return GB_UNDEFINED; break;
        }
      }
      else{
        printf("cannot parse input: unknown GTS irrep %s\n",gts_irrep);
        return GB_UNDEFINED;
      } /* GTS irrep */
    }
    else{
      printf("invalid number of strange quarks %i for symmetry %s and GTS irrep %s\n",
        ns,sym_label,gts_irrep);
      return GB_UNDEFINED;
    } /* Strangeness */
  } /* Antisymmetric */

  //else if(strcmp("M0",sym_label) == 0 || strcmp("M0*",sym_label) == 0){
  //  /* only ns=1 accepted for I=0 */
  //  if (ns == 1){
  //    if (strcmp("8",gts_irrep) == 0){
  //      switch (cls) {
  //        case 2:  return GB_8_M_I0_C2;   break;
  //        case 3:  return GB_8_M_I0_C3;   break;
  //        case 41: // allow for both
  //        case 42: return GB_8_M_I0_C4;   break;
  //        case 5:  return GB_8_M_I0_C5;   break;
  //        case 61: return GB_8_M_I0_C6_1; break;
  //        case 62: return GB_8_M_I0_C6_2; break;
  //        default:
  //          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
  //            cls,sym_label,gts_irrep);
  //          return GB_UNDEFINED; break;
  //      }
  //    }
  //    else if (strcmp("8'",gts_irrep) == 0){
  //      switch (cls)
  //      {
  //        case 41: // allow for both
  //        case 42: return GB_8P_M_I0_C4; break;
  //        default:
  //          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
  //            cls,sym_label,gts_irrep);
  //          return GB_UNDEFINED; break;
  //      }
  //    }
  //    else if (strcmp("16+",gts_irrep) == 0){
  //      switch (cls)
  //      {
  //        case 2:  return GB_16U_M_I0_C2;   break;
  //        case 3:  return GB_16U_M_I0_C3;   break;
  //        case 41: return GB_16U_M_I0_C4_1; break;
  //        case 42: return GB_16U_M_I0_C4_2; break;
  //        case 61: return GB_16U_M_I0_C6_1; break;
  //        case 62: return GB_16U_M_I0_C6_2; break;
  //        case 7:  return GB_16U_M_I0_C7;   break;
  //        default:
  //          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
  //            cls,sym_label,gts_irrep);
  //          return GB_UNDEFINED; break;
  //      }
  //    }
  //    else if (strcmp("16-",gts_irrep) == 0){
  //      switch (cls)
  //      {
  //        case 2:  return GB_16D_M_I0_C2;   break;
  //        case 3:  return GB_16D_M_I0_C3;   break;
  //        case 41: return GB_16D_M_I0_C4_1; break;
  //        case 42: return GB_16D_M_I0_C4_2; break;
  //        case 61: return GB_16D_M_I0_C6_1; break;
  //        case 62: return GB_16D_M_I0_C6_2; break;
  //        case 7:  return GB_16D_M_I0_C7;   break;
  //        default:
  //          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
  //            cls,sym_label,gts_irrep);
  //          return GB_UNDEFINED; break;
  //      }
  //    }
  //    else{
  //      printf("cannot parse input: unknown GTS irrep %s\n",gts_irrep);
  //      return GB_UNDEFINED;
  //    } /* GTS irrep */
  //  }
  //  else{
  //    printf("invalid number of strange quarks %i for symmetry %s and GTS irrep %s\n",
  //      ns,sym_label,gts_irrep);
  //    return GB_UNDEFINED;
  //  } /* Strangeness */
  //} /* Mixed symmetric I_3 = 0*/

  else if( ((strcmp("M1/2",sym_label) == 0 || strcmp("M1/2*",sym_label) == 0))
  //else if( ((strcmp("M1/2",sym_label) == 0 || strcmp("M1/2*",sym_label) == 0) && (ns == 0))
        //|| ((strcmp("M1",sym_label)   == 0 || strcmp("M1*",sym_label)   == 0) && ns == 1 )
         ){
    if (strcmp("8",gts_irrep) == 0){
      switch (cls)
      {
        case 2:  return GB_8_M_I1_C2;   break;
        case 3:  return GB_8_M_I1_C3;   break;
        case 41: // allow for both
        case 42: return GB_8_M_I1_C4;   break;
        case 5:  return GB_8_M_I1_C5;   break;
        case 61: return GB_8_M_I1_C6_1; break;
        case 62: return GB_8_M_I1_C6_2; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("8'",gts_irrep) == 0){
      switch (cls)
      {
        case 41: // allow for both
        case 42: return GB_8P_M_I1_C4; break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("16+",gts_irrep) == 0){
      switch (cls)
      {
        case 2:  return GB_16U_M_I1_C2;   break;
        case 3:  return GB_16U_M_I1_C3;   break;
        case 41: return GB_16U_M_I1_C4_1; break;
        case 42: return GB_16U_M_I1_C4_2; break;
        case 61: return GB_16U_M_I1_C6_1; break;
        case 62: return GB_16U_M_I1_C6_2; break;
        case 7:  return GB_16U_M_I1_C7;   break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else if (strcmp("16-",gts_irrep) == 0){
      switch (cls)
      {
        case 2:  return GB_16D_M_I1_C2;   break;
        case 3:  return GB_16D_M_I1_C3;   break;
        case 41: return GB_16D_M_I1_C4_1; break;
        case 42: return GB_16D_M_I1_C4_2; break;
        case 61: return GB_16D_M_I1_C6_1; break;
        case 62: return GB_16D_M_I1_C6_2; break;
        case 7:  return GB_16D_M_I1_C7;   break;
        default:
          printf("invalid class %i for symmetry %s and GTS irrep %s\n",
            cls,sym_label,gts_irrep);
          return GB_UNDEFINED; break;
      }
    }
    else{
      printf("cannot parse input: unknown GTS irrep %s\n",gts_irrep);
      return GB_UNDEFINED;
    } /* GTS irrep */
  }

  else{
    printf("invalid number of strange quarks %i for symmetry %s and GTS irrep %s\n",
      ns,sym_label,gts_irrep);
    return GB_UNDEFINED;
  }

}

#endif //GB_BARYON
