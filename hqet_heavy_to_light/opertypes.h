/****************************** opertypes.h ********************************/
/* Definitions of form factor components calculated in this code */

/* 
   The first part of the name refers to the final state meson, the
   second part, the current vertex operator.  To allow for operator
   improvement, we provide for tensor corrections to vector operators
   and vice versa, and for pseudoscalar corrections to axial vector
   operators.

   Some improvement schemes also include SW rotations of the light
   quark fields.  Thus each operator listed below is calculated
   twice: once with and once without the SW rotation.  Hence only
   half of the computed operator types and names are listed below.

   This list is deliberately inclusive, since the computational cost
   of additional operators is relatively small.
   
   Thus, depending on the momentum of the current and the velocity of
   the heavy-light meson, some operators will be related by symmetry
   and some will be zero. 
   */

/* Modifications
   C. DeTar 5/01/97 Original version
   C. DeTar 5/24/97 Cut 3 pt list in 1/2 and added "_ROT" suffix
   */

enum threept_opertype { 
  PSEUDO_VECTOR_X ,
  PSEUDO_VECTOR_Y ,
  PSEUDO_VECTOR_Z , 
  PSEUDO_VECTOR_T , 
  PSEUDO_SIGMA_XY ,
  PSEUDO_SIGMA_XZ ,
  PSEUDO_SIGMA_XT ,
  PSEUDO_SIGMA_YZ ,
  PSEUDO_SIGMA_YT ,
  PSEUDO_SIGMA_ZT ,

  RHOX_SIGMA_XY ,
  RHOX_SIGMA_XZ ,
  RHOX_SIGMA_XT ,
  RHOX_SIGMA_YZ ,
  RHOX_SIGMA_YT ,
  RHOX_SIGMA_ZT ,
  RHOX_VECTOR_X , 
  RHOX_VECTOR_Y , 
  RHOX_VECTOR_Z , 
  RHOX_VECTOR_T , 
  RHOX_AXIAL_VECTOR_X , 
  RHOX_AXIAL_VECTOR_Y , 
  RHOX_AXIAL_VECTOR_Z , 
  RHOX_AXIAL_VECTOR_T , 
  RHOX_PSEUDO ,

  RHOY_SIGMA_XY ,
  RHOY_SIGMA_XZ ,
  RHOY_SIGMA_XT ,
  RHOY_SIGMA_YZ ,
  RHOY_SIGMA_YT ,
  RHOY_SIGMA_ZT ,
  RHOY_VECTOR_X , 
  RHOY_VECTOR_Y , 
  RHOY_VECTOR_Z , 
  RHOY_VECTOR_T , 
  RHOY_AXIAL_VECTOR_X , 
  RHOY_AXIAL_VECTOR_Y , 
  RHOY_AXIAL_VECTOR_Z , 
  RHOY_AXIAL_VECTOR_T , 
  RHOY_PSEUDO ,

  RHOZ_SIGMA_XY ,
  RHOZ_SIGMA_XZ ,
  RHOZ_SIGMA_XT ,
  RHOZ_SIGMA_YZ ,
  RHOZ_SIGMA_YT ,
  RHOZ_SIGMA_ZT ,
  RHOZ_VECTOR_X , 
  RHOZ_VECTOR_Y , 
  RHOZ_VECTOR_Z , 
  RHOZ_VECTOR_T , 
  RHOZ_AXIAL_VECTOR_X , 
  RHOZ_AXIAL_VECTOR_Y , 
  RHOZ_AXIAL_VECTOR_Z , 
  RHOZ_AXIAL_VECTOR_T , 
  RHOZ_PSEUDO ,

  HALF_NO_THREEPT_OPERS  } ;  /* Required last entry */

#define NO_THREEPT_OPERS HALF_NO_THREEPT_OPERS*2;

/* Maximum length of correlator names.  Must leave room for 
   appending suffix to denote SW rotation */

#define MAXNAME 32
static char SW_suffix[5] = "_ROT";


/* The next list must be in 1-to-1 correspondence with the preceding list */

static char name_3pt[ HALF_NO_THREEPT_OPERS ][ MAXNAME ] = {
  "PSEUDO_VECTOR_X" ,            /* B -> pi + leptons */
  "PSEUDO_VECTOR_Y" ,
  "PSEUDO_VECTOR_Z" , 
  "PSEUDO_VECTOR_T" , 
  "PSEUDO_SIGMA_XY" ,
  "PSEUDO_SIGMA_XZ" ,
  "PSEUDO_SIGMA_XT" ,
  "PSEUDO_SIGMA_YZ" ,
  "PSEUDO_SIGMA_YT" ,
  "PSEUDO_SIGMA_ZT" ,

  "RHOX_SIGMA_XY" ,                /* B -> K^* + gamma */
  "RHOX_SIGMA_XZ" ,                /* and B -> rho + leptons */
  "RHOX_SIGMA_XT" ,
  "RHOX_SIGMA_YZ" ,
  "RHOX_SIGMA_YT" ,
  "RHOX_SIGMA_ZT" ,
  "RHOX_VECTOR_X" , 
  "RHOX_VECTOR_Y" , 
  "RHOX_VECTOR_Z" , 
  "RHOX_VECTOR_T" , 
  "RHOX_AXIAL_VECTOR_X" , 
  "RHOX_AXIAL_VECTOR_Y" , 
  "RHOX_AXIAL_VECTOR_Z" , 
  "RHOX_AXIAL_VECTOR_T" , 
  "RHOX_PSEUDO" ,

  "RHOY_SIGMA_XY" ,
  "RHOY_SIGMA_XZ" ,
  "RHOY_SIGMA_XT" ,
  "RHOY_SIGMA_YZ" ,
  "RHOY_SIGMA_YT" ,
  "RHOY_SIGMA_ZT" ,
  "RHOY_VECTOR_X" , 
  "RHOY_VECTOR_Y" , 
  "RHOY_VECTOR_Z" , 
  "RHOY_VECTOR_T" , 
  "RHOY_AXIAL_VECTOR_X" , 
  "RHOY_AXIAL_VECTOR_Y" , 
  "RHOY_AXIAL_VECTOR_Z" , 
  "RHOY_AXIAL_VECTOR_T" , 
  "RHOY_PSEUDO" ,

  "RHOZ_SIGMA_XY" ,
  "RHOZ_SIGMA_XZ" ,
  "RHOZ_SIGMA_XT" ,
  "RHOZ_SIGMA_YZ" ,
  "RHOZ_SIGMA_YT" ,
  "RHOZ_SIGMA_ZT" ,
  "RHOZ_VECTOR_X" , 
  "RHOZ_VECTOR_Y" , 
  "RHOZ_VECTOR_Z" , 
  "RHOZ_VECTOR_T" , 
  "RHOZ_AXIAL_VECTOR_X" , 
  "RHOZ_AXIAL_VECTOR_Y" , 
  "RHOZ_AXIAL_VECTOR_Z" , 
  "RHOZ_AXIAL_VECTOR_T" , 
  "RHOZ_PSEUDO"
  } ;

/* Definitions of two-point functions calculated in this code */

enum twopt_opertype { 
  G5_G5 , 
  G5G4_G5 , 
  G5_G5G4 , 
  G1_G1 , 
  G2_G2 , 
  G3_G3 , 
  NO_TWOPT_OPERS } ; /* Required last entry */

/* The next list must be in 1-to-1 correspondence with the preceding list */

static char twopt_name[NO_TWOPT_OPERS][MAXNAME] = {
  "G5_G5" , 
  "G5G4_G5" , 
  "G5_G5G4" , 
  "G1_G1" ,
  "G2_G2" ,
  "G3_G3"
} ;


