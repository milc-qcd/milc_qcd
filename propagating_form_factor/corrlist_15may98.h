/****************************** corrlist.h ********************************/
/* Two-point correlators computed in this code */
/* May 15, 1998.  Used in generating beta = 5.50 am = 0.025 form factors 
   through Sept 1998 */
#define MAX_TWOPT 12
static gamma_corr two_pt[MAX_TWOPT] = {
  G5,   G5,   G5_G5,
  G5T,  G5,   G5G4_G5,
  G5,   G5T,  G5_G5G4,
  GX,   GX,   G1_G1,
  GY,   GX,   G2_G1,
  GZ,   GX,   G3_G1,
  GX,   GY,   G1_G2,
  GY,   GY,   G2_G2,
  GZ,   GY,   G3_G2,
  GX,   GZ,   G1_G3,
  GY,   GZ,   G2_G3,
  GZ,   GZ,   G3_G3
};

/* Three point correlators computed in this code */

#define MAX_THREEPT 55
static gamma_corr three_pt[MAX_THREEPT] = {
  GX,   G5,   PSEUDO_VECTOR_X,   /* B --> pion + leptons */
  GY,   G5,   PSEUDO_VECTOR_Y,
  GZ,   G5,   PSEUDO_VECTOR_Z,
  GT,   G5,   PSEUDO_VECTOR_T,

  GXY,  G5,   PSEUDO_SIGMA_XY,   /* B --> pion V-T operator mixing */
  GYZ,  G5,   PSEUDO_SIGMA_YZ,
  GZX,  G5,   PSEUDO_SIGMA_ZX,
  GXT,  G5,   PSEUDO_SIGMA_XT,
  GYT,  G5,   PSEUDO_SIGMA_YT,
  GZT,  G5,   PSEUDO_SIGMA_ZT,

  GXY,  GX,   RHOX_SIGMA_XY,     /* B --> K* + photon */
  GYZ,  GX,   RHOX_SIGMA_YZ,     /* and B --> rho operator mixing */
  GZX,  GX,   RHOX_SIGMA_ZX,
  GXT,  GX,   RHOX_SIGMA_XT,
  GYT,  GX,   RHOX_SIGMA_YT,
  GZT,  GX,   RHOX_SIGMA_ZT,

  GX,   GX,   RHOX_VECTOR_X,     /* B --> rho + leptons */
  GY,   GX,   RHOX_VECTOR_Y,
  GZ,   GX,   RHOX_VECTOR_Z,
  GT,   GX,   RHOX_VECTOR_T,

  G5X,  GX,   RHOX_AXIAL_VECTOR_X,     /* B --> rho + leptons */
  G5Y,  GX,   RHOX_AXIAL_VECTOR_Y,
  G5Z,  GX,   RHOX_AXIAL_VECTOR_Z,
  G5T,  GX,   RHOX_AXIAL_VECTOR_T,

  G5,   GX,   RHOX_PSEUDO,     /* B --> rho + leptons */

  GXY,  GY,   RHOY_SIGMA_XY,     /* B --> K* + photon */
  GYZ,  GY,   RHOY_SIGMA_YZ,     /* and B --> rho operator mixing */
  GZX,  GY,   RHOY_SIGMA_ZX,
  GXT,  GY,   RHOY_SIGMA_XT,
  GYT,  GY,   RHOY_SIGMA_YT,
  GZT,  GY,   RHOY_SIGMA_ZT,

  GX,   GY,   RHOY_VECTOR_X,     /* B --> rho + leptons */
  GY,   GY,   RHOY_VECTOR_Y,
  GZ,   GY,   RHOY_VECTOR_Z,
  GT,   GY,   RHOY_VECTOR_T,

  G5X,  GY,   RHOY_AXIAL_VECTOR_X,     /* B --> rho + leptons */
  G5Y,  GY,   RHOY_AXIAL_VECTOR_Y,
  G5Z,  GY,   RHOY_AXIAL_VECTOR_Z,
  G5T,  GY,   RHOY_AXIAL_VECTOR_T,

  G5,   GY,   RHOY_PSEUDO,     /* B --> rho + leptons */

  GXY,  GZ,   RHOZ_SIGMA_XY,     /* B --> K* + photon */
  GYZ,  GZ,   RHOZ_SIGMA_YZ,     /* and B --> rho operator mixing */
  GZX,  GZ,   RHOZ_SIGMA_ZX,
  GXT,  GZ,   RHOZ_SIGMA_XT,
  GYT,  GZ,   RHOZ_SIGMA_YT,
  GZT,  GZ,   RHOZ_SIGMA_ZT,

  GX,   GZ,   RHOZ_VECTOR_X,     /* B --> rho + leptons */
  GY,   GZ,   RHOZ_VECTOR_Y,
  GZ,   GZ,   RHOZ_VECTOR_Z,
  GT,   GZ,   RHOZ_VECTOR_T,

  G5X,  GZ,   RHOZ_AXIAL_VECTOR_X,     /* B --> rho + leptons */
  G5Y,  GZ,   RHOZ_AXIAL_VECTOR_Y,
  G5Z,  GZ,   RHOZ_AXIAL_VECTOR_Z,
  G5T,  GZ,   RHOZ_AXIAL_VECTOR_T,

  G5,   GZ,   RHOZ_PSEUDO      /* B --> rho + leptons */
};





