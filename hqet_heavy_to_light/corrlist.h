/****************************** corrlist.h ********************************/
/* Two-point correlators computed in this code */
#define MAX_TWOPT 6
static gamma_corr two_pt[MAX_TWOPT] = {
  { G5,   G5,   G5_G5   }, 
  { G5T,  G5,   G5G4_G5 },
  { G5,   G5T,  G5_G5G4 },
  { GX,   GX,   G1_G1   },
  { GY,   GY,   G2_G2   },
  { GZ,   GZ,   G3_G3   }
};

/* Three point correlators computed in this code */

#define MAX_THREEPT 14
static gamma_corr hqet_to_light[MAX_THREEPT] = {
  { GX,   G5,   PSEUDO_VECTOR_X     },   /* B --> pion + leptons */
  { GY,   G5,   PSEUDO_VECTOR_Y     },
  { GZ,   G5,   PSEUDO_VECTOR_Z     },
  { GT,   G5,   PSEUDO_VECTOR_T     },
  { GXY,  GX,   RHOX_SIGMA_XY       },     /* B --> K* + photon */
  { GXT,  GX,   RHOX_SIGMA_XT       },
  { GYT,  GX,   RHOX_SIGMA_YT       },
  { GYZ,  GX,   RHOX_SIGMA_YZ       },
  { GT,   GX,   RHOX_VECTOR_T       },     /* B --> rho + leptons */
  { GX,   GX,   RHOX_VECTOR_X       },
  { GY,   GX,   RHOX_VECTOR_Y       },
  { G5T,  GX,   RHOX_AXIAL_VECTOR_T },
  { G5X,  GX,   RHOX_AXIAL_VECTOR_X },
  { G5Y,  GX,   RHOX_AXIAL_VECTOR_Y }
};





