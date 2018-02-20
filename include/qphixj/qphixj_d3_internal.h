#ifndef QPHIXJ_D3_INTERNAL_H
#define QPHIXJ_D3_INTERNAL_H

#include <qphix/geometry.h>

/* Define the single-precision opaque objects */
/* This file is included in qphixj_internal.h */

typedef typename QPhiX::Geometry<double,VECLEN_DP,QPHIX_SOALEN,COMPRESS> Geom_D;
typedef typename QPhiX::Geometry<double,VECLEN_DP,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor_D;
typedef typename QPhiX::Geometry<double,VECLEN_DP,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock Gauge_D;
typedef typename QPhiX::Geometry<double,VECLEN_DP,QPHIX_SOALEN,COMPRESS>::CloverBlock Clover_D;

struct QPHIXJ_D3_DiracFermion_struct {
  QPHIXJ_evenodd_t parity;
  Spinor_D *p_even;
  Spinor_D *p_odd;
};

  /* Wilson datatypes */

struct QPHIXJ_D3_FermionLinksWilson_struct {
  QPHIXJ_evenodd_t parity;
  Clover_D *A_cb0;
  Clover_D *A_inv_cb1;
  Gauge_D *packed_gauge_cb0;
  Gauge_D *packed_gauge_cb1;
};


#endif /* QPHIXJ_D3_INTERNAL_H */
