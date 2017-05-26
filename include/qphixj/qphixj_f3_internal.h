#ifndef QPHIXJ_F3_INTERNAL_H
#define QPHIXJ_F3_INTERNAL_H

#include <qphix/geometry.h>

/* Define the single-precision opaque objects */
/* This file is included in qphixj_internal.h */

typedef typename QPhiX::Geometry<float,VECLEN_SP,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor_F;
typedef typename QPhiX::Geometry<float,VECLEN_SP,QPHIX_SOALEN,COMPRESS> Geom_F;
typedef typename QPhiX::Geometry<float,VECLEN_SP,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock Gauge_F;
typedef typename QPhiX::Geometry<float,VECLEN_SP,QPHIX_SOALEN,COMPRESS>::CloverBlock Clover_F;

/* Wilson datatypes */

struct QPHIXJ_F3_DiracFermion_struct {
  QPHIXJ_evenodd_t parity;
  Spinor_F *p_even;
  Spinor_F *p_odd;
};


struct QPHIXJ_F3_FermionLinksWilson_struct {
  QPHIXJ_evenodd_t parity;
  Clover_F *A_cb0;
  Clover_F *A_inv_cb1;
  Gauge_F *packed_gauge_cb0;
  Gauge_F *packed_gauge_cb1;
};


#endif /* QPHIXJ_F3_INTERNAL_H */
