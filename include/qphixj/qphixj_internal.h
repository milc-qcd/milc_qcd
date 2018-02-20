#ifndef QPHIXJ_INTERNAL_H
#define QPHIXJ_INTERNAL_H

#include <qphix/qphix_config.h>
#include "../include/qphixj/qphixj_arch.h"
#include "../include/qphixj/qphixj_int.h"
#include <qphix/geometry.h>

template<typename FT, int V>
struct QPHIXJ_DiracFermion_struct {
  QPHIXJ_evenodd_t parity;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock *p_even;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock *p_odd;
};

template<typename FT, int V>
struct QPHIXJ_FermionLinksWilson_struct {
  QPHIXJ_evenodd_t parity;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::CloverBlock *A_cb0;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::CloverBlock *A_inv_cb1;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock *packed_gauge_cb0;
  typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock *packed_gauge_cb1;
};

typedef QPHIXJ_DiracFermion_struct<float, VECLEN_SP> QPHIXJ_F3_DiracFermion_struct;
typedef QPHIXJ_DiracFermion_struct<double, VECLEN_DP> QPHIXJ_D3_DiracFermion_struct;

typedef QPHIXJ_FermionLinksWilson_struct<float, VECLEN_SP> QPHIXJ_F3_FermionLinksWilson_struct;
typedef QPHIXJ_FermionLinksWilson_struct<double, VECLEN_DP> QPHIXJ_D3_FermionLinksWilson_struct;

#endif // QPHIXJ_INTERNAL_H
