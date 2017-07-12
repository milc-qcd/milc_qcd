#ifndef _MGRID_INTERNAL_H
#define _MGRID_INTERNAL_H

#include <Grid/Grid.h>

using namespace Grid;
using namespace Grid::QCD;

// The color-vector and gauge field objects used in the interface

template<typename ImprovedStaggeredFermion>
struct GRID_ColorVector_struct {
  typename ImprovedStaggeredFermion::FermionField *cv;
};

template<typename ImprovedStaggeredFermion5D>
struct GRID_ColorVector5D_struct {
  typename ImprovedStaggeredFermion5D::FermionField *cv;
};

template<typename LatticeGaugeField>
struct GRID_FermionLinksAsqtad_struct {
  LatticeGaugeField *thnlinks;
  LatticeGaugeField *fatlinks;
  LatticeGaugeField *lnglinks;
};

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionF>  GRID_F3_ColorVector_struct;
typedef struct GRID_ColorVector5D_struct<ImprovedStaggeredFermion5DF>  GRID_F3_ColorVector5D_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldF> GRID_F3_FermionLinksAsqtad_struct;

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionD> GRID_D3_ColorVector_struct;
typedef struct GRID_ColorVector5D_struct<ImprovedStaggeredFermion5DD> GRID_D3_ColorVector5D_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldD> GRID_D3_FermionLinksAsqtad_struct;

#endif /* _MGRID_INTERNAL_H */

