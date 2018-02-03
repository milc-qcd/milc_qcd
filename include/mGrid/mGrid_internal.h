#ifndef _MGRID_INTERNAL_H
#define _MGRID_INTERNAL_H

#include <Grid/Grid.h>

using namespace Grid;
using namespace Grid::QCD;

// Containers for grid definitions

struct GRID_4Dgrid_struct { GridCartesian *grid; };
struct GRID_4DRBgrid_struct { GridRedBlackCartesian *grid; };
struct GRID_5Dgrid_struct { GridCartesian *grid; };
struct GRID_5DRBgrid_struct { GridRedBlackCartesian *grid; };

// The color-vector and gauge field objects used in the interface

template<typename ISF>
struct GRID_ColorVector_struct {
  typename ISF::FermionField *cv;
};

template<typename ISF5D>
struct GRID_ColorVectorBlock_struct {
  typename ISF5D::FermionField *cv;
};

template<typename LGF>
struct GRID_FermionLinksAsqtad_struct {
  LGF *thnlinks;
  LGF *fatlinks;
  LGF *lnglinks;
};

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionF>  GRID_F3_ColorVector_struct;
typedef struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DF>  GRID_F3_ColorVectorBlock_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldF> GRID_F3_FermionLinksAsqtad_struct;

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionD> GRID_D3_ColorVector_struct;
typedef struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DD> GRID_D3_ColorVectorBlock_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldD> GRID_D3_FermionLinksAsqtad_struct;

#endif /* _MGRID_INTERNAL_H */

