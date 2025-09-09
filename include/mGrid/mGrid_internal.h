#ifndef _MGRID_INTERNAL_H
#define _MGRID_INTERNAL_H

#include <Grid/Grid.h>

using namespace Grid;
//using namespace Grid::QCD;

// Containers for grid definitions

struct GRID_4Dgrid_struct { GridCartesian *gridF; GridCartesian *gridD; };
struct GRID_4DRBgrid_struct { GridRedBlackCartesian *gridF; GridRedBlackCartesian *gridD; };
struct GRID_5Dgrid_struct { GridCartesian *gridF; GridCartesian *gridD; };
struct GRID_5DRBgrid_struct { GridRedBlackCartesian *gridF; GridRedBlackCartesian *gridD; };

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
struct GRID_ColorMatrix_struct {
  LGF *links;
};

template<typename LGF>
struct GRID_FermionLinksAsqtad_struct {
  LGF *thnlinks;
  LGF *fatlinks;
  LGF *lnglinks;
};

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionF>  GRID_F3_ColorVector_struct;
typedef struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DF>  GRID_F3_ColorVectorBlock_struct;
typedef struct GRID_ColorMatrix_struct<LatticeGaugeFieldF> GRID_F3_ColorMatrix_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldF> GRID_F3_FermionLinksAsqtad_struct;

typedef struct GRID_ColorVector_struct<ImprovedStaggeredFermionD> GRID_D3_ColorVector_struct;
typedef struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DD> GRID_D3_ColorVectorBlock_struct;
typedef struct GRID_ColorMatrix_struct<LatticeGaugeFieldD> GRID_D3_ColorMatrix_struct;
typedef struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldD> GRID_D3_FermionLinksAsqtad_struct;

template< typename ISF >
struct GRID_ColorVectorArray_struct
{
  std::vector< typename ISF::FermionField > * cv;
  int N;
};
typedef struct GRID_ColorVectorArray_struct<ImprovedStaggeredFermionF> GRID_F3_ColorVectorArray_struct;
typedef struct GRID_ColorVectorArray_struct<ImprovedStaggeredFermionD> GRID_D3_ColorVectorArray_struct;


#endif /* _MGRID_INTERNAL_H */

