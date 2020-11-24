#ifndef CLOV_MAP_H
#define CLOV_MAP_H

#include <qphix/qphix_config.h>
#include "../include/qphixj/qphixj_arch.h"
#include "../include/qphixj/qphixj_internal.h"
#include "../include/qphixj/qphixj.h"

extern "C" {
#include "generic_includes.h"
}

class QPHIXJClovMap {

 public:

 QPHIXJClovMap( const int *geometry, QPHIXJ_vec_layout_t *vl ) :
            lX(nx/geometry[0]), lY(ny/geometry[1]), 
	    lZ(nz/geometry[2]), lT(nt/geometry[3]), 
	    By(vl->By), Bz(vl->Bz), NCores(vl->NCores), Sy(vl->Sy), Sz(vl->Sz), PadXY(vl->PadXY), 
	    PadXYZ(vl->PadXYZ), MinCt(vl->MinCt), S(QPHIX_SOALEN) {
	    subLattSize[0] = lX; subLattSize[1] = lY; subLattSize[2] = lZ; subLattSize[3] = lT;}
	  
template<typename FT, int V>
  QPHIXJ_DiracFermion_struct<FT, V> *
  create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd);

template<typename FT, int V>
  void 
  extract_D_to_wvec( wilson_vector *dest, QPHIXJ_DiracFermion_struct<FT, V> *qdf, 
		     QPHIXJ_evenodd_t evenodd );

template<typename FT, int V>
  void 
  destroy_D( QPHIXJ_DiracFermion_struct<FT, V> *qdf );

template<typename FT, int V>
  QPHIXJ_FermionLinksWilson_struct<FT, V>  *
  wilson_create_L_from_MILC( su3_matrix *backlinks, clover *clov, Real kappa, QPHIXJ_evenodd_t evenodd );

template<typename FT, int V>
  void  
  wilson_destroy_L( QPHIXJ_FermionLinksWilson_struct<FT, V> *ql );

private:

  int subLattSize[4];
  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;
  const int lX;
  const int lY;
  const int lZ;
  const int lT;
  const int S;
};

#endif // CLOV_MAP_H
