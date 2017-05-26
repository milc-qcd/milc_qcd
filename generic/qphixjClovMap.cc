// Methods for the clovMap class for mapping between MILC and QPHIXJ types
#include "qphixjClovMap.h"
#include <omp.h>

#include <qphix/clover_dslash_def.h>
#include <qphix/clover_dslash_body.h>
#include <qphix/clover.h>
#include <qphix/invcg.h>
#include <qphix/invbicgstab.h>

#include <cstdlib>
#include <assert.h>

using namespace std;
using namespace QPhiX;

extern QPHIXJ_vec_layout_t vecLayout;

template<typename FT, int V>
QPHIXJ_DiracFermion_struct<FT, V> *
QPHIXJClovMap::create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd){
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS> Geom;
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor;
  Geom geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  // Allocate data for the spinors
  Spinor *p_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor *p_odd=(Spinor*)geom.allocCBFourSpinor();

  // Point to the second block of the array. Now there is padding on both ends.
  Spinor *psi_s[2] = { p_even, p_odd };
  
  masterPrintf("Filling Input spinor: ");
  
  double start=omp_get_wtime();
  
#pragma omp parallel for collapse(4)    
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) { 
	  for(int spin=0; spin < 4; spin++) { 
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) { 
		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		for(int cb=0; cb < 2; cb++){
		  int x_coord;
		  if((y+z+t)%2 == 0)
		    x_coord = 2*(x + s*S) + cb;
		  else
		    x_coord = 2*(x + s*S) + 1 - cb;
		  int indMILC = node_index(x_coord, y, z, t);
		  if( (cb == 0) && (evenodd != QPHIXJ_ODD)
		   || (cb == 1) && (evenodd != QPHIXJ_EVEN))
		    {
		      psi_s[cb][ind][col][spin][0][x] = src[indMILC].d[spin].c[col].real;
		      psi_s[cb][ind][col][spin][1][x] = src[indMILC].d[spin].c[col].imag;
		    } 
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);

  // Pack the DiracFermion structure
  QPHIXJ_DiracFermion_struct<FT, V> *qdf = (QPHIXJ_DiracFermion_struct<FT, V> *)malloc(sizeof(QPHIXJ_DiracFermion_struct<FT, V>));
  qdf->parity = QPHIXJ_EVEN;  // Hard coded for now
  qdf->p_even = p_even;
  qdf->p_odd  = p_odd;

  return qdf;
}

// Map a Dirac vector field from QPhiX layout to MILC layout
template<typename FT, int V>
void QPHIXJClovMap::extract_D_to_wvec( wilson_vector *dest, QPHIXJ_DiracFermion_struct<FT, V> *qdf, 
				       QPHIXJ_evenodd_t evenodd ){
  
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS> Geom;
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor;
  Geom geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  // Unpack the DiracFermion structure

  //  assert( evenodd == qdf->parity );
  Spinor *p_even = qdf->p_even;
  Spinor *p_odd = qdf->p_odd;

  // Point to the second block of the array. Now there is padding on both ends.
  Spinor *psi_s[2] = { p_even, p_odd };
  
  masterPrintf("Filling output spinor: ");
  
  double start=omp_get_wtime();

#pragma omp parallel for collapse(4)    
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) { 
	  for(int spin=0; spin < 4; spin++) { 
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) { 
		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		for(int cb=0; cb < 2; cb++){
		  int x_coord;
		  if((y+z+t)%2 == 0)
		    x_coord = 2*(x + s*S) + cb;
		  else
		    x_coord = 2*(x + s*S) + 1 - cb;
		  int indMILC = node_index(x_coord, y, z, t);
		  if( (cb == 0) && (evenodd != QPHIXJ_ODD)
		   || (cb == 1) && (evenodd != QPHIXJ_EVEN))
		    {
		      dest[indMILC].d[spin].c[col].real = psi_s[cb][ind][col][spin][0][x];
		      dest[indMILC].d[spin].c[col].imag = psi_s[cb][ind][col][spin][1][x];
		    } 
		}
	      }
	    }
	  }
	}
      }
    }
  }

  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);
}

// free Dirac vector
template<typename FT, int V>
void QPHIXJClovMap::destroy_D( QPHIXJ_DiracFermion_struct<FT, V> *qdf ){

  if ( qdf == NULL ) return;

  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS> Geom;
  Geom geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

  geom.free(qdf->p_even);
  geom.free(qdf->p_odd);
  free(qdf);
}

// Create wilson fermion links from MILC fields
// Clover from the MILC clover structure
// Forward gauge links from the global MILC lattice[i].link site structure
// Backward gauge links from the raw4 layout
// Precision conversion takes place in the copies if need be

template<typename FT, int V>
QPHIXJ_FermionLinksWilson_struct<FT, V>  *
QPHIXJClovMap::wilson_create_L_from_MILC( su3_matrix *backlinks, clover *clov, Real kappa,
					  QPHIXJ_evenodd_t evenodd ){

  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS> Geom;
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor;
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::CloverBlock Clover;

  Geom geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  bool compress = vecLayout.compress12;
  bool verbose = false;

  Clover *A_cb0 = (Clover*)geom.allocCBClov();
  Clover *A_inv_cb1 = (Clover*)geom.allocCBClov();

  assert( evenodd != QPHIXJ_EVENODD );  //Can't support dual parity
  assert( evenodd == QPHIXJ_EVEN );  // Support only even for now
  assert( compress == COMPRESS );

  masterPrintf("Filling the Clover Term\n");

  double start=omp_get_wtime();
  double scale = 4.*kappa*kappa;

  // See ~/QPhiX_Joo/package-qphix-9-26-16/src/qphix/include/qphix/qdp_packer_parscalar.h
  
#pragma omp parallel for collapse(4)
  for(int t = 0; t < lT; t++) {
    for(int z = 0; z < lZ; z++) {
      for(int y = 0; y < lY; y++) {
	for(int s = 0; s < nvecs; s++) {
	  for(int x = 0; x < S; x++) {
	    int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
	    // This will work out to be between 0 and veclen
	    int xx = (y%nyg)*S+x;
	    for(int cb=0; cb < 2; cb++){
	      int x_coord;
	      if((y+z+t)%2 == 0)
		x_coord = 2*(x + S*s) + cb;
	      else
		x_coord = 2*(x + S*s) + 1-cb;
	      int indMILC = node_index(x_coord, y, z, t);
	      // MILC clov must have the inverse on the other checkerboard, scaled by 4*kappa^2
	      for(int i=0; i < 6; i++){ 
		if( (cb == 0) && (evenodd == QPHIXJ_EVEN)
		    || (cb == 1) && (evenodd == QPHIXJ_ODD))
		  {
		    A_cb0[block].diag1[i][xx] = clov->clov_diag[indMILC].di[0][i];
		    A_cb0[block].diag2[i][xx] = clov->clov_diag[indMILC].di[1][i];
		  } 
		else 
		  {
		    A_inv_cb1[block].diag1[i][xx] = clov->clov_diag[indMILC].di[0][i]*scale;
		    A_inv_cb1[block].diag2[i][xx] = clov->clov_diag[indMILC].di[1][i]*scale;
		  } 
	      }
	      
	      for(int i=0; i < 15; i++){ 
		if( (cb == 0) && (evenodd == QPHIXJ_EVEN)
		    || (cb == 1) && (evenodd == QPHIXJ_ODD))
		  {
		    A_cb0[block].off_diag1[i][0][xx] = clov->clov[indMILC].tr[0][i].real;
		    A_cb0[block].off_diag1[i][1][xx] = clov->clov[indMILC].tr[0][i].imag;
		    A_cb0[block].off_diag2[i][0][xx] = clov->clov[indMILC].tr[1][i].real;
		    A_cb0[block].off_diag2[i][1][xx] = clov->clov[indMILC].tr[1][i].imag;
		  } 
		else 
		  {
		    A_inv_cb1[block].off_diag1[i][0][xx] = clov->clov[indMILC].tr[0][i].real*scale;
		    A_inv_cb1[block].off_diag1[i][1][xx] = clov->clov[indMILC].tr[0][i].imag*scale;
		    A_inv_cb1[block].off_diag2[i][0][xx] = clov->clov[indMILC].tr[1][i].real*scale;
		    A_inv_cb1[block].off_diag2[i][1][xx] = clov->clov[indMILC].tr[1][i].imag*scale;
		  } 
	      } // i
	    } // cb
	  } // x
	}
      }
    }
  }

  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);

  masterPrintf("Filling the Gauge links\n");

  // Allocate data for the gauges
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();

  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  
  start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for(int t = 0; t < lT; t++) {
    for(int z = 0; z < lZ; z++) {
      for(int y = 0; y < lY; y++) {
	for(int s = 0; s < nvecs; s++) {
	  for(int mu = 0; mu < 8; mu++) {
	    for(int c = 0; c < (compress ? 2 : 3) ; c++) {
	      for(int c2 = 0; c2 < 3; c2++) {
		for(int x = 0; x < S; x++) {
		  int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
		  // This will work out to be between 0 and veclen
		  int xx = (y%nyg)*S+x;
		  for(int cb=0; cb < 2; cb++){
		    int x_coord;
		    if((y+z+t)%2 == 0)
		      x_coord = 2*(x + S*s) + cb;
		    else
		      x_coord = 2*(x + S*s) + 1-cb;
		    int indMILC = node_index(x_coord, y, z, t);
		    int dir = mu/2;
		    if(mu%2 == 0){
		      // Row and column indices are reversed from MILC
		      // Adjoint of MILC
		      u_packed[cb][block][mu][c][c2][RE][xx] =  backlinks[4*indMILC+dir].e[c][c2].real;
		      u_packed[cb][block][mu][c][c2][IM][xx] = -backlinks[4*indMILC+dir].e[c][c2].imag;
		    } else {
		      u_packed[cb][block][mu][c2][c][RE][xx] =  lattice[indMILC].link[dir].e[c][c2].real;
		      u_packed[cb][block][mu][c2][c][IM][xx] =  lattice[indMILC].link[dir].e[c][c2].imag;
		    }
		  }
		}
	      }
	    } // row
	  }
	}
      }
    }
  }

  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);

  // Pack and return the wilson fermion links structure

  QPHIXJ_FermionLinksWilson_struct<FT, V> *ql = 
    (QPHIXJ_FermionLinksWilson_struct<FT, V> *)malloc(sizeof(QPHIXJ_FermionLinksWilson_struct<FT, V>));
  ql->parity = QPHIXJ_EVEN;  // Hard-coded for now
  ql->A_cb0 = A_cb0;
  ql->A_inv_cb1 = A_inv_cb1;
  ql->packed_gauge_cb0 = packed_gauge_cb0;
  ql->packed_gauge_cb1 = packed_gauge_cb1;

  return ql;
}

// free wilson fermion links
template<typename FT, int V>
void  
QPHIXJClovMap::wilson_destroy_L( QPHIXJ_FermionLinksWilson_struct<FT, V> *ql ){

  if ( ql == NULL ) return;

  typedef typename QPhiX::Geometry<FT,V,QPHIX_SOALEN,COMPRESS> Geom;
  Geom geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

  geom.free(ql->A_cb0);
  geom.free(ql->A_inv_cb1);
  geom.free(ql->packed_gauge_cb0);
  geom.free(ql->packed_gauge_cb1);

  free(ql);
}

//====================================================================//
// The QPHIXJ C API for mapping between MILC and QPHIXJ types

extern QPHIXJ_vec_layout_t vecLayout;

// Map a Dirac vector field from MILC layout to QPhiX layout
QPHIXJ_F3_DiracFermion *
QPHIXJ_F3_create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.create_D_from_wvec<float, VECLEN_SP>( src, evenodd );
}
  
// Map a Dirac vector field from MILC layout to QPhiX layout
QPHIXJ_D3_DiracFermion *
QPHIXJ_D3_create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.create_D_from_wvec<double, VECLEN_DP>( src, evenodd );
}
  
// Map a Dirac vector field from QPhiX layout to MILC layout
void 
QPHIXJ_F3_extract_D_to_wvec( wilson_vector *dest, QPHIXJ_F3_DiracFermion *qdf, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  CM.extract_D_to_wvec<float, VECLEN_SP>( dest, qdf, evenodd );
}

// Map a Dirac vector field from QPhiX layout to MILC layout
void 
QPHIXJ_D3_extract_D_to_wvec( wilson_vector *dest, QPHIXJ_D3_DiracFermion *qdf, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  CM.extract_D_to_wvec<double, VECLEN_DP>( dest, qdf, evenodd );
}

// free Dirac vector
void  
QPHIXJ_F3_destroy_D( QPHIXJ_F3_DiracFermion *qdf ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.destroy_D<float, VECLEN_SP>( qdf );
}

// free Dirac vector
void  
QPHIXJ_D3_destroy_D( QPHIXJ_D3_DiracFermion *qdf ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.destroy_D<double, VECLEN_DP>( qdf );
}

// create wilson fermion links from MILC
QPHIXJ_F3_FermionLinksWilson  *
QPHIXJ_F3_wilson_create_L_from_MILC( su3_matrix *backlinks, clover *clov, Real kappa, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.wilson_create_L_from_MILC<float, VECLEN_SP>( backlinks, clov, kappa, evenodd );
}

// create wilson fermion links from MILC
QPHIXJ_D3_FermionLinksWilson  *
QPHIXJ_D3_wilson_create_L_from_MILC( su3_matrix *backlinks, clover *clov, Real kappa, QPHIXJ_evenodd_t evenodd ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.wilson_create_L_from_MILC<double, VECLEN_DP>( backlinks, clov, kappa, evenodd );
}

// free wilson fermion links
void  
QPHIXJ_F3_wilson_destroy_L( QPHIXJ_F3_FermionLinksWilson *ql ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.wilson_destroy_L<float, VECLEN_SP>( ql );
}

// free wilson fermion links
void  
QPHIXJ_D3_wilson_destroy_L( QPHIXJ_D3_FermionLinksWilson *ql ){
  QPHIXJClovMap CM( get_logical_dimensions(), &vecLayout);
  return CM.wilson_destroy_L<double, VECLEN_DP>( ql );
}


