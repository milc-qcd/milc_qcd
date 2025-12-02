/*
Copyright (c) 1998-2018 MILC Developers and Contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice, this permission notice and the following
disclaimers shall be included in all copies or substantial portions of
the Software.

Neither the name of the MILC collaboration, nor the names of its
developers or contributors may be used to endorse or promote products
derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
/*  END LEGAL */

/**
 * @brief Wrappers for Grid staggered link fattening
 * @author Carleton DeTar & Curtis Taylor Peterson
 * @details
 * This wraps the Grid HISQ smearing routines for use in MILC. Please see 
 * documentation within HighlyImprovedStaggeredFermionImpl.h for details on the
 * HISQ smearing implementation.
 *
 * Acknowledgements:
 *   Curtis Taylor Peterson would like to thank James Osborn for developing/testing
 *   the implementations of HISQ in QEX and QOPQDP, from which the "fast" option for 
 *   the fat7/asqtad smearing and derivativative are based and has been tested against.
 *   
 *   This material is based upon work supported by the U.S. Department of Energy, 
 *   Office of Science, Office of Advanced Scientific Computing Research, Scientific 
 *   Discovery through Advanced Computing (SciDAC) program.
 */

#include <omp.h>

#undef HMC

#include "../include/macros.h"
extern "C" {
#include "../include/fermion_links.h"
}

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
#include "../include/mGrid/mGrid_assert.h"
#include "../include/mGrid/mGrid_timing.h"

#include "../generic/gridMap.h"

#include <Grid/Grid.h>
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>

using namespace Grid;

//
// C++ wrappers for Grid HISQ API
//

template<typename LatticeGaugeField, typename Gimpl, typename Complex> 
static void hisqLinks(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* fat,
  su3_matrix* lng,
  su3_matrix* in,
  GridCartesian* CGrid,
  bool reunitarizeFatLinks = false
) {
  auto timer = GridHISQTimer();

  LatticeGaugeField U(CGrid), V(CGrid), W(CGrid), X(CGrid);

  GRID_ASSERT(&X != NULL, GRID_MEM_ERROR);
  GRID_ASSERT(&WWW != NULL, GRID_MEM_ERROR);

  HighlyImprovedStaggeredFermionImpl<Gimpl> hisq(&CGrid, false); 

  // smear
  timer.tic();
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(in, &U);
  hisq.smear(V, U);
  hisq.project(W, V);
  if (lng != NULL) {
    LatticeGaugeField WWW(CGrid);
    hisq.smear(X, WWW, W);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &WWW);
  } else { hisq.smear(X, X, W, false); }
  gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &X);
  timer.toc("Grid smearing");

  timer.epoch("Generate fat and long links");
}

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void hisqAuxLinks(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* U, 
  su3_matrix* V, 
  su3_matrix* W,
  GridCartesian* CGrid
) {
  auto timer = GridHISQTimer();

  LatticeGaugeField gV(CGrid), gW(CGrid);

  hisqLinks<LatticeGaugeField, Gimpl, Complex>(info, path_coeff, V, NULL, U, CGrid);

  // smearing + additional projection of "auxiliary" link
  timer.tic();
  HighlyImprovedStaggeredFermionImpl<Gimpl> hisq(&CGrid, false); 
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(V, &gV);
  hisq.project(gW, gV);
  timer.toc("Grid auxiliary link projection");
  
  gridToMilcGaugeField<LatticeGaugeField, Complex>(V, &gV);
  gridToMilcGaugeField<LatticeGaugeField, Complex>(W, &gW);

  timer.epoch("Generate HISQ aux links");
}

//
// the Grid C API for link fattening
//

void GRID_F3_hisq_links(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* fat,
  su3_matrix* lng,
  su3_matrix* in,
  GRID_4Dgrid* grid_full
) {
  //  std::cout << "GRID_F3_hisq_links is not supported yet" << std::endl;
  //  assert(0);
  // hisqLinks<LatticeGaugeFieldF, StaggeredImplF, ComplexF>(
  //   info, path_coeff, fat, lng, in, grid_full->gridF
  // );
}

void GRID_D3_hisq_links(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* fat,
  su3_matrix* lng,
  su3_matrix* in,
  GRID_4Dgrid* grid_full
) {
  hisqLinks<LatticeGaugeFieldD, StaggeredImplD, ComplexD>(
    info, path_coeff, fat, lng, in, grid_full->gridD
  );
}

void GRID_F3_hisq_aux_links(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* U, 
  su3_matrix* V, 
  su3_matrix* W,
  GRID_4Dgrid* grid_full
) {
  std::cout << "GRID_F3_hisq_aux_links" << std::endl;
  assert(0);
  // hisqAuxLinks<LatticeGaugeFieldF, StaggeredImplF, ComplexF>(
  //   info, path_coeff, U, V, W, grid_full->gridF
  //);
}

void GRID_D3_hisq_aux_links(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* U, 
  su3_matrix* V, 
  su3_matrix* W,
  GRID_4Dgrid* grid_full
) {
  hisqAuxLinks<LatticeGaugeFieldD, StaggeredImplD, ComplexD>(
    info, path_coeff, U, V, W, grid_full->gridD
  );
}