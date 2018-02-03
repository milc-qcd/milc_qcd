#ifndef CLOV_INVERT
#define CLOV_INVERT

#include "../include/qphixj/qphixj_internal.h"
#include "../include/qphixj/qphixj.h"
#include "../include/qphixj/qphixj_arch.h"

extern "C" {
#include "generic_clover_includes.h"
}

class QPHIXJClovInvert { 

public: 

QPHIXJClovInvert( QPHIXJ_vec_layout_t *vl, bool do_dslash_, bool do_m_, bool do_bicgstab_) : 
    By(vl->By), Bz(vl->Bz), NCores(vl->NCores), Sy(vl->Sy), Sz(vl->Sz), PadXY(vl->PadXY), 
    PadXYZ(vl->PadXYZ), MinCt(vl->MinCt), N_simt(vl->Sy*vl->Sz), compress12(vl->compress12), 
    do_dslash(do_dslash_), do_m(do_m_), do_bicgstab(do_bicgstab_) {}

template<typename FT, int V>
  void run(QPHIXJ_info_t *info,
	   QPHIXJ_FermionLinksWilson_struct<FT, V> *ql, QPHIXJ_invert_arg_t *inv_arg,
	   QPHIXJ_resid_arg_t *res_arg, FT kappa,
	   QPHIXJ_DiracFermion_struct<FT, V> *qdf_out, 
	   QPHIXJ_DiracFermion_struct<FT, V> *qdf_in); 
 
private:

template<typename FT, int V>
  void runClov(QPHIXJ_info_t *info,
	       QPHIXJ_FermionLinksWilson_struct<FT, V> *ql, QPHIXJ_invert_arg_t *inv_arg,
	       QPHIXJ_resid_arg_t *res_arg, FT kappa,
	       QPHIXJ_DiracFermion_struct<FT, V> *qdf_out, 
	       QPHIXJ_DiracFermion_struct<FT, V> *qdf_in);

  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;
  const int N_simt;
  const bool compress12;
  const bool do_dslash;
  const bool do_m;
  const bool do_bicgstab;
};

#endif // CLOV_INVERT
