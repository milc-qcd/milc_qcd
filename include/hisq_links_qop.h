#ifndef _HISQ_LINKS_QOP_H
#define _HISQ_LINKS_QOP_H
/************************ hisq_links_qop.h ***************************/
/* MILC version 7 */

#include "../include/complex.h"
#include "../include/precision.h"
#include "../include/su3.h"
#include "../include/link_phase_info.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/fn_links.h"
#include "../include/fn_links_qop.h"
#include <qop.h>

/* The hisq_links_qop_t "class" */

typedef struct {
  link_phase_info_t *phase;
  QOP_F3_FermionLinksHisq *hl_F;
  QOP_D3_FermionLinksHisq *hl_D;
} hisq_links_qop_t;

hisq_links_qop_t *create_hisq_links_qop(QOP_hisq_coeffs_t *hc, int precision,
					su3_matrix *links, int want_deps, int want_aux);

void destroy_hisq_links_qop(hisq_links_qop_t *hl);

QOP_F3_FermionLinksHisq *get_F_hisq_links_qop(hisq_links_qop_t *hl);
QOP_D3_FermionLinksHisq *get_D_hisq_links_qop(hisq_links_qop_t *hl);

void set_asqtad_links_from_hisq(fn_links_qop_t *fn, hisq_links_qop_t *hl, int i);
void unset_asqtad_links_from_hisq(fn_links_qop_t *fn);

void set_asqtad_deps_links_from_hisq(fn_links_qop_t *fn, hisq_links_qop_t *hl);
void unset_asqtad_deps_links_from_hisq(fn_links_qop_t *fn);

#endif /* _HISQ_LINKS_QOP_H */
