#ifndef _KS_ACTION_COEFFS_QOP_H
#define _KS_ACTION_COEFFS_QOP_H

#include "../include/complex.h"
#include "../include/ks_action_paths.h"
#include <qop.h>

/* ks_action_coeffs_qop.c */
QOP_asqtad_coeffs_t *create_asqtad_coeffs_qop(ks_action_paths *ap);
void destroy_asqtad_coeffs_qop(QOP_asqtad_coeffs_t *ac);
char *get_ac_string_asqtad_qop(QOP_asqtad_coeffs_t *ac);

QOP_hisq_coeffs_t *create_hisq_coeffs_qop(ks_action_paths_hisq *ap);
void destroy_hisq_coeffs_qop(QOP_hisq_coeffs_t *ac);
int get_n_naiks_qop(QOP_hisq_coeffs_t *ac);
double *get_eps_naik_qop(QOP_hisq_coeffs_t *ac);
int get_ugroup_qop(QOP_hisq_coeffs_t *ac);
char *get_ac_string_hisq_qop(QOP_hisq_coeffs_t *ac);

#endif /* _KS_ACTION_COEFFS_QOP_H */
