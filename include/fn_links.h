#ifndef _FN_LINKS_H
#define _FN_LINKS_H
/* MILC version 7 */

#include "../include/complex.h"
#include "../include/precision.h"
#include "../include/su3.h"
#include "../include/link_phase_info.h"

/* The fn_links_t "class" */

typedef struct {
  link_phase_info_t *phase;
  su3_matrix *fat;
  su3_matrix *lng;
  su3_matrix *fatback;  // NULL if unused
  su3_matrix *lngback;  // NULL if unused
} fn_links_t;

su3_matrix *create_lnglinks(void);
void destroy_lnglinks(su3_matrix *lng);
su3_matrix *create_fatlinks(void);
void destroy_fatlinks(su3_matrix *fat);
void load_fn_backlinks(fn_links_t *fn);
void destroy_fn_backlinks(fn_links_t *fn);

fn_links_t *create_fn_links(void);
void destroy_fn_links(fn_links_t *fn);

void copy_fn(fn_links_t *fn_src, fn_links_t *fn_dst);
void scalar_mult_fn(fn_links_t *fnsrc, Real s, fn_links_t *fndst);
void add_fn(fn_links_t *fnA, fn_links_t *fnB, fn_links_t *fnC);

#endif /* _FN_LINKS_H */
