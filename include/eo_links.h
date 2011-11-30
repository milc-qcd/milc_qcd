#ifndef _EO_LINKS_H
#define _EO_LINKS_H

/* A simplified links structure for generic actions */

/* We do not support boundary twists for these actions */

#include "../include/ks_action_paths.h"

typedef struct {
  ks_action_paths *ap;
} eo_links_t;

#endif /* _EO_LINKS_H */

eo_links_t *create_eo_links(void);
void destroy_eo_links(eo_links_t *fn);
