#ifndef _LINK_PHASE_INFO_H
#define _LINK_PHASE_INFO_H
/*************** link_phase_info.h *************************************/
/* MILC version 7 */

/* Phase information */

typedef struct {
  int twist_in;       /* Whether boundary-twist phases are in */
  Real bdry_phase[4]; /* The boundary-twist phases in current use */
  int r0[4];          /* The effective lattice origin */
} link_phase_info_t;

link_phase_info_t *create_link_phase_info(void);
void destroy_link_phase_info(link_phase_info_t *lp);
void copy_link_phase_info(link_phase_info_t *dest, link_phase_info_t *src);


#endif /* _LINK_PHASE_INFO_H */
