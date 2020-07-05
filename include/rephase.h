#ifndef _REPHASE_H
#define _REPHASE_H

/* rephase.c */
void apply_apbc( su3_matrix *links, int r0t );
void phaseset(void);
void rephase( int flag );
void rephase_field_offset( su3_matrix *internal_links, int flag, 
			   int *status_now, int r0[] );
void rephase_offset( int flag, int r0[] );

#endif /* _REPHASE_H */
