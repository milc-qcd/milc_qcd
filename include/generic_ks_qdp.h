#include "../include/generic_ks.h"
#include <qdp.h>

void dslash_fn_qdp(QDP_ColorVector *src, QDP_ColorVector *dest,
		   QDP_Subset parity);
void dslash_fn_special_qdp(QDP_ColorVector *src, QDP_ColorVector *dest,
			   QDP_Subset parity, QDP_ColorVector *temp[]);
void dslash_fn_special2_qdp(QDP_ColorVector *src, QDP_ColorVector *dest,
			    QDP_Subset parity, QDP_ColorVector *temp[]);
int ks_congrad_qdp(QDP_ColorVector *src, QDP_ColorVector *dest, QLA_Real mass,
		   int niter, QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);
int ks_multicg_qdp(QDP_ColorVector *src, QDP_ColorVector **dest,
		   QLA_Real *masses, int num_masses, int niter,
		   QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);
