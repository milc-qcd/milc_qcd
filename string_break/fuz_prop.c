/************************** fuz_prop.c *****************************/
/* MIMD version 6 */

/* Fuzz a quark_propagator (type su2_vector) at the sink.
   The fuzzed propagator overwrites the unfuzzed one;
   use xxx, ttt, cg_p and resid as temporaries. */

#include "string_break_includes.h"

void fuz_prop(field_offset fprop, int r0)
{
register int i;
register site *s; 

int  dir, k;

msg_tag *tag0, *tag1;

    /* Save unfuzzed propagator in xxx */
    copy_latvec( fprop, F_OFFSET(xxx), EVENANDODD);
    /* Give central value weight two */
    scalar_mult_latvec( F_OFFSET(xxx), 2.0, fprop, EVENANDODD);

    if (r0 > 0){

	for(dir=XUP;dir<=ZUP;dir++){
	    /* Start gathering for 'backward' link-product */
	    tag0 = start_gather_site(F_OFFSET(xxx), sizeof(su3_vector),
		dir, EVENANDODD, gen_pt[0]);

	    /* Start 'forward' link-product */
	    FORALLSITES(i,s) {
		mult_adj_su3_mat_vec(&(s->link[dir]), &(s->xxx), &(s->ttt));
	    }
	    tag1 = start_gather_site(F_OFFSET(ttt), sizeof(su3_vector),
		OPP_DIR(dir), EVENANDODD, gen_pt[1]);

	    for(k=1;k<r0;k++) {
		wait_gather(tag0);
/*
		copy_latvec( (field_offset)gen_pt[0], F_OFFSET(resid),
		    EVENANDODD);
*/
		FORALLSITES(i,s) {
		    su3vec_copy((su3_vector *)gen_pt[0][i], &(s->resid));
		}
		FORALLSITES(i,s) {
		    mult_su3_mat_vec(&(s->link[dir]), &(s->resid), &(s->cg_p));
		}
		if(k==1) {
		    cleanup_gather(tag0);
		    tag0 = start_gather_site(F_OFFSET(cg_p), sizeof(su3_vector),
			dir, EVENANDODD, gen_pt[0]);
		}
		else {
		    restart_gather_site(F_OFFSET(cg_p), sizeof(su3_vector),
			dir, EVENANDODD, gen_pt[0], tag0);
		}

		wait_gather(tag1);
/*
		copy_latvec( (field_offset)gen_pt[1], F_OFFSET(resid),
		    EVENANDODD);
*/
		FORALLSITES(i,s) {
		    su3vec_copy((su3_vector *)gen_pt[1][i], &(s->resid));
		}
		FORALLSITES(i,s) {
		    mult_adj_su3_mat_vec(&(s->link[dir]), &(s->resid),
			&(s->ttt));
		}
		restart_gather_site(F_OFFSET(ttt), sizeof(su3_vector),
		    OPP_DIR(dir), EVENANDODD, gen_pt[1], tag1);

	    } /* k<r0 */

	    wait_gather(tag0);
	    FORALLSITES(i,s) {
		mult_su3_mat_vec_sum(&(s->link[dir]),
		    (su3_vector *)(gen_pt[0][i]), (su3_vector *)F_PT(s,fprop));
	    }
	    cleanup_gather(tag0);

	    wait_gather(tag1);
	    FORALLSITES(i,s) {
		add_su3_vector((su3_vector *)(gen_pt[1][i]),
		    (su3_vector *)F_PT(s,fprop), (su3_vector *)F_PT(s,fprop));
	    }
	    cleanup_gather(tag1);

	} /* directions */

    } /* r0 > 0 */

} /* fuz_prop */

