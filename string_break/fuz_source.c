/************************** fuz_source.c *******************************/
/* MIMD version 6 */
/* version of 8/20/98 by UMH */

/* Make a gaussian random source, spread out in the on-axis directions
   by distance r0 and make gauge co-variant with smeared link fields.
   Gaussian source is in grd, fuzzed source in fsrc.
   ttt, cg_p and resid are use as temporaries. */

#include "string_break_includes.h"

/** #define DBG **/

void fuz_source(field_offset grd, field_offset fsrc, int r0, int t0, int step) {
register int i, dir, j, k;
int x, y, z, t;
int t00, t_step;
register site *s;
msg_tag *tag0, *tag1;
Real ftmp;

    clear_latvec( grd, EVENANDODD);
    clear_latvec( F_OFFSET(ttt), EVENANDODD);

    if( t0 >= nt ){
	t00 = 0;
	t_step = step;
    }
    else{
	t00 = t0;
	t_step = nt;
    }

#ifndef DBG
    /* Make a gaussian random vector */
    for(t=t00;t<nt;t+=t_step){
	for(x=0;x<nx;x+=step)for(y=0;y<ny;y+=step)for(z=0;z<nz;z+=step){
	    if( node_number(x,y,z,t) != mynode() ) continue;
	    i = node_index(x,y,z,t);
	    s = &(lattice[i]);

#ifdef RAN_GAUGE
	    for(j=0;j<3;j++){
#ifdef SITERAND
		s->ttt.c[j].real = gaussian_rand_no(&(s->site_prn));
		s->ttt.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
		s->ttt.c[j].real = gaussian_rand_no(&node_prn);
		s->ttt.c[j].imag = gaussian_rand_no(&node_prn);
#endif
	    }
	    mult_adj_su3_mat_vec(&(s->rgt), &(s->ttt),
		(su3_vector *)F_PT(s,grd));
#else	/* RAN_GAUGE */
	    for(j=0;j<3;j++){
#ifdef SITERAND
		((su3_vector *)F_PT(s,grd))->c[j].real =
		    gaussian_rand_no(&(s->site_prn));
		((su3_vector *)F_PT(s,grd))->c[j].imag =
		    gaussian_rand_no(&(s->site_prn));
#else
		((su3_vector *)F_PT(s,grd))->c[j].real =
		    gaussian_rand_no(&node_prn);
		((su3_vector *)F_PT(s,grd))->c[j].imag =
		    gaussian_rand_no(&node_prn);
#endif
	    }
#endif	/* RAN_GAUGE */
	}
    }

#else	/* ifdef DBG */
    j = off_axis_flag % 3;
    k = off_axis_flag / 3;
    t00 += (nt - t_step);
    switch(k){
      case 0:
	x = 0;
	y = 0;
	z = 0;
        ftmp = 1.0;
	break;
      case 1:
	x = 1;
	y = 0;
	z = 0;
        t00 = (t00 + 1) % nt;
        ftmp = 2.0;
	break;
      case 2:
	x = 1;
	y = 1;
	z = 0;
        t00 = (t00 + 2) % nt;
        ftmp = 4.0;
	break;
      case 3:
	x = 1;
	y = 1;
	z = 1;
        t00 = (t00 + 1) % nt;
        ftmp = 8.0;
	break;
      case 4:
	x = 2;
	y = 1;
	z = 0;
        t00 = (t00 + 1) % nt;
        ftmp = 16.0;
	break;
      case 5:
	x = 0;
	y = 1;
	z = 0;
        t00 = (t00 + 1) % nt;
        ftmp = 32.0;
	break;
      case 6:
	x = 1;
	y = 2;
	z = 0;
        t00 = (t00 + 1) % nt;
        ftmp = 64.0;
	break;
      case 7:
	x = 0;
	y = 0;
	z = 1;
        t00 = (t00 + 1) % nt;
        ftmp = 128.0;
	break;
    }
    if( node_number(x,y,z,t00) == mynode() )
    {
	i = node_index(x,y,z,t00);
	s = &(lattice[i]);
#ifdef RAN_GAUGE
        lattice[i].ttt.c[j].real = ftmp;
	mult_adj_su3_mat_vec(&(s->rgt), &(s->ttt),
	    (su3_vector *)F_PT(s,grd));
#else
	((su3_vector *)F_PT(s,grd))->c[j].real = ftmp;
#endif
    }
#endif	/* ifdef DBG */

    /* Give central value weight two */
    scalar_mult_latvec( grd, 2.0, fsrc, EVENANDODD);

    if (r0 > 0){

	for(dir=XUP;dir<=ZUP;dir++){
	    /* Start gathering for 'backward' link-product */
	    tag0 = start_gather_site(grd, sizeof(su3_vector),
		dir, EVENANDODD, gen_pt[0]);

	    /* Start 'forward' link-product */
	    FORALLSITES(i,s){
		mult_adj_su3_mat_vec(&(s->link[dir]),
		    (su3_vector *)F_PT(s,grd), &(s->ttt));
	    }
	    tag1 = start_gather_site(F_OFFSET(ttt), sizeof(su3_vector),
		OPP_DIR(dir), EVENANDODD, gen_pt[1]);

	    for(k=1;k<r0;k++){
		wait_gather(tag0);
/*
		copy_latvec( (field_offset)gen_pt[0], F_OFFSET(resid),
		    EVENANDODD);
*/
		FORALLSITES(i,s){
		    su3vec_copy((su3_vector *)gen_pt[0][i], &(s->resid));
		}
		FORALLSITES(i,s){
		    mult_su3_mat_vec(&(s->link[dir]), &(s->resid), &(s->cg_p));
		}
		if(k==1){
		    cleanup_gather(tag0);
		    tag0 = start_gather_site(F_OFFSET(cg_p), sizeof(su3_vector),
			dir, EVENANDODD, gen_pt[0]);
		}
		else{
		    restart_gather_site(F_OFFSET(cg_p), sizeof(su3_vector),
			dir, EVENANDODD, gen_pt[0], tag0);
		}

		wait_gather(tag1);
/*
		copy_latvec( (field_offset)gen_pt[1], F_OFFSET(resid),
		    EVENANDODD);
*/
		FORALLSITES(i,s){
		    su3vec_copy((su3_vector *)gen_pt[1][i], &(s->resid));
		}
		FORALLSITES(i,s){
		    mult_adj_su3_mat_vec(&(s->link[dir]), &(s->resid),
			&(s->ttt));
		}
		restart_gather_site(F_OFFSET(ttt), sizeof(su3_vector),
		    OPP_DIR(dir), EVENANDODD, gen_pt[1], tag1);
	    } /* k<r0 */

	    wait_gather(tag0);
	    FORALLSITES(i,s){
		mult_su3_mat_vec_sum(&(s->link[dir]),
		    (su3_vector *)(gen_pt[0][i]), (su3_vector *)F_PT(s,fsrc));
	    }
	    cleanup_gather(tag0);

	    wait_gather(tag1);
	    FORALLSITES(i,s){
		add_su3_vector((su3_vector *)(gen_pt[1][i]),
		    (su3_vector *)F_PT(s,fsrc), (su3_vector *)F_PT(s,fsrc));
	    }
	    cleanup_gather(tag1);

	} /* directions */

    } /* r0 > 0 */

} /* fuz_source */
