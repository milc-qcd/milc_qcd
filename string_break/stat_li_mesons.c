/************************** stat_li_mesons.c *****************************/
/* MIMD version 6 */

/* Compute the static-light meson propagotors from the light
   quark propagators qprop[j] with gaussian random sources g_rand[j],
   and time-like (but in axial gauge!) products of gauge
   fields representing the static (anti-) quark propagators. */

#include "string_break_includes.h"

void stat_li_mesons(int tot_smear, int step)
{
register int i;
register site *s; 
int nth, t, j;
complex *sl_mes, *sl_mes_t, cc;
msg_tag *tag0;

    nth = nt/2+1;
    sl_mes = (complex *)malloc(nth*sizeof(complex));
    sl_mes_t = (complex *)malloc(nth*sizeof(complex));

    for(t=0; t<nth; t++)
	sl_mes_t[t] = cmplx(0.0,0.0);

    for(j=0; j<num_src; j++){

	/* Start gather from next time slice */
	tag0 = start_gather_site( F_OFFSET(qprop[j]), sizeof(su3_vector),
	    TUP, EVENANDODD, gen_pt[0]);

	for(t=0; t<nth; t++)
	    sl_mes[t] = cmplx(0.0,0.0);

	FORALLSITES(i,s){
	    cc = su3_dot(&(s->g_rand[j]), &(s->qprop[j]));
	    sl_mes[0].real += cc.real;
	    sl_mes[0].imag += cc.imag;
	}

	for(t=1; t<nth; t++){
	    wait_gather(tag0);
/*
	    copy_latvec( (field_offset)gen_pt[0], F_OFFSET(resid), EVENANDODD);
*/
	    FORALLSITES(i,s){
		su3vec_copy((su3_vector *)gen_pt[0][i], &(s->resid));
	    }
	    copy_latvec( F_OFFSET(resid), F_OFFSET(ttt), EVENANDODD);

	    if(t==1){
		cleanup_gather(tag0);
		tag0 = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector),
		    TUP, EVENANDODD, gen_pt[0]);
	    }
	    else if(t<(nth-1)){
		restart_gather_site( F_OFFSET(ttt), sizeof(su3_vector),
		    TUP, EVENANDODD, gen_pt[0], tag0);
	    }
	    else
		cleanup_gather(tag0);

	    FORALLSITES(i,s){
		if( ((s->t)+t)>=nt ){
		    mult_su3_mat_vec(&(s->link[TUP]), &(s->resid), &(s->cg_p));
		    cc = su3_dot(&(s->g_rand[j]), &(s->cg_p));
		}
		else
		    cc = su3_dot(&(s->g_rand[j]), &(s->resid));
		sl_mes[t].real += cc.real;
		sl_mes[t].imag += cc.imag;
	    }

	} /* t < nth */

	g_sync();

	for(t=0; t<nth; t++){
	    g_floatsum( &sl_mes[t].real );
	    sl_mes[t].real /= (Real)volume;
	    g_floatsum( &sl_mes[t].imag );
	    sl_mes[t].imag /= (Real)volume;
	    if(this_node == 0)
		printf("SL%d_MESON_%d %d %d %e %e\n", step, tot_smear, j, t,
		     (double)sl_mes[t].real, (double)sl_mes[t].imag);
	    sl_mes_t[t].real += sl_mes[t].real;
	    sl_mes_t[t].imag += sl_mes[t].imag;
	}

    } /* j < num_src */

    for(t=0; t<nth; t++){
	sl_mes_t[t].real /= (Real)num_src;
	sl_mes_t[t].imag /= (Real)num_src;
	if(this_node == 0)
	    printf("SL%d_MESON_T%d %d %e %e\n", step, tot_smear, t,
		(double)sl_mes_t[t].real, (double)sl_mes_t[t].imag);
    }

    free(sl_mes);
    free(sl_mes_t);

} /* stat_li_mesons */

