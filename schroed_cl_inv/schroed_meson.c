/*********** schroed_meson.c *************/
/* MIMD version 7 */
/* UMH April 96 */

#include "schroed_cl_includes.h"

void schroed_meson(
    field_offset src1,	/* src1 is type spin_wilson_vector */
    field_offset src2,	/* src1 is type spin_wilson_vector */
    complex *prop[10],
    int max_prop)
{

int gamma_in[4],gamma_out[4];
int n_in,n_out;

/* gamma_in = source Dirac matrix, gamma_out = sink Dirac matrix */

/* note (gamma_5 gamma_mu)^\dagger = gamma_5 gamma_mu */

    /* PION */
    n_in=1;n_out=1;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= GAMMAFIVE;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[0]);


    /* PS505 */
    n_in=1;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= TUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[1]);


    if(max_prop == 2) return;
    /* RHO33 (we use also the other polarizations) */
    n_in=1;n_out=1;
    gamma_in[0]= ZUP;
    gamma_out[0]= ZUP;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[2]);

    gamma_in[0]= XUP;
    gamma_out[0]= XUP;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[2]);

    gamma_in[0]= YUP;
    gamma_out[0]= YUP;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[2]);


    if(max_prop == 3) return;
    /* SCALAR */
    n_in=0;n_out=0;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[3]);


    if(max_prop == 4) return;
    /* PV35 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_in[1]= ZUP;
    gamma_out[0]= ZUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[1]= XUP;
    gamma_out[0]= XUP;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[1]= YUP;
    gamma_out[0]= YUP;
    meson_cont_site(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

} /* w_meson */

/************************** w_source.c *****************************/

void w_source_sf_site(field_offset src, my_quark_source *wqs)
{
  register int i,j;
  register site *s; 
  msg_tag *tag0;

  int color, spin, source_type;
  Real Kappa;

  /* Unpack structure */
  color = wqs->color;
  spin = wqs->spin;
  source_type = wqs->type;
  Kappa = wqs->kappa;

    /*printf("WSOURCE: source = %d\n",source); */
	
    /* zero src to be safe */
    FORALLSITES(i,s) {
	clear_wvec((wilson_vector *)F_PT(s,src)); 
    }

    if(source_type == PLUS) {
	/* Source on first time slice: 2*Kappa*U^\dagger(0,TUP)*P_- */
	/* In MILC convention P_+ and P_- appear interchanged */

	/* get link[TUP] from direction -TUP */
	tag0 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
	    OPP_DIR(TUP), EVENANDODD, gen_pt[0]);

	wait_gather(tag0);
	/* Here spin stands actually for spin+2 */
	if(spin == 0){
	    FORALLSITES(i,s) {
		if(s->t != 1)continue;	/* only do this if t==1 */

		for(j=0;j<3;j++){
		    ((wilson_vector *)F_PT(s,src))->d[0].c[j].real =
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].real;
		    ((wilson_vector *)F_PT(s,src))->d[0].c[j].imag = -
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].imag;
		    ((wilson_vector *)F_PT(s,src))->d[2].c[j].real = -
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].real;
		    ((wilson_vector *)F_PT(s,src))->d[2].c[j].imag =
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].imag;
		}
	    }
	}
	else if(spin == 1){
	    FORALLSITES(i,s) {
		if(s->t != 1)continue;	/* only do this if t==1 */

		for(j=0;j<3;j++){
		    ((wilson_vector *)F_PT(s,src))->d[1].c[j].real =
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].real;
		    ((wilson_vector *)F_PT(s,src))->d[1].c[j].imag = -
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].imag;
		    ((wilson_vector *)F_PT(s,src))->d[3].c[j].real = -
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].real;
		    ((wilson_vector *)F_PT(s,src))->d[3].c[j].imag =
			Kappa*((su3_matrix *)gen_pt[0][i])->e[color][j].imag;
		}
	    }
	}
	else{
	    printf("W_SOURCE: spin must be 0 or 1, not %d\n", spin);
	    terminate(1);
	}
	cleanup_gather(tag0);
    }
    else if(source_type == MINUS) {
	/* Source on time slice nt-1: 2*Kappa*U(nt-1,TUP)*P_- */
	/* In MILC convention P_+ and P_- appear interchanged */

	if(spin == 0){
	    FORALLSITES(i,s) {
		if(s->t != (nt-1))continue;	/* only do this if t==(nt-1) */

		for(j=0;j<3;j++){
		    ((wilson_vector *)F_PT(s,src))->d[0].c[j].real =
			Kappa*s->link[TUP].e[j][color].real;
		    ((wilson_vector *)F_PT(s,src))->d[0].c[j].imag =
			Kappa*s->link[TUP].e[j][color].imag;
		    ((wilson_vector *)F_PT(s,src))->d[2].c[j].real =
			Kappa*s->link[TUP].e[j][color].real;
		    ((wilson_vector *)F_PT(s,src))->d[2].c[j].imag =
			Kappa*s->link[TUP].e[j][color].imag;
		}
	    }
	}
	else if(spin == 1){
	    FORALLSITES(i,s) {
		if(s->t != (nt-1))continue;	/* only do this if t==(nt-1) */

		for(j=0;j<3;j++){
		    ((wilson_vector *)F_PT(s,src))->d[1].c[j].real =
			Kappa*s->link[TUP].e[j][color].real;
		    ((wilson_vector *)F_PT(s,src))->d[1].c[j].imag =
			Kappa*s->link[TUP].e[j][color].imag;
		    ((wilson_vector *)F_PT(s,src))->d[3].c[j].real =
			Kappa*s->link[TUP].e[j][color].real;
		    ((wilson_vector *)F_PT(s,src))->d[3].c[j].imag =
			Kappa*s->link[TUP].e[j][color].imag;
		}
	    }
	}
	else{
	    printf("W_SOURCE: spin must be 0 or 1, not %d\n", spin);
	    terminate(1);
	}
    }

} /* w_source */

