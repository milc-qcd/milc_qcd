/******* d_congrad5_fn_tmp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* 12/2006 C. D. NEEDS UPDATING IF WE ARE GOING TO USE IT */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.
	This is different than our old definition of the stopping
	criterion.  To convert an old stopping residual to the new
	one, multiply the old one by sqrt( (2/3)/(8+2*m) )
        This is because the source is obtained from
        a random vector with average squared magnitude 3 on each site.
        Then, on 1/2 the sites, we gather and sum the eight neighboring
        random vectors and add 2*m times the local vector.
            source = M_adjoint*R, on even sites
   reinitialize after niters iterations and try once more.
   parity=EVEN = do only even sites, parity=ODD = do odd sites,
   parity=EVENANDODD = do all sites
*/
#define _GNU_SOURCE
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include <lattice_qdp.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

extern int QDP_Debug;

void
print_mem(void)
{
  if(QDP_this_node==0) {
    int pid;
    char *s=NULL;

    //sleep(10);
    pid = getpid();
    printf("vsize = ");
    fflush(stdout);
    asprintf(&s, "ps -o vsz= -p %i", pid);
    system(s);
    free(s);  s = NULL;
  }
}

static int congrad_setup=0;
static QDP_ColorVector *ttt, *tttt, *resid, *cg_p;
static QDP_ColorVector *temp1[16], *temp2[16];
extern QDP_ColorMatrix *bcklink[8];

void setup_dslash(void);
void unset_dslash(void);

static void
setup_congrad(void)
{
  if(!congrad_setup) {
    int i;
    setup_dslash();
    ttt = QDP_create_V();
    tttt = QDP_create_V();
    resid = QDP_create_V();
    cg_p = QDP_create_V();
    for(i=0; i<16; i++) {
      temp1[i] = QDP_create_V();
      temp2[i] = QDP_create_V();
    }
    congrad_setup = 1;
  }
}

static void
unset_congrad(void)
{
  if(congrad_setup) {
    int i;
    unset_dslash();
    QDP_destroy_V(ttt);       ttt = NULL;
    QDP_destroy_V(tttt);      tttt = NULL;
    QDP_destroy_V(resid);     resid = NULL;
    QDP_destroy_V(cg_p);      cg_p = NULL;
    for(i=0; i<16; i++) {  
      QDP_destroy_V(temp1[i]); temp1[i] = NULL;
      QDP_destroy_V(temp2[i]); temp2[i] = NULL;
    }
    congrad_setup = 0;
  }
}

int
ks_congrad_qdp(QDP_ColorVector *src, QDP_ColorVector *dest, QLA_Real mass,
	       int niter, int nrestart, QLA_Real rsqmin, QDP_Subset parity,
	       QLA_Real *final_rsq_ptr)
{
  QLA_Real a,b;		 /* Sugar's a,b,resid**2,last resid*2 */
  QLA_Real rsq,oldrsq,pkp; /* pkp = cg_p.K.cg_p */
  QLA_Real msq_x4;	 /* 4*mass*mass */
  QLA_Real source_norm;	 /* squared magnitude of source vector */
  QLA_Real rsqstop;	 /* stopping residual normalized by source norm */
  QDP_Subset q_parity, q_otherparity;
  int iteration;	 /* counter for iterations */

/* Timing */
#ifdef CGTIME
  double dtimec;
  double nflop;
  nflop = 1187;
  if(parity==QDP_all) nflop *= 2;
#endif

  setup_congrad();
  //print_mem();

  if(parity==QDP_odd) {
    q_parity = QDP_odd;
    q_otherparity = QDP_even;
  } else if(parity==QDP_even) {
    q_parity = QDP_even;
    q_otherparity = QDP_odd;
  } else {
    q_parity = QDP_all;
    q_otherparity = QDP_all;
  }
  msq_x4 = 4.0*mass*mass;
  iteration = 0;

  load_ferm_links();
  set4_M_from_field(fatlinks, t_fatlink, EVENANDODD);
  set4_M_from_field(longlinks, t_longlink, EVENANDODD);

  //#if 0
  {
    QDP_ColorMatrix *tcm;
    int i;
    //QDP_Debug=9;
    tcm = QDP_create_M();
    for(i=0; i<8; ++i) {
      QDP_M_eq_sM(tcm, implinks[i], shiftdirs[i], QDP_backward, QDP_all);
      QDP_M_eqm_Ma(bcklink[i], tcm, QDP_all);
    }
    QDP_destroy_M(tcm); tcm = NULL;
  }
  //#endif

#if 0
  {
    int i;
    for(i=0; i<8; ++i) {
      QDP_M_eq_sM(bcklink[i], implinks[i], shiftdirs[i], QDP_backward, QDP_all);
    }
  }
#endif
  //QDP_set_block_size(256);

#if 0
    //QDP_Debug=0;
    for(i=0; i<8; ++i) {
      QDP_destroy_M(implinks[i]); implinks[i] = NULL;
    }
    for(i=0; i<8; ++i) {
      implinks[i] = QDP_create_M();
    }
    set4_M_from_temp(fatlinks, t_fatlink, EVENANDODD);
    set4_M_from_temp(longlinks, t_longlink, EVENANDODD);
  }
  //print_mem();
#endif

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  /* source_norm = |src|^2 */
  QDP_r_eq_norm2_V(&source_norm, src, q_parity);
  rsqstop = rsqmin * source_norm;

  do {
    /* initialization process
       ttt <- (-1)*M_adjoint*M*dest
       resid <- src + ttt
       rsq = |resid|^2              */

    QDP_V_eq_V(cg_p, dest, q_parity);

    dslash_qdp_fn_special2(cg_p, tttt, q_otherparity, temp1);
    dslash_qdp_fn_special2(tttt, ttt, q_parity, temp2);
    QDP_V_meq_r_times_V(ttt, &msq_x4, cg_p, q_parity);
    iteration++;    /* iteration counts multiplications by M_adjoint*M */

    QDP_V_eq_V_plus_V(resid, src, ttt, q_parity);
    QDP_r_eq_norm2_V(&rsq, resid, q_parity);

    if(rsq>rsqstop) {

      QDP_V_eq_V(cg_p, resid, q_parity);

      /* main loop - do until convergence or time to restart */
      /*
	oldrsq <- rsq
	ttt <- (-1)*M_adjoint*M*cg_p
	pkp <- (-1)*cg_p.M_adjoint*M.cg_p
	a <- -rsq/pkp
	dest <- dest + a*cg_p
	resid <- resid + a*ttt
	rsq <- |resid|^2
	b <- rsq/oldrsq
	cg_p <- resid + b*cg_p
      */
      //print_mem();
      while(1) {
	oldrsq = rsq;

	dslash_qdp_fn_special2(cg_p, tttt, q_otherparity, temp1);
	dslash_qdp_fn_special2(tttt, ttt, q_parity, temp2);
	QDP_V_meq_r_times_V(ttt, &msq_x4, cg_p, q_parity);
	iteration++;

	/* pkp <- cg_p . ttt */
	QDP_r_eq_re_V_dot_V(&pkp, cg_p, ttt, q_parity);

	a = rsq / pkp;

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
	QDP_V_meq_r_times_V(dest, &a, cg_p, q_parity);
	QDP_V_meq_r_times_V(resid, &a, ttt, q_parity);
	QDP_r_eq_norm2_V(&rsq, resid, q_parity);

	if( (rsq<rsqstop) || (iteration%niter==0) ) break;

	b = rsq / oldrsq;
	/* cg_p  <- resid + b*cg_p */
	QDP_V_eq_r_times_V_plus_V(cg_p, &b, cg_p, resid, q_parity);
      }
    }

  } while( (rsq>rsqstop) && (iteration<nrestart*niter) );

  if( rsq <= rsqstop ) {
#ifdef CGTIME
    dtimec += dclock();
    if(QDP_this_node==0) {
      printf("CONGRAD5: time = %e (fn_alt1_qdp) masses = 1 iters = %d mflops = %e\n",
	     dtimec, iteration,
	     (double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);
    }
#endif
  } else {
    if(QDP_this_node==0) {
      printf("ks_congrad_qdp: CG not converged after %d iterations, res. = %e wanted %e\n",
	     iteration, rsq, rsqstop);
      fflush(stdout);
    }
  }
  *final_rsq_ptr = rsq;
  total_iters += iteration;
  //print_mem();
  //unset_congrad();
  return(iteration);
}


/* prec argument is ignored */
int
ks_congrad(field_offset f_src, field_offset f_dest, Real mass,
  int niter, int nrestart, Real rsqmin, int prec, int parity, 
  Real *final_rsq_ptr)
{
  QLA_Real qmass, qrsqmin, qfinal_rsq_ptr;
  QDP_ColorVector *src, *dest;
  QDP_Subset q_parity;
  int iteration;

  print_mem();

  switch(parity) {
  case(EVEN):  q_parity = QDP_even; break;
  case(ODD):   q_parity = QDP_odd;  break;
  default:     q_parity = QDP_all;  break;
  }

  src = QDP_create_V();
  dest = QDP_create_V();

  set_V_from_site(src, f_src,EVENANDODD);
  set_V_from_site(dest, f_dest,EVENANDODD);

  qmass = (QLA_Real) mass;
  qrsqmin = (QLA_Real) rsqmin;
  iteration = ks_congrad_qdp(src, dest, qmass, niter, nrestart, 
    qrsqmin, q_parity, &qfinal_rsq_ptr);
  *final_rsq_ptr = (Real) qfinal_rsq_ptr;

  set_site_from_V(f_dest, dest, parity);

  QDP_destroy_V(dest); dest = NULL;
  QDP_destroy_V(src);  src  = NULL;

  print_mem();
  if(QDP_this_node==0) printf("end\n");

  return(iteration);
}

/***************************/
/* original MILC functions */
/***************************/

#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

/* clear an su3_vector in the lattice */
void clear_latvec(field_offset v,int parity){
register int i,j;
register site *s;
register su3_vector *vv;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case ODD: FORODDSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
    } 
}

/* copy an su3_vector in the lattice */
void copy_latvec(field_offset src,field_offset dest,int parity){
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case ODD: FORODDSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
    } 
}

/* scalar multiply and add an SU3 vector in the lattice */
void scalar_mult_add_latvec(field_offset src1,field_offset src2,
			    Real scalar,field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
               spt1 = (su3_vector *)F_PT(s,src1);
                spt2 = (su3_vector *)F_PT(s,src2);
                dpt = (su3_vector *)F_PT(s,dest);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV( (su3_vector *)F_PT((s+FETCH_UP),src1),
				(su3_vector *)F_PT((s+FETCH_UP),src2),
				(su3_vector *)F_PT((s+FETCH_UP),dest) );
		}
                scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt);
       } END_LOOP
}

void scalar2_mult_add_su3_vector(su3_vector *a, Real s1, su3_vector *b, 
				 Real s2, su3_vector *c){
register int i;
    for(i=0;i<3;i++){
        c->c[i].real = s1*a->c[i].real + s2*b->c[i].real;
        c->c[i].imag = s1*a->c[i].imag + s2*b->c[i].imag;
    }
}

/* scalar multiply two SU3 vector and add through the lattice */
void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt  = (su3_vector *)F_PT(s,dest);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV((su3_vector *)F_PT((s+FETCH_UP),src1),
			       (su3_vector *)F_PT((s+FETCH_UP),src2),
			       (su3_vector *)F_PT((s+FETCH_UP),dest) );
		}
		scalar2_mult_add_su3_vector( spt1, scalar1, spt2, scalar2, dpt);
       } END_LOOP
}

/* scalar multiply an SU3 vector in the lattice */
void scalar_mult_latvec(field_offset src,Real scalar,
			field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case ODD: FORODDSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
    } 
}
