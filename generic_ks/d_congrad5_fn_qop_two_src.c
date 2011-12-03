/******* d_congrad5_fn_qop_two_src.c - conjugate gradient for SU3/fermions **/
/* MIMD version 7 */

/* This is the two-source MILC wrapper for the SciDAC Level 3 QOP inverter 
   using the raw interface */
/* 2/2005 D. Renner and C. Jung */
/* 5/2005 C. DeTar two source version eliminates one remapping */
/* 9/2005 C. DeTar converted to C code */
/* 12/2005 C. DeTar upgraded to new Level 3 API */
/* 12/07/2006 C. DeTar mixed precision conversion */

#include "generic_ks_includes.h"
#include "../include/generic_ks_qop.h"

/* The standard MILC interface for two sources */


static int 
ks_congrad_two_src_F(	/* Return value is number of iterations taken */
    field_offset src1,     /* source vector (type su3_vector) */
    field_offset src2,
    field_offset sol1,	/* solution vectors */
    field_offset sol2,
    quark_invert_control *qic,
    Real mass1,
    Real mass2,
    imp_ferm_links_t *fn 
    )
{
  int iterations_used;
  static float t_mass1;
  static float t_mass2;
  float *masses[2];
  int nmass[2], nsrc;
  su3_vector *srcs[2], *sols0[1], *sols1[1], **sols[2];

  /* Map masses and source fields from site structure to temporary fields */

  nsrc = 2;
  nmass[0] = 1;
  nmass[1] = 1;
  t_mass1 = mass1;
  t_mass2 = mass2;  
  masses[0] = &t_mass1;
  masses[1] = &t_mass2;
  srcs[0] = create_v_field_from_site_member(src1);
  srcs[1] = create_v_field_from_site_member(src2);

  /* Make room for temporary solution fields */

  sols0[0] = create_v_field();
  sols1[0] = create_v_field();
  sols[0] = sols0;
  sols[1] = sols1;

  iterations_used = 
    ks_congrad_qop_F_field2field( qic, masses, nmass, srcs,
				  sols, nsrc, fn );

  /* Copy solutions to site structure */

  copy_site_member_from_v_field(sol1, sols0[0]);
  copy_site_member_from_v_field(sol2, sols1[0]);

  /* Cleanup */

  destroy_v_field(sols0[0]);
  destroy_v_field(sols1[0]);
  destroy_v_field(srcs[1]);
  destroy_v_field(srcs[0]);

  return iterations_used;
}

static int 
ks_congrad_two_src_D(	/* Return value is number of iterations taken */
    field_offset src1,     /* source vector (type su3_vector) */
    field_offset src2,
    field_offset sol1,	/* solution vectors */
    field_offset sol2,
    quark_invert_control *qic,
    Real mass1,
    Real mass2,
    imp_ferm_links_t *fn 
    )
{
  int iterations_used;
  static double t_mass1;
  static double t_mass2;
  double *masses[2];
  int nmass[2], nsrc;
  su3_vector *srcs[2], *sols0[1], *sols1[1], **sols[2];

  /* Map masses and source fields from site structure to temporary fields */

  nsrc = 2;
  nmass[0] = 1;
  nmass[1] = 1;
  t_mass1 = mass1;
  t_mass2 = mass2;  
  masses[0] = &t_mass1;
  masses[1] = &t_mass2;
  srcs[0] = create_v_field_from_site_member(src1);
  srcs[1] = create_v_field_from_site_member(src2);

  /* Make room for temporary solution fields */

  sols0[0] = create_v_field();
  sols1[0] = create_v_field();
  sols[0] = sols0;
  sols[1] = sols1;

  iterations_used = 
    ks_congrad_qop_D_field2field( qic, masses, nmass, srcs,
				  sols, nsrc, fn );

  /* Copy solutions to site structure */

  copy_site_member_from_v_field(sol1, sols0[0]);
  copy_site_member_from_v_field(sol2, sols1[0]);

  /* Cleanup */

  destroy_v_field(sols0[0]);
  destroy_v_field(sols1[0]);
  destroy_v_field(srcs[1]);
  destroy_v_field(srcs[0]);

  return iterations_used;
}

int
ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type su3_vector) */
    field_offset src2,
    field_offset dest1,	/* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,		/* maximal number of CG interations */
    int nrestart,       /* maximal number of CG restarts */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for the inversion (ignored) */
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr, /* final residue squared */
    imp_ferm_links_t *fn       /* Storage for fermion links */
    )
{
  int iterations_used;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;  /* Currently ignored */
  qic.min       = 0;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.parity    = parity;
  qic.start_flag = 0;
  qic.nsrc      = 1;
  qic.resid     = sqrt(rsqmin);
  qic.relresid  = 0;     /* Suppresses this test */

  if(PRECISION == 1)
    iterations_used = ks_congrad_two_src_F(src1, src2, dest1, 
					   dest2, &qic, mass1, mass2, fn);
  else
    iterations_used = ks_congrad_two_src_D(src1, src2, dest1, 
					   dest2, &qic, mass1, mass2, fn);
  report_status(&qic);
  report_status(&qic);  // So regression tests see the same number of lines
  return iterations_used;
}
