/******************* smearing.c ****************************************/

/* MIMD version 7 */
/* original version of 5/23/94 by UMH */
/* 2/19/98 Version 5 port CD */
/* 11/25/01 Removed redundant code. Changed to call ape_smear_dir CD */

/* Perform one iteration of "APE smearing" on the space-like links.
   For the projection back to SU(3) use multiples of 3 hits over the
   3 different SU(2) subgroups until the procedure converges.
   Note: only this complete projection makes the smearing truly
   gauge invariant. */

#include "../generic/generic_includes.h"
#include "hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#define TOL 1e-5
#define MAXCOUNT 50

void smearing( void ) {
  register int mu,i;
  register site *s;
  su3_matrix **apelinks=NULL;
  su3_matrix *src=NULL,*diag=NULL;
  int idir;
#ifdef STAP_APE
  int stap[4]={0,0,0,0};
  int ndirs=0;

#ifdef APE_4D_SMEARING
//  node0_printf("4D-");
  for ( mu = XUP; mu <= TUP; mu++ ) stap[xc[mu]] = 1;
#endif
#if defined APE_3D_SMEARING 
//  node0_printf("3D-");
  for ( mu = XUP; mu <= TUP; mu++ )  stap[xc[mu]] = ( cor_dir != xc[mu] ? 1 : 0 );
#endif
#ifdef APE_2D_SMEARING
//  node0_printf("2D-");
  for ( mu = XUP; mu <= TUP; mu++ ) stap[xc[mu]] = ( maxc[mu] > 0 && cor_dir != xc[mu] ? 1 : 0 );
#endif
#ifdef APE_1D_SMEARING
//  node0_printf("1D-");
  for ( mu = XUP; mu <= TUP; mu++ ) stap[xc[mu]] = ( xc[mu] == stap_dir && maxc[stap_dir] == 0 && cor_dir != xc[mu] ? 1 : 0 );
#endif
#ifdef APE_1D2_SMEARING 
//  node0_printf("Spatial and temporal 1D2-");
  for ( mu = XUP; mu <= TUP; mu++ ) stap[xc[mu]] = ( xc[mu] == stap_dir && maxc[stap_dir] == 0 ? 1 : 0 );
#endif
  for ( mu = XUP; mu <= TUP; mu++ ) ndirs += stap[xc[mu]];
#ifdef APE_1D2_SMEARING 
  ndirs++;
#endif
#else
#ifdef APE_4D_SMEARING
//  node0_printf("4D-APE smear\n");
  int ndirs = 4;
#endif
#ifdef APE_3D_SMEARING
//  node0_printf("Spatial 3D-APE smear\n");
  int ndirs = 3;
#endif
#endif

#if ( (!(defined COULOMB) && (defined APE_4D_SMEARING || defined APE_1D2_SMEARING ) ) )
  if ( fixflag == AXIAL_GAUGE_FIX ) { 
#define SET2SU3UNITMAT( mat ) \
{ clear_su3mat(&(mat)); \
mat.e[0][0].real=1.; \
mat.e[1][1].real=1.; \
mat.e[2][2].real=1.; }
  /* if in axial gauge, then first replace dummy temporal links by unit matrix */
    FORALLSITES_OMP(i,s, default(shared) ){ if ( site_coord(s,xc[TUP]) < nc[TUP]-1) { SET2SU3UNITMAT( s->link[xc[TUP]] ); } } END_LOOP_OMP;
  }
#endif

  /* Allocate space for temporary storage of results */
  src  = create_G_from_site();
  diag = hqp_alloc_su3mat_buffer(4);
  apelinks = (su3_matrix **)malloc(ndirs*sizeof(su3_matrix*));
  assert( ( apelinks !=NULL ) );
  for ( idir = 0; idir < ndirs; idir++ )  
    apelinks[idir] =hqp_alloc_su3mat_buffer(1); 

  /* Loop over the space directions and APE-smear the links.
     The results will temporarily be stored in x_link, y_link
     and z_link for mu=XUP, YUP and ZUP and copied do link[mu]
     at the end. */
  
#if (defined APE_4D_SMEARING || defined APE_1D2_SMEARING)
  for( idir = 0, mu = XUP; mu <= TUP; mu++ ){ 
#else
  for( idir = 0, mu = XUP; mu <= ZUP; mu++ ){
#endif
    if ( idir >= ndirs || maxc[mu] == 0 ) continue;
    /* Do the APE smearing and SU(3) projection for all links in mu */
#ifdef STAP_APE
    if ( stap_dir == xc[mu] ) continue;
    ape_smear_field_dir_stap(src,xc[mu],diag,
		  1./smear_fac,1., stap ,3*MAXCOUNT,TOL);
#else
#ifdef APE_4D_SMEARING
    ape_smear_field_dir(src,xc[mu],diag,
		  1./smear_fac,1., 0 ,3*MAXCOUNT,TOL);
#elif (defined APE_3D_SMEARING)
    ape_smear_field_dir(src,xc[mu],diag,
		  1./smear_fac,1., 1 ,3*MAXCOUNT,TOL);
#endif
#endif
    
    /* Storage management: Temporarily store the new link */
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(diag[xc[mu]+4*i]), &(apelinks[idir][i]) ); } END_LOOP_OMP;
    
    idir++; 
  } /* end loop over mu */
  
  /* Now copy the smeared links, overwriting the previous ones */
#if (defined APE_4D_SMEARING || defined APE_1D2_SMEARING)
  for( idir = 0, mu = XUP; mu <= TUP; mu++ ){ 
#else
  for( idir = 0, mu = XUP; mu <= ZUP; mu++ ){
#endif
    if ( idir >= ndirs || maxc[mu] == 0 ) continue;
#ifdef STAP_APE
    if ( stap_dir == xc[mu] ) continue;
#endif
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(apelinks[idir][i]), &(s->link[xc[mu]]) ); } END_LOOP_OMP;
    idir++; 
  } /* end loop over mu */

  hqp_free_su3mat_buffer( src );
  hqp_free_su3mat_buffer( diag );
  for ( idir = 0; idir < ndirs; idir++ ) { hqp_free_su3mat_buffer(apelinks[idir]); }
  free(apelinks);

#if ( (defined AX_GAUGE && (defined APE_4D_SMEARING || defined APE_1D2_SMEARING ) ) )
  /* restore axial gauge if needed */
  if ( fixflag == AXIAL_GAUGE_FIX ) {
    //node0_printf("Recover axial gauge buffers\n");
    ax_gauge(); 
  }
#endif

} /* smearing */
