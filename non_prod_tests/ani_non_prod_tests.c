/********** non_prod_tests.c *********************************************/
/* MIMD version 7 */

/************************************************************
 * This library contains non-production level tests for verifying the 
 * correct functionality of the MILC code. They usually assume unit 
 * links or specific field configurations designed for the tests.
 *
 * List of contents:
 * void free_KS_ani_test( ... )
 *
 * The test automatically fails wth longlink actions for nx,ny,nz<6 
 * (due to mixing between 1-link and 3-link terms into the same sites) 
 * and for nt <8 (due to antiperiodic boundary conditions).
 *
 ***********************************************************/

#include "ani_non_prod_tests_includes.h"
#include <assert.h>

#ifdef ANISOTROPY 

#ifdef FREE_KS_ANI_TEST
void free_KS_ani_test( su3_vector **src,
                       quark_source *my_ksqs,
                       ks_param my_ksp[],
                       int num_prop,
                       su3_vector **dst ) { 
  site *s; 
  int i,j;
  int isosites=0, anisites=0;;
  int r0[4],dr[4],parity;
  double isosum=0.,anisum=0.;
  imp_ferm_links_t **fn = NULL;
  imp_ferm_links_t **fn_multi = NULL;

  /* Restore fermion links */
  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
  fn = get_fm_links(fn_links);
  /* Copy pointers for fermion links, based on Naik epsilon indices */
  fn_multi = (imp_ferm_links_t **) malloc(sizeof(imp_ferm_links_t *)*num_prop);
  for(j = 0; j < num_prop; j++) { fn_multi[j] = fn[my_ksp[j].naik_term_epsilon_index]; }

  r0[XUP] = my_ksqs->x0;
  r0[YUP] = my_ksqs->y0;
  r0[ZUP] = my_ksqs->z0;
  r0[TUP] = my_ksqs->t0;
  parity=(r0[XUP]+r0[YUP]+r0[ZUP]+r0[TUP])%2;

  for(j = 0; j < num_prop; j++ ){
    if ( my_ksqs->type == POINT ) { 
node0_printf("###################################################\n\
Free Field Theory Anisotropy Test for set %d\n\
###################################################\n",j);
      fflush(stdout);

      ks_dirac_op( src[j], dst[j], (j+1.), parity , fn_multi[j]);
//      ks_dirac_op( src[j], dst[j], (j+1.), EVENANDODD , fn_multi[j]);

      FORALLSITES_OMP(i,s, private(dr) reduction(+:isosum,anisum) ){
        dr[XUP]= ( s->x-r0[XUP] <= nx/2 && (s->x-r0[XUP]) >= -nx/2 ? (s->x-r0[XUP]) : 
                   s->x-r0[XUP] >= nx/2 ?  (s->x-r0[XUP])-nx : -nx -(s->x-r0[XUP]) );
        dr[YUP]= ( s->y-r0[YUP] <= ny/2 && (s->y-r0[YUP]) >= -ny/2 ? (s->y-r0[YUP]) : 
                   s->y-r0[YUP] >= ny/2 ?  (s->y-r0[YUP])-ny : -ny -(s->y-r0[YUP]) );
        dr[ZUP]= ( s->z-r0[ZUP] <= nz/2 && (s->z-r0[ZUP]) >= -nz/2 ? (s->z-r0[ZUP]) : 
                   s->z-r0[ZUP] >= nz/2 ?  (s->z-r0[ZUP])-nz : -nz -(s->z-r0[ZUP]) );
        dr[TUP]= ( s->t-r0[TUP] <= nt/2 && (s->t-r0[TUP]) >= -nt/2 ? (s->t-r0[TUP]) : 
                   s->t-r0[TUP] >= nt/2 ?  (s->t-r0[TUP])-nt : -nt -(s->t-r0[TUP]) );
        if ( dr[XUP]*dr[XUP] + dr[YUP]*dr[YUP] + dr[ZUP]*dr[ZUP] + dr[TUP]*dr[TUP] == 1 
        || ( dr[XUP]*dr[XUP] == 9 && dr[YUP]*dr[YUP] + dr[ZUP]*dr[ZUP] + dr[TUP]*dr[TUP] == 0 ) 
        || ( dr[YUP]*dr[YUP] == 9 && dr[XUP]*dr[XUP] + dr[ZUP]*dr[ZUP] + dr[TUP]*dr[TUP] == 0 ) 
        || ( dr[ZUP]*dr[ZUP] == 9 && dr[XUP]*dr[XUP] + dr[YUP]*dr[YUP] + dr[TUP]*dr[TUP] == 0 ) 
        || ( dr[TUP]*dr[TUP] == 9 && dr[XUP]*dr[XUP] + dr[YUP]*dr[YUP] + dr[ZUP]*dr[ZUP] == 0 ) 
           ) {
#ifdef ANI_VERBOSE
          printf("\n###################################################\n\
                  \nnode %d, dx dy dz dt %d %d %d %d mass %d\n\n",mynode(),
                  -dr[XUP],-dr[YUP],-dr[ZUP],-dr[TUP],(j+1)); fflush(stdout);
#endif
          if ( dr[ani_dir]*dr[ani_dir] > 0 ) { 
#ifndef ONEDIM_ANISO_TEST
#ifdef ANI_VERBOSE
            printf("Rescale back with inverse anisotropy 1/xi=%.3f\n\n",1./ani_xiq);
#endif
            scalar_mult_su3_vector( &(dst[j][i]),1./ani_xiq,&(dst[j][i]) );
            if ( (double)magsq_su3vec( &(dst[j][i]) )>0 ) {
#else 
#ifdef ANI_VERBOSE
            printf("Rescale with opposite anisotropic factor xi=%.3f\n\n",iso_xiq);
#endif
            scalar_mult_su3_vector( &(dst[j][i]),iso_xiq,&(dst[j][i]) );
            if ( (double)magsq_su3vec( &(dst[j][i]) )>0 || iso_xiq==0) {
#endif
              anisum+=(double)magsq_su3vec( &(dst[j][i]) );
              anisites++;
            }
          }
          else { 
#ifndef ONEDIM_ANISO_TEST
            if ( (double)magsq_su3vec( &(dst[j][i]) )>0 ) {
#else
#ifdef ANI_VERBOSE
            printf("Rescale with opposite anisotropic factor xi=%.3f\n\n",ani_xiq);
#endif
            if ( (double)magsq_su3vec( &(dst[j][i]) )>0 || iso_xiq==0 ) {
              scalar_mult_su3_vector( &(dst[j][i]),ani_xiq,&(dst[j][i]) );
#endif
              isosum+=(double)magsq_su3vec( &(dst[j][i]) );
              isosites++;
            }
          }
#ifdef ANI_VERBOSE
          dumpvec(&(dst[j][i])); fflush(stdout);
          printf("\n###################################################\n"); fflush(stdout);
#endif
        }
      } END_LOOP_OMP   
      g_doublesum( &isosum);  
      g_doublesum( &anisum);  
#ifndef ONEDIM_ANISO_TEST
      node0_printf("Avg of the unresc mag at %2d 3d-isotropic sites =%.12f\n",isosites,isosum/isosites);
      node0_printf("Avg of 1./xi-resc mag at %2d  anisotropic sites =%.12f\n",anisites,anisum/anisites);
#else
      node0_printf("Avg of opp fact resc  mag at %2d 3d-isotropic sites =%.12f\n",isosites,isosum/isosites);
      node0_printf("Avg of opp fact resc  mag at %2d  anisotropic sites =%.12f\n",anisites,anisum/anisites);
#endif
      node0_printf("###################################################\n"); fflush(stdout);
      assert( (isosum/3.-anisum)*(isosum/3.-anisum) < 1.e-24 );
      node0_printf("The anisotropic Dirac operator behaves correctly for free fields with a point source.\n\
###################################################\n\n"); fflush(stdout);
    }   
  }   

  if(fn_multi != NULL)free(fn_multi);
}
#endif

#ifdef FUNNYLINKS
void funnylinks( Real umu[4] ){
  int mu, i;
  site *s=NULL;

  node0_printf("umu %.2f %.2f %.2f %.2f\n",umu[XUP],umu[YUP],umu[ZUP],umu[TUP]);
  
  FORALLSITES_OMP(i,s, private(mu) ){
    for (mu=XUP; mu<=TUP; mu++) scalar_mult_su3_matrix( &(s->link[mu]), umu[mu],  &(s->link[mu]) );
  } END_LOOP_OMP;
}

void noreunit( void ) { };
#endif

#endif
