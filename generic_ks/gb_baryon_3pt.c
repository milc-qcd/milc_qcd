#include "generic_ks_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>
/*
 Unit test -
  Source should be point-split point source
  Sink should be only single point sources in the positive spatial octant, not point split!
  Only ties up a single site at the sink
  Used with local scalar-scalar spin-taste current, any baryon
  Tests gauge invariance
  --> Both coulomb and no gauge fix should give same answer
 Methodology -
  Only use one link direction at sink to avoid multi-site source objects (i.e. cross terms)
  Use sources defined on single sites to construct an object which has all the correct links
  We don't care about what the correlator actually is - we only want a gauge invariant object!
  This tests gauge invariance of source/sink objects
  Use only GAUGE_INV_SRCSNK with local current to test gauge invariance of source/sink
  Use both SRCSNK and CUR to test gauge invariance of current insertion
  When using GAUGE_INV_CUR, use the NEGATIVE_SNK_*_LINK to reverse all the links
    which are used in point-splitting the current insertion
    (i.e. a G1-GX current should use NEGATIVE_SNK_X_LINK only)
 IMPORTANT - if using mmap, need to uncomment same defines in gb_baryon_mmap.c
 Example tests:
  16+ S2-S2, with G1-G1  current
  16+ S2-S2, with G1-GX  current
  16+ S2-S3, with G1-GXY current
  16+ S2-S3, with G3-G5T current
*/
//#define UNIT_TEST_GAUGE_INV_SRCSNK
//#define UNIT_TEST_GAUGE_INV_CUR
//#define UNIT_TEST_NEGATIVE_SNK_X_LINK
//#define UNIT_TEST_NEGATIVE_SNK_Y_LINK
//#define UNIT_TEST_NEGATIVE_SNK_Z_LINK

/* Unit test - only sum over current insertion origin */
//#define UNIT_TEST_CURRENT_ORIGIN

/* Unit test - move the correlator pieces around within unit cube */
//#define UNIT_TEST_TRANSLATE_SOURCE 6
//#define UNIT_TEST_TRANSLATE_SINK 0
//#define UNIT_TEST_TRANSLATE_CURRENT 6

//#define NO_SINK_LINKS

//#define UNIT_TEST_ORIGIN_TRANSLATE 7

// objects mapped to memory
//#ifdef GB_SPEC_MAP
//static su3_vector *spec_map = NULL; // hold spectator contribution
//static short      *spec_fill= NULL; // hold bit to indicate whether slot is filled
//#endif
//#ifdef GB_INT_MAP
//#define MAX_INT_MAP_CURR 16
//static su3_vector **int_map [MAX_INT_MAP_CURR]= {NULL};// hold interacting contribution
//static int          int_fill[MAX_INT_MAP_CURR]= {0};// spin-taste index identifier for each buffer
//static int      num_int_buff = 0; // integer to specify number of buffers
//#endif

#ifdef GB_BARYON
/*------------------------------------------------------------------*/
static void conj_v_field(su3_vector *dest,su3_vector *src){
  int i,c;
  site *s;
  #ifdef OMP
    #pragma omp parallel for \
      shared(src, dest, lattice, sites_on_node) \
      private(i,s,c) \
      schedule(static)
  #endif // OMP
  for(i=0;i<sites_on_node;i++) {
   s = lattice + i;
   for(c=0;c<3;c++){
    CONJG(src[i].c[c],dest[i].c[c]);
   }
  }
}

/*------------------------------------------------------------------*/
/**
   Symmetric shift a quark at the sink.
   Normal operating direction is for doBW = 0x0, but doBW = 0x1 is
   needed to link with the interacting quark. Links are applied
   to the spectator pair rather than the interacting quark to prevent
   the need for a factor of 8 more inversions.
 */
static void
sym_shift_3pt(int dir, short doBW, su3_vector *dest, su3_vector *src, su3_matrix *links)
{
  register int i;
  msg_tag *tag[2] = {NULL};
  su3_vector *cvec1;

  // shifting in dir moves +dir to 0
  tag[0] = start_gather_field(src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0]);
  #ifdef NO_SINK_LINKS
    cvec1 = src;
  #else
    cvec1 = create_v_field();
    FORALLFIELDSITES_OMP(i,) { 
      mult_adj_su3_mat_vec( links+4*i+dir, src+i, cvec1+i ); 
    } END_LOOP_OMP
  #endif // NO_SINK_LINKS

  #ifndef ONE_SIDED_SHIFT_GB
    // shifting in opp_dir moves 0 to +dir
    tag[1] = start_gather_field(cvec1, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1]);
  #endif
  wait_gather(tag[0]);
  #ifndef ONE_SIDED_SHIFT_GB
    wait_gather(tag[1]);
  #endif // ONE_SIDED_SHIFT_GB

  FORALLFIELDSITES_OMP(i,) {
    #ifdef NO_SINK_LINKS 
      su3vec_copy( (su3_vector *)gen_pt[0][i], dest+i); 
    #else
      mult_su3_mat_vec( links+4*i+dir, (su3_vector*)gen_pt[0][i], dest+i );
    #endif // NO_SINK_LINKS
    #ifndef ONE_SIDED_SHIFT_GB
      add_su3_vector(dest+i, (su3_vector*)gen_pt[1][i], dest+i ); 
      scalar_mult_su3_vector( dest+i, .5, dest+i ); 
    #endif // ONE_SIDED_SHIFT_GB
  } END_LOOP_OMP

  cleanup_gather(tag[0]);
  #ifndef ONE_SIDED_SHIFT_GB
    cleanup_gather(tag[1]);
  #endif // ONE_SIDED_SHIFT_GB

  #ifndef NO_SINK_LINKS
    destroy_v_field(cvec1);
  #endif
}

/*------------------------------------------------------------------*/
static void
apply_sym_shift_3pt(int n, int *d, int *r0, short doBW, ks_prop_field *dest,
                    ks_prop_field *src, su3_matrix *links){
  /* Apply the symmetric shifts listed in d             *
   * Similar to zeta_shift_field in spin_taste_ops.c    *
   * but with different sign factors                    */
  int i,c;
  ks_prop_field *tmp0 = create_ksp_field(3);

  for(i=0;i<n;i++)
  {
    /* Do the shift in d[i] */
    if(i==0)
    {
      /* first time from source */
      for(c=0;c<3;c++) sym_shift_3pt(d[i], doBW, dest->v[c], src->v[c], links);
    }
    else
    {
      /* other times from dest */
      for(c=0;c<3;c++) sym_shift_3pt(d[i], doBW, tmp0->v[c], dest->v[c], links);
      copy_ksp_field(dest,tmp0);
    }
  }
  destroy_ksp_field(tmp0);
}

/*------------------------------------------------------------------*/
void
apply_par_xport_3pt(ks_prop_field *dest, ks_prop_field *src,
                    int n, int dir[], int r0[], short doBW, su3_matrix *links){
  double start = dclock();
  site *s;
  int i,j,c;
  // the specific case we care about has no sink links and n = 0 or 1
  #ifdef NO_SINK_LINKS
  assert(n <= 1);
  #else
  ks_prop_field *tvec0 = create_ksp_field(3);
  ks_prop_field *tvec1 = create_ksp_field(3);
  #endif
  int d[6][3] =
    {{XUP,YUP,ZUP},
     {YUP,ZUP,XUP},
     {ZUP,XUP,YUP},
     {XUP,ZUP,YUP},
     {YUP,XUP,ZUP},
     {ZUP,YUP,XUP}};

  if(n == 0){
    copy_ksp_field(dest,src);
  }
  #ifdef NO_SINK_LINKS
  else {
    apply_sym_shift_3pt(n,dir,r0,doBW,dest,src,links);
  }
  #else
  else if(n == 1){
    /* one link */
    d[0][0] = dir[0];
    apply_sym_shift_3pt(n,d[0],r0,doBW,dest,src,links);
  }
  else if (n == 2){
    /* two link */
    d[0][0] = dir[0]; d[0][1] = dir[1];
    d[1][1] = dir[0]; d[1][0] = dir[1];
    apply_sym_shift_3pt(n,d[0],r0,doBW,tvec0,src,links);
    apply_sym_shift_3pt(n,d[1],r0,doBW,tvec1,src,links);

    #pragma omp parallel for private(c,i) collapse(2)
    for(c=0;c<3;c++){
      for(i=0;i<sites_on_node;i++) {
        add_su3_vector( &tvec0->v[c][i], &tvec1->v[c][i], &dest->v[c][i] );
        scalar_mult_su3_vector( &dest->v[c][i], 0.5, &dest->v[c][i] );
      }
    } // END OMP LOOPS
  }
  else if (n == 3){
    /* three link */
    /* use the d given */  
    for(j=0;j<6;j++){
      apply_sym_shift_3pt(n,d[j],r0,doBW,tvec0,src,links);
      if(j==0) {
        copy_ksp_field(tvec1,tvec0);
      } else {
        #pragma omp parallel for private(c,i) collapse(2)
        for(c=0;c<3;c++){
          for(i=0;i<sites_on_node;i++) {
            add_su3_vector(&tvec1->v[c][i],&tvec0->v[c][i],&tvec1->v[c][i]);
          }
        } // END OMP LOOPS
      }
    }
    #pragma omp parallel for private(c,i) collapse(2)
    for(c=0;c<3;c++){
      for(i=0;i<sites_on_node;i++) {
        scalar_mult_su3_vector( &tvec1->v[c][i], 1./6., &dest->v[c][i] );
      }
    } // END OMP LOOPS
  } // END c==3
  destroy_ksp_field(tvec1);
  destroy_ksp_field(tvec0);
  #endif
  double end = dclock();
  }

#endif // GB_BARYON
