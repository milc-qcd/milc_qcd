#ifdef GB_BARYON
#include "generic_ks_includes.h"
#include "../include/generic_ks.h"

#define GB_DIRECTORY "/tmp"
//#define NO_SINK_LINKS

//#define DO_MSYNC

//#define UNIT_TEST_GAUGE_INV_SRCSNK
//#define UNIT_TEST_GAUGE_INV_CUR
//#define UNIT_TEST_NEGATIVE_SNK_X_LINK
//#define UNIT_TEST_NEGATIVE_SNK_Y_LINK
//#define UNIT_TEST_NEGATIVE_SNK_Z_LINK

static mmap_cache **gb_qk_cache = NULL;
static mmap_cache *gb_src_cache = NULL;

/** Solve_ksprop has no info about source number, so keep track of it manually */
static int gb_src_iter = 0;

void set_mmap_src_number(int num){
  gb_src_iter = num;
}

static void
sym_shift_mmap(int dir, su3_vector *dest, su3_vector *src, su3_matrix *links)
{
  register int i;
  msg_tag *tag[2] = {NULL};
  su3_vector *tvec  = create_v_field();

/** Looks like spaghetti code, but written this way so that
    the unit tests and the generic code are using the same
    set of functions, rather than an independent set (or as
    close as possible). That way if changes are made,
    the changes can be directly checked.

    SRCSNK unit test defaults always to the same path,
    CUR unit test chooses whether to take the default or the
     opposite path based on which NEGATIVE_SNK_*_LINKs are defined,
    and no unit test definitions always take both directions. */
// TODO: redo this, combine with other three copies in gb baryon files
#ifdef UNIT_TEST_GAUGE_INV_SRCSNK
#ifdef UNIT_TEST_GAUGE_INV_CUR
  switch(dir){
#ifdef UNIT_TEST_NEGATIVE_SNK_X_LINK
  case XUP:
#endif
#ifdef UNIT_TEST_NEGATIVE_SNK_Y_LINK
  case YUP:
#endif
#ifdef UNIT_TEST_NEGATIVE_SNK_Z_LINK
  case ZUP:
#endif
#ifdef NO_SINK_LINKS
    copy_v_field(tvec,src);
#else
    FORALLFIELDSITES(i) { mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec+i ); }
#endif // NO_SINK_LINKS
    tag[1] = start_gather_field( tvec, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1] );
    break;
  default:
#endif
    tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
#else
    tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );

#ifdef NO_SINK_LINKS
    copy_v_field(tvec,src);
#else
    FORALLFIELDSITES(i) { mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec+i ); }
#endif // NO_SINK_LINKS
    tag[1] = start_gather_field( tvec, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1] );
#endif // UNIT_TEST_GAUGE_INV_SRCSNK

#ifdef UNIT_TEST_GAUGE_INV_SRCSNK
#ifdef UNIT_TEST_GAUGE_INV_CUR
    break; }
  switch(dir){
#ifdef UNIT_TEST_NEGATIVE_SNK_X_LINK
  case XUP:
#endif
#ifdef UNIT_TEST_NEGATIVE_SNK_Y_LINK
  case YUP:
#endif
#ifdef UNIT_TEST_NEGATIVE_SNK_Z_LINK
  case ZUP:
#endif
    wait_gather(tag[1]);
    FORALLFIELDSITES(i) { su3vec_copy( (su3_vector*)gen_pt[1][i], dest+i ); }
    cleanup_gather(tag[1]);
    break;
  default:
#endif

    wait_gather(tag[0]);
#ifdef NO_SINK_LINKS
    FORALLFIELDSITES(i){ su3vec_copy( (su3_vector *)gen_pt[0][i], dest+i ); }
#else
    FORALLFIELDSITES(i){ mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i ); }
#endif // NO_SINK_LINKS
    cleanup_gather(tag[0]);
#ifdef UNIT_TEST_GAUGE_INV_CUR
    break; }
#endif
#else
    wait_gather(tag[0]);
#ifdef NO_SINK_LINKS
    FORALLFIELDSITES(i){ su3vec_copy( (su3_vector *)gen_pt[0][i], dest+i ); }
#else
    FORALLFIELDSITES(i){ mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i ); }
#endif // NO_SINK_LINKS
    cleanup_gather(tag[0]);

    wait_gather(tag[1]);
    FORALLFIELDSITES(i) { add_su3_vector(dest+i, (su3_vector*)gen_pt[1][i], dest+i ); }
    FORALLFIELDSITES(i) { scalar_mult_su3_vector( dest+i, .5, dest+i ) ; }
    cleanup_gather(tag[1]);
#endif // UNIT_TEST_GAUGE_INV_SRCSNK

  destroy_v_field(tvec);
}

/*------------------------------------------------------------------*/
static void
apply_sym_shift_mmap(int n, int *d, int *r0, ks_prop_field *dest,
                    ks_prop_field *src, su3_matrix *links){
  int i,c;
  ks_prop_field *tmp0 = create_ksp_field(3);

  for(i=0;i<n;i++)
  {
    /* Do the shift in d[i] */
    if(i==0)
    {
      /* first time from source */
      for(c=0;c<3;c++) sym_shift_mmap(d[i], dest->v[c], src->v[c], links);
    }
    else
    {
      /* other times from dest */
      for(c=0;c<3;c++) sym_shift_mmap(d[i], tmp0->v[c], dest->v[c], links);
      copy_ksp_field(dest,tmp0);
    }
  }
  destroy_ksp_field(tmp0);
}

/*------------------------------------------------------------------*/
void
apply_par_xport_mmap(ks_prop_field *dest, ks_prop_field *src,
                    int n, int dir[], int r0[], su3_matrix *links){
  site *s;
  int i,j,c;
  ks_prop_field *tvec0 = create_ksp_field(3);
  ks_prop_field *tvec1 = create_ksp_field(3);
  ks_prop_field *tsrc  = create_ksp_field(3);
  int d[6][3] =
    {{XUP,YUP,ZUP},
     {YUP,ZUP,XUP},
     {ZUP,XUP,YUP},
     {XUP,ZUP,YUP},
     {YUP,XUP,ZUP},
     {ZUP,YUP,XUP}};
  copy_ksp_field(tsrc,src);

  if(n == 0){
    copy_ksp_field(dest,src);
  }
  else if(n == 1){
    /* one link */
    d[0][0] = dir[0];
    apply_sym_shift_mmap(n,d[0],r0,dest,tsrc,links);
  }
  else if (n == 2){
    /* two link */
    d[0][0] = dir[0]; d[0][1] = dir[1];
    d[1][1] = dir[0]; d[1][0] = dir[1];
    apply_sym_shift_mmap(n,d[0],r0,tvec0,tsrc,links);
    apply_sym_shift_mmap(n,d[1],r0,tvec1,tsrc,links);

    for(c=0;c<3;c++){
      FORALLSITES(i,s){
        add_su3_vector( &tvec0->v[c][i], &tvec1->v[c][i], &dest->v[c][i] );
        scalar_mult_su3_vector( &dest->v[c][i], 0.5, &dest->v[c][i] );
      }
    }
  }
  else if (n == 3){
    /* three link */
    /* use the d given */
    for(j=0;j<6;j++){
      apply_sym_shift_mmap(n,d[j],r0,tvec0,tsrc,links);
      if(j==0) copy_ksp_field(tvec1,tvec0);
      else for(c=0;c<3;c++){
        FORALLFIELDSITES(i){
          add_su3_vector(&tvec1->v[c][i],&tvec0->v[c][i],&tvec1->v[c][i]);
        }
      }
    }
    for(c=0;c<3;c++){FORALLSITES(i,s){
      scalar_mult_su3_vector( &tvec1->v[c][i], 1./6., &dest->v[c][i] );
    }}
  }
  destroy_ksp_field(tvec1);
  destroy_ksp_field(tvec0);
  destroy_ksp_field(tsrc);
}

/** need special routine to set ks_prop_field so it doesn't try to reallocate space */
void fetch_ksp_from_cache(ks_prop_field **dest,int qknum,int scIdx,int skIdx){
  int c;
  (*dest) = (ks_prop_field *)malloc(sizeof(ks_prop_field));
  if((*dest) == NULL){
    node0_printf("fetch_ksp_from_cache(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  (*dest)->nc = 3;
  (*dest)->v = (su3_vector **) malloc(3*sizeof(su3_vector *));
  if((*dest)->v == NULL){
    node0_printf("fetch_ksp_from_cache(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  for(c=0;c<3;c++){
    (*dest)->v[c] = (su3_vector *) gb_qk_cache[qknum]->buffer[24*scIdx+3*skIdx+c];
  }
}

/** special routine to fetch a source su3_vector, no need to toss when done */
void fetch_su3v_from_cache(su3_vector **dest,int qknum,int skIdx,int disp,int color){
  *dest = (su3_vector *) gb_qk_cache[qknum]->buffer[24*skIdx+3*disp+color];
}

/** special routine to fetch a source su3_vector, no need to toss when done */
void fetch_src_from_cache(su3_vector **dest,int srcnum,int color){
  *dest = (su3_vector *) gb_src_cache->buffer[3*srcnum+color];
}

/** ensure that a su3v buffer is flushed before continuing */
void msync_su3v_from_cache(int qknum, int skIdx, int disp, int color){
  int ecode;
  node0_printf("Syncronizing su3_vector buffer (%d,%d,%d,%d)\n",qknum,skIdx,disp,color);
  if((ecode = msync_mmap_buffer(gb_qk_cache[qknum],24*skIdx+3*disp+color))){
    printf("Could not msync su3_vector buffer %d! Error %d on node %d\n",
     24*skIdx+3*disp+color,ecode,this_node);
    terminate(1);
  }
}

/** ensure that a su3v buffer is flushed before continuing */
void msync_src_from_cache(int srcnum, int color){
  int ecode;
  node0_printf("Syncronizing source buffer (%d,%d)\n",srcnum,color);
  if((ecode = msync_mmap_buffer(gb_src_cache,3*srcnum+color))){
    printf("Could not msync source buffer %d! Error %d on node %d\n",
     3*srcnum+color,ecode,this_node);
    terminate(1);
  }
}

/** ensure that a ksp buffer is flushed before continuing */
void msync_ksp_from_cache(int qknum, int skIdx, int disp){
  int c;
  //node0_printf("Syncronizing ksp buffer (%d,%d,%d)\n",qknum,skIdx,disp);
  for (c=0;c<3;c++){
    msync_su3v_from_cache(qknum,skIdx,disp,c);
  }
}

/** need special routine to free ks_prop_field so it doesn't try to free mmap */
void toss_ksp_from_cache(ks_prop_field **ksp){
  free((*ksp)->v);
  free((*ksp));
  *ksp = NULL;
}

/** fill a quark octet memory map with all 8 parallel transports for all 8 quarks*/
void populate_qk_oct_point_split(ks_prop_field **qko,int qknum,int r0[],su3_matrix *links){
  int scIdx,skIdx;
  int n [8];
  int dir[8][3] = {0};
  ks_prop_field * kspc [8];
  ks_prop_field *ksp = create_ksp_field(3);

  for(scIdx=1;scIdx<7;scIdx++){
   node0_printf("Creating mmap cache %d for quark object %d...\n",qknum,scIdx);
   for(skIdx=1;skIdx<7;skIdx++){
    n[skIdx] = singlet_index_to_disp(skIdx);
    singlet_index_to_dir(skIdx,dir[skIdx]);
    fetch_ksp_from_cache(kspc+skIdx,qknum,scIdx,skIdx); // get pointer
    // check to see if there is a previous partial solution
    // if some links are already smeared, we can use this as starting point
    int oldIdx = -1;
    if (n[skIdx] > 1) {
      oldIdx = skIdx - 1;
      for (; oldIdx >= 1; oldIdx --) {
        if (n[oldIdx] == n[skIdx] - 1) {
          bool dirsMatch = true;
          for (int i = 0; i < n[oldIdx]; i ++)
            dirsMatch &= (dir[oldIdx][i] == dir[skIdx][i]);
          if (dirsMatch) break;
        }
      }
    }
    if (oldIdx != -1) {
      int tempN = 1;
      int tempDir [3] = {dir[skIdx][n[skIdx]-1], 0, 0};
      apply_par_xport_3pt(kspc[skIdx],kspc[oldIdx],tempN,tempDir,r0,0x0,links);
    } else {
      apply_par_xport_3pt(kspc[skIdx],qko[scIdx],n[skIdx],dir[skIdx],r0,0x0,links);
    }

#ifdef DO_MSYNC
    msync_ksp_from_cache(qknum,scIdx,skIdx);
#endif
   }
   for(skIdx=1;skIdx<7;skIdx++)
    toss_ksp_from_cache(kspc+skIdx); // toss pointer
  }
  destroy_ksp_field(ksp);
  node0_printf("Done with mmap cache %d for quark object\n",qknum);
}

/** fill a source octet memory map with all 8 parallel transports for all 8 sources*/
void populate_sink_oct_point_split(su3_vector **sko,int qknum,int r0[],su3_matrix *links){
  int skIdx,disp;
  int n,c;
  int dir[3] = {0};
  ks_prop_field *kspc = NULL;
  ks_prop_field *ksp0 = create_ksp_field(3);
  ks_prop_field *ksp1 = create_ksp_field(3);

  for(skIdx=0;skIdx<8;skIdx++){
   node0_printf("Creating mmap cache for sink object %d...\n",skIdx);
   for(disp=0;disp<8;disp++){
    for(c=0;c<3;c++){
     copy_v_field(ksp0->v[c],sko[3*skIdx+c]);
    }
    n = singlet_index_to_disp(disp);
    singlet_index_to_dir(disp,dir);
    apply_par_xport_mmap(ksp1,ksp0,n,dir,r0,links);

    fetch_ksp_from_cache(&kspc,qknum,skIdx,disp); // get pointer, ksp pointers are okay
    copy_ksp_field(kspc,ksp1); // save once
#ifdef DO_MSYNC
    msync_ksp_from_cache(qknum,skIdx,disp);
#endif
    toss_ksp_from_cache(&kspc); // toss pointer
   }
  }
  destroy_ksp_field(ksp0);
  destroy_ksp_field(ksp1);
}

/** save the next source using the static iterator available */
void populate_next_mmap_src(su3_vector *src,int color){
  su3_vector *vec0 = NULL;
  fetch_src_from_cache(&vec0,gb_src_iter,color); // no check against going over max
  //print_unit_cube8("inSrc",src);
  copy_v_field(vec0,src);
#ifdef DO_MSYNC
  msync_src_from_cache(gb_src_iter,color);
#endif
  node0_printf("Saving source %d, color %d to mmap\n",gb_src_iter,color);
}

/** create a cache to store a quark octet*/
void create_qk_oct_cache(ks_prop_field **qko,int qknum,
                         int r0[],su3_matrix *links){
  /** create a cache for the octet
      there need to be three source color vector fields per prop field, all allocated separately
      8 sources x 8 point-splittings x 3 su3_vector fields
       = 192 su3_vector fields
       = 64 ks_prop_fields
   */

  /* malloc a cache so pointers can be copied rather than copying the object itself */
  gb_qk_cache[qknum] = (mmap_cache*)malloc(sizeof(mmap_cache));
  if(alloc_mmap_cache(gb_qk_cache[qknum],192,sites_on_node*sizeof(su3_vector),GB_DIRECTORY)){
    node0_printf("Failed to create gb baryon quark mmap cache!\n");
    terminate(1);
  }

  populate_qk_oct_point_split(qko,qknum,r0,links);
}

/** create a cache to store a source octet*/
void create_sink_oct_cache(su3_vector **sko,int qknum,
                         int r0[],su3_matrix *links){
  /* malloc a cache so pointers can be copied rather than copying the object itself */
  gb_qk_cache[qknum] = (mmap_cache*)malloc(sizeof(mmap_cache));
  if(alloc_mmap_cache(gb_qk_cache[qknum],192,sites_on_node*sizeof(su3_vector),GB_DIRECTORY)){
    node0_printf("Failed to create gb baryon sink mmap cache!\n");
    terminate(1);
  }

  populate_sink_oct_point_split(sko,qknum,r0,links);
}

void destroy_qk_oct_cache(int qknum){
  /** free the memory containing cache */
  if(!gb_qk_cache[qknum]){
   node0_printf("Warning: cache missing, cannot destroy\n");
   return;
  }
  free_mmap_cache(gb_qk_cache[qknum]);
  free(gb_qk_cache[qknum]);
}

void unmap_qk_oct_cache(int qknum){
  /** keep the memory but unmap the pointer, or memory has already been wiped */
  gb_qk_cache[qknum] = NULL;
}

mmap_cache* get_qk_cache_pointer(int qknum){
  return gb_qk_cache[qknum];
}

void assign_qk_cache_pointer(mmap_cache* qcache, int qknum){
  gb_qk_cache[qknum] = qcache;
}

/** degenerate with destroy_qk_oct_cache */
void destroy_sink_oct_cache(int qknum){
  destroy_qk_oct_cache(qknum);
}

/** copy the pointer rather than copying the entire cache */
void copy_qk_oct_cache(int qkcpy, int qkone){
  node0_printf("copying quark cache %d to cache %d\n",qkone,qkcpy);
  gb_qk_cache[qkcpy] = gb_qk_cache[qkone];
}

/** create a pointer array to store cache pointers*/
void create_gb_qk_cache(int numqk){
  /* malloc cache pointers so we can duplicate the pointers easily */
  gb_qk_cache = (mmap_cache**)malloc(numqk*sizeof(mmap_cache*));
  if(gb_qk_cache == NULL){
    printf("create_gb_qk_cache(%d): No room for cache\n",this_node);
    terminate(1);
  }
}

void destroy_gb_qk_cache(){
  free(gb_qk_cache);
}

/** create a cache to store source objects for recall later*/
void create_gb_src_cache(int numsrc){
  gb_src_cache = (mmap_cache*)malloc(sizeof(mmap_cache));
  if(gb_src_cache == NULL){
    printf("create_gb_src_cache(%d): No room for source cache\n",this_node);
    terminate(1);
  }
  if(alloc_mmap_cache(gb_src_cache,3*numsrc,sites_on_node*sizeof(su3_vector),GB_DIRECTORY)){
    node0_printf("Failed to create gb baryon source mmap cache!\n");
    terminate(1);
  }
}

void destroy_gb_src_cache(){
  free_mmap_cache(gb_src_cache);
  free(gb_src_cache);
}

#ifdef GB_INT_MAP
/*---------------------------------------------------------------------*/
/**
  Map indices for interacting quarks
  */
int
int_map_indices(int ci,int ki,int t,int cidx,int momi){
  /* organized by order of access */
  return (((offset_singlet_index(ci,ki)*8 +ci)*nt +t)*num_int_mom +momi)*3 +cidx;
  //return (((ci*8 +ki)*nt +t)*3 +cidx)*num_int_mom +momi;
}

/*---------------------------------------------------------------------*/
/**
  Check for existing interacting quark buffer for chosen spin-taste
  */
void
int_map_add_buffer(int stIdx){

  node0_printf("adding gb interaction buffer number %d for spin-taste %d\n",
    num_int_buf,stIdx);
  int_fill[num_int_buf] = stIdx;
  num_int_buf++;
}

/*---------------------------------------------------------------------*/
/**
  Check for existing interacting quark buffer for chosen spin-taste
  */
int
int_map_buffer_num(int stIdx){
  int i;
  for(i=0;i<num_int_buf;i++){
    if(stIdx == int_fill[i]) { return i;}
  }
  return -1;
}

/** special routine to fetch a source su3_vector, no need to toss when done */
void fetch_int_from_cache(su3_vector **dest,int stIdx){
  if(this_node == 0){
    *dest = (su3_vector *) gb_int_cache->buffer[int_map_buffer_num(stIdx)];
  }
}

void create_gb_int_cache(int numcur,int nummom){
  /* only need on head node! */
  /* be careful about putting this_node==0 statements around code blocks*/
  node0_printf("creating gb interaction cache\n");
  int i;
  if(this_node == 0){
    gb_int_cache = (mmap_cache*)malloc(sizeof(mmap_cache));
    int_fill = (int*)malloc(numcur*sizeof(int));
    if(gb_int_cache == NULL || int_fill == NULL){
      printf("create_gb_int_cache(%d): No room for interaction cache\n",this_node);
      terminate(1);
    }
    if(alloc_mmap_cache(gb_int_cache,numcur,nummom*nt*192*sizeof(su3_vector),GB_DIRECTORY)){
      node0_printf("Failed to create gb baryon interaction mmap cache!\n");
      terminate(1);
    }
    num_int_cur = numcur;
    num_int_mom = nummom;
    for(i=0;i<numcur;i++){
      int_fill[i] = 0;
    }
  } // this_node == 0
  else {
    /* need int_fill so other nodes know not to recreate buffers */
    int_fill = (int*)malloc(numcur*sizeof(int));
    gb_int_cache = NULL;
  }
  num_int_buf = 0;
  node0_printf("done creating interaction cache\n");
}

void destroy_gb_int_cache(){
  if(this_node == 0){
    free_mmap_cache(gb_int_cache);
    free(gb_int_cache);
    free(int_fill);
  }
  else {
    free(int_fill);
  }
  num_int_buf = 0;
}
#endif

#ifdef GB_SPEC_MAP
/**
  Map indices for spectator quarks
  */
int
spec_map_indices(int iqkn,int ci,int cj,int ki,int kj,int kk,int disp,int moms,int sc0){
  if(iqkn>2||iqkn<0){
    node0_printf("spec_map_indices given invalid quark index: %d\n",iqkn);
    terminate(1);
  }
  if(ci>7||cj>7||ci<0||cj<0){
    node0_printf("spec_map_indices given invalid source index set: %d, %d\n",ci,cj);
    terminate(1);
  }
  if(ki>7||kj>7||kk>7||ki<0||kj<0||kk<0){
    node0_printf("spec_map_indices given invalid sink index set: %d, %d, %d\n",ki,kj,kk);
    terminate(1);
  }
  if(disp>7||disp<0){
    node0_printf("spec_map_indices given invalid sink displacement: %d\n",disp);
    terminate(1);
  }

  return (iqkn*262144 +disp*32768 +ci*4096 +cj*512 +ki*64 +kj*8 +kk)*8*num_spec_mom +8*moms +sc0;
}

/**
  Check if a spectator map slot is filled
  */
short
get_spec_fill(int i){
  /* factor of 8 for neglected unit cube corner index -
     all cube corners assigned at same time, so redundant */
  return spec_fill[i/8];
}

/**
  Set a spectator map slot to some value
  */
void
set_spec_fill(int i,short val){
  /* factor of 8 for neglected unit cube corner index */
  spec_fill[i/8] = val;
}

/** special routine to fetch spectator su3_vector map, no need to toss when done */
void fetch_spec_from_cache(su3_vector **dest){
  if(this_node == 0){
    *dest = (su3_vector *) gb_spec_cache->buffer[0]; // only one buffer
  }
}

/*---------------------------------------------------------------------*/
void
create_gb_spec_cache(int nummom){
  node0_printf("creating gb spectator cache\n");
  if(this_node == 0){
    if(gb_spec_cache == NULL){
      gb_spec_cache = (mmap_cache*)malloc(sizeof(mmap_cache));
      /* one buffer, contains all data - otherwise, too many buffers */
      /* 3 src color x 3 interacting qk x 8^7 point splitting */
      if(gb_spec_cache == NULL){
        node0_printf("create_gb_spec_cache: no room for spectator map!\n");
        terminate(1);
      }
      if(alloc_mmap_cache(gb_spec_cache,1,18874368*nummom*sizeof(su3_vector),GB_DIRECTORY)){
        node0_printf("Failed to create gb baryon source mmap cache!\n");
        terminate(1);
      }
      num_spec_mom = nummom;
      /* 3 interacting qk x 8^6 point splitting (all epsilons and momenta done simultaneously) */
      spec_fill = (short *)calloc(786432*nummom,sizeof(short)); // calloc to zero memory
    } else {
      node0_printf("create_gb_spec_cache: spectator map already exists!\n");
      terminate(1);
    }
  } // this_node == 0
  else {
    /* still need spec_fill to keep track of missing slots on other nodes */
    if(spec_fill == NULL){
      /* 3 interacting qk x 8^6 point splitting */
      spec_fill = (short *)calloc(786432*nummom,sizeof(short)); // calloc to zero memory
      num_spec_mom = nummom;
      if(spec_fill == NULL){
        printf("node(%d): create_gb_spec_cache: no room for spectator map!\n",this_node);
        terminate(1);
      }
    } else {
      printf("node(%d): create_gb_spec_cache: spectator map already exists!\n",this_node);
      terminate(1);
    }
  }
}

void
destroy_gb_spec_cache(){
  if(this_node == 0){
    free_mmap_cache(gb_spec_cache);
    free(gb_spec_cache);
    free(spec_fill);
    gb_spec_cache = NULL;
    spec_fill = NULL;
  } // this_node == 0
  else {
    free(spec_fill);
    spec_fill = NULL;
  }
}
/*---------------------------------------------------------------------*/
#endif
#endif /* GB baryon */
