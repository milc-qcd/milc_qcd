/******************* milc_to_qphix_utilities.c ************************/
/* For the QPhiX interface */
/* MIMD version 7 */

/* 11/28/15 Created by Dhiraj Khalamkar */

#include "../include/generic_qphix.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

#define CG_DEBUG 1
#define UNOPTIMIZED_PACK_UNPACK 0

static int is_qphix_env_setup = 0;

/* The MILC layout routines list the coordinates and assume 4D */

static int milc_node_number(const int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int milc_node_index(const int coords[]){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

QPHIX_evenodd_t milc2qphix_parity(int milc_parity){
  switch(milc_parity){
  case(EVEN):       return QPHIX_EVEN;
  case(ODD ):       return QPHIX_ODD;
  case(EVENANDODD): return QPHIX_EVENODD;
  default:
    printf("milc2qphix_parity: Bad MILC parity %d\n", milc_parity);
    terminate(1);
  }
  return (QPHIX_evenodd_t)-999;
}

int qphix2milc_parity(QPHIX_evenodd_t qphix_parity){
  switch(qphix_parity){
  case(QPHIX_EVEN):     return EVEN;
  case(QPHIX_ODD ):     return ODD;
  case(QPHIX_EVENODD ): return EVENANDODD;
  default:
    printf("qphix2milc_parity: Bad QPHIX parity %d\n", qphix_parity);
    terminate(1);
  }
  return -999;
}

QPHIX_status_t 
initialize_qphix(int precision){
  /* \fixme - pass the qphix params via setup */
  /* create mbench */

  int minCt = 1;
  int threads_per_core = 1;
  int numThreads = omp_get_max_threads();
  if(getenv("MINCT")) minCt = atoi(getenv("MINCT"));
  if(getenv("THREADS_PER_CORE")) threads_per_core = atoi(getenv("THREADS_PER_CORE"));
  int numCores = numThreads / threads_per_core;

  int status = 0;

  static int latsize[4];
  static QPHIX_layout_t layout;
  
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;

  if(is_qphix_env_setup > 0){
    if(is_qphix_env_setup == precision)
      return status;
    else 
      /* Finalize so we can then initialize with the new precision */
      finalize_qphix();
  }

  layout.node_number = milc_node_number;
  layout.node_index = milc_node_index;
  layout.latdim = 4;
  layout.latsize = latsize;
  layout.machdim = 4;
  layout.machsize = (int *)get_logical_dimensions();
  layout.this_node = this_node;
  layout.even_sites_on_node = even_sites_on_node;
  layout.sites_on_node = sites_on_node;

  node0_printf("Initializing QPhiX for precision %d\n", precision);
  node0_printf("NumCores = %d, ThreadsPerCore = %d, minCt = %d\n", numCores, threads_per_core, minCt);

  status = QPHIX_init(&layout, precision);

  if(status){
    node0_printf("Error initializing QPhiX\n");
    fflush(stdout);
    terminate(1);
  }

  fflush(stdout);

  is_qphix_env_setup = precision;
  return status;
}

void
finalize_qphix(void){
  QPHIX_finalize();
}

#if 0
/*! \brief Constructor to create the single global qphix_env_t object. */ 
void
create_qphix_env ( int nx
                 , int ny
                 , int nz
                 , int nt
                 , int ncores
                 , int threads_per_core
                 , int min_ct
                 //, int soalen
                 //, int by
                 //, int bz
                 //, bool use_compressed12
              )
{
    void *gfl_e, *gfl_o, *gll_e, *gll_o;
    int neigh_ranks[8];
    int latt[4] = { nx, ny, nz, nt};
    //int *myCoord = get_logical_coordinates();
    int myCoord[4]; // = get_logical_coordinates();
    int *geom = get_logical_dimensions();
    
    int tmp = this_node;

    myCoord[0] = tmp % geom[0];
    tmp = tmp / geom[0];
    myCoord[1] = tmp % geom[1];
    tmp = tmp / geom[1];
    myCoord[2] = tmp % geom[2];
    tmp = tmp / geom[2];
    myCoord[3] = tmp % geom[3];

		if(this_node !=  myCoord[0] + geom[0]*( myCoord[1] + geom[1]*( myCoord[2] + geom[2]*( myCoord[3] )))) {
			printf("my Coord are wrong, this_node=%d\n", this_node);
			exit(1);
		}
    qphix_env_t *env = (qphix_env_t*)malloc(sizeof(qphix_env_t));
    if(env == NULL) {
        fprintf( stderr, "ERROR: Could not alloc mbench_t. Exiting." );
        exit(1);
    }
		for(int i = 0; i < 4; i++) {
			int neigh_coord[4] = {myCoord[0], myCoord[1], myCoord[2], myCoord[3]};
			neigh_coord[i] = (myCoord[i] == 0 ? geom[i] - 1 : myCoord[i] - 1);
			neigh_ranks[2*i] = neigh_coord[0] + geom[0]*( neigh_coord[1] + geom[1]*( neigh_coord[2] + geom[2]*( neigh_coord[3] )));;

			neigh_coord[i] = (myCoord[i] == geom[i] - 1 ? 0 : myCoord[i] + 1);
			neigh_ranks[2*i+1] = neigh_coord[0] + geom[0]*( neigh_coord[1] + geom[1]*( neigh_coord[2] + geom[2]*( neigh_coord[3] )));;
		}

    setup_mbench(latt, geom, myCoord,  neigh_ranks, this_node, ncores, threads_per_core, min_ct);
    /* Allocate the data structures for mbench */
    gfl_e = allocGauge18();
    gfl_o = allocGauge18();
    gll_e = allocGauge();
    gll_o = allocGauge();

    env->ks_src1   = allocKS();
    env->ks_dest1  = allocKS();
    env->ks_src2   = allocKS();
    env->ks_dest2  = allocKS();

    env->gll[0]   = gll_e;
    env->gll[1]   = gll_o;
    env->gfl[0]   = gfl_e;
    env->gfl[1]   = gfl_o;    
    
    /* Set the global values */
    g_qphix_env_obj = env;    
    is_qphix_env_setup = true;
}              

/*! \brief Destroy the global qphix_env object */
void
destroy_qphix_env (void)
{
    freeKS(g_qphix_env_obj->ks_src1);
    freeKS(g_qphix_env_obj->ks_dest1);
    freeKS(g_qphix_env_obj->ks_src2);
    freeKS(g_qphix_env_obj->ks_dest2);
    freeGauge(g_qphix_env_obj->gll[0]);
    freeGauge(g_qphix_env_obj->gll[1]);
    freeGauge18(g_qphix_env_obj->gfl[0]);
    freeGauge18(g_qphix_env_obj->gfl[1]);
    
    free(g_qphix_env_obj);
    is_qphix_env_setup = false;
}

/*!
 * Copy of load_longbacklinks (fermion_links_helpers.c), without the adjoint.
 * We do adjoint in the generated code for the dslash kernels for the back
 * links.
 */
su3_matrix *
load_longbacklinks_without_adjoint(fn_links_t *fn)
{
    su3_matrix *t_lbl = NULL;
    su3_matrix *t_ll = get_lnglinks(fn);
    register int i;
    register site *s;
    int dir;
    su3_matrix *tempmat1 = NULL;
    msg_tag *tag[4];
    char myname[] = "load_longbacklinks";

    /* Allocate space for t_lbl if NULL */
    t_lbl = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
    if(t_lbl==NULL){
        printf("%s(%d): no room for t_lbl\n",myname,this_node);
        terminate(1);
    }
        
    tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
        printf("%s: Can't malloc temporary\n",myname);
        terminate(1);
    }
    
    /* gather backwards longlinks */
    for( dir=XUP; dir<=TUP; dir ++){
    #pragma omp parallel for 
        for (i = 0; i < sites_on_node; i++) {
            tempmat1[i] = t_ll[dir+4*i];
        }
        tag[dir] = start_gather_field( tempmat1
                                     , sizeof(su3_matrix)
                                     , OPP_3_DIR(DIR3(dir))
                                     , EVENANDODD
                                     , gen_pt[dir] );
        wait_gather( tag[dir] );
        #pragma omp parallel for 
        for (i = 0; i < sites_on_node; i++) {
            su3_matrix * temp_ = (t_lbl + dir + 4*i);
            *temp_ = *((su3_matrix *)gen_pt[dir][i]);
        }
        cleanup_gather( tag[dir] );
    }
    
    free(tempmat1); 
    tempmat1 = NULL;
    
    return t_lbl;
}

/*!
 * Copy of load_fatbacklinks (fermion_links_helpers.c), without the adjoint.
 * We do adjoint in the generated code for the dslash kernels for the back
 * links.
 */
su3_matrix *
load_fatbacklinks_without_adjoint(fn_links_t *fn)
{
    su3_matrix *t_fbl = NULL;
    su3_matrix *t_fl = get_fatlinks(fn);
    register int i;
    register site *s;
    int dir;
    su3_matrix *tempmat1 = NULL;
    msg_tag *tag[4];
    char myname[] = "load_fatbacklinks";

    /* Allocate space for t_fbl if NULL */
    t_fbl = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
    if(t_fbl==NULL){
        printf("%s(%d): no room for t_fbl\n",myname,this_node);
        terminate(1);
    }
    
    tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
        printf("%s: Can't malloc temporary\n",myname);
        terminate(1);
    }
    
    /* gather backwards fatlinks */
    for( dir=XUP; dir<=TUP; dir ++){
        #pragma omp parallel for 
        for (i = 0; i < sites_on_node; i++) {
            tempmat1[i] = t_fl[dir+4*i];
        }
        tag[dir] = start_gather_field( tempmat1
                                       , sizeof(su3_matrix)
                                       , OPP_DIR(dir)
                                       , EVENANDODD
                                       , gen_pt[dir] 
                                       );
        wait_gather( tag[dir] );
        #pragma omp parallel for 
        for (i = 0; i < sites_on_node; i++) {
            su3_matrix* temp_ = (t_fbl + dir + 4*i);
            *temp_ = *((su3_matrix *)gen_pt[dir][i]);
        }
        cleanup_gather( tag[dir] );
    }
    
    free(tempmat1); 
    tempmat1 = NULL;
    
    return t_fbl;
}

/*! 
 * Convenience function to create and array of su3_vectors from the lattice 
 * site array. The argument defines if we want to copy over the source or the
 * destination.
 */

/* SEE su3_vector *create_v_field_from_site_member(field_offset sv) */

void
get_links_from_lattice (void* fgauge[2], void *lgauge[2], fn_links_t *fn)
{
    int i, dir, p;
    site *s;
    su3_matrix* t_fatlink = get_fatlinks(fn);
    su3_matrix* t_lnglink = get_lnglinks(fn);
    double t0 = __rdtsc();

    msg_tag *ftag[4], *ltag[4], *tag;

    su3_matrix *ftempmat = NULL, *ltempmat = NULL;
    ind_t *arr = (ind_t *)malloc(8*sites_on_node*sizeof(ind_t));
    ftempmat = (su3_matrix *)malloc(8*sites_on_node*sizeof(su3_matrix));
    ltempmat = (su3_matrix *)malloc(8*sites_on_node*sizeof(su3_matrix));
   
    /* gather backwards fatlinks */
    for( dir=XUP; dir<=TUP; dir++){
      ftag[dir] = start_gather_field( &t_fatlink[dir]
                                       , 4*sizeof(su3_matrix)
                                       , OPP_DIR(dir)
                                       , EVENANDODD
                                       , gen_pt[dir] 
                                       );
      ltag[dir] = start_gather_field( &t_lnglink[dir]
                                       , 4*sizeof(su3_matrix)
                                       , OPP_3_DIR(DIR3(dir))
                                       , EVENANDODD
                                       , gen_pt[4+dir] 
                                       );
    }
    for( dir=XUP; dir<=TUP; dir++){
        wait_gather( ftag[dir] );
        wait_gather( ltag[dir] );
    }
 
    double t1 = __rdtsc();
    node0_printf("BackLnkGather = %10.4g\n", t1-t0);

    //for(p = 0; p < 2; p++) 
    {
    //    int parity = (p == 0 ? EVEN : ODD);
        #pragma omp parallel for private(s) default(shared)
    //    FORSOMEPARITY(i,s,parity) {
        for(i=0;i<sites_on_node;i++) {
            site *s = &(lattice[i]);
            su3_matrix *fat4, *lng4;
            fat4     = &(t_fatlink[4*i]);
            lng4     = &(t_lnglink[4*i]);
            //int ii = (((s->t *nz+ s->z)*ny+s->y)*nx/2+s->x/2) + (parity == ODD ? even_sites_on_node : 0);
	    //if(ii != i) printf("YYY ii=%d not same as  i=%d\n", ii, i); 
            /* plug the 4 su3_matrices in thee forward dirs into the qphix array */
            for( int dir=XUP; dir <= TUP; dir++) {
                ftempmat[8*i+MILC2MBENCH[OPP_DIR(dir)]] = *((su3_matrix*)gen_pt[dir][i]);
                ftempmat[8*i+MILC2MBENCH[dir]] = fat4[dir];
                ltempmat[8*i+MILC2MBENCH[OPP_DIR(dir)]] = *((su3_matrix*)gen_pt[4+dir][i]);
                ltempmat[8*i+MILC2MBENCH[dir]] = lng4[dir];
            }
        }    
    }
    for( dir=XUP; dir<=TUP; dir++){
        cleanup_gather( ftag[dir] );
        cleanup_gather( ltag[dir] );
    }
    double t2 = __rdtsc();
    node0_printf("BackLnkCopy = %10.4g\n", t2-t1);

#if 1
    setFullGauge18(fgauge[0], &ftempmat[0], even_sites_on_node*8);
    setFullGauge18(fgauge[1], &ftempmat[even_sites_on_node*8], 8*even_sites_on_node);
    setFullGauge  (lgauge[0], &ltempmat[0], even_sites_on_node*8);
    setFullGauge  (lgauge[1], &ltempmat[even_sites_on_node*8], even_sites_on_node*8);
#endif
   double t3 = __rdtsc();
    node0_printf("BackLnkConv = %10.4g\n", t3-t2);
    free(arr);
    free(ftempmat);
    free(ltempmat);
}

void
get_fatlinks_from_lattice (void* gauge18s, int parity, fn_links_t *fn)
{
    int i;
    site *s;
    su3_matrix* t_fatlink = get_fatlinks(fn);
    su3_matrix* t_fatback = load_fatbacklinks_without_adjoint(fn);
    su3_matrix *fat4, *backfat4;

#if _OPENMP
	#pragma omp parallel for private(s) default(shared)
#endif
    FORSOMEPARITY(i,s,parity) {
        fat4     = &(t_fatlink[4*i]);
        backfat4 = &(t_fatback[4*i]);
        /* plug the 4 su3_matrices in thee forward dirs into the qphix array */
        for( int dir=XUP; dir <= TUP; dir++) {
            setGauge18( gauge18s
                      , &fat4[dir]
                      , MILC2MBENCH[dir]
                      , s->x/2
                      , s->y
                      , s->z
                      , s->t
                      );
            setGauge18( gauge18s
                      , &backfat4[dir]
                      , MILC2MBENCH[OPP_DIR(dir)]
                      , s->x/2
                      , s->y
                      , s->z
                      , s->t
                      );
        }
    }    
    free(t_fatback);
}

void
get_longlinks_from_lattice (void* gauges, int parity, fn_links_t *fn)
{
    int i;
    site *s;
    su3_matrix* t_longlink = get_lnglinks(fn);
    su3_matrix* t_longback = load_longbacklinks_without_adjoint(fn);
    su3_matrix *lng4, *backlng4;

#if _OPENMP
	#pragma omp parallel for private(s) default(shared)
#endif
    FORSOMEPARITY(i,s,parity) {
        lng4     = &(t_longlink[4*i]);
        backlng4 = &(t_longback[4*i]);

        for( int dir=XUP; dir <= TUP; dir++) {
            setGauge( gauges
                    , &lng4[dir]
                    , MILC2MBENCH[dir]
                    , s->x/2
                    , s->y
                    , s->z
                    , s->t
                    );
            
            setGauge( gauges
                    , &backlng4[dir]
                    , MILC2MBENCH[OPP_DIR(dir)]
                    , s->x/2
                    , s->y
                    , s->z
                    , s->t
                    );
        }
    }   
    free(t_longback);
}


/*! 
 * Get the su3_vectors from the site major lattice and convert them into an
 * array of ks_spinors for qphix
 */
void
get_ks_spinors_from_lattice (su3_vector *ks, void* ks_spinor, int parity)
{
    int i;
    site *s;

#if _OPENMP
    //#pragma omp parallel for private(s) default(shared)
#warning using omp    
	int loopend= (parity)==EVEN ? even_sites_on_node : sites_on_node ; 
	#pragma omp parallel for private(s) default(shared)
	for( int i=((parity)==ODD ? even_sites_on_node : 0 ); i<loopend; i++){
		s = &(lattice[i]);
#else
    FORSOMEPARITY(i,s,parity) {
#endif
        setKS( ks_spinor
               , ks + i
               , s->x/2, s->y, s->z, s->t 
             ); 
    }
}

void
set_ks_spinors_into_lattice (su3_vector *ks, void* ks_spinor, int parity)
{
    int i;
    site *s;

#if 1
#if _OPENMP
    #pragma omp parallel for private(s) default(shared)
#endif
    FORSOMEPARITY(i,s,parity) {
      //        getKS(ks_spinor, (su3_vector *)F_PT(s, ks_off)
        getKS(ks_spinor, ks + i
              , s->x/2, s->y, s->z, s->t 
              ); 
    }
#else
    // Get the spinors into a separate array. We pass the whole array to mbench
    //    su3_vector *spinors = get_su3_vectors_from_lattice(ks_off);
    su3_vector *spinors = ks;
#endif
}

#endif

