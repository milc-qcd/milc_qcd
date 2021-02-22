/****** new_hvy_pot.c  -- ******************/

/***********************************************************************
* Heavy quark potential
* MIMD version 7
*
************************************************************************
*
************************************************************************
* Old algorithm (newer versions):
************************************************************************
* 1) uses general gather to do xc[XUP]-xc[YUP] shifts all at once
* 2) assigns additional su3_mat valued field buffer that  
*    contains xc[XUP]-xc[YUP shifted Wilson lines; 
*    xc[XUP]-shifts performed as general gathers (except mlat[XUP]==1),
*    xc[YUP]-shifts performed sequentially as single step gathers, 
*    after xc[ZUP]-shifts: xc[YUP]-shifted field recovered from buffer
*
* Typical speedup is nearly 40% less execution time compared to 
* the original code in generic/hvy_pot.c. JHW 2/24/15
*
*
*
************************************************************************
* New algorithm (revised versions):
************************************************************************
* Gather entire local lattice blocks of the fields and 
* interpret their entries by adding to the corresponding 
* site indices.
* 
* PLEASE BE AWARE THAT THERE MIGHT BE A CONFLICT!
* Introduced modified multi-distance gathers in the scheme
* of the 3n gathers for the naik term (-> setup.c). This 
* is done by replacing the 3n gathers for the naik term with
* fixed nzyxt for transporting the full local lattice.
*
* Evaluate in different spatial directions, to check rotational
* invariance.  Gauge fix to Coulomb gauge to do spatial segments.
*
* Modified these things in "setup.c":
* 1) Get the size of the local lattices
* 2) Make gathers for displacements in local lattice
*
* These things are done in this file (new_hvy_pot.c):
* 3) Initialise Wilson lines in buffer 
*
* OUTER LOOP OF LOCAL LATTICES
* 4) Decide and request which sublattice is gathered next 
*    via function: static next_gather find_gather( ... ) 
*   INNER LOOP OF LOCAL LATTICES
*   5) Do proper all-to-all contractions on local lattice 
*      The inner loop is a double loop. The inner inner loop is OMP parallelized.
* GO ON WITH OUTER LOOP 
* 6) Wait for whatever was gathered before 
* 7) Global sums once all local lattices are done 
* 8) Copy original Wilson line into shifted Wilson line 
* 9) Check loop end condition and control next shift 
*
* ALL field loops are OMP parallelized, except 
* THE OUTER LOOP of the innermost double loop
*
*
************************************************************************/

#include "../generic/generic_includes.h"	/* definitions files and prototypes */
#include "./hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

//#define DEBUG
#ifdef NEW_HVY_POT 

#ifdef DEBUG
// show only first row of su3_matrix for debugging
static void show_mat_row( su3_matrix *m ){

  register int i;
  for(i=0;i<3;i++)
    printf("(%.0f,%.0f)\t", m->e[0][i].real,m->e[0][i].imag);
  printf("adress %ld\n",(long)m);
}

// set special field values for debugging
static void set_special_row( su3_matrix *m ){

  register int i;
  memset(m,'\0',sites_on_node*sizeof(su3_matrix));
  site *s;
  FORALLSITES_OMP(i,s, default(shared) ){
    m[i].e[0][0].real = (site_coord(s,xc[XUP]));
    m[i].e[0][0].imag = (site_coord(s,xc[YUP]));
    m[i].e[0][1].real = (site_coord(s,xc[ZUP]));
    m[i].e[0][1].imag = (site_coord(s,xc[TUP]));
    m[i].e[0][2].real = this_node;
    m[i].e[0][2].imag = i;
  } END_LOOP_OMP;
}
#endif

/************************************************************************/
/* static variables, functions and macros needed in all cases */
/************************************************************************/

#define HQPTIME
#ifdef HQPTIME
enum { HQPTIME_COMM=0, HQPTIME_CONTRACT=1, HQPTIME_OUTPUT=2, HQPTIME_TOTAL=3 };
static double hqptime[HQPTIME_TOTAL]={0.,0.,0.};
#endif
static void hqp_print_timings(char myalg[]);

#ifdef PLOOPCOR_MEAS 
/* Indicates the status of the Polyakov loop calculation */
#define P_NOTDONE 0
#define P_INDOING 1
#define P_NOTWRITTEN 2
#define P_WRITTEN 3
#endif

#ifdef PLOOPCOR_MEAS 
#define WANTCONTRACTIONS( mlt ) \
  ( mlt <= maxc[TUP] || mlt == nc[TUP] ) 
#else 
#define WANTCONTRACTIONS( mlt ) \
  ( mlt <= maxc[TUP] ) 
#endif 

static int buflvl=0;
static int llat[4]={0,0,0,0};
static hqp_geom *geom=NULL;

/* static variables that are used in determining the gathers,
 * useful only for choosing a better communication paradigm */
static int nstale, ngather; 

/* buffers for correlators */
#ifdef DEBUG
static double *counter=NULL;
#endif
#ifdef COULOMB
static double *wlc=NULL;
#endif
#ifdef DIQUARK
static double *dqc=NULL;
#endif
#ifdef PLOOPCOR_MEAS
static int ploop_done=P_NOTDONE;
static double *plc=NULL,*plcr=NULL,*plci=NULL;
#endif

/* buffers for fields */
static su3_matrix *owline=NULL,*swline=NULL,*buffer=NULL;

/* static functions used in all algorithms */
#ifdef PLOOPCOR_MEAS
static inline double_complex add_to_ploop( int i, double_complex *ploop );
#endif
static inline void clear_buffers( int mysize );
static inline void contractions( int gi, int i, int j, double_complex *ctr1);
static inline void disp_from_index( int gi, int *disp );
static inline void disp_from_sites( int mlat[], site *s, site *t, int *disp );
static inline void global_sums( int mysize );

static inline int hqp_buffer_dim (void );
static void hqp_free_buffers( void );
static void hqp_setup_buffers( su3_matrix *links, int mysize );

static inline int index_from_disp ( int disp[] );
static inline void output_all( char smtag[], int mi, int mlt);

#define NOTDONE 0
#define DONE 1

/* Indicates the overall type of the next gather operation */
#define GENGATHER 1
#define SIMGATHER 0
#define NOGATHER -1

/* Indicates if some points are inside of the 3-sphere with radius max_r */
#define OUT 0
#define IN  1

//#define MAXCLINE 256
#define MAXCNAME 16
#define MAXSMTAG 4

/************************************************************************/
/* static variables, functions and macros needed only for algorithm 1 */
/************************************************************************/

#define BUFMAT2

#ifdef BUFMAT0
#define GENSHIFT
#endif
#ifdef BUFMAT1
#define GENSHIFT
#define BUFMAT
#endif
#ifdef BUFMAT2
static su3_matrix *buffer2=NULL;;
#define GENSHIFT
#define BUFMAT
#endif
#ifdef BUFMAT
static su3_matrix *buffer1=NULL;
#endif

static inline void contractions_old( int gi, int i, su3_matrix *sfield, double_complex *ctr1);

static void shiftmat( su3_matrix *src, su3_matrix *dest, int dir );
#ifdef GENSHIFT
static void general_shiftmat( su3_matrix *src, su3_matrix *dest, int *disp );
#endif

/************************************************************************/
/* static variables, functions and macros needed only for algorithm 2 */
/************************************************************************/

#define HASHLOOP
#ifdef HASHLOOP
#define HASH(ct,cs,ncs) hash[ct*ncs+cs]
static inline int hash_ct_from_site( site *s, int ncs );
static int *setup_coord_hash( int *ncs );
static void free_coord_hash( int *hash );
#endif

/* FLUTE and SNAKE are two different schemes for traversing the lattice 
 * that dramatically improve performance on large lattices for massively 
 * parallelized calculations, i.e. if the local volume per node can be 
 * made very small; one has to be used; for more details see 
 * the detailed description at the end of this file */ 
#define FLUTE
#ifdef SNAKE
static int raising_y=1;
static int raising_z=1;
#endif
#if ( !( (defined FLUTE) ^ (defined SNAKE) ) )
BOMB THE COMPILE
#endif
 
/* new data type that contains all required information about the next gather */
typedef struct {
  int gengather;
  int raise;
  int insphere;
  int lastraise;
  int raisedir[3];
} next_gather;

static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) ;

/* define macros that are used to judge necessity 
 * for further spatial shifts of local lattices */
#define SPATIAL(dir) ((dir >= XUP && dir <TUP) || (dir > TDOWN && dir <=XDOWN))

#ifndef ANISOTROPY
#define INSPHERE(base,llat,max_r2) \
  ( ( +base[xc[XUP]]*base[xc[XUP]]+base[xc[YUP]]*base[xc[YUP]]+base[xc[ZUP]]*base[xc[ZUP]] <= max_r2 \
   || (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1)+base[xc[YUP]]*base[xc[YUP]]+base[xc[ZUP]]*base[xc[ZUP]] <= max_r2 \
   || +base[xc[XUP]]*base[xc[XUP]]+(base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1)+base[xc[ZUP]]*base[xc[ZUP]] <= max_r2 \
   || +base[xc[XUP]]*base[xc[XUP]]+base[xc[YUP]]*base[xc[YUP]]+(base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) <= max_r2 \
   || (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1)+base[xc[YUP]]*base[xc[YUP]]+(base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) <= max_r2 \
   || (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1)+(base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1)+base[xc[ZUP]]*base[xc[ZUP]] <= max_r2 \
   || +base[xc[XUP]]*base[xc[XUP]]+(base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1)+(base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) <= max_r2 \
   || (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1)+(base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) \
     +(base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) <= max_r2 )? IN : OUT )
#else 
#ifdef ONEDIM_ANISO_TEST
#define INSPHERE(base,llat,max_r2) \
  ( ( (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * base[xc[XUP]]*base[xc[XUP]] \
      +(iso_xiq*iso_xiq) * ( base[xc[YUP]]*base[xc[YUP]] + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +(iso_xiq*iso_xiq) * ( base[xc[YUP]]*base[xc[YUP]] + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * base[xc[XUP]]*base[xc[XUP]] \
      +(iso_xiq*iso_xiq) * ( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * base[xc[XUP]]*base[xc[XUP]] \
      +(iso_xiq*iso_xiq) * ( base[xc[YUP]]*base[xc[YUP]] + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +(iso_xiq*iso_xiq) * ( base[xc[YUP]]*base[xc[YUP]] + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +(iso_xiq*iso_xiq) * ( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * base[xc[XUP]]*base[xc[XUP]] \
      +(iso_xiq*iso_xiq) * ( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +(iso_xiq*iso_xiq) * ( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 )? IN : OUT )
#else
#define INSPHERE(base,llat,max_r2) \
  ( ( (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * base[xc[XUP]]*base[xc[XUP]] \
      +( base[xc[YUP]]*base[xc[YUP]] + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +( base[xc[YUP]]*base[xc[YUP]] + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * base[xc[XUP]]*base[xc[XUP]] \
      +( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * base[xc[XUP]]*base[xc[XUP]] \
      +( base[xc[YUP]]*base[xc[YUP]] + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +( base[xc[YUP]]*base[xc[YUP]] + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + base[xc[ZUP]]*base[xc[ZUP]] ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * base[xc[XUP]]*base[xc[XUP]] \
      +( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 \
   || (cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * (base[xc[XUP]]-llat[xc[XUP]]+1)*(base[xc[XUP]]-llat[xc[XUP]]+1) \
      +( (base[xc[YUP]]-llat[xc[YUP]]+1)*(base[xc[YUP]]-llat[xc[YUP]]+1) + (base[xc[ZUP]]-llat[xc[ZUP]]+1)*(base[xc[ZUP]]-llat[xc[ZUP]]+1) ) <= max_r2 )? IN : OUT )
#endif
#endif
/************************************************************************/
/* externally accessible function for cycling spatial components of maxc */
/************************************************************************/

void hqp_cycle_spatial( su3_matrix *links, int hqp_alg ) {
  char myname[] = "hqp_cycle_spatial";

  int mmxc[4];
  memcpy(mmxc,maxc,4*sizeof(int));
  int any_dirs_are_equal=(maxc[XUP]==maxc[YUP])
                        +(maxc[YUP]==maxc[ZUP])
                        +(maxc[ZUP]==maxc[XUP]);
  if ( any_dirs_are_equal>1 ) {
    hqp_switch(links, hqp_alg );
  } else {
    for (int pxc=1;pxc<=6;pxc++) {
      if (pxc>3) {
        if ( any_dirs_are_equal==0 ) {  
          maxc[XUP] = mmxc[(YUP+pxc)%3];
          maxc[YUP] = mmxc[(XUP+pxc)%3];
        } else 
          continue;
      } else {
        maxc[XUP] = mmxc[(XUP+pxc)%3];
        maxc[YUP] = mmxc[(YUP+pxc)%3];
      }    
      maxc[ZUP] = mmxc[(ZUP+pxc)%3];
      hqp_switch(links, hqp_alg );
    }    
  }    
  memcpy(maxc,mmxc,4*sizeof(int));
}

/************************************************************************/
/* externally accessible function for switching between both algorithms */
/************************************************************************/

void hqp_switch( su3_matrix *links, int hqp_alg ) {
  local_lattice_size(llat);
  int plls=llat[xc[XUP]]*llat[xc[YUP]]*llat[xc[ZUP]];
  char myname[] = "hqp_switch";

  switch(hqp_alg) {
    case(HQPALG_ONE): 
    {
      hvy_pot_alg_old( links );
      break;
    }
    case(HQPALG_TWO): 
    {
      hvy_pot_alg_new( links );
      break;
    }
    case(HQPALG_THREE): 
    {
      llat[xc[ZUP]]=1;
      hvy_pot_alg_new( links );
      break;
    }
    case(HQPALG_FOUR): 
    {
      llat[xc[ZUP]]=1;
      llat[xc[YUP]]=1;
      hvy_pot_alg_new( links );
      break;
    }
    case(HQPALG_FIVE): 
    {
      llat[xc[ZUP]]=1;
      llat[xc[YUP]]=1;
      llat[xc[ZUP]]=1;
      hvy_pot_alg_new( links );
      break;
    }
    case(HQPALG_TIME): 
    {
      char myalg[16];

      sprintf(myalg,"%s","alg1");
      node0_printf("%s -- optimized, traditional MILC algorithm\n",myalg);
      double alg1time = -dclock();
      hvy_pot_alg_old( links );
      alg1time += dclock();
      hqp_print_timings(myalg);

      sprintf(myalg,"%s","alg2");
      node0_printf("%s -- multistep communication in 3 directions\n",myalg);
      double alg2time = -dclock();
      hvy_pot_alg_new( links );
      alg2time += dclock();
      hqp_print_timings(myalg);

      sprintf(myalg,"%s","alg3");
      node0_printf("%s -- multistep communication in 2 directions\n",myalg);
      llat[xc[ZUP]]=1;
      double alg3time = -dclock();
      hvy_pot_alg_new( links );
      alg3time += dclock();
      hqp_print_timings(myalg);

      sprintf(myalg,"%s","alg4");
      node0_printf("%s -- multistep communication in 1 directions\n",myalg);
      llat[xc[YUP]]=1;
      double alg4time = -dclock();
      hvy_pot_alg_new( links );
      alg4time += dclock();
      hqp_print_timings(myalg);

      sprintf(myalg,"%s","alg5");
      node0_printf("%s -- multistep communication in 0 directions\n",myalg);
      llat[xc[XUP]]=1;
      double alg5time = -dclock();
      hvy_pot_alg_new( links );
      alg5time += dclock();
      hqp_print_timings(myalg);

      local_lattice_size(llat);
      plls=llat[xc[XUP]]*llat[xc[YUP]]*llat[xc[ZUP]];
      node0_printf("%s: local lattice [%d,%d,%d,%d] -- local spatial volume %d (<%d): total times\n\
alg1 hqp timed total %.4e\n\
alg2 hqp timed total %.4e\n\
alg3 hqp timed total %.4e\n\
alg4 hqp timed total %.4e\n\
alg5 hqp timed total %.4e\n\n",
        myname,llat[xc[XUP]],llat[xc[YUP]],llat[xc[ZUP]],llat[xc[TUP]], plls, LLS_THRESH, 
        alg1time, alg2time, alg3time, alg4time, alg5time);

      break;
    }
    case(HQPALG_AUTO): 
    default:
    {
      /* LLS_THRESH needs to be defined as the threshold for switching between both algorithms */
      if ( plls > LLS_THRESH ) { 
        plls /= llat[xc[ZUP]];
        llat[xc[ZUP]]=1;
        if ( plls > LLS_THRESH ) { 
          plls /= llat[xc[YUP]];
          llat[xc[YUP]]=1;
          if ( plls > LLS_THRESH ) { 
            plls /= llat[xc[XUP]];
            llat[xc[XUP]]=1;
          }
        }
      }
      node0_printf("%s: auto-select communication scheme in steps of [%d,%d,%d] \
with double loop over local volume %d\n",myname,
      llat[xc[XUP]],llat[xc[YUP]],llat[xc[ZUP]], plls);
      hvy_pot_alg_new( links );
      break;
    }   
  }   
}

/************************************************************************/
/* externally accessible function of old algorithm */
/************************************************************************/

void hvy_pot_alg_old( su3_matrix *links ) { 

  register int i; 
  register int mu;
  int genv[4]={0,0,0,0};
  int disp[4]={0,0,0,0};
  int mlat[4]={0,0,0,0};
  su3_matrix *oldmat=NULL, *newmat=NULL, *tt=NULL;
  msg_tag *mtag0;
  char myname[] = "hvy_pot_alg_old";

  char smtag[MAXSMTAG]="";
#ifdef SMEARING
  sprintf(smtag,"_%d",tot_smear);
#endif

  double_complex ctr1, ploop={0,0};
  // Trigger allocation of additional buffers
#ifdef BUFMAT
  buflvl=1;
#  ifdef BUFMAT2
  buflvl=2;
#  endif
#endif

  geom = hqp_geometry( myname );
  int mi = hqp_buffer_dim();
  hqp_setup_buffers(links,mi);
#ifdef BUFMAT
  su3_matrix *bufmat1=buffer1; 
#  ifdef BUFMAT2
  su3_matrix *bufmat2=buffer2; 
#  endif
#endif

#ifdef HQPTIME
  memset(hqptime,'\0',HQPTIME_TOTAL*sizeof(double));
#endif

  ngather=0;
  nstale=0;
  /* Use owline to construct t-direction path from each point */
#ifdef PLOOPCOR_MEAS
  for ( mlat[TUP]=1; (mlat[TUP]<=nc[TUP]); mlat[TUP]++ )
#else
  for ( mlat[TUP]=1; (mlat[TUP]<=maxc[TUP]); mlat[TUP]++ )
#endif
  {
#ifdef HQPTIME
    hqptime[HQPTIME_COMM] -= dclock();
#endif
    if ( mlat[TUP]>1 ) {
      mtag0 = start_gather_field( owline, sizeof(su3_matrix),
              xc[TUP], EVENANDODD, gen_pt[0] );
      wait_gather(mtag0);
      FORALLFIELDSITES(i) 
        mult_su3_nn( links+4*i+xc[TUP], (su3_matrix *)gen_pt[0][i], buffer+i ); 
      cleanup_gather( mtag0 ); 
      ngather++;
      memcpy(owline,buffer,sites_on_node*sizeof(su3_matrix));
    }
    /* now owline is path of length mlat[TUP] in xc[TUP] direction */
    oldmat = swline; //swline;
    newmat = buffer; //staple;        /* will switch these two */
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] += dclock();
#endif

#ifdef BUFMAT
  bufmat1=buffer1; 
#  ifdef BUFMAT2
  bufmat2=buffer2; 
#  endif
#endif

    for( mlat[XUP]=0; mlat[XUP]<=maxc[XUP]; mlat[XUP]++ ){
      for( mlat[YUP]=0; mlat[YUP]<=maxc[YUP]; mlat[YUP]+= 1 ){

        /* now gather from spatial dirs, compute products of paths */
#ifdef HQPTIME
        hqptime[HQPTIME_COMM] -= dclock();
#endif
#ifndef GENSHIFT
        memcpy(swline,owline,sites_on_node*sizeof(su3_matrix));
        for(i=0;i<mlat[XUP],maxc[XUP]>0;i++){
          /* shift oldmat by one further step in mlat[XUP],
           * then swap pointers */
          shiftmat( oldmat, newmat, xc[XUP] ); 
          ngather++;
          tt=oldmat; oldmat=newmat; newmat=tt;
        }
        for(i=0;i<mlat[YUP],maxc[YUP]>0;i++){
          /* shift oldmat by one further step in mlat[YUP],
           * then swap pointers */
          shiftmat( oldmat, newmat, xc[YUP] ); 
          ngather++;
          tt=oldmat; oldmat=newmat; newmat=tt;
        }
#else
#  ifdef BUFMAT0
        // uses general gather 
        memcpy(oldmat,owline,sites_on_node*sizeof(su3_matrix));
        genv[XUP]=mlat[XUP]; 
        genv[YUP]=mlat[YUP];
        general_shiftmat( oldmat, newmat, genv ); 
        ngather++;
        tt=oldmat; oldmat=newmat; newmat=tt;
#  endif //BUFMAT0

        /* use minimal shifts at the expense of extra memory */
#  ifdef BUFMAT
        // shift matrix by mlat[XUP] steps only 
        // if mlat[YUP] loop starts from zero, otherwise keep buffered
        if ( mlat[YUP]==0 ) {
          /* first put xc[XUP]-shifted wline into bufmat
           * then, after xc[YUP] shift continue shifting 
           * oldmat by mlat[XUP] steps in neg. xc[XUP] direction */
          genv[XUP]=mlat[XUP];
#    ifdef BUFMAT2
          if ( mlat[XUP] >= 1) {
            /* shift bufmat2 by one further step in mlat[XUP],
             * then swap pointers */
            shiftmat( bufmat2, bufmat1, xc[XUP] ); 
            ngather++;
            tt=bufmat2; bufmat2=bufmat1; bufmat1=tt;
          } else 
            // put wline into bufmat2 
            memcpy(bufmat2,owline,sites_on_node*sizeof(su3_matrix));
          // copy shifted wline from bufmat2 to bufmat1
          memcpy(bufmat1,bufmat2,sites_on_node*sizeof(su3_matrix));
#    else
          // put wline into oldmat
          memcpy(oldmat,owline,sites_on_node*sizeof(su3_matrix));
          /* shift oldmat by one or more steps in mlat[XUP],
           * then swap pointers */
          if ( mlat[XUP] >= 1 && mlat[XUP] <= maxc[XUP] ) {
            if ( mlat[XUP] > 1) {
              // shift oldmat by mlat[XUP] steps 
              general_shiftmat( oldmat, newmat, genv ); 
            } else 
              // shift oldmat by just one step in mlat[XUP]
              shiftmat( oldmat, newmat, xc[XUP] ); 
            ngather++;
            tt=oldmat; oldmat=newmat; newmat=tt;
          }
          // copy shifted wline from oldmat to bufmat1
          memcpy(bufmat1,oldmat,sites_on_node*sizeof(su3_matrix));
#    endif
        } else 
          if ( mlat[YUP] > 0 && mlat[YUP] <= maxc[YUP] ) {
          // shift x-shifted buffer one further step in mlat[YUP], 
          // keep xy-shifted matrix buffered
            shiftmat( bufmat1, newmat, xc[YUP] ); 
            ngather++;
            tt=bufmat1; bufmat1=newmat; newmat=tt;
          }
#  endif //BUFMAT
#endif //ifndef GENSHIFT
#ifdef HQPTIME
        hqptime[HQPTIME_COMM] += dclock();
#endif

        for( mlat[ZUP]=0; mlat[ZUP]<=maxc[ZUP]; mlat[ZUP]++ ){

#ifdef HQPTIME
          hqptime[HQPTIME_CONTRACT] -= dclock();
#endif

          int stale=1; /* stale == 1 -> all displacements out of the range ? */
          /* 5a)  Decide if we want to contractions */
          if (WANTCONTRACTIONS(mlat[TUP])) 
          {   

            /* 5b) Do proper all-to-all contractions on local lattice */
            site *s=NULL;
            FORALLSITES_OMP(i,s, shared(stale) ) {

              memcpy(disp,mlat,4*sizeof(int));
#ifdef PLOOPCOR_MEAS
              ctr1 = add_to_ploop( i, ( mlat[XUP]==0 && mlat[YUP]==0 && mlat[ZUP]==0
                                     && mlat[TUP] == nc[TUP] ?&ploop:NULL) );
#endif

              su3_matrix *sfield = oldmat;
#ifdef BUFMAT
              if (mlat[ZUP]==0)  
                sfield = bufmat1;
#endif //ifndef BUFMAT

              /* Compute displacement disp and squared radius. If okay, add to correlator buffers? */
              if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) {
                contractions_old ( index_from_disp(disp), i, sfield, &ctr1 );
                stale=0;
              } 

            } END_LOOP_OMP;
          } // END if (WANTCONTRACTIONS(mlat[TUP])) 

#ifdef HQPTIME
          hqptime[HQPTIME_CONTRACT] += dclock();
#endif

#ifdef PLOOPCOR_MEAS
          if ( ploop_done == P_INDOING ) 
            ploop_done = P_NOTWRITTEN; 
#endif

          if ( mlat[ZUP] < maxc[ZUP] ) {
          /* before we increment mlat[ZUP], shift in xc[ZUP] direction */
#ifdef HQPTIME
            hqptime[HQPTIME_COMM] -= dclock();
#endif
#ifndef BUFMAT
            shiftmat( oldmat, newmat, xc[ZUP] ); 
#else
            if (mlat[ZUP]>0) {
              shiftmat( oldmat, newmat, xc[ZUP] ); 
            } else {
              shiftmat( bufmat1, newmat, xc[ZUP] ); 
            }
#endif //ifndef BUFMAT
            ngather++;
            nstale +=stale;
            tt=oldmat; oldmat=newmat; newmat=tt;
#ifdef HQPTIME
            hqptime[HQPTIME_COMM] += dclock();
#endif
          } /* if ( mlat[ZUP] < maxc[ZUP] ) */
        } /* mlat[ZUP] distance */
      } /* mlat[YUP] distance */
    } /* mlat[XUP] distance */

#ifdef HQPTIME
    hqptime[HQPTIME_OUTPUT] -= dclock();
#endif
  /* 7a)  Decide if we actually did the contractions and want to output */
    if (WANTCONTRACTIONS(mlat[TUP])) 
    {
      /* 7b) Global sums once all local lattices are done */
#ifdef PLOOPCOR_MEAS
      if ( ploop_done == P_NOTWRITTEN && mlat[TUP] == nc[TUP] ) { 
        g_dcomplexsum(&ploop); 
        node0_printf("P_LOOP%s:\t%.6e\t%.6e\n", smtag, ploop.real/volume, ploop.imag/volume );
        ploop_done = P_WRITTEN;
      }
#endif
      global_sums( mi );

      /* output only data we are interested in (and have actually computed) */
      output_all( smtag, mi , mlat[TUP] );
    } // END OF DECISION WHETHER WE DID CONTRACTIONS
    fflush(stdout);
#ifdef HQPTIME
    hqptime[HQPTIME_OUTPUT] += dclock();
#endif
    /* Reset correlator buffers */
    clear_buffers( mi );

  } /* mlat[TUP] */

  hqp_free_buffers();
  hqp_free_geometry( geom );

#ifdef HQPTIME
  node0_printf("Number of gathers = %d(%d stale)\n",ngather,nstale);
#endif
} /* void hvy_pot_alg_old( su3_matrix *links ) */

/************************************************************************/
/* externally accessible function of new algorithm */
/************************************************************************/

void hvy_pot_alg_new( su3_matrix *links ) { 

  register int i; 
  register int mu;
  int disp[4]={0,0,0,1};
  int mlat[4]={0,0,0,0};
  int base[4]={0,0,0,0},nlat[4]={0,0,0,0};

#ifdef HASHLOOP
  int *hash=NULL;
#endif

  next_gather next;
#ifdef SNAKE
  next.lastraise=TUP;
  for ( mu=XUP; mu<=ZUP; mu++ ) { next.raisedir[mu]=xc[mu]; }
#endif
  msg_tag *mtag0=NULL,*mtag1=NULL;

  double_complex ctr1, ploop={0,0};
  char myname[] = "hvy_pot_alg_new";

  char smtag[MAXSMTAG]="";
#ifdef SMEARING
  sprintf(smtag,"_%d",tot_smear);
#endif

  geom = hqp_geometry( myname );
  int mi = hqp_buffer_dim();
  hqp_setup_buffers(links,mi);

  if (llat[XUP]==0) 
    local_lattice_size(llat);
  int ncs = llat[xc[XUP]]*llat[xc[YUP]]*llat[xc[ZUP]];
#ifdef HASHLOOP
  hash = setup_coord_hash( &ncs );   
#endif

#ifdef SNAKE
  raising_y=1;
  raising_z=1;
#endif

  /****************************************************************
   * 
   * Proceed to loop over local sublattices
   * using either the FLUTE or SNAKE paths for the local lattice shifts
   * and correlator contractions for all requested observables
   * 
   ***************************************************************/

#ifdef HQPTIME
  memset(hqptime,'\0',HQPTIME_TOTAL*sizeof(double));
#endif

  ngather = 0;
  nstale = 0;
  mlat[TUP] = 1;
  do { /* BEGIN while ( mlat[TUP] <= maxc[TUP] ); */
#ifdef HQPTIME
    hqptime[HQPTIME_COMM] -= dclock();
#endif
    for ( mu=XUP; mu<TUP; mu++ ) 
      nlat[mu] = mlat[mu]; 
    next.raise = NODIR; 

    /* 4) Work on sublattice mlat with contractions; decide and request which sublattice nlat is gathered next */
#if ((defined FLUTE) ^ (defined SNAKE) )
    next = find_gather ( mlat, nlat, base, next ) ;
#  ifdef DEBUG
    node0_printf("next: rd %d %d %d gg %d lr %d raise %d\n",next.raisedir[XUP],next.raisedir[YUP],next.raisedir[ZUP],next.gengather,xc[next.lastraise], next.raise);
    node0_printf("mlat %d %d %d %d max %d %d %d %d maxr2 %d raise %d insphere %d gengather %d\n",
      mlat[XUP],mlat[YUP],mlat[ZUP],mlat[TUP],maxc[XUP], maxc[YUP], maxc[ZUP], maxc[TUP], geom->max_r2,next.raise,next.insphere,next.gengather);
#  endif
#endif // END #if ((defined FLUTE) ^ (defined SNAKE) )

    /* initiate next gather */
    if ( next.gengather == GENGATHER ) { 
      mtag1 = start_general_gather_field( (void *)owline, sizeof(su3_matrix), base, EVENANDODD, gen_pt[1] );
    } else { // ELSE OF if (gengather)
      switch(next.raise) {
        case(NODIR):{ 
          break; 
        }
        case(TUP)  :{ 
          mtag0 = start_gather_field( owline, sizeof(su3_matrix),xc[TUP], EVENANDODD, gen_pt[0] ); 
          break; 
        }
        case(TDOWN):
        case(XDOWN):{ 
          node0_printf("%s(%d): No suitable gather direction, raise=%d\n",myname,this_node,next.raise); terminate(1); 
        }
        default    :{ 
          mtag1 = start_gather_field( swline, sizeof(su3_matrix),
                  ( llat[(next.raise<=TUP?xc[next.raise]:xc[OPP_DIR(next.raise)])]>1?DIR_NLL(xc[next.raise]):(xc[next.raise]) ), 
                  EVENANDODD, gen_pt[1] ); 
          break; 
        }
      }
    } // END if ( gengather == GENGATHER )
#ifdef HQPTIME
    hqptime[HQPTIME_COMM] += dclock();
#endif

#ifdef HQPTIME
    hqptime[HQPTIME_CONTRACT] -= dclock();
#endif
    int stale=1; /* stale == 1 -> all displacements out of the range ? */
    /* 5a)  Decide if we want to contractions */
    if (WANTCONTRACTIONS(mlat[TUP])) 
    {   
      double dtime=start_timing();

      /* 5b) Do proper all-to-all contractions on local lattice */
      register int cs;
      site *u=NULL;
#ifdef HASHLOOP
      for ( cs=0; cs<ncs; cs++) 
#else
      FORALLSITES(cs,u) 
#endif
      {
        site *s=NULL;
        FORALLSITES_OMP(i,s, shared(cs,ncs,u,stale) ) {

          /* Obtain displacement disp and shifted wilson line */
          register int j;
          site *t=NULL;
          if (ncs>1) {
#ifdef HASHLOOP
            j = HASH(hash_ct_from_site( s,ncs ),cs,ncs); 
            t = &(lattice[j]);
#else
            if ( site_coord(s,xc[TUP]) != site_coord(u,xc[TUP]) ) 
              continue;
            j = cs;
            t = u;
#endif
            disp_from_sites (mlat, s, t, disp );
          } else {
            j = i;
            t = s;
            memcpy(disp,mlat,4*sizeof(int));
          }
#ifdef PLOOPCOR_MEAS
          ctr1 = add_to_ploop( i, ( cs==0 && mlat[TUP] == nc[TUP] ?&ploop:NULL) );
#endif

          /* If squared radius is okay, add to correlator buffers */
          if ( ncs == 1 || hqp_disp_rsq_ok( disp, geom ) == 1 ) {
            contractions ( index_from_disp(disp), i, j, &ctr1 );
            stale=0;
          } 

        } END_LOOP_OMP;
        if ( ncs == 1 ) 
          break;
      }

    } // END if (WANTCONTRACTIONS(mlat[TUP])) 
#ifdef HQPTIME
    hqptime[HQPTIME_CONTRACT] += dclock();
#endif


#ifdef PLOOPCOR_MEAS
    if ( ploop_done == P_INDOING ) 
      ploop_done = P_NOTWRITTEN; 
#endif

    /* 6) Wait for whatever was gathered before */
    if ( SPATIAL(next.raise) ) {
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] -= dclock();
#endif
      if ( next.gengather == GENGATHER ) {
        wait_general_gather( mtag1 );
        FORALLFIELDSITES_OMP(i,default(shared)){ su3mat_copy( (su3_matrix *)gen_pt[1][i] , buffer+i ); } END_LOOP_OMP;
        memcpy(swline,buffer,sites_on_node*sizeof(su3_matrix));
        cleanup_general_gather( mtag1 );
      } else { 
        wait_gather(mtag1);
        FORALLFIELDSITES_OMP(i,default(shared)){ su3mat_copy( (su3_matrix *)gen_pt[1][i] , buffer+i ); } END_LOOP_OMP;
        memcpy(swline,buffer,sites_on_node*sizeof(su3_matrix));
        cleanup_gather( mtag1 );
      }
      ngather++;
      nstale += stale;
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] += dclock();
#endif
    } else { // if next.raise == TUP
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] -= dclock();
#endif
      if ( next.raise == TUP ) {
        wait_gather(mtag0);
        FORALLFIELDSITES_OMP(i,default(shared)){ mult_su3_nn( links+4*i+xc[TUP], (su3_matrix *)gen_pt[0][i], buffer+i ); } END_LOOP_OMP;
        memcpy(owline,buffer,sites_on_node*sizeof(su3_matrix));

        cleanup_gather( mtag0 );
        ngather++;
        nstale += stale;
      }
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] += dclock();
#endif

#ifdef HQPTIME
      hqptime[HQPTIME_OUTPUT] -= dclock();
#endif
      /* 7a)  Decide if we actually did the contractions and want to output */
      if (WANTCONTRACTIONS(mlat[TUP])) 
      {
    /* 7b) Global sums once all local lattices are done */
#ifdef PLOOPCOR_MEAS
        if ( ploop_done == P_NOTWRITTEN && mlat[TUP] == nc[TUP] ) { 
          g_dcomplexsum(&ploop); 
          node0_printf("P_LOOP%s:\t%.6e\t%.6e\n", smtag, ploop.real/volume, ploop.imag/volume );
          ploop_done = P_WRITTEN;
        }
#endif
        global_sums( mi );

        /* output only data we are interested in (and have actually computed) */
        output_all( smtag, mi, mlat[TUP] );
      } // END OF DECISION WHETHER WE DID CONTRACTIONS
fflush(stdout);
#ifdef HQPTIME
      hqptime[HQPTIME_OUTPUT] += dclock();
#endif

#ifdef HQPTIME
      hqptime[HQPTIME_COMM] -= dclock();
#endif
      /* Reset correlator buffers */
      clear_buffers( mi );

      /* 8) Copy original Wilson line into shifted Wilson line */
#ifdef PLOOPCOR_MEAS
      if ( mlat[TUP] <= maxc[TUP] || mlat[TUP] == nc[TUP]-1 )
#else
      if ( mlat[TUP] <= maxc[TUP] )
#endif
        memcpy(swline,owline,sites_on_node*sizeof(su3_matrix));
#ifdef HQPTIME
      hqptime[HQPTIME_COMM] += dclock();
#endif
      } // END ELSE BRANCH if SPATIAL(next.raise)

    /* 9) Check loop end condition and control next shift */
    for ( mu=XUP; mu<TUP; mu++ ) 
      mlat[mu] = nlat[mu]; 
    if ( next.raise == TUP ) 
      disp[TUP] = ++mlat[TUP];
#ifdef SNAKE
    next.lastraise = next.raise; 
#endif 
    if ( next.raise == NODIR ) 
      break;  
#ifdef DEBUG
    node0_printf("mlat %d %d %d %d max %d %d %d %d maxr2 %d raise %d insphere %d gengather %d\n",
      mlat[XUP],mlat[YUP],mlat[ZUP],mlat[TUP],maxc[XUP], maxc[YUP], maxc[ZUP], maxc[TUP], geom->max_r2,next.raise,next.insphere,next.gengather);
#endif
  } while ( mlat[TUP] <= geom->maxlen );

#ifdef HASHLOOP
  free_coord_hash(hash); 
#endif
  hqp_free_buffers();
  hqp_free_geometry( geom );

#ifdef HQPTIME
  node0_printf("Number of gathers = %d(%d stale)\n",ngather,nstale);
#endif
} /* void hvy_pot_alg_new( su3_matrix *links ) */



/************************************************************************/
/* internally accessible static function for all algorithms */
/************************************************************************/

#ifdef PLOOPCOR_MEAS
static inline double_complex add_to_ploop( int i, double_complex *ploop ) {

  double_complex ctr1;
#if (MILC_PRECISION == 2)
  ctr1 = trace_su3( owline+i );
#else   
  complex ctr = trace_su3( owline+i );
  ctr1.real = (double)(ctr.real);
  ctr1.imag = (double)(ctr.imag);
#endif  
  /* Calculate (smeared) Polaykov loops once it is needed */
  if ( ploop != NULL && ploop_done < P_NOTWRITTEN ) { 
    CSUM( *ploop,ctr1 );
    ploop_done=P_INDOING; 
  }
  return ( ctr1 );
} 
#endif

static inline void clear_buffers( int mysize ) {

#ifdef COULOMB 
  memset (wlc,0,mysize*sizeof(double));
#endif
#ifdef DIQUARK 
  memset (dqc,0,mysize*sizeof(double));
#endif
#ifdef PLOOPCOR_MEAS
  memset (plc,0,mysize*sizeof(double));
  memset (plcr,0,mysize*sizeof(double));
  memset (plci,0,mysize*sizeof(double));
#endif

}

static inline void contractions( int gi, int i, int j, double_complex *ctr1) {

#ifdef DEBUG
  counter[gi] +=1;
#endif

#ifdef COULOMB
  wlc[gi] += (double)realtrace_su3( owline+i,swline+j );
#  ifdef DIQUARK
  su3_matrix dqm;
  su3_adjoint(swline+j, &dqm );
  dqc[gi] += (double)realtrace_su3( owline+i,&dqm );
#  endif 
#endif
#ifdef PLOOPCOR_MEAS
  double_complex ctr2, cprod;
#  if (MILC_PRECISION == 2)
  ctr2 = trace_su3( swline+j ); /* ctr1 already done */
#  else
  complex ctr = trace_su3( swline+j ); /* ctr1 already done */
  ctr2.real = (double)(ctr.real);
  ctr2.imag = (double)(ctr.imag);
#  endif
  CMUL_J( *ctr1, ctr2, cprod );
  plc[gi] += cprod.real;
  plcr[gi]+=ctr1->real*ctr2.real;
  plci[gi]+=ctr1->imag*ctr2.imag;
#endif 

}

static inline void disp_from_index( int gi, int *disp ) {

  disp[XUP] = gi/((maxc[YUP]+1)*(maxc[ZUP]+1));
  disp[YUP] = (gi-disp[XUP]*(maxc[YUP]+1)*(maxc[ZUP]+1))/(maxc[ZUP]+1);
  disp[ZUP] =  (gi-(disp[XUP]*(maxc[YUP]+1)+disp[YUP])*(maxc[ZUP]+1));
}

static inline void disp_from_sites( int mlat[], site *s, site *t, int *disp ) {

  register int mu;
  for ( mu=XUP; mu<=ZUP; mu++) 
    disp[mu]=((mlat[mu]*llat[xc[mu]]+site_coord(t,xc[mu])-site_coord(s,xc[mu]))
              +nc[mu])%nc[mu];
}

static inline void global_sums( int mysize ) {

#ifdef DEBUG
  g_vecdoublesum(counter,mysize);
#endif 

#ifdef COULOMB
  g_vecdoublesum(wlc,mysize); 
#endif
#ifdef DIQUARK
  g_vecdoublesum(dqc,mysize); 
#endif
#ifdef PLOOPCOR_MEAS
  g_vecdoublesum(plc,mysize);
  g_vecdoublesum(plcr,mysize);
  g_vecdoublesum(plci,mysize);
#endif

}

static inline int hqp_buffer_dim (void ) {
  return ((maxc[XUP]+1)*(maxc[YUP]+1)*(maxc[ZUP]+1));
}

static void hqp_free_buffers( void ) { 

  hqp_free_su3mat_buffer(buffer);
  hqp_free_su3mat_buffer(swline);
  hqp_free_su3mat_buffer(owline);

#ifdef BUFMAT2
  if (buflvl==2) {
    hqp_free_su3mat_buffer(buffer2);
    buflvl--;
  }
#endif
#ifdef BUFMAT
  if (buflvl>=1) {
    hqp_free_su3mat_buffer(buffer1);
    buflvl--;
  }
#endif

#ifdef DEBUG
  hqp_free_dble_buffer(counter);
#endif

#ifdef COULOMB 
  hqp_free_dble_buffer(wlc); 
#endif
#ifdef DIQUARK
  hqp_free_dble_buffer(dqc); 
#endif
#ifdef PLOOPCOR_MEAS
  hqp_free_dble_buffer(plc); 
  hqp_free_dble_buffer(plcr); 
  hqp_free_dble_buffer(plci); 
#endif

}

static void hqp_setup_buffers( su3_matrix *links, int mysize ) { 

  /****************************************************************
   * 
   * Setup of buffers and prefabrication of 
   * elements to insert in static FORCE calculation (not active)
   * 
   ***************************************************************/

#ifdef PLOOPCOR_MEAS 
  ploop_done = P_NOTDONE;
#endif
  register int i; 

  owline = hqp_alloc_su3mat_buffer( 1 );
  swline = hqp_alloc_su3mat_buffer( 1 );
  buffer = hqp_alloc_su3mat_buffer( 1 );

#ifdef BUFMAT
  if (buflvl>=1) 
    buffer1 = hqp_alloc_su3mat_buffer( 1 );
#endif
#ifdef BUFMAT2
  if (buflvl==2) 
    buffer2 = hqp_alloc_su3mat_buffer( 1 );
#endif

#ifdef DEBUG
  counter = hqp_alloc_dble_buffer( mysize );
#endif

#ifdef COULOMB
  wlc = hqp_alloc_dble_buffer( mysize );
#endif
#ifdef DIQUARK
  dqc = hqp_alloc_dble_buffer( mysize );
#endif
#ifdef PLOOPCOR_MEAS
  plc  = hqp_alloc_dble_buffer( mysize );
  plcr = hqp_alloc_dble_buffer( mysize );
  plci = hqp_alloc_dble_buffer( mysize );
#endif

  /*S: 3) Initialise Wilson lines in buffer */
  FORALLFIELDSITES_OMP(i, default(shared) ){
    su3mat_copy( links+4*i+xc[TUP], owline+i );
  } END_LOOP_OMP;

// set special values for the links (coordinates, then node, then local idx in the first row) to help with debugging
#ifdef DEBUG
  set_special_row( owline );
#endif

  memcpy(swline,owline,sites_on_node*sizeof(su3_matrix)); 

}

static inline int index_from_disp ( int disp[] ) {

   return (disp[XUP]*((maxc[YUP]+1)*(maxc[ZUP]+1)) + disp[YUP]*(maxc[ZUP]+1) + disp[ZUP]);
}

static inline void output_all( char smtag[], int mi , int mlt ) {

  int disp[4];
  disp[TUP] = mlt;
  for ( register int gi=0; gi<mi; gi++ ) {
    disp_from_index( gi, disp );
    if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) {

#ifdef DEBUG
      hqp_output_corr("CCOUNTER",smtag, disp, counter[gi]);
#endif

#ifdef COULOMB 
      hqp_output_corr("POT_LOOP",smtag, disp, wlc[gi]);
#  ifdef DIQUARK
      hqp_output_corr("DIQUARK",smtag, disp, dqc[gi]);
#  endif 
#endif 
#ifdef PLOOPCOR_MEAS
      hqp_output_corr("POL_CORR",smtag, disp, plc[gi]);
      hqp_output_corr("POL_C_RE",smtag, disp, plcr[gi]);
      hqp_output_corr("POL_C_IM",smtag, disp, plci[gi]);
#endif 

    }
  } 

}

static void hqp_print_timings(char myalg[]) {

#ifdef HQPTIME
  node0_printf("Resolved timings of %s: \n\
%s hqp timed communication %.4e\n\
%s hqp timed contractions  %.4e\n\
%s hqp timed output        %.4e\n",
  myalg,
  myalg,hqptime[HQPTIME_COMM],
  myalg,hqptime[HQPTIME_CONTRACT],
  myalg,hqptime[HQPTIME_OUTPUT]);
#endif
  fflush(stdout);

}

/************************************************************************/
/* static variables, functions and macros needed only for algorithm 1 */
/************************************************************************/

static inline void contractions_old( int gi, int i, su3_matrix *sfield, double_complex *ctr1) {

#ifdef DEBUG
  counter[gi]+=1;
#endif

#ifdef COULOMB
  wlc[gi] += (double)realtrace_su3( owline+i,sfield+i );
#  ifdef DIQUARK
  su3_matrix dqm;
  su3_adjoint(sfield+i, &dqm );
  dqc[gi] += (double)realtrace_su3( owline+i,&dqm );
#  endif 
#endif
#ifdef PLOOPCOR_MEAS
  double_complex ctr2, cprod;
#  if (MILC_PRECISION == 2)
  ctr2 = trace_su3( sfield+i ); /* ctr1 already done */
#  else
  complex ctr = trace_su3( sfield+i ); /* ctr1 already done */
  ctr2.real = (double)(ctr.real);
  ctr2.imag = (double)(ctr.imag);
#  endif
  CMUL_J( *ctr1, ctr2, cprod );
  plc[gi] += cprod.real;
  plcr[gi]+=ctr1->real*ctr2.real;
  plci[gi]+=ctr1->imag*ctr2.imag;
#endif 

}

/* shift, without parallel transport, a matrix from direction "dir" */
static void shiftmat( su3_matrix *src, su3_matrix *dest, int mu ){
  register int i;
  msg_tag *mtag;
  mtag = start_gather_field( src, sizeof(su3_matrix),
      mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag);
  FORALLFIELDSITES(i){
      su3mat_copy( (su3_matrix *)gen_pt[0][i], dest+i );
  }
  cleanup_gather( mtag );
}

#ifdef GENSHIFT
/* shift, without parallel transport, a matrix from displacement "disp" */
static void general_shiftmat( su3_matrix *src, su3_matrix *dest, int *disp ){
  register int i;
  msg_tag *mtag;
  mtag = start_general_gather_field( (void *)src, sizeof(su3_matrix),
      disp, EVENANDODD, gen_pt[0] );
  wait_general_gather(mtag);
  FORALLFIELDSITES(i){
      su3mat_copy( (su3_matrix *)gen_pt[0][i], dest+i );
  }
  cleanup_general_gather( mtag );
}
#endif

/************************************************************************/
/* internally accessible static function for algorithm 2 */
/************************************************************************/

#ifdef HASHLOOP
static inline int hash_ct_from_site( site *s, int ncs ) {

  int local[4]={site_coord(s,xc[XUP])%geom->llat[xc[XUP]],
                site_coord(s,xc[YUP])%geom->llat[xc[YUP]],
                site_coord(s,xc[ZUP])%geom->llat[xc[ZUP]],
                site_coord(s,xc[TUP])%geom->llat[xc[TUP]]};
  char myname[]="hash_ct_from_site"; 

  if ( ncs == 1 ) 
    return((int)(s-&(lattice[0])));
  if ( ncs == geom->llat[xc[XUP]] )
    return(local[YUP]+geom->llat[xc[YUP]]*
          (local[ZUP]+geom->llat[xc[ZUP]]*
          (local[TUP])));
  if ( ncs == geom->llat[xc[XUP]]*geom->llat[xc[YUP]] ) 
    return(local[ZUP]+geom->llat[xc[ZUP]]*
          (local[TUP]));
  if ( ncs == geom->llat[xc[XUP]]*geom->llat[xc[YUP]]*geom->llat[xc[ZUP]] ) 
    return(local[TUP]);

  node0_printf("%s: HASH FAIL -- illegal value of ncs = %d\n",myname,ncs);
  terminate(1);
  return(-1); // to evade a needless compiler warning

}

static int * setup_coord_hash( int *ncs ) {

  /* idea: split the index such that the fast index runs over the directions 
           where llat[xc[mu]]>1 and communication is done in local lattice chunks 
           and the slow index runs over the directions where llat[xc[mu]]==1 (or 
           the cor_dir), where communication is done in steps of exactly one.
           when resolving the hash table, we premit nontrivial displacements 
           in the fast directions and require equality in the slow directions. */
  register int i,mu;
  int cs,ct; 
  site *s=NULL;

  int *hash=(int*)(malloc(sites_on_node*sizeof(int))); 
  assert( hash!=NULL );
  *ncs=sites_on_node/llat[xc[TUP]];

  for (mu=ZUP;mu>=XUP;mu--) 
    if ( llat[xc[mu]] == 1 ) {
    *ncs /= geom->llat[xc[mu]];
    } else 
      break;
  int nct = sites_on_node/(*ncs);

  if (*ncs > 1) {
    FORALLSITES_OMP(i,s, private(ct,cs,mu) ) {
      int svol=1, tvol=1;
      for ( cs=0,ct=0,mu=XUP;mu<=ZUP;mu++) {
        if (llat[xc[mu]]>1) { 
          cs+=(site_coord(s,xc[mu])%geom->llat[xc[mu]])*svol;
          svol *= geom->llat[xc[mu]];
        } else {
          ct+=(site_coord(s,xc[mu])%geom->llat[xc[mu]])*tvol;
          tvol *= geom->llat[xc[mu]];
        } 
      }
      ct+=(site_coord(s,xc[TUP])%geom->llat[xc[TUP]])*tvol;
      HASH(ct,cs,*ncs)=i;

      if ( ct < 0 || ct >= nct
        || cs < 0 || cs >= *ncs 
        || HASH(ct,cs,*ncs)<0 || HASH(ct,cs,*ncs)>=sites_on_node) { 
        printf("HASH FAIL on node %d site %d %d %d %d ct %d(%d) cs %d(%d) i %d hash %d (of %ld)\n",
               this_node,
               site_coord(s,xc[XUP]),site_coord(s,xc[YUP]),site_coord(s,xc[ZUP]),site_coord(s,xc[TUP]),
               ct,nct,cs,*ncs,i,HASH(ct,cs,*ncs),sites_on_node); terminate(1); 
      }
    } END_LOOP_OMP;
  } else {
    FORALLSITES_OMP(i,s, private(ct,cs,mu) ) {
      HASH(i,0,1)=i;
    } END_LOOP_OMP;
  }
  fflush(stdout); g_sync();

  return(hash);
}

static void free_coord_hash( int *hash ) {
  free(hash);
  hash=NULL;
}
#endif

/*****************************************************************
 * 
 * This code required a full overhaul to be useful for large lattices.
 * In this ovehaul, we replaced the previous communication structure 
 * with a completely new setup, that is very different from the previous 
 * MILC v6 code. Namely, complete one-node local-lattices are gathered 
 * asynchronously. Then all possible contractions between two local 
 * lattices are performed and accumulated in different locations of the 
 * correlator arrays.
 * 
 * We have two concurring communcation paradigms, using overlapping 
 * general and standard (single direction) gathers.
 * 1) the pure SNAKE paradigm, i.e. walking in even x-slices as
 *  >------+
 *  ^------<
 *  >------^
 * and in odd x-slices
 *  v------<
 *  >------v
 *  +------<
 * which is optimal for walking the full cube, using minimal 
 * numbers of exclusively next neighbor NLL gathers, but 
 * wasteful for just the sphere, using in 8/( 4pi/3 ) ~ 6/pi 
 * times too many gathers (in the limit of single-site local 
 * lattices) and many wasteful loops where nothing happens.
 * Note: this has been fixed now, by reverting to a FLUTE-type step, 
 * since it produced errors in some weird non-symmetric geometric
 * geometries. Still, the algorithm is prone to errors, so you better 
 * make a test run first or use FLUTE.
 * 2) the FLUTE paradigm, i.e. walking in each x-slice in 
 * the same way along the z-direction, i.e.
 *  g
 *  >-g
 *  >--g
 *  >---g 
 * which needs general gathers wherever we have g in the 
 * sketch. As a consequence, this scheme is more expensive 
 * in terms of communication for the case of the full cube.
 *
 ****************************************************************/

#ifdef FLUTE
#define FLUTEZSHIFT(mlat,llat,maxc) \
  ( mlat[ZUP]*llat[xc[ZUP]] < maxc[ZUP] ) 
#define FLUTEYSHIFT(mlat,llat,max_x) \
  ( mlat[YUP]*llat[xc[YUP]] < maxc[YUP] ) 
#define FLUTEXSHIFT(mlat,llat,max_x) \
  ( mlat[XUP]*llat[xc[XUP]] < maxc[XUP] )

static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) {
  int mu;
  //node0_printf("###111\n\%d gathers (%d stale) done \nmlat %d %d %d %d parsing\n",
  //  ngather,nstale, mlat[XUP],mlat[YUP],mlat[ZUP],mlat[TUP]);
  for ( this.insphere=OUT, this.gengather=SIMGATHER; this.insphere<IN && this.raise<TUP; ) {
    if ( llat[xc[ZUP]] < nc[ZUP] && FLUTEZSHIFT(mlat,llat,maxc) ) { 
      /* gather next xc[ZUP] local lattice */
      this.raise = ZUP;
      this.gengather = SIMGATHER;
    } else { if ( llat[xc[YUP]] < nc[YUP] && FLUTEYSHIFT(mlat,llat,maxc) ) { 
      /* gather next xc[YUP] local lattice */
      this.raise = YUP;
      nlat[ZUP] = 0;
      this.gengather = GENGATHER;
    } else { if ( llat[xc[XUP]] < nc[XUP] && FLUTEXSHIFT(mlat,llat,maxc) ) { 
      /* gather next xc[XUP] local lattice */
      this.raise = XUP;
      nlat[ZUP] = 0;
      nlat[YUP] = 0;
      this.gengather=GENGATHER;
    } else { 
      if ( mlat[TUP] < geom->maxlen ) 
        /* gather next wilson line segment from xc[TUP] cor-dir */
      {
        this.raise=TUP;
        nlat[ZUP] = 0;
        nlat[YUP] = 0;
        nlat[XUP] = 0;
        this.gengather = SIMGATHER;
      }
    }}}
    if ( this.raise > NODIR && this.raise < TUP ) 
      nlat[this.raise]++; 
    for ( mu=XUP; mu<TUP; mu++ ) 
      base[xc[mu]] = nlat[mu]*llat[xc[mu]]; 

    this.insphere=INSPHERE(base,llat,geom->max_r2); /* within sphere? */
    //node0_printf("###222\nnlat %d %d %d %d raise %d insphere %d gengather %d\n",
    //  nlat[XUP],nlat[YUP],nlat[ZUP],nlat[TUP],this.raise,this.insphere,this.gengather);
    while ( this.insphere == OUT && this.raise > NODIR && this.raise < TUP ) {
      nlat[this.raise] = 0;
      base[xc[this.raise]] = 0;
      if ( this.raise > XUP && this.raise < TUP ) {
        nlat[--this.raise]++;
        base[xc[this.raise]] += llat[xc[this.raise]];
    //node0_printf("###333\nnlat %d %d %d %d raise %d insphere %d gengather %d\n",
    //  nlat[XUP],nlat[YUP],nlat[ZUP],nlat[TUP],this.raise,this.insphere,this.gengather);
      } else {
        this.raise = TUP;
        nlat[ZUP] = 0;
        nlat[YUP] = 0;
        nlat[XUP] = 0;
        this.gengather = SIMGATHER;
        break;
      }
      if (base[xc[this.raise]] <=maxc[this.raise]) { 
        this.insphere = INSPHERE(base,llat,geom->max_r2); /* within sphere? */
        this.gengather = GENGATHER;
      } else 
        if (nlat[this.raise] == 1 && llat[xc[this.raise]] < nc[this.raise]) {
          this.insphere = IN;
          this.gengather = SIMGATHER;
        }
    }
  }
  //node0_printf("###444\nnlat %d %d %d %d raise %d insphere %d gengather %d\n",
  //  nlat[XUP],nlat[YUP],nlat[ZUP],nlat[TUP],this.raise,this.insphere,this.gengather);
  return (this);
}
#endif

#ifdef SNAKE
#define SNAKEZSHIFT(mlat,llat,maxc) \
  ( ( raising_z == 1 && mlat[ZUP]*llat[xc[ZUP]] < maxc[ZUP] ) \
 || ( raising_z == 0 && mlat[ZUP]*llat[xc[ZUP]] > 0 ) ) 
#define SNAKEYSHIFT(mlat,llat,maxc) \
  ( ( raising_y == 1 && mlat[YUP]*llat[xc[YUP]] < maxc[YUP] ) \
 || ( raising_y == 0 && mlat[YUP]*llat[xc[YUP]] > 0 ) ) 
#define SNAKEXSHIFT(mlat,llat,maxc) \
  ( mlat[XUP]*llat[xc[XUP]] < maxc[XUP] )

static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) {
  int mu; 
  for ( this.insphere=OUT, this.gengather=SIMGATHER; this.insphere<IN && this.raise < TUP; ) {
    //node0_printf("###111\n\%d gathers (%d stale) done lastraise %d raising (y,z) ? (%d,%d) \nmlat %d %d %d %d parsing\n",
    //  ngather,nstale, xc[this.lastraise], raising_y,raising_z, mlat[XUP],mlat[YUP],mlat[ZUP],mlat[TUP]);
    if ( llat[xc[ZUP]] < nc[ZUP] && SNAKEZSHIFT(mlat,llat,maxc) ) { 
      /* gather next z local lattice */
      if ( mlat[ZUP] == 0 ) {
        raising_z=1;
      } else 
        if ( (mlat[ZUP]*llat[xc[ZUP]] >= maxc[ZUP]) ) 
          raising_z=0;
      this.raisedir[ZUP] = ( raising_z >= 1 ? ZUP : ZDOWN); 
      this.raise = this.raisedir[ZUP];
    } else { 
      if ( llat[xc[YUP]] < nc[YUP] && SNAKEYSHIFT(mlat,llat,maxc) ) { 
        /* gather next y local lattice */
        if ( mlat[YUP] == 0 ) {
          raising_y=1;
        } else 
          if ( mlat[YUP]*llat[xc[YUP]] >= maxc[YUP] ) 
            raising_y=0;
        this.raisedir[YUP] = ( raising_y >= 1 ? YUP : YDOWN); 
        this.raise = this.raisedir[YUP];
        raising_z=1-raising_z;

      } else { 
        if ( llat[xc[XUP]] < nc[XUP] && SNAKEXSHIFT(mlat,llat,maxc) ) { 
          /* gather next x local lattice from xc[XUP] dir */ 
          this.raise = XUP; 
          if ( mlat[YUP] == 0 ) 
            raising_y=1;
          if ( mlat[YUP]*llat[xc[YUP]] >= maxc[YUP] ) 
            raising_y=0;
          if ( mlat[ZUP] == 0 ) 
            raising_z=1;
          if ( mlat[ZUP]*llat[xc[ZUP]] >= maxc[ZUP] ) 
            raising_z=0;
        } else {
          if ( mlat[TUP] < geom->maxlen ) 
          {    
            this.raise=TUP; 
            raising_z=1;
            raising_y=1;
            nlat[TUP]++; 
            nlat[ZUP]=0; 
            nlat[YUP]=0; 
            nlat[XUP]=0; 
          }
        }
      }
    }

    if ( SPATIAL(this.raise) ) { 
      if ( this.raise<TUP ) { 
        nlat[this.raise]++; 
      } else 
        nlat[OPP_DIR(this.raise)]--; 
    } 
    for ( mu=XUP; mu<TUP; mu++ ) 
      base[xc[mu]]=nlat[mu]*llat[xc[mu]]; 
    this.insphere=INSPHERE(base,llat,geom->max_r2); // within sphere? 
    //node0_printf("###222\nnlat %d %d %d %d raise %d insphere %d gengather %d\n",
    //             nlat[XUP],nlat[YUP],nlat[ZUP],nlat[TUP],xc[this.raise],this.insphere,this.gengather);

    while ( this.insphere == OUT && SPATIAL(this.raise) ) {
      if ( this.raise < TUP ) { 
        if ( this.raise == XUP ) {
          if ( nlat[XUP]*llat[xc[XUP]] < maxc[XUP] ) {
            raising_z=1;
            raising_y=1;
            nlat[ZUP]=0; 
            nlat[YUP]=0; 
          } 
        }
        if ( this.raise == YUP ) {
          if ( nlat[ZUP] > 0 ) { 
            raising_z=1;
            nlat[ZUP]=0; 
          } else 
            if ( mlat[YUP]*llat[xc[YUP]] >= maxc[YUP] ) {
              raising_z=0;
              raising_y=1;
              nlat[ZUP]=0; 
              nlat[YUP]=0; 
            }
        }
        if ( this.raise == ZUP ) {
          nlat[ZUP]=0;
          raising_z=1;
          if ( mlat[YUP]*llat[xc[YUP]] >= maxc[YUP] ) {
            raising_y=1;
            nlat[YUP]=0; 
            if ( nlat[XUP]*llat[xc[XUP]] < maxc[XUP] ) {
              this.raise = XUP;
              nlat[XUP]++;
            } else {
              if ( mlat[TUP] < geom->maxlen ) 
              {    
                this.raise=TUP; 
                raising_z=1;
                raising_y=1;
                nlat[TUP]++; 
                nlat[ZUP]=0; 
                nlat[YUP]=0; 
                nlat[XUP]=0; 
              } else 
                this.raise = NODIR;
            }
          } else {
            this.raise = YUP;
            nlat[YUP]++;
            for ( mu=XUP; mu<TUP; mu++ ) 
              base[xc[mu]]=nlat[mu]*llat[xc[mu]]; 
            this.insphere=INSPHERE(base,llat,geom->max_r2); // within sphere? 
            if ( this.insphere == OUT ) {
              raising_y=1;
              nlat[YUP]=0; 
              if ( nlat[XUP]*llat[xc[XUP]] < maxc[XUP] ) {
                this.raise = XUP;
                nlat[XUP]++;
              } else {
                if ( mlat[TUP] < geom->maxlen ) 
                {    
                  this.raise=TUP; 
                  raising_z=1;
                  raising_y=1;
                  nlat[TUP]++; 
                  nlat[ZUP]=0; 
                  nlat[YUP]=0; 
                  nlat[XUP]=0; 
                } else 
                  this.raise = NODIR;
              }
            }
          }
        }
      }
      for ( mu=XUP; mu<TUP; mu++ ) 
        base[xc[mu]]=nlat[mu]*llat[xc[mu]]; 
      this.insphere=INSPHERE(base,llat,geom->max_r2); // within sphere? 
      if ( SPATIAL(this.raise) ) 
        this.gengather = GENGATHER;
    }

    if ( this.raise == NODIR ) { break; }
  }
  //node0_printf("###444\nnlat %d %d %d %d raise %d insphere %d gengather %d\n",
  //             nlat[XUP],nlat[YUP],nlat[ZUP],nlat[TUP],xc[this.raise],this.insphere,this.gengather);
  return(this);
}
#endif
#if ( !( (defined FLUTE) ^ (defined SNAKE) ) )
static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) {
  terminate(1);
  return(NULL);
}
#endif

#endif
