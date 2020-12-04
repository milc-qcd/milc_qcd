/****** new_hvy_pot.c  -- ******************/

/***********************************************************************
* Heavy quark potential
* MIMD version 7
*
************************************************************************
*
* New Strategy:
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
************************************************************************/

/***********************************************************************
*
* Iteration 1)
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
* Iteration 2)
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

#define HASHLOOP
#ifdef HASHLOOP
#define HASH(ct,cs) hash[ct*ncs+cs]
#endif

#ifdef NEW_HVY_POT

#define NOTDONE 0
#define DONE 1

/* FLUTE and SNAKE are two different schemes for traversing the lattice
 * that dramatically improve performance on large lattices for massively
 * parallelized calculations; one has to be used; for more details see
 * the detailed description at the end of this file */
#define SNAKE
#if ( !( (defined FLUTE) ^ (defined SNAKE) ) )
BOMB THE COMPILE
#endif

/* Indicates the overall type of the next gather operationr */
#define GENGATHER 1
#define SIMGATHER 0
#define NOGATHER -1

/* Indicates if some points are inside of the 3-sphere with radius max_r */
#define OUT 0
#define IN  1

/* Indicates the status of the Polyakov loop calculation */
#define P_NOTDONE 0
#define P_INDOING 1
#define P_NOTWRITTEN 2
#define P_WRITTEN 3

#ifdef FORCE
static void free_gen(su3_matrix *genT) { free(genT); }

static su3_matrix *set_gen(void) ;
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

/* static variables that are used extensively in determining the gathers,
 * which are not updated anymore after their initialization */
static int llat[4];
static int nll[3];
static int max_r2;
static int nstale, ngather; /* useful only for choosing the better communication paradigm */

/* main new version of the routine */
void new_hvy_pot( su3_matrix *links ) {
    register int i,j,gi;
    register int mu, dok;
    register int mi=(maxc[XUP]+1)*(maxc[YUP]+1)*(maxc[ZUP]+1);
    int disp[3]={0,0,0};
    int mlat[4]={0,0,0,0};
    int base[4]={0,0,0,0},nlat[4]={0,0,0,0};
#ifdef HASHLOOP
    int ct,cs,ncs;
    int *hash=NULL;
#endif
    next_gather next;
#ifdef SNAKE
    next.lastraise=TUP;
    for ( mu=XUP; mu<=ZUP; mu++ ) { next.raisedir[mu]=xc[mu]; }
#endif
    Real drad=0.;
    int stale; /* useful only for choosing the better communication paradigm */
    int ploop_done=P_NOTDONE;
    site *s=NULL,*t=NULL;
    msg_tag *mtag0=NULL,*mtag1=NULL;
    su3_matrix *owline=NULL,*swline=NULL,*buffer=NULL;
    complex ctr;
    double_complex ctr1, ploop={0,0};
    complex ctr2, cprod;
#ifdef COULOMB
    double *wlc=NULL;
#endif
#ifdef DIQUARK
    su3_matrix dqm;
    double *dqc=NULL;
#endif
#ifdef PLOOPCOR_MEAS
    double *plc=NULL,*plcr=NULL,*plci=NULL;
#endif
    char myname[] = "new_hvy_pot";

    /*******************************************************
     *
     * Start of header:
     *
     * Text output concerning the geometry and its consistency.
     * In particular, the correlation directions 'xc[mu]' are
     * permuted from the original lattice directions 'mu'
     * to permit static quark correlations with arbitrary
     * correlation time direction and arbitrary anisotropic
     * direction.
     *
     ******************************************************/
    node0_printf("%s with %d times smeared links\n",myname,tot_smear);
#ifdef DEBUG
    node0_printf(" [min,max](CT)[%d] = [%d,%d], \
                \n max(C0)[%x] = %d, max(C1)[%x] = %d, max(C2)[%x] = %d, \
                \n max_r = %.3g\n",cor_dir, min_ct,maxc[TUP], xc[XUP],maxc[XUP] ,xc[YUP],maxc[YUP] ,xc[ZUP],maxc[ZUP] ,max_r);
#endif
    assert ( maxc[XUP] <= nc[XUP]/2 && maxc[YUP] <= nc[YUP]/2 && maxc[ZUP] <= nc[ZUP]/2 );
    assert ( min_ct <= maxc[TUP] && maxc[TUP] <= nc[TUP] );
    assert ( max_r >= 0 );
    if ( max_r == 0 ) {
#ifndef ANISOTROPY
      max_r=sqrt(maxc[XUP]*maxc[XUP] + maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP] ) + 1e-6;
#else
#ifndef ONEDIM_ANISO_TEST
      max_r=sqrt( ( ani_dir != cor_dir ? ani_xiq*ani_xiq : 1 ) * maxc[XUP]*maxc[XUP] + maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP] ) + 1e-6;
#else
      max_r=sqrt( ( ani_dir != cor_dir ? ani_xiq*ani_xiq : iso_xiq ) * maxc[XUP]*maxc[XUP] + (iso_xiq*iso_xiq)*(maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP]) ) + 1e-6;
#endif
#endif
      node0_printf("max_r being automatically  reset to = %g\n",max_r);
    }
    for (mu=XUP; mu <TUP; mu++) {
#ifndef ANISOTROPY
      if ( max_r < maxc[mu] ) {
        maxc[mu] = floor(max_r);
#else
      if ( mu==XUP && max_r < ( ani_dir != cor_dir ? ani_xiq : 1 ) * maxc[mu] ) {
        maxc[mu] = floor(max_r/ani_xiq);
#ifndef ONEDIM_ANISO_TEST
      if ( max_r < ( mu == XUP && ani_dir != cor_dir ? ani_xiq : iso_xiq ) * maxc[mu] ) {
        maxc[mu] = floor(max_r/( mu == XUP && ani_dir != cor_dir ? ani_xiq : iso_xiq ));
#endif
#endif
        node0_printf("Decrease to max(C%d) = %d\n",mu,maxc[mu]);
      }
    }
    max_r2=floor(max_r*max_r);
    local_lattice_size(llat);
    nll[XUP]=nc[XUP]/llat[xc[XUP]]; nll[YUP]=nc[YUP]/llat[xc[YUP]]; nll[ZUP]=nc[ZUP]/llat[xc[ZUP]];
#ifdef DEBUG
    if ( numnodes() > 1 ) {
      node0_printf("####################################\
                    \n Using a distributed grid of %dx%dx%d(x%d) local lattices\n",
                    nll[XUP],nll[YUP],nll[ZUP],nc[TUP]/llat[xc[TUP]]);
    }
    node0_printf("####################################\n\n");
#endif
fflush(stdout);

    /****************************************************************
     *
     * End of header
     * Setup of buffers and prefabrication of
     * elements to insert in static FORCE calculation (not active)
     *
     ***************************************************************/

#ifdef HASHLOOP
    hash=(int*)(malloc(sites_on_node*sizeof(int)));
    assert( hash!=NULL );
    ncs=sites_on_node/llat[xc[TUP]];
    FORALLSITES_OMP(i,s, private(ct,cs,mu) ) {
      for (cs=(site_coord(s,xc[XUP])%llat[xc[XUP]]),mu=YUP;mu<=ZUP;mu++) {
        cs*=llat[xc[mu]];
        cs+=(site_coord(s,xc[mu])%llat[xc[mu]]);
      }
      ct=site_coord(s,xc[TUP])%llat[xc[TUP]];
      HASH(ct,cs)=i;
      if ( ct<0 || ct >= nc[TUP]
        || cs <0 || cs >= ncs
        || HASH(ct,cs)<0 || HASH(ct,cs)>=sites_on_node) {
        printf("HASH FAIL on node %d site %d %d %d %d ct %d cs %d i %d hash %d (of %ld)\n",
               this_node,
               site_coord(s,xc[XUP]),site_coord(s,xc[YUP]),site_coord(s,xc[ZUP]),site_coord(s,xc[TUP]),
               ct,cs,i,HASH(ct,cs),sites_on_node); terminate(1);
      }
    } END_LOOP_OMP;
    g_sync();
node0_printf("HASH TABLE ASSEMBLED\n");
#endif

    owline = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    swline = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    buffer = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    assert( (owline!=NULL)&&(swline!=NULL)&&(buffer!=NULL) );
    memset (buffer,0,sites_on_node*sizeof(su3_matrix));
    memset (owline,0,sites_on_node*sizeof(su3_matrix));
    memset (swline,0,sites_on_node*sizeof(su3_matrix));

#ifdef COULOMB
    wlc = (double *)malloc(mi*sizeof(double));
    assert( (wlc != NULL ) );
    memset (wlc,0,mi*sizeof(double));
#endif

#ifdef DIQUARK
    dqc = (double *)malloc(mi*sizeof(double));
    assert( (dqc != NULL ) );
    memset (dqc,0,mi*sizeof(double));
#endif

#ifdef PLOOPCOR_MEAS
    plc = (double *)malloc(mi*sizeof(double));
    plcr= (double *)malloc(mi*sizeof(double));
    plci= (double *)malloc(mi*sizeof(double));
    assert( (plc != NULL ) && (plcr != NULL ) && (plci != NULL ) );
    memset (plc,0,mi*sizeof(double));
    memset (plcr,0,mi*sizeof(double));
    memset (plci,0,mi*sizeof(double));
#endif

    /* 3) Initialise Wilson lines in buffer */
    FORALLFIELDSITES_OMP(i, default(shared) ){
      su3mat_copy( links+4*i+xc[TUP], owline+i );
      su3mat_copy( owline+i, swline+i );
    } END_LOOP_OMP;

    /****************************************************************
     *
     * End of setup
     * Proceed to loop over local sublattices
     * using either the FLUTE or SNAKE paths for the local lattice shifts
     * and correlator contractions for all requested observables
     *
     ***************************************************************/
    nstale = 0;
    ngather = 0;
    mlat[xc[TUP]] = 1;
    do { /* BEGIN while ( mlat[xc[TUP]] <= maxc[TUP] ); */
      for ( mu=XUP; mu<ZUP; mu++ ) { nlat[xc[mu]] = mlat[xc[mu]]; }
      next.raise = NODIR;
      /* 4) Work on sublattice mlat with contractions; decide and request which sublattice nlat is gathered next */
#ifdef FLUTE
      next = find_gather ( mlat, nlat, base, next ) ;
#endif // END #ifdef FLUTE

#ifdef SNAKE
      next = find_gather ( mlat, nlat, base, next ) ;
#endif // END #ifdef SNAKE

      if ( next.gengather == GENGATHER ) {
        mtag1 = start_general_gather_field( (void *)owline, sizeof(su3_matrix), base, EVENANDODD, gen_pt[1] );
      } else { // ELSE OF if (gengather)
        switch(next.raise) {
          case(NODIR):{ break; }
          case(TUP)  :{
            mtag0 = start_gather_field( owline, sizeof(su3_matrix),xc[TUP], EVENANDODD, gen_pt[0] );
            break;
          }
          case(TDOWN):
          case(XDOWN):{
            node0_printf("%s(%d): No suitable gather direction, raise=%d\n",myname,this_node,next.raise); terminate(1);
          }
          default    :{
            mtag1 = start_gather_field( swline, sizeof(su3_matrix),DIR_NLL(xc[next.raise]), EVENANDODD, gen_pt[1] );
            break;
          }
        }
      } // END if ( gengather == GENGATHER )
      memset (buffer,0,sites_on_node*sizeof(su3_matrix));

      stale=1; /* stale == 1 -> all displacements out of the range ? */

      /* 5a)  Decide if we want to contractions */
#ifdef PLOOPCOR_MEAS
      if ( mlat[xc[TUP]] <= maxc[TUP] || mlat[xc[TUP]] == nc[TUP] )
#else
      if ( mlat[xc[TUP]] <= maxc[TUP] )
#endif
      {
        double dtime=start_timing();
        /* 5b) Do proper all-to-all contractions on local lattice */
        FORALLSITES(i,s) {
#ifdef PLOOPCOR_MEAS
#if (MILC_PRECISION == 2)
          ctr1 = trace_su3( owline+i );
#else
          ctr = trace_su3( owline+i );
          ctr1.real = (double)(ctr.real);
          ctr1.imag = (double)(ctr.imag);
#endif
          /* Calculate (smeared) Polaykov loops once it is needed */
          if ( ploop_done < P_NOTWRITTEN && mlat[xc[TUP]] == nc[TUP] ) {
            CSUM(ploop,ctr1);
            ploop_done=P_INDOING;
          }
#endif

#ifdef HASHLOOP
          ct=(site_coord(s,xc[TUP])%llat[xc[TUP]]);
#ifdef OMP
          {
#pragma omp parallel for private(cs,j,t, disp,dok,gi,mu, ctr,ctr2,cprod) shared(ct,ctr1) schedule(guided)
          for ( cs=0; cs<ncs; cs++) {
#else
          for ( cs=0; cs<ncs; cs++) {
#endif //ifdef OMP
            j=HASH(ct,cs);
            t=&(lattice[j]);
            if ( ct<0 || ct >= nc[TUP]
              || cs <0 || cs >= ncs
              || t < &(lattice[0]) || t > &(lattice[sites_on_node-1])
              || j <0 || j >= sites_on_node
              || HASH(ct,cs)<0 || HASH(ct,cs)>=sites_on_node) {
              printf("HASH FAIL on node %d site %d %d %d %d ct %d cs %d j %d hash %d (of %ld)\n",
                     this_node,
                     site_coord(t,xc[XUP]),site_coord(t,xc[YUP]),site_coord(t,xc[ZUP]),site_coord(t,xc[TUP]),
                     ct,cs,j,HASH(ct,cs),sites_on_node); terminate(1);
            }
#else
          FORALLSITES_OMP(j,t,private(disp,dok,gi,mu, ctr,ctr2,cprod ) shared(ctr1) ) {
            if ( site_coord(s,xc[TUP]) != site_coord(t,xc[TUP]) ) continue;
#endif // HASHLOOP

            dok=1; drad=0.;
          /* Compute displacement disp, and index gi of the correlator buffers */
            for ( mu=XUP; mu<=ZUP; mu++) {
              disp[mu]=((mlat[xc[mu]]*llat[xc[mu]]+site_coord(t,xc[mu])-site_coord(s,xc[mu]))+nc[mu])%nc[mu];
              if ( disp[mu] > maxc[mu] || disp[mu] < 0 ) {
                dok = 0;
              } else { /* add up to displacement radius squared */
                drad +=
#ifdef ANISOTROPY
#ifdef ONEDIM_ANISO_TEST
                        ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) *
#else
                        ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) *
#endif
#endif
                        disp[mu]*disp[mu];
              }
            }

            /* Are interested in this displacement? Then fill correlator buffers */
            if ( dok == 1 && drad <= max_r2 ) {
              gi = disp[XUP]*((maxc[YUP]+1)*(maxc[ZUP]+1)) + disp[YUP]*(maxc[ZUP]+1) + disp[ZUP];
              stale=0;

#ifdef COULOMB
              wlc[gi] += (double)realtrace_su3( owline+i,swline+j );
#endif

#ifdef DIQUARK
              su3_adjoint(swline+j, &dqm );
              dqc[gi] += (double)realtrace_su3( owline+i,&dqm );
#endif


#ifdef PLOOPCOR_MEAS
#if (MILC_PRECISION == 2)
              ctr2 = trace_su3( swline+j ); /* ctr1 already done */
#else
              ctr = trace_su3( swline+j ); /* ctr1 already done */
              ctr2.real = (double)(ctr.real);
              ctr2.imag = (double)(ctr.imag);
#endif
              CMUL_J( ctr1, ctr2, cprod );

              plc[gi] += cprod.real;
              plcr[gi]+=ctr1.real*ctr2.real;
              plci[gi]+=ctr1.imag*ctr2.imag;
#endif

            } // END if ( dok > 0 && drad <= max_r2 )
#ifdef HASHLOOP
          }
#ifdef OMP
          }
#endif
#else
          } END_LOOP_OMP;  // END FORALLSITES_OMP(j,t,private(...) )
#endif // HASHLOOP
        } // END FORALLSITES(i,s)
#ifdef DEBUG
        char a2atag[64];
        sprintf(a2atag,"all-to-all correlations on shift %d,%d,%d,%d",
                mlat[xc[XUP]],mlat[xc[YUP]],mlat[xc[ZUP]],mlat[xc[TUP]]);
        print_timing(dtime, a2atag);
#endif
      } // END OF DECISION WHETHER TO DO CONTRACTIONS


#ifdef PLOOPCOR_MEAS
      if ( ploop_done == P_INDOING ) { ploop_done = P_NOTWRITTEN; }
#endif

      /* 6) Wait for whatever was gathered before */
      if ( SPATIAL(next.raise) ) {
        if ( next.gengather == GENGATHER ) {
          wait_general_gather( mtag1 );
          FORALLFIELDSITES(i){ su3mat_copy( (su3_matrix *)gen_pt[1][i] , buffer+i ); }
          FORALLFIELDSITES(i){ su3mat_copy( buffer+i, swline+i ); }
          cleanup_general_gather( mtag1 );
        } else {
          wait_gather(mtag1);
          FORALLFIELDSITES(i){ su3mat_copy( (su3_matrix *)gen_pt[1][i] , buffer+i ); }
          FORALLFIELDSITES(i){ su3mat_copy( buffer+i, swline+i ); }
          cleanup_gather( mtag1 );
        }
        ngather++;
        nstale += stale;
      } else { // if next.raise == TUP
        if ( next.raise == TUP ) {
          wait_gather(mtag0);
          FORALLFIELDSITES(i){ mult_su3_nn( links+4*i+xc[TUP], (su3_matrix *)gen_pt[0][i], buffer+i ); }
          FORALLFIELDSITES(i){ su3mat_copy( buffer+i, owline+i ); }

          cleanup_gather( mtag0 );
          ngather++;
          nstale += stale;
        }

      /* 7a)  Decide if we did the contractions */
#ifdef PLOOPCOR_MEAS
        if ( mlat[xc[TUP]] <= maxc[TUP] || mlat[xc[TUP]] == nc[TUP] )
#else
        if ( mlat[xc[TUP]] <= maxc[TUP] )
#endif
        {
      /* 7b) Global sums once all local lattices are done */
#ifdef PLOOPCOR_MEAS
          if ( ploop_done == P_NOTWRITTEN && mlat[xc[TUP]] == nc[TUP] ) {
            g_dcomplexsum(&ploop);
#ifdef SMEARING
            node0_printf("P_LOOP_%d:\t%.6e\t%.6e\n", tot_smear, ploop.real/volume, ploop.imag/volume );
#else
            node0_printf("P_LOOP:\t%.6e\t%.6e\n", ploop.real/volume, ploop.imag/volume );
#endif
            ploop_done = P_WRITTEN;
          }
#endif

#ifdef COULOMB
          g_vecdoublesum(wlc,mi);
#endif

#ifdef DIQUARK
          g_vecdoublesum(dqc,mi);
#endif

#ifdef PLOOPCOR_MEAS
          g_vecdoublesum(plc,mi);
          g_vecdoublesum(plcr,mi);
          g_vecdoublesum(plci,mi);
#endif

          for ( gi=0; gi<mi; gi++ ) {
            disp[XUP] = gi/((maxc[YUP]+1)*(maxc[ZUP]+1));
            disp[YUP] = (gi-disp[XUP]*(maxc[YUP]+1)*(maxc[ZUP]+1))/(maxc[ZUP]+1);
            disp[ZUP] =  (gi-(disp[XUP]*(maxc[YUP]+1)+disp[YUP])*(maxc[ZUP]+1));
            dok = 1;
            drad = 0.;
          /* output only data we are interested in (and have actually computed) */
            for ( mu=XUP; mu<=ZUP; mu++ ) {
              if (disp[mu] > maxc[mu] ) { dok = 0;
              } else {
                drad +=
#ifdef ANISOTROPY
#ifdef ONEDIM_ANISO_TEST
                        ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) *
#else
                        ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) *
#endif
#endif
                        disp[mu]*disp[mu];
              }
            }
            if ( dok == 1 && drad <= max_r2 ) {

#ifdef COULOMB
#ifdef SMEARING
              node0_printf("POT_LOOP_%d:   %d  %d  %d   %d \t%.6e\n", tot_smear,
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], wlc[gi]/volume );
#else
              node0_printf("POT_LOOP:  %d  %d  %d   %d \t%.6e\n",
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], wlc[gi]/volume );
#endif
#endif

#ifdef DIQUARK
#ifdef SMEARING
              node0_printf("DIQUARK_%d:    %d  %d  %d   %d \t%.6e\n", tot_smear,
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], dqc[gi]/volume );
#else
              node0_printf("DIQUARK:    %d  %d  %d   %d \t%.6e\n",
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], dqc[gi]/volume );
#endif
#endif

#ifdef PLOOPCOR_MEAS
#ifdef SMEARING
              node0_printf("POL_CORR_%d:   %d  %d  %d   %d \t%.6e\n", tot_smear,
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plc[gi]/volume );
              node0_printf("POL_C_RE_%d:   %d  %d  %d   %d \t%.6e\n", tot_smear,
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plcr[gi]/volume );
              node0_printf("POL_C_IM_%d:   %d  %d  %d   %d \t%.6e\n", tot_smear,
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plci[gi]/volume );
#else
              node0_printf("POL_CORR:   %d  %d  %d   %d \t%.6e\n",
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plc[gi]/volume );
              node0_printf("POL_C_RE:   %d  %d  %d   %d \t%.6e\n",
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plcr[gi]/volume );
              node0_printf("POL_C_IM:   %d  %d  %d   %d \t%.6e\n",
                disp[XUP], disp[YUP], disp[ZUP], mlat[xc[TUP]], plci[gi]/volume );
#endif
#endif

            } // END if ( dok == 1 && drad <= max_r2 )
          } // END for ( gi=0; gi<mi; gi++ )
        } // END OF DECISION WHETHER WE DID CONTRACTIONS
fflush(stdout);

        /* Reset correlator buffers */
#ifdef COULOMB
        memset (wlc,0,mi*sizeof(double));
#endif

#ifdef DIQUARK
        memset (dqc,0,mi*sizeof(double));
#endif

#ifdef PLOOPCOR_MEAS
        memset (plc,0,mi*sizeof(double));
        memset (plcr,0,mi*sizeof(double));
        memset (plci,0,mi*sizeof(double));
#endif

      /* 8) Copy original Wilson line into shifted Wilson line */
#ifdef PLOOPCOR_MEAS
        if ( mlat[xc[TUP]] <= maxc[TUP] || mlat[xc[TUP]] == nc[TUP]-1 )
#else
        if ( mlat[xc[TUP]] <= maxc[TUP] )
#endif
        {
          FORALLFIELDSITES_OMP(i,default(shared)){ su3mat_copy( owline+i,swline+i ); } END_LOOP_OMP;
        }
      } // END ELSE BRANCH if SPATIAL(next.raise)

      /* 9) Check loop end condition and control next shift */
      for ( mu=XUP; mu<TUP; mu++ ) { mlat[xc[mu]] = nlat[xc[mu]]; }
      if ( next.raise == TUP ) { mlat[xc[next.raise]]++; }
#ifdef SNAKE
      next.lastraise = next.raise;
#endif
     if ( next.raise == NODIR ) { break; }
    }
#ifdef PLOOPCOR_MEAS
    while ( mlat[xc[TUP]] <= nc[TUP] );
#else
    while ( mlat[xc[TUP]] <= maxc[TUP] );
#endif

#ifdef COULOMB
    free(wlc);
#endif

#ifdef DIQUARK
    free(dqc);
#endif

#ifdef PLOOPCOR_MEAS
    free(plc);
    free(plcr);
    free(plci);
#endif

    free(buffer);
    free(swline);
    free(owline);

#ifdef HASHLOOP
    free(hash); hash=NULL;
#endif

#ifdef DEBUG
    node0_printf("####################################\n");
    node0_printf("Number of gathers = %d(%d stale)\n",ngather,nstale);
    node0_printf("####################################\n");
#endif
} /* new_hvy_pot() */

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
  ( mlat[xc[ZUP]]*llat[xc[ZUP]] < maxc[ZUP] )
#define FLUTEYSHIFT(mlat,llat,max_x) \
  ( mlat[xc[YUP]]*llat[xc[YUP]] < maxc[YUP] )
#define FLUTEXSHIFT(mlat,llat,max_x) \
  ( mlat[xc[XUP]]*llat[xc[XUP]] < maxc[XUP] )

static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) {
  int mu;
  for ( this.insphere=OUT, this.gengather=SIMGATHER; this.insphere<IN && this.raise<TUP; ) {
    if ( llat[xc[ZUP]] < nc[ZUP] && FLUTEZSHIFT(mlat,llat,maxc) ) { /* gather next xc[ZUP] local lattice */
      this.raise = ZUP;
      this.gengather = SIMGATHER;
    } else { if ( llat[xc[YUP]] < nc[YUP] && FLUTEYSHIFT(mlat,llat,maxc) ) { /* gather next xc[YUP] local lattice */
      this.raise = YUP;
      nlat[xc[ZUP]] = 0;
      this.gengather = GENGATHER;
    } else { if ( llat[xc[XUP]] < nc[XUP] && FLUTEXSHIFT(mlat,llat,maxc) ) { /* gather next xc[XUP] local lattice */
      this.raise = XUP;
      nlat[xc[ZUP]] = 0;
      nlat[xc[YUP]] = 0;
      this.gengather=GENGATHER;
    } else {
#ifdef PLOOPCOR_MEAS
      if ( mlat[xc[TUP]] < nc[TUP] ) /* gather next wilson line segment from xc[TUP] cor-dir */
#else
      if ( mlat[xc[TUP]] < maxc[TUP] ) /* gather next wilson line segment from xc[TUP] cor-dir */
#endif
      {
        this.raise=TUP;
        nlat[xc[ZUP]] = 0;
        nlat[xc[YUP]] = 0;
        nlat[xc[XUP]] = 0;
        this.gengather = SIMGATHER;
      }
    }}}
    if ( this.raise > NODIR && this.raise < TUP ) { nlat[xc[this.raise]]++; }
    for ( mu=XUP; mu<TUP; mu++ ) { base[xc[mu]] = nlat[xc[mu]]*llat[xc[mu]]; }
    this.insphere=INSPHERE(base,llat,max_r2); /* within sphere? */
    while ( this.insphere == OUT && this.raise > NODIR && this.raise < TUP ) {
      nlat[xc[this.raise]] = 0;
      base[xc[this.raise]] = 0;
      if ( this.raise > XUP && this.raise < TUP ) {
        nlat[xc[--this.raise]]++;
        base[xc[this.raise]] += llat[xc[this.raise]];
      } else {
        this.raise = TUP;
        nlat[xc[ZUP]] = 0;
        nlat[xc[YUP]] = 0;
        nlat[xc[XUP]] = 0;
        this.gengather = SIMGATHER;
        break;
      }
      this.insphere = INSPHERE(base,llat,max_r2); /* within sphere? */
      this.gengather = GENGATHER;
    }
  }
  return (this);
}
#endif

#ifdef SNAKE
#define SNAKEZSHIFT(mlat,llat,maxc) \
  ( ( !((mlat[xc[XUP]]+mlat[xc[YUP]])%2) && mlat[xc[ZUP]]*llat[xc[ZUP]] < maxc[ZUP] ) \
   ||( ((mlat[xc[XUP]]+mlat[xc[YUP]])%2) && mlat[xc[ZUP]]*llat[xc[ZUP]] > 0 ) )
#define SNAKEYSHIFT(mlat,llat,maxc) \
  ( ( !(mlat[xc[XUP]]%2) && mlat[xc[YUP]]*llat[xc[YUP]] < maxc[YUP] ) \
   ||( (mlat[xc[XUP]]%2) && mlat[xc[YUP]]*llat[xc[YUP]] > 0 ) )
#define SNAKEXSHIFT(mlat,llat,maxc) \
  ( mlat[xc[XUP]]*llat[xc[XUP]] < maxc[XUP] )

static next_gather find_gather( int mlat[], int nlat[], int base[], next_gather this ) {
  int mu;
  for ( this.insphere=OUT, this.gengather=SIMGATHER; this.insphere<IN && this.raise < TUP; ) {
    if ( llat[xc[ZUP]] < nc[ZUP] && SNAKEZSHIFT(mlat,llat,maxc) ) { /* gather next z local lattice */
      if ((mlat[xc[XUP]]+mlat[xc[YUP]])%2 ) { this.raise = ZDOWN; /* shift from xc[ZDOWN] dir */
      } else { this.raise=ZUP; /* shift from xc[ZUP] dir */ }
      if ( mlat[xc[ZUP]] == 0 ) { this.raisedir[ZUP] = ZUP;
      } else { if ( !( this.lastraise == ZUP || this.lastraise == ZDOWN ) ) { this.raisedir[ZUP] = OPP_DIR(this.raisedir[ZUP]); }; };
      this.raise = this.raisedir[ZUP];
    } else {
      if ( llat[xc[YUP]] < nc[YUP] && SNAKEYSHIFT(mlat,llat,maxc) ) { /* gather next y local lattice */
        if ( mlat[xc[XUP]]%2 ) { this.raise=YDOWN; /* shift from xc[YDOWN] dir */
        } else { this.raise=YUP; /* shift from ixc[YUP] dir */ }
        if ( mlat[xc[YUP]] == 0 ) { this.raisedir[YUP]=YUP; };
        this.raise = this.raisedir[YUP];
      } else {
        if ( llat[xc[XUP]] < nc[XUP] && SNAKEXSHIFT(mlat,llat,maxc) ) { /* gather next x local lattice from xc[XUP] dir */
          this.raise = XUP; /* the value of this.raisedir[XUP] never matters */
        } else {
#ifdef PLOOPCOR_MEAS
          if ( mlat[xc[TUP]] < nc[TUP] ) /* gather next wilson line segment from xc[TUP] cor-dir */
#else
          if ( mlat[xc[TUP]] < maxc[TUP] ) /* gather next wilson line segment from xc[TUP] cor-dir */
#endif
          {
            this.raise=TUP;
            nlat[xc[ZUP]]=0;
            nlat[xc[YUP]]=0;
            nlat[xc[XUP]]=0;
          }
        }
      }
    }
    if ( SPATIAL(this.raise) ) {
      if ( this.raise<TUP ) { nlat[xc[this.raise]]++;
      } else  { nlat[xc[OPP_DIR(this.raise)]]--; }
    }
    for ( mu=XUP; mu<TUP; mu++ ) { base[xc[mu]]=nlat[xc[mu]]*llat[xc[mu]]; }
    this.insphere=INSPHERE(base,llat,max_r2); /* within sphere? */
    while ( this.insphere == OUT && SPATIAL(this.raise) ) {
      if ( this.raise < TUP ) {
        if ( this.raise == XUP ) {
          if ( nlat[xc[ZUP]] > 0 ) {
            nlat[xc[ZUP]]--;
            base[xc[ZUP]] -= llat[xc[ZUP]];
          } else {
            this.raisedir[ZUP]=ZDOWN;
            if ( nlat[xc[YUP]] > 0 ) {
              nlat[xc[YUP]]--;
              base[xc[YUP]] -= llat[xc[YUP]];
            } else {
              this.raisedir[YUP]=YDOWN;
              if ( mlat[xc[TUP]] < maxc[TUP] ) {
                this.raise=TUP;
                nlat[xc[ZUP]] = 0;
                nlat[xc[YUP]] = 0;
                nlat[xc[XUP]] = 0;
                this.gengather = SIMGATHER;
                break;
              } else this.raise=NODIR;
            }
          }
        } else {
          if ( nlat[xc[this.raise]] > 0 ) {
            nlat[xc[this.raise]]--;
            base[xc[this.raise]] -= llat[xc[this.raise]];
          }
        }
        if ( this.raise > XUP ) {
          this.raise = this.raisedir[this.raise-1];
          if ( this.raise < TUP ) {
            nlat[xc[this.raise]]++;
            base[xc[this.raise]] += llat[xc[this.raise]];
          } else {
            if ( this.raise < XDOWN ) {
              if ( nlat[xc[OPP_DIR(this.raise)]] > 0 ) {
                nlat[xc[OPP_DIR(this.raise)]]--;
                base[xc[OPP_DIR(this.raise)]] -= llat[xc[OPP_DIR(this.raise)]];
              }
            } else {
              if ( mlat[xc[TUP]] < maxc[TUP] ) {
                this.raise=TUP;
                nlat[xc[ZUP]] = 0;
                nlat[xc[YUP]] = 0;
                nlat[xc[XUP]] = 0;
                this.gengather = SIMGATHER;
                break;
              } else this.raise = NODIR;
            }
          }
        }
      }
      for ( mu=XUP; mu<TUP; mu++) { base[xc[mu]] = nlat[xc[mu]]*llat[xc[mu]]; }
      this.insphere = INSPHERE(base,llat,max_r2); /* within sphere? */
      this.gengather = GENGATHER;
    }
    switch(this.raise){
      case(XUP):{
        this.raisedir[YUP] = OPP_DIR(this.raisedir[YUP]);
        break;
      }
      case(TUP):{
        nlat[xc[ZUP]] = 0;
        nlat[xc[YUP]] = 0;
        nlat[xc[XUP]] = 0;
        break;
      }
    }
    if ( this.raise == NODIR ) { break; }
  }
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
