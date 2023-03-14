/****** gauge_stuff.c  -- ******************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file */

/**#define GFTIME**/ /* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"
#if defined (HAVE_QUDA)
#include "../include/generic_quda.h"
#endif

/* I don't understand the advantage of following SG 5/25/17 will comment out */
/*
#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        // define it to be nothing
#endif
*/
/* I don't understand the advantage of preceeding SG 5/25/17 will comment out */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void printpath( int *path, int length );

#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART1

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

char gauge_action_description[128];
int  gauge_action_nloops=NLOOP;
int  gauge_action_nreps=NREPS;
static int loop_length[NLOOP];	/* lengths of various kinds of loops */
static int loop_num[NLOOP];	/* number of rotations/reflections  for each kind */

    /* table of directions, for each rotation and reflection of each kind of
	loop.  tabulated with "canonical" starting point and direction. */
static int*** loop_table;
    /* table of coefficients in action, for various "representations" (actually,
	powers of the trace) */
static Real **loop_coeff;
    /* for each rotation/reflection, an integer distinct for each starting
	point, or each cyclic permutation of the links */
int loop_char[MAX_NUM];

#ifdef ANISOTROPY
    /* for each rotation/reflection, an integer indicating if the path
       is spatial (=0) or temporal (=1) */
int **loop_st;
#endif

static void char_num( int *dig, int *chr, int length);

/* Make table of loops in action */
void make_loop_table() {

    int perm[8],pp[8],ir[4];
    int length,iloop,i,j,chr;
    int vec[MAX_LENGTH];
    int count,flag;
    int total_dyn_flavors;
    char myname[] = "make_loop_table";

    total_dyn_flavors = 0;
    for(i = 0; i < n_dyn_masses; i++){
      total_dyn_flavors += dyn_flavors[i];
    }

    /* Allocate as if loop_table[NLOOP][MAX_NUM][MAX_LENGTH] */

    loop_table = (int ***)malloc(sizeof(int **)*NLOOP);
    if(loop_table == NULL){
      printf("%s(%d): No room for loop_table\n",myname,this_node);
      terminate(1);
    }

    for(iloop = 0; iloop < NLOOP; iloop++){
      loop_table[iloop] = (int **)malloc(sizeof(int *)*MAX_NUM);
      if(loop_table[iloop] == NULL){
	printf("%s(%d): No room for loop_table\n",myname,this_node);
	terminate(1);
      }

      for(count = 0; count < MAX_NUM; count++){
	loop_table[iloop][count] = (int *)malloc(sizeof(int)*MAX_LENGTH);
	if(loop_table[iloop][count] == NULL){
	  printf("%s(%d): No room for loop_table\n",myname,this_node);
	  terminate(1);
	}
      }
    }

    /* Allocate as if loop_coeff[NLOOP][NREPS] */

    loop_coeff = (Real **)malloc(sizeof(Real *)*NLOOP);
    if(loop_coeff == NULL){
      printf("%s(%d): No room for loop_coeff\n",myname,this_node);
      terminate(1);
    }

    for(iloop = 0; iloop < NLOOP; iloop++){
      loop_coeff[iloop] = (Real *)malloc(sizeof(Real)*NREPS);
      if(loop_coeff[iloop] == NULL){
	printf("%s(%d): No room for loop_coeff\n",myname,this_node);
	terminate(1);
      }
    }

#define GAUGE_ACTION_PART2
/* defines all loops and their coefficients */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART2

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	count=0;
	/* permutations */
	for(perm[0]=0;perm[0]<4;perm[0]++)
	for(perm[1]=0;perm[1]<4;perm[1]++)
	for(perm[2]=0;perm[2]<4;perm[2]++)
	for(perm[3]=0;perm[3]<4;perm[3]++){
	    if(perm[0] != perm[1] && perm[0] != perm[2] 
		&& perm[0] != perm[3] && perm[1] != perm[2]
	    	&& perm[1] != perm[3] && perm[2] != perm[3] ) {
	        /* reflections*/
	 	for(ir[0]=0;ir[0]<2;ir[0]++)
		for(ir[1]=0;ir[1]<2;ir[1]++)
		for(ir[2]=0;ir[2]<2;ir[2]++)
		for(ir[3]=0;ir[3]<2;ir[3]++){
		    for(j=0;j<4;j++){
			pp[j]=perm[j];

			if(ir[j] == 1) pp[j]=7-pp[j];
			pp[7-j]=7-pp[j];
		    }
		    /* create new vector*/
		    for(j=0;j<length;j++) vec[j]=pp[loop_ind[iloop][j]];

         	    char_num(vec,&chr,length);
         	    flag=0;
		    /* check if it's a new set: */
		    for(j=0;j<count;j++) if(chr == loop_char[j])flag=1;
		    if(flag == 0 ){
			loop_char[count]=chr;
			for(j=0;j<length;j++)
			    loop_table[iloop][count][j]=vec[j];
			count++;
/**node0_printf("ADD LOOP: "); printpath( vec, length );**/
		    }
		    if(count>MAX_NUM){
			node0_printf("OOPS: MAX_NUM too small\n");
			exit(0);
		    }
		    loop_num[iloop]=count;

		} /* end reflection*/
	    } /* end permutation if block */
	} /* end permutation */
    } /* end iloop */

    /* print out the loop coefficients */
    node0_printf("gauge_action: total_dyn_flavors = %d\n",total_dyn_flavors);
    node0_printf("loop coefficients: nloop rep loop_coeff  multiplicity\n");
    for(i=0;i<NREPS;i++) for(j=0;j<NLOOP;j++) {
	node0_printf("                    %d %d      %e     %d\n",
	    j,i,loop_coeff[j][i],loop_num[j]);
    }

} /* make_loop_table */

int get_max_length(){
  return MAX_LENGTH;
}

int get_nloop(){
  return NLOOP;
}

int get_nreps(){
  return NREPS;
}

int *get_loop_length(){
  return loop_length;
}

int *get_loop_num(){
  return loop_num;
}

int ***get_loop_table(){
  return loop_table;
}

Real **get_loop_coeff(){
  return loop_coeff;
}

/* find a number uniquely identifying the cyclic permutation of a path,
   or the starting point on the path.  Backwards paths are considered
   equivalent here, so scan those too. */
static void char_num( int *dig, int *chr, int length){
    int j;
    int bdig[MAX_LENGTH],tenl,newv,old;
    /* "dig" is array of directions.  "bdig" is array of directions for
	backwards path. */

  tenl=1;
  for(j=0;j<length-1;j++) tenl=tenl*10;

  *chr=dig[length-1];
  for(j=length-2;j>=0;j--) *chr= *chr*10+dig[j];

  /* forward*/
  old=*chr;
  for(j=length-1;j>=1;j--){
       newv=old-tenl*dig[j];
       newv=newv*10+dig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

   /* backward*/
   for(j=0;j<length;j++)bdig[j]=7-dig[length-j-1];
   old=bdig[length-1];
   for(j=length-2;j>=0;j--) old=old*10+bdig[j];
   if(old < *chr ) *chr=old;
   for(j=length-1;j>=1;j--){
       newv=old-tenl*bdig[j];
       newv=newv*10+bdig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

} /* char_num */

/* Measure gauge observables:
    Loops in action (time and space directions treated differently)
    Polyakov loop
*/

#if defined (HAVE_QUDA) && defined(USE_GA_GPU) && !defined(ANISOTROPY) && !defined(BPCORR) && NREPS == 1
void g_measure_gpu( ) {
    complex p_loop;
    register int i;
    register site *s;
    double ss_plaquette, st_plaquette;
    complex trace;
    double average,action,act2,total_action;
    double this_total_action; /* need for loop over sitest */
    int length;
    su3_matrix *tempmat1;
    su3_matrix *links;
    /* these are for loop_table  */
    int ln,iloop,rep;

    /* Count total number of loops */
    int num_paths = 0;
    for (iloop = 0; iloop < NLOOP; iloop++)
        for (ln = 0; ln < loop_num[iloop]; ln++)
            num_paths++;

    /* Max length */
    int max_length = get_max_length();

    /* Storage for traces */
    double *traces = (double*)malloc(2 * num_paths * sizeof(double));

    /* Storage for input paths */
    int **input_path_buf = (int**)malloc(num_paths * sizeof(int*));
    for (i = 0; i < num_paths; i++)
        input_path_buf[i] = (int*)malloc(max_length * sizeof(int));

    /* Storage for path lengths */
    int *path_length = (int*)malloc(num_paths * sizeof(int));

    /* Storage for loop coefficients */
    double *loop_coeff = (double*)malloc(num_paths * sizeof(double));

    /* Overall scaling factor */
    double factor = 1. / volume;

    num_paths = 0;
    for (iloop = 0; iloop < NLOOP; iloop++) {
        length = loop_length[iloop];
        for (ln = 0; ln < loop_num[iloop]; ln++) {
            path_length[num_paths] = length; /* path length */
            loop_coeff[num_paths] = 1.0; /* due to the "3. - [...]" convention below, we'll wait to scale then */
            for (i = 0; i < length; i++)
                input_path_buf[num_paths][i] = loop_table[iloop][ln][i];
            num_paths++;
        }
    }

    Real **loop_coeff_milc = get_loop_coeff();
    double plaq_array[3];
    double ploop_array[2];

    initialize_quda();

    QudaMILCSiteArg_t arg = newQudaMILCSiteArg();

    /* Fused kernel that computes the plaquette, temporal Polyakov loop, and gauge loop traces */
    qudaGaugeMeasurementsPhased(MILC_PRECISION, plaq_array, ploop_array, 3, traces, input_path_buf, path_length,
                                loop_coeff, num_paths, max_length, factor, &arg, phases_in);


    ss_plaquette = 3.0 * plaq_array[1];
    st_plaquette = 3.0 * plaq_array[2];

#if (MILC_PRECISION==1)
    node0_printf("PLAQ:\t%f\t%f\n", ss_plaquette, st_plaquette );
#else
    node0_printf("PLAQ:\t%.16f\t%.16f\n", ss_plaquette, st_plaquette );
#endif

    node0_printf("P_LOOP:\t%e\t%e\n", ploop_array[0], ploop_array[1] );

    /* Accumulate the actions out of the gauge loop traces */
    num_paths = 0;
    total_action = 0.0;
    for (iloop = 0; iloop < NLOOP; iloop++) {
        int length = loop_length[iloop];
        /* loop over rotations and reflections */
        for (ln = 0; ln < loop_num[iloop]; ln++) {
            this_total_action = 0.;
            average = traces[2 * num_paths]; // extract real part
            action = 3.0 - traces[2 * num_paths];
            this_total_action = (double)loop_coeff_milc[iloop][0] * action;
            total_action += this_total_action;
            
            /* dump the loop */
            node0_printf("G_LOOP:  %d  %d  %d   ", iloop, ln, length);
#if (MILC_PRECISION==1)
            node0_printf("\t%e", average);
#else
            node0_printf("\t%.16e", average);
#endif
            node0_printf("\t( ");
            for (i = 0; i < length; i++) node0_printf("%d ", loop_table[iloop][ln][i]);
            node0_printf(" )\n");

            num_paths++;
        } /* ln */
    } /* iloop */

    node0_printf("GACTION: %e\n", total_action);
    /**node0_printf("CHECK:   %e   %e\n",total_action,imp_gauge_action_gpu() );**/

    if(this_node==0)fflush(stdout);

    free(loop_coeff);
    free(path_length);
    for (i = 0; i < num_paths; i++)
        free(input_path_buf[i]);
    free(input_path_buf);
    free(traces);

} /* g_measure_gpu */
#endif

void g_measure( ){
    double ss_plaquette, st_plaquette;
    complex p_loop;
    register int i;
    register site *s;
    complex trace;
    double average[NREPS],action,act2,total_action;
    double this_total_action; /* need for loop over sitest */
    int length;
    su3_matrix *tempmat1;
    su3_matrix *links;
    /* these are for loop_table  */
    int ln,iloop,rep;

    tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
      printf("g_measure: Can't malloc temporary\n");
      terminate(1);
    }
    links = create_G_from_site();

    /* KS and BC minus signs should be out for this routine */
    d_plaquette( &ss_plaquette, &st_plaquette );
#if (MILC_PRECISION==1)
    if(this_node==0)printf("PLAQ:\t%f\t%f\n", ss_plaquette, st_plaquette );
#else
    if(this_node==0)printf("PLAQ:\t%.16f\t%.16f\n", ss_plaquette, st_plaquette );
#endif
    p_loop = ploop();
    if(this_node==0)printf("P_LOOP:\t%e\t%e\n", p_loop.real, p_loop.imag );

    /* gauge action, all loops that contribute */
    total_action=0.0;
    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product_fields(links, loop_table[iloop][ln] , length, tempmat1 );

	    for(rep=0;rep<NREPS;rep++)average[rep] = 0.0;
    	    this_total_action=0.;
	    //	    FORALLSITES_OMP(i,s,private(trace,action,act2,rep) reduction(+:this_total_action,average)){
	    FORALLFIELDSITES(i){
		trace=trace_su3( &tempmat1[i] );
		average[0] += (double)trace.real;
		action =  3.0 - (double)trace.real;
#ifndef ANISOTROPY
		this_total_action += (double)loop_coeff[iloop][0]*action;
		/* need the "3 -" for higher characters */
#else
		/* NOTE: in the total action calculation
                   for anisotropic case every loop is multiplied by
		   the corresponding spatial (beta[0]) or temporal (beta[1])
		   coupling, while in the isotropic case all loops are
		   added together and are NOT multiplied by beta in
		   this function */
		this_total_action += (double)loop_coeff[iloop][0]*action
		 		*beta[loop_st[iloop][ln]];
		/* loop_st[iloop][ln] is either 0 or 1 */
#endif
        	act2=action;
		for(rep=1;rep<NREPS;rep++){
		    act2 *= action;
		    average[rep] += act2;
		    this_total_action += (double)loop_coeff[iloop][rep]*act2;
		} /* reps */
		//	    } END_LOOP_OMP; /* sites */
	    } /* sites */
	    g_vecdoublesum( average, NREPS );
	    total_action += this_total_action;
	    /* dump the loop */
	    node0_printf("G_LOOP:  %d  %d  %d   ",iloop,ln,length);
#if (MILC_PRECISION==1)
	    for(rep=0;rep<NREPS;rep++)node0_printf("\t%e",average[rep]/volume);
#else
	    for(rep=0;rep<NREPS;rep++)node0_printf("\t%.16e",average[rep]/volume);
#endif
	    node0_printf("\t( ");
	    for(i=0;i<length;i++)node0_printf("%d ",loop_table[iloop][ln][i]);
	    node0_printf(" )\n");
	} /* ln */
    } /* iloop */
    g_doublesum( &total_action );
    node0_printf("GACTION: %e\n",total_action/volume);
    /**node0_printf("CHECK:   %e   %e\n",total_action,imp_gauge_action() );**/

    if(this_node==0)fflush(stdout);
    destroy_G(links);
    special_free(tempmat1);
}

void printpath( int *path, int length ){
    register int i;
    node0_printf("\t( ");
    for(i=0;i<length;i++)node0_printf("%d ",path[i]);
    node0_printf(",  L = %d )\n", length );
}

#ifdef ANISOTROPY
/* Auxilliary function that goes through all possible paths rotations
   and reflections and records if the path is spatial or temporal.
   The results are stored in loop_st[NLOOP][MAX_NUM] array */
void path_determine_st() {

    int iloop, count, i;
    char myname[] = "path_determine_st";

    /* Allocate as if loop_st[NLOOP][MAX_NUM] */

    loop_st = (int **)malloc(sizeof(int *)*NLOOP);
    if(loop_st == NULL){
      printf("%s(%d): No room for loop_st\n",myname,this_node);
      terminate(1);
    }

    for(iloop = 0; iloop < NLOOP; iloop++){
      loop_st[iloop] = (int *)malloc(sizeof(int)*MAX_NUM);
      if(loop_st[iloop] == NULL){
        printf("%s(%d): No room for loop_st\n",myname,this_node);
        terminate(1);
      }

      /* loop over rotations/reflections */
      for(count = 0; count < loop_num[iloop]; count++){

        /* set initially as a spatial path */
        loop_st[iloop][count] = 0;

        /* loop over directions in the path */
        for(i = 0; i < loop_length[iloop]; i++){
          if( loop_table[iloop][count][i]==TUP ||
              loop_table[iloop][count][i]==TDOWN ) {
            loop_st[iloop][count] = 1;
            break;
          }
        }
      }
    }
}
#endif

#ifdef N_SUBL32
/*** code from symanzik_sl32/dsdu_qhb.c  -- compute the staple ***/
/* This is a version for extended actions where 32 sublattices are
   needed to make the links independent. */
/* U.M. Heller August 1997 */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */
#include <assert.h>

/* NOTE: the staple is returned in the site strcture for use in monte
   and relax in symanzik_sl32! */

void dsdu_qhb_subl(int dir, int subl)
{
register site *st;
register int i;
int iloop, ln, k, j;
int dirs[MAX_LENGTH], length;
int path_dir[MAX_LENGTH], path_length;
su3_matrix tmat1, *tempmat1;
int fsubl;
#ifdef ANISOTROPY
int is_temporal; /* to decide what kind of staple we have:
                    0 - spatial, 1 - temporal */
#endif

 assert(NREPS==1);   /* This procedure designed only for NREPS = 1 */

 tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
 if(tempmat1 == NULL){
   printf("dsdu_qhb_subl: Can't malloc temporary\n");
   terminate(1);
 }

    FORSOMESUBLATTICE_OMP(i,st,subl,default(shared)) {
	clear_su3mat(&(st->staple));
#ifdef ANISOTROPY
	clear_su3mat(&(st->staple_a[0]));
	clear_su3mat(&(st->staple_a[1]));
#endif
    } END_LOOP_OMP;

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	for(ln=0;ln<loop_num[iloop];ln++){
	    /* set up dirs.  we are looking at loop starting in "XUP"
	       direction, rotate so it starts in "dir" direction. */
#ifdef ANISOTROPY
            /* initialize staple flag as spatial */
            is_temporal = 0;
#endif
	    for(k=0;k<length;k++){
		if( GOES_FORWARDS(loop_table[iloop][ln][k]) ){
		    dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
		}
		else {
		    dirs[k]=OPP_DIR(
			(dir+OPP_DIR(loop_table[iloop][ln][k]))%4 );
		}
#ifdef ANISOTROPY
		/* flip the flag if a temporal link is encountered */
		if( is_temporal==0 && ( dirs[k]==TUP || dirs[k]==TDOWN ) )
		    is_temporal=1;
#endif
	    }

	    path_length = length-1;	/* generalized "staple" */
	    /* The path starts at the forward end of the link */
	    fsubl = neighsubl[subl][dir];

	    /* check for links in direction of link to be updated.
	       Note the direction of the path - opposite the link. */
	    for(k=0;k<length;k++)if( dirs[k]==dir||dirs[k]==OPP_DIR(dir)) {
		if( GOES_FORWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[j] = dirs[(k+j+1)%length];
		}
		if( GOES_BACKWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[path_length-1-j] =
			OPP_DIR(dirs[(k+j+1)%length]);
		}
		path_prod_subl(path_dir, path_length, fsubl, tempmat1);

		/* We took the path in the other direction from our old
		   convention in order to get it to end up "at our site".
		   So now take adjoint */
		FORSOMESUBLATTICE_OMP(i,st,subl,private(tmat1)) {
		    su3_adjoint( &tempmat1[i], &tmat1 );
#ifndef ANISOTROPY
		    scalar_mult_add_su3_matrix(&(st->staple), &tmat1,
			loop_coeff[iloop][0], &(st->staple) );
#else
		    scalar_mult_add_su3_matrix(&(st->staple_a[is_temporal]),
			&tmat1, loop_coeff[iloop][0],
			&(st->staple_a[is_temporal]) );
#endif
		} END_LOOP_OMP;
	    } /* k (location in path) */
	} /* ln */
    } /* iloop */

#ifdef ANISOTROPY
    /* Add spatial and temporal staples weighted by betas to the
       "staple" variable */
    FORSOMESUBLATTICE(i,st,subl) {
	scalar_mult_add_su3_matrix(&(st->staple), &(st->staple_a[0]),
	    beta[0], &(st->staple) );
	scalar_mult_add_su3_matrix(&(st->staple), &(st->staple_a[1]),
	    beta[1], &(st->staple) );
    }
#endif
    special_free(tempmat1);
    g_sync();

} /* dsdu_qhb */

#endif /* N_SUBL32 */
