/****** gauge_stuff.c  -- ******************/
/* MIMD version 6 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file */

/**#define GFTIME**/ /* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */

#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        /* define it to be nothing */
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void printpath( int *path, int length );

#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART1

char gauge_action_description[128];
int  gauge_action_nloops=NLOOP;
int  gauge_action_nreps=NREPS;
int loop_length[NLOOP];	/* lengths of various kinds of loops */
int loop_num[NLOOP];	/* number of rotations/reflections  for each kind */

    /* table of directions, 1 for each kind of loop */
int loop_ind[NLOOP][MAX_LENGTH];
    /* table of directions, for each rotation and reflection of each kind of
	loop.  tabulated with "canonical" starting point and direction. */
int loop_table[NLOOP][MAX_NUM][MAX_LENGTH];
    /* table of coefficients in action, for various "representations" (actually,
	powers of the trace) */
Real loop_coeff[NLOOP][NREPS];
    /* for each rotation/reflection, an integer distinct for each starting
	point, or each cyclic permutation of the links */
int loop_char[MAX_NUM];
    /* for each kind of loop for each rotation/reflection, the expectation
	value of the loop */
double loop_expect[NLOOP][NREPS][MAX_NUM];


/* Make table of loops in action */
void make_loop_table() {

    int perm[8],pp[8],ir[4];
    int length,iloop,i,j,chr;
    int vec[MAX_LENGTH];
    int count,flag;
    void char_num( int *dig, int *chr, int length);

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
    node0_printf("loop coefficients: nloop rep loop_coeff  multiplicity\n");
    for(i=0;i<NREPS;i++) for(j=0;j<NLOOP;j++) {
	node0_printf("                    %d %d      %e     %d\n",
	    j,i,loop_coeff[j][i],loop_num[j]);
    }

} /* make_loop_table */


/* find a number uniquely identifying the cyclic permutation of a path,
   or the starting point on the path.  Backwards paths are considered
   equivalent here, so scan those too. */
void char_num( int *dig, int *chr, int length){
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

double imp_gauge_action() {
    register int i;
    int rep;
    register site *s;
    complex trace;
    double g_action;
    double action,act2,total_action;
    int length;

    /* these are for loop_table  */
    int ln,iloop;

    g_action=0.0;

    /* gauge action */
    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product( loop_table[iloop][ln] , length );

	    FORALLSITES(i,s){
		trace=trace_su3( &s->tempmat1 );
		action =  3.0 - (double)trace.real;
		/* need the "3 -" for higher characters */
        	total_action= (double)loop_coeff[iloop][0]*action;
        	act2=action;
		for(rep=1;rep<NREPS;rep++){
		    act2 *= action;
		    total_action += (double)loop_coeff[iloop][rep]*act2;
		}

        	g_action  += total_action;

	    } END_LOOP /* sites */
	} /* ln */
    } /* iloop */

    g_doublesum( &g_action );
    return( g_action );
} /* imp_gauge_action */


/* update the momenta with the gauge force */
void imp_gauge_force( Real eps, field_offset mom_off ){
    register int i,dir;
    register site *st;
    su3_matrix tmat1,tmat2;
    register Real eb3;
    register anti_hermitmat* momentum;
    int nflop = 153004;  /* For Symanzik1 action */
#ifdef GFTIME
    double dtime;
#endif

    int j,k;
    int dirs[MAX_LENGTH],length;
    int path_dir[MAX_LENGTH],path_length;

    int ln,iloop;
    Real action,act2,new_term;

    int ncount;

#ifdef GFTIME
dtime=-dclock();
#endif
    eb3 = eps*beta/3.0;

    /* Loop over directions, update mom[dir] */
    for(dir=XUP; dir<=TUP; dir++){

	FORALLSITES(i,st)for(j=0;j<3;j++)for(k=0;k<3;k++){
			st->staple.e[j][k]=cmplx(0.0,0.0);
	} END_LOOP

	ncount=0;
	for(iloop=0;iloop<NLOOP;iloop++){
	    length=loop_length[iloop];
	    for(ln=0;ln<loop_num[iloop];ln++){
/**printf("UPD:  "); printpath( loop_table[iloop][ln], length );**/
		/* set up dirs.  we are looking at loop starting in "XUP"
		   direction, rotate so it starts in "dir" direction. */
		for(k=0;k<length;k++){
                    if( GOES_FORWARDS(loop_table[iloop][ln][k]) ){
                	dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
		    }
            	    else {
                        dirs[k]=OPP_DIR(
			    (dir+OPP_DIR(loop_table[iloop][ln][k]))%4 );
		    }
		}

		path_length= length-1;  /* generalized "staple" */

		/* check for links in direction of momentum to be
		   updated, each such link gives a contribution. Note
		   the direction of the path - opposite the link. */
		for(k=0;k<length;k++)if( dirs[k]==dir||dirs[k]==OPP_DIR(dir)) {
		    if( GOES_FORWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
			path_dir[j] = dirs[(k+j+1)%length];
		    }
		    if( GOES_BACKWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
			path_dir[path_length-1-j] =
			    OPP_DIR(dirs[(k+j+1)%length]);
		    }
/**if(dir==XUP)printf("X_UPDATE PATH: "); printpath( path_dir, path_length );**/
		    path_product(path_dir,path_length);

		    /* We took the path in the other direction from our
			old convention in order to get it to end up
			"at our site", so now take adjoint */
		    /* then compute "single_action" contribution to
			staple */
		    FORALLSITES(i,st){
			su3_adjoint( &(st->tempmat1), &tmat1 );
			/* first we compute the fundamental term */
			new_term = loop_coeff[iloop][0];

			/* now we add in the higher representations */
			if(NREPS > 1){
node0_printf("WARNING: THIS CODE IS NOT TESTED\n"); exit(0);
			    act2=1.0;
			    action = 3.0 - realtrace_su3(&(st->link[dir]),
				&tmat1 ); 

			    for(j=1;j<NREPS;j++){
				act2 *= action;
				new_term +=
				    loop_coeff[iloop][j]*act2*(Real)(j+1);
			    }
			}  /* end if NREPS > 1 */

			scalar_mult_add_su3_matrix( &(st->staple), &tmat1,
				new_term, &(st->staple) );

		    } END_LOOP

		    ncount++;

		} /* k (location in path) */
	    } /* ln */
	} /* iloop */

	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	    mult_su3_na( &(st->link[dir]), &(st->staple), &tmat1 );
	    momentum = (anti_hermitmat *)F_PT(st,mom_off);
	    uncompress_anti_hermitian( &momentum[dir], &tmat2 );
	    scalar_mult_sub_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &momentum[dir] );
	} END_LOOP
    } /* dir loop */
#ifdef GFTIME
dtime+=dclock();
node0_printf("GFTIME:   time = %e (Symanzik1) mflops = %e\n",dtime,
	     nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif
} /* imp_gauge_force.c */

/* Measure gauge observables:
    Loops in action (time and space directions treated differently)
    Polyakov loop

*/
void g_measure( ){
    double ss_plaquette, st_plaquette;
    complex p_loop;
    register int i;
    register site *s;
    complex trace;
    double average[NREPS],action,act2,total_action;
    int length;
    /* these are for loop_table  */
    int ln,iloop,rep;

    /* KS and BC minus signs should be out for this routine */
    d_plaquette( &ss_plaquette, &st_plaquette );
    if(this_node==0)printf("PLAQ:\t%f\t%f\n", ss_plaquette, st_plaquette );

    p_loop = ploop();
    if(this_node==0)printf("P_LOOP:\t%e\t%e\n", p_loop.real, p_loop.imag );

    /* gauge action, all loops that contribute */
    total_action=0.0;
    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product( loop_table[iloop][ln] , length );

	    for(rep=0;rep<NREPS;rep++)average[rep] = 0.0;
	    FORALLSITES(i,s){
		trace=trace_su3( &s->tempmat1 );
		average[0] += (double)trace.real;
		action =  3.0 - (double)trace.real;
		total_action += (double)loop_coeff[iloop][0]*action;
		/* need the "3 -" for higher characters */
        	act2=action;
		for(rep=1;rep<NREPS;rep++){
		    act2 *= action;
		    average[rep] += act2;
		    total_action += (double)loop_coeff[iloop][rep]*act2;
		} /* reps */
	    } END_LOOP /* sites */
	    g_vecdoublesum( average, NREPS );
	    /* dump the loop */
	    node0_printf("G_LOOP:  %d  %d  %d   ",iloop,ln,length);
	    for(rep=0;rep<NREPS;rep++)node0_printf("\t%e",average[rep]/volume);
	    node0_printf("\t( ");
	    for(i=0;i<length;i++)node0_printf("%d ",loop_table[iloop][ln][i]);
	    node0_printf(" )\n");
	} /* ln */
    } /* iloop */
    g_doublesum( &total_action );
    node0_printf("GACTION: %e\n",total_action/volume);
    /**node0_printf("CHECK:   %e   %e\n",total_action,imp_gauge_action() );**/

    if(this_node==0)fflush(stdout);

}

void printpath( int *path, int length ){
    register int i;
    node0_printf("\t( ");
    for(i=0;i<length;i++)node0_printf("%d ",path[i]);
    node0_printf(",  L = %d )\n", length );
}

#ifdef N_SUBL32
/*** code from symanzik_sl32/dsdu_qhb.c  -- compute the staple ***/
/* This is a version for extended actions where 32 sublattices are
   needed to make the links independent. */
/* U.M. Heller August 1997 */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */
#include <assert.h>
void dsdu_qhb_subl(int dir, int subl)
{
register site *st;
register int i;
int iloop, ln, k, j;
int dirs[MAX_LENGTH], length;
int path_dir[MAX_LENGTH], path_length;
su3_matrix tmat1;
int fsubl;

 assert(NREPS==1);   /* This procedure designed only for NREPS = 1 */

    FORSOMESUBLATTICE(i,st,subl) {
	clear_su3mat(&(st->staple));
    }

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	for(ln=0;ln<loop_num[iloop];ln++){
	    /* set up dirs.  we are looking at loop starting in "XUP"
	       direction, rotate so it starts in "dir" direction. */
	    for(k=0;k<length;k++){
		if( GOES_FORWARDS(loop_table[iloop][ln][k]) ){
		    dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
		}
		else {
		    dirs[k]=OPP_DIR(
			(dir+OPP_DIR(loop_table[iloop][ln][k]))%4 );
		}
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
		path_prod_subl(path_dir, path_length, fsubl);

		/* We took the path in the other direction from our old
		   convention in order to get it to end up "at our site".
		   So now take adjoint */
		FORSOMESUBLATTICE(i,st,subl) {
		    su3_adjoint(&(st->tempmat1), &tmat1 );
		    scalar_mult_add_su3_matrix(&(st->staple), &tmat1,
			loop_coeff[iloop][0], &(st->staple) );
		}
	    } /* k (location in path) */
	} /* ln */
    } /* iloop */

    g_sync();

} /* dsdu_qhb */

#endif /* N_SUBL32 */
