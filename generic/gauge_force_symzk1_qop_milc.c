/***************** gauge_force_symzk1_qop_milc.c  -- ***********************/
/* MIMD version 7 */
/* MILC implementation of QOP gauge force routine for testing the interface */
/* This routine is intended to be specific to the Symanzik 1 loop gauge action
 * but it uses the general improved action algorithm.
 *
 * T.D. and A.H. general gauge action updating code
 * D.T. modified  5/97
 * D.T. modified 12/97, optimized gauge_force a little
 * D.T. modified 3/99, gauge action in include file
 * C.D. split from gauge_stuff.c 10/06 */

#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/qop_milc.h"

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

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/* update the momenta with the gauge force */
void QOP_symanzik_1loop_gauge_force(QOP_info_t *info, QOP_GaugeField *gauge, 
		    QOP_Force *force, QOP_gauge_coeffs_t *coeffs, Real eps)
{
    register int i,dir;
    register site *st;
    su3_matrix tmat1;
    register Real eb3;    /* Note: eps now includes eps*beta */
    register su3_matrix* momentum;
    su3_matrix *staple, *tempmat1;

    /* lengths of various kinds of loops */
    int *loop_length = get_loop_length();
    /* number of rotations/reflections  for each kind */
    int *loop_num = get_loop_num();
    /* table of directions, 1 for each kind of loop */
    int ***loop_table = get_loop_table();
    /* table of coefficients in action, for various "representations"
	(actually, powers of the trace) */
    Real **loop_coeff = get_loop_coeff(); /* We make our own */
    int max_length = get_max_length(); /* For Symanzik 1 loop! */
    int nloop = get_nloop();
    int nreps = get_nreps();
    su3_matrix *forwardlink[4];
    su3_matrix *tmpmom[4];

    int nflop = 153004;  /* For Symanzik1 action */
    Real final_flop;
    double dtime;
    int j,k;
    int *dirs,length;
    int *path_dir,path_length;

    int ln,iloop;
    Real action,act2,new_term;

    int ncount;
    char myname[] = "imp_gauge_force";

    dtime=-dclock();

    info->status = QOP_FAIL;

    /* Parity requirements */
    if(gauge->evenodd != QOP_EVENODD ||
       force->evenodd != QOP_EVENODD
       )
      {
	printf("QOP_asqtad_force: Bad parity gauge %d force %d\n",
	       gauge->evenodd, force->evenodd);
	return;
      }

    /* Map field pointers to local static pointers */
    
    FORALLUPDIR(dir){
      forwardlink[dir] = gauge->g + dir*sites_on_node;
      tmpmom[dir]  = force->f + dir*sites_on_node;
    }
    /* Check loop coefficients */

    if(coeffs->plaquette != loop_coeff[0][0] ||
       coeffs->rectangle != loop_coeff[1][0] ||
       coeffs->parallelogram != loop_coeff[2][0])
      {
	printf("%s(%d): Path coeffs don't match\n",myname,this_node);
	return;
      }

    /* Allocate arrays according to action */
    dirs = (int *)malloc(max_length*sizeof(int));
    if(dirs == NULL){
      printf("%s(%d): Can't malloc dirs\n",myname,this_node);
      return;
    }

    path_dir = (int *)malloc(max_length*sizeof(int));
    if(path_dir == NULL){
      printf("%s(%d): Can't malloc path_dir\n",myname,this_node);
      return;
    }
    staple = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
    if(staple == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      return;
    }

    tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      return;
    }

    eb3 = eps/3.0;

    /* Loop over directions, update mom[dir] */
    for(dir=XUP; dir<=TUP; dir++){

	FORALLSITES(i,st)for(j=0;j<3;j++)for(k=0;k<3;k++){
			staple[i].e[j][k]=cmplx(0.0,0.0);
	} END_LOOP

	ncount=0;
	for(iloop=0;iloop<nloop;iloop++){
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
		    path_product(path_dir,path_length, tempmat1);

		    /* We took the path in the other direction from our
			old convention in order to get it to end up
			"at our site", so now take adjoint */
		    /* then compute "single_action" contribution to
			staple */
		    FORALLSITES(i,st){
			su3_adjoint( &(tempmat1[i]), &tmat1 );
			/* first we compute the fundamental term */
			new_term = loop_coeff[iloop][0];

			/* now we add in the higher representations */
			if(nreps > 1){
node0_printf("WARNING: THIS CODE IS NOT TESTED\n"); exit(0);
			    act2=1.0;
			    action = 3.0 - realtrace_su3(forwardlink[dir]+i,
			      &tmat1 ); 

			    for(j=1;j<nreps;j++){
				act2 *= action;
				new_term +=
				    loop_coeff[iloop][j]*act2*(Real)(j+1);
			    }
			}  /* end if nreps > 1 */

			scalar_mult_add_su3_matrix( &(staple[i]), &tmat1,
				new_term, &(staple[i]) );

		    } END_LOOP

		    ncount++;

		} /* k (location in path) */
	    } /* ln */
	} /* iloop */

	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	    mult_su3_na( forwardlink[dir]+i, &(staple[i]), &tmat1 );
	    momentum = tmpmom[dir] + i;
	    scalar_mult_sub_su3_matrix( momentum, &tmat1,
		eb3, momentum );
	} END_LOOP
    } /* dir loop */


  free(dirs);
  free(path_dir);
  free(staple); 
  free(tempmat1); 

  dtime+=dclock();
  final_flop = (Real)nflop*volume/numnodes();
  info->final_sec = dtime; 
  info->final_flop = final_flop;
  info->status = QOP_SUCCESS; 
} /* QOP_symanzik_1loop_gauge_force.c */

