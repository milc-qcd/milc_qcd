/****** quark_stuff.c  -- ******************/
/* MIMD version 6 */
/* quark action stuff for improved action
* D.T. 1/28/98, starting from gauge_stuff.c
* K.O. 3/99 Added optimized fattening for Asq actions
* D.T. 4/99 Combine force calculations for both mass quarks
* K.O. 4/99 Optimized force for Asq action
* S.G. 7/01, modified to use t_longlink and t_fatlink
* C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
* T.B. 11/01, Added d(M)/d(u0) dependencies for equation of state calc's with
*             Asqtad action - ATTN: site structure needs 'dfatlink_du0' in
*             addition to 'fatlink': #define DM_DU0
*
* J.O. 3/04 Rearranged loops for optimization
* J.O. C.D. 3/04 Copied forward links for optimization and 
*                kept mtags open for restart_gather
*                Worked with pointers where possible to avoid copying.
* C.D. 3/05 Moved fermion force and dslash_eo to separate files.

* This code combines quark_stuff.c and quark_stuff_tmp.c

* In this directory, assume all paths connect even to odd sites, etc.
* Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
* as separate parity transforms of the fundamental paths.  They will
* generally need a negative sign in Dslash.  See bottom for a long
* comment on sign conventions.
*/

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fatlinks:       61632 for load_fatlinks
 */

/**#define FFTIME**/
/**#define LLTIME**/

#include "generic_ks_includes.h"	/* definitions files and prototypes */


#define NULL_FP -1 /* NULL field_offset to be used in the optimized version *
                    * of the load_fatlinks subroutine */

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom. */
#include <quark_action.h>
    /* Include file specifies the basic paths */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

void printpath( int *path, int length );

#ifdef DM_DU0

#ifndef TADPOLE_IMPROVE
BOMB THE COMPILE
#endif

void compute_gen_staple(field_offset staple, int mu, int nu,
                        field_offset link, Real coef, Real coef2 ) ;
#else /* not DM_DU0 */
void compute_gen_staple(field_offset staple, int mu, int nu,
                        field_offset link, Real coef ) ;
#endif

int path_num[MAX_BASIC_PATHS];	/* number of rotations/reflections for each 
					kind */
static Real act_path_coeff[MAX_BASIC_PATHS]; /* actual path coefficient     *
                                               * it is equal to path_coeff   *
                                               * if not tadpole improvement  *
                                               * is specified                *
                                               * or path_coeff*u_0^(L-1) when*
                                               * tadpole improvement is      *
                                               * specified                   */
#ifdef DM_DU0
static Real act_path_coeff2[MAX_BASIC_PATHS]; /* coefficient for
						  d(Dslash)/d(u0) */
#endif

/* Array of structures, for each rotation and reflection of each kind of
	path.  */
static Q_path q_paths[MAX_NUM];
static int num_q_paths;	/* number of paths in dslash */
int num_basic_paths;	/* number of paths before rotation/reflection */

int is_path_equal( int *path1, int* path2, int length );
int add_basic_path( int *vec, int length, Real coeff );

#ifdef QSINLINE
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) {\
\
  Real _a0r,_a0i,_a1r,_a1i,_a2r,_a2i;\
  Real _b0r,_b0i,_b1r,_b1i,_b2r,_b2i;\
  \
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
}

#else /* External versions */
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) mult_su3_mat_hwvec( mat, src, dest ) 
#endif

/********************************************************************/
/* Make table of paths in action */
/********************************************************************/
void make_path_table() {

    int i,j;
#ifdef TADPOLE_IMPROVE
    int k;
#endif

    /* table of directions, 1 for each kind of path */
    /**int path_ind[MAX_BASIC_PATHS][MAX_LENGTH];**/
    /* table of coefficients in action, for each path */

    node0_printf("%s\n",quark_action_description);
    num_q_paths = 0;
    num_basic_paths = 0;
    if(MAX_LENGTH > MAX_PATH_LENGTH){
      printf("Path length for this action is too long.  Recompile.\n");
      terminate(1);
    }

    /* add rots. and reflects to table, print out the path coefficients */
    node0_printf("path coefficients: npath  path_coeff  multiplicity\n");
    for(j=0;j<quark_action_npaths;j++) {
	Real this_coeff;
	this_coeff = path_coeff[j];
#ifdef TADPOLE_IMPROVE
	for(k=1;k< path_length_in[j];k++)this_coeff /= u0;
#endif
	act_path_coeff[j] = this_coeff ;
#ifdef DM_DU0
	act_path_coeff2[j] = this_coeff*(1-path_length_in[j])/u0;
#endif
	i = add_basic_path( path_ind[j], path_length_in[j],
	    this_coeff );
	node0_printf("                    %d      %e     %d\n",
	    j,this_coeff,i);
    }
}

/* Accessors for path table */
int get_num_q_paths(){
  return num_q_paths;
}

Q_path *get_q_paths(){
  return q_paths;
}


/* Accessor for quark action information */
Real *get_quark_path_coeff(){
  return act_path_coeff;
}

#ifdef DM_DU0
/* Accessor for quark action information */
Real *get_quark_path_coeff2(){
  return act_path_coeff2;
}
#endif

/********************************************************************/
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */
/********************************************************************/
int add_basic_path( int *basic_vec, int length, Real coeff ) {

    int perm[8],pp[8],ir[4];
    int j,path_num;
    int vec[MAX_LENGTH];
    int flag;

    path_num = 0;
    /* now fill the long table with all rotations and reflections
	of the fundamental path.  The path presented to us is for
        the positive x component of dslash, so if the x coordinate
        is reflected it will appear with a negative sign. */
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

	      if(ir[j] == 1) pp[j]=OPP_DIR(pp[j]);
	      pp[OPP_DIR(j)]=OPP_DIR(pp[j]);
	    }
	    /* create new vector*/
	    for(j=0;j<length;j++) vec[j]=pp[basic_vec[j]];
	    for(j=length;j<MAX_LENGTH;j++) vec[j]=NODIR;

            flag=0;
	    /* check if it's a new set: */
	    for(j=0;j<num_q_paths;j++){
	      flag = is_path_equal( vec, q_paths[j].dir, MAX_LENGTH );
	      if(flag==1)break;
	    }
	    if(flag == 0 ){
	      if(num_q_paths>=MAX_NUM){
		node0_printf("OOPS: MAX_NUM too small\n");
		exit(0);
	      }
	      q_paths[num_q_paths].length=length;
	      for(j=0;j<MAX_LENGTH;j++) q_paths[num_q_paths].dir[j]=vec[j];
		/* remember to copy NODIR's, or comparison will get confused */
	      if(ir[0]==0){
		q_paths[num_q_paths].coeff =  coeff;
#ifdef DM_DU0
		q_paths[num_q_paths].coeff2 = coeff*(1-length)/u0;
#endif
		q_paths[num_q_paths].forwback =  +1;
	      }
	      else{
		q_paths[num_q_paths].coeff = -coeff;
#ifdef DM_DU0
		q_paths[num_q_paths].coeff2 = -coeff*(1-length)/u0;
#endif
		q_paths[num_q_paths].forwback = -1;
	      }
	      num_q_paths++;
	      path_num++;
	      /**node0_printf("ADD PATH %d:  rx=%d ",num_q_paths-1,ir[0]);
		 printpath( vec, length );**/
	    }

	  } /* end reflection*/
        } /* end permutation if block */
      } /* end permutation */
    num_basic_paths++;
    return(path_num);
} /* add_basic_path */

/********************************************************************/
/* parallel transport a vector to the current site along a path.
   For example, if the path is "XUP", bring in the vector from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   vectors rather than stuff in lattice.h */
/********************************************************************/
void path_transport( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    su3_vector *tmp_pt; /* scratch */
    int tmp_parity, tmp_otherparity; /* parity for this step */

  if( length > 0 ){
    tmp_src = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_dest = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_work = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(su3_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_from_temp( tmp_src, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_vec( &(s->link[dir[j]]),
		    (su3_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_vec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_from_temp( tmp_work, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(su3_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
    }
  }
} /* path_transport */

/********************************************************************/
/* Path transport a half_wilson_vector */
/********************************************************************/
void path_transport_hwv( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    half_wilson_vector *tmp_src,*tmp_dest,*tmp_work; /*source, dest and workspace*/
    half_wilson_vector *tmp_pt; /* scratch */
    int tmp_parity, tmp_otherparity; /* parity for this step */

  if( length > 0 ){
    tmp_src = (half_wilson_vector *)malloc(
      sites_on_node*sizeof(half_wilson_vector) );
    tmp_dest = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );
    tmp_work = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );

    for( j=length-1; j>=0; j-- ){
	/* figure out parities for this step */
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case EVEN: tmp_otherparity=ODD; break;
		case ODD: tmp_otherparity=EVEN; break;
		case EVENANDODD: tmp_otherparity=EVENANDODD; break;
	    }
	}
	else { /* odd # step */
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case EVEN: tmp_parity=ODD; break;
		case ODD: tmp_parity=EVEN; break;
		case EVENANDODD: tmp_parity=EVENANDODD; break;
	    }
	}

	if( j==length-1 ){
	    FORSOMEPARITY(i,s,tmp_otherparity){
	        tmp_src[i] = *(half_wilson_vector *)F_PT(s,src);
	    }
	}

	if( GOES_FORWARDS(dir[j]) ) {
	    mtag0 = start_gather_from_temp( tmp_src,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		mult_su3_mat_hwvec_for_inline( &(s->link[dir[j]]),
		    (half_wilson_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{ /* GOES_BACKWARDS(dir[j]) */
	    FORSOMEPARITY(i,s,tmp_otherparity){
		mult_adj_su3_mat_hwvec( &(s->link[OPP_DIR(dir[j])]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_from_temp( tmp_work,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    FORSOMEPARITY(i,s,tmp_parity){
		 tmp_dest[i] = *(half_wilson_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	/* src for next step is dest for this one. */
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }  /* j=link in path */
    /* done, copy result into real dest. (tmp_src now points to result) */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  } /* end if(length>0) */
  else if( src != dest ){ /* for length=0 */
    FORSOMEPARITY(i,s,parity){
        *(half_wilson_vector *)F_PT(s,dest) =
	  *(half_wilson_vector *)F_PT(s,src);
    }
  }
} /* path_transport_hwv */

#ifdef FN
/********************************************************************/
/* Sum over paths connecting to nearest neighbor point (fat link) and to third
   nearest neighbor (longlinks) */
/********************************************************************/
/* Doug Toussaint 2/4/98 */
/* modified to use t_longlinks, S. Gottlieb 7/13/01 */
/* long link calculating routine */
/* path_product() follows the path starting at step 0, and
   leaves the answer at the end of the path.  We want the answer
   at the site where the path begins.  So we look for paths with
   the opposite displacement from the displacement of the point
   that we want to transport to this site, and take the adjoint
   of the matrix at the end. clear? */
/* KS phases and APBC must be in the links. See long comment at bottom*/
void load_longlinks() {
  register int i;
  register site *s;
  int ipath,dir;
  int disp[4];
  int nflop = 1804;
  register su3_matrix *long1;

#ifdef LLTIME
double dtime;
dtime=-dclock();
#endif
  if( phases_in != 1){
    node0_printf("BOTCH: load_longlinks needs phases in\n");
    terminate(0);
  }
  for (dir=XUP; dir<=TUP; dir++){ /* loop over longlink directions */
    /* set longlink to zero */
    FORALLSITES(i,s){
      long1 = &(t_longlink[4*i+dir]);
      clear_su3mat( long1 );
    }

    /* loop over paths, checking for ones with total displacement 3*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=3,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("ipath = %d, found a path:  ",ipath);
for(j=0;j<q_paths[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product( q_paths[ipath].dir, q_paths[ipath].length );
	FORALLSITES(i,s){
	  su3_adjoint( &(s->tempmat1), &(s->staple) );
	  long1 = &(t_longlink[4*i+dir]);
          scalar_mult_add_su3_matrix( long1,
	    &(s->staple), -q_paths[ipath].coeff, long1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */

  } /* loop over directions */

  valid_longlinks = 1;
#ifdef LLTIME
dtime += dclock();
node0_printf("LLTIME(long): time =  %e (Naik) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_longlinks() */

/* KS phases and APBC must be in the links. See long comment at bottom*/
void load_fatlinks() {
  register int i;
  register site *s;
  int dir;
  register su3_matrix *fat1;
#ifdef DM_DU0
  register su3_matrix *fat2;
#endif
#ifdef ASQ_OPTIMIZED_FATTENING
  int  nu,rho,sig ;
  Real one_link,one_link2; /* needed to fix the problem with the Lepage
			       term */
#else
  int ipath;
  int disp[4];
#endif

#ifdef LLTIME
  int nflop = 61632;
double dtime;
dtime=-dclock();
#endif
  if( phases_in != 1){
    node0_printf("BOTCH: load_fatlinks needs phases in\n");
    terminate(0);
  }

#ifndef  ASQ_OPTIMIZED_FATTENING   /* general case code */
  for (dir=XUP; dir<=TUP; dir++){ /* loop over fatlink directions */
    /* set fatlink to zero */
    FORALLSITES(i,s){
      fat1 = &(t_fatlink[4*i+dir]);
      clear_su3mat( fat1 );
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+dir]);
      clear_su3mat( fat2 );
#endif
    }
    
    /* loop over paths, checking for ones with total displacement 1*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=1,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("dir = %d, found a path:  ",dir);
for(j=0;j<q_paths.[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product( q_paths[ipath].dir, q_paths[ipath].length );
	FORALLSITES(i,s){
	  su3_adjoint( &(s->tempmat1), &(s->staple) );
	  fat1 = &(t_fatlink[4*i+dir]);
          scalar_mult_add_su3_matrix( fat1,
	    &(s->staple), -q_paths[ipath].coeff, fat1 );
		/* minus sign in coeff. because we used backward path*/
#ifdef DM_DU0
	  fat2 = &(t_dfatlink_du0[4*i+dir]);
          scalar_mult_add_su3_matrix( fat2,
	    &(s->staple), -q_paths[ipath].coeff2, fat2 );
		/* minus sign in coeff. because we used backward path*/
#endif
	}
    } /* ipath */
  } /* loop over directions */
#else	/* ASQ_OPTIMIZED_FATTENING, for Asq and Asqtad actions */
/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
 
 /* to fix up the Lepage term, included by a trick below */
 one_link = (act_path_coeff[0] - 6.0*act_path_coeff[5]);
#ifdef DM_DU0
 one_link2 = (act_path_coeff2[0] - 6.0*act_path_coeff2[5]);
#endif
 
 for (dir=XUP; dir<=TUP; dir++){
   FORALLSITES(i,s) /* Intialize fat links with c_1*U_\mu(x) */
     {
       fat1 = &(t_fatlink[4*i+dir]);
       scalar_mult_su3_matrix(&(s->link[dir]), one_link ,
			      fat1 );
#ifdef DM_DU0
       fat2 = &(t_dfatlink_du0[4*i+dir]);
       scalar_mult_su3_matrix(&(s->link[dir]), one_link2 ,
			      fat2 );
#endif
     }
   for(nu=XUP; nu<=TUP; nu++) if(nu!=dir)
     {
#ifdef DM_DU0
       compute_gen_staple(F_OFFSET(staple),dir,nu,F_OFFSET(link[dir]),
			  act_path_coeff[2],act_path_coeff2[2]);
       /* The Lepage term */
       /* Note this also involves modifying c_1 (above) */
       compute_gen_staple(NULL_FP,dir,nu,F_OFFSET(staple),act_path_coeff[5],
			  act_path_coeff2[5]);
       for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   compute_gen_staple(F_OFFSET(tempmat1),dir,rho,F_OFFSET(staple),
			      act_path_coeff[3],act_path_coeff2[3]);
	   for(sig=XUP; sig<=TUP; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 compute_gen_staple(NULL_FP,dir,sig,
				    F_OFFSET(tempmat1),
				    act_path_coeff[4],act_path_coeff2[4]);
	       } /* sig */
	 } /* rho */
#else
       compute_gen_staple(F_OFFSET(staple),dir,nu,F_OFFSET(link[dir]),
			  act_path_coeff[2]);
       /* The Lepage term */
       /* Note this also involves modifying c_1 (above) */
       compute_gen_staple(NULL_FP,dir,nu,F_OFFSET(staple),act_path_coeff[5]);
       for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   compute_gen_staple(F_OFFSET(tempmat1),dir,rho,F_OFFSET(staple),
			      act_path_coeff[3]);
	   for(sig=XUP; sig<=TUP; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 compute_gen_staple(NULL_FP,dir,sig,
				    F_OFFSET(tempmat1),
				    act_path_coeff[4]);
	       } /* sig */
	 } /* rho */
#endif
     } /* nu */
 }/* dir */  
#endif

  valid_fatlinks = 1;
#ifdef LLTIME
dtime += dclock();
 node0_printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,
	      (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_fatlinks() */
#endif /* ifdef FN */


/********************************************************************/
/* compare two paths, return 1 if equal, else zero */
/********************************************************************/
int is_path_equal( int *path1, int* path2, int length ){
   register int i;
   for(i=0;i<length;i++)if(path1[i]!=path2[i])return(0);
   return(1);
}


#ifdef  ASQ_OPTIMIZED_FATTENING   /* Asqtad action only, "_fn" executables */
#ifndef FN
BOMB THE COMPILE
#endif
#ifdef DM_DU0
void compute_gen_staple(field_offset staple, int mu, int nu, 
			field_offset link, Real coef, Real coef2) {
#else
void compute_gen_staple(field_offset staple, int mu, int nu, 
			field_offset link, Real coef) {
#endif
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat ;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;
#ifdef DM_DU0
  register su3_matrix *fat2;
#endif

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any su3_matrix. The result is saved in staple.
    if staple==NULL_FP then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather( link, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL_FP){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, (su3_matrix *)F_PT(s,staple) );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &tmat2 );
      fat1 = &(t_fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix(fat1, &tmat2, coef,
				 fat1) ;
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+mu]);
      scalar_mult_add_su3_matrix(fat2, &tmat2, coef2,
				 fat2) ;
#endif
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather( F_OFFSET(link[nu]),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( &(s->link[nu]),(su3_matrix *)F_PT(s,link), &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_from_temp( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL_FP){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( (su3_matrix *)F_PT(s,staple),(su3_matrix *)gen_pt[0][i], 
		      (su3_matrix *)F_PT(s,staple) );
      fat1 = &(t_fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)F_PT(s,staple), coef, 
				 fat1 );
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+mu]);
      scalar_mult_add_su3_matrix( fat2,
				 (su3_matrix *)F_PT(s,staple), coef2, 
				 fat2 );
#endif
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(t_fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+mu]);
      scalar_mult_add_su3_matrix( fat2,
				 (su3_matrix *)gen_pt[0][i], coef2, 
				 fat2 );
#endif
    }
  }

  free(tempmat);
  cleanup_gather(mtag0);
}
#endif  /* ASQ_OPTIMIZED_FATTENING   */
/* quark_stuff.c */
