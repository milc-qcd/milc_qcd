/****** quark_stuff.c  -- ******************/
/* MIMD version 7 */
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
*                kept mtags open for restart_gather_site
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

#include "generic_ks_includes.h"	/* definitions files and prototypes */


    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash_site().  Rotations
       and reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at end of
       fermion_force_general.c. */
#include <quark_action.h>
    /* Include file specifies the basic paths */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

void printpath( int *path, int length );

#ifdef DM_DU0

#ifndef TADPOLE_IMPROVE
BOMB THE COMPILE
#endif

void compute_gen_staple_site(su3_matrix *staple, int mu, int nu,
			     field_offset link, Real coef, Real coef2 ) ;
void compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			      su3_matrix *link, Real coef, Real coef2);
#else /* not DM_DU0 */
void compute_gen_staple_site(su3_matrix *staple, int mu, int nu,
			     field_offset link, Real coef ) ;
void compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			      su3_matrix *link, Real coef);
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
void compute_gen_staple_site(su3_matrix *staple, int mu, int nu, 
			field_offset link, Real coef, Real coef2) {
#else
void compute_gen_staple_site(su3_matrix *staple, int mu, int nu, 
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
    if staple==NULL then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_site( link, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &staple[i] );
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
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( &(s->link[nu]),(su3_matrix *)F_PT(s,link), &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( &staple[i],(su3_matrix *)gen_pt[0][i], 
		      &staple[i] );
      fat1 = &(t_fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 &staple[i], coef, 
				 fat1 );
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+mu]);
      scalar_mult_add_su3_matrix( fat2,
				 &staple[i], coef2, 
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
} /* compute_gen_staple_site */
#endif  /* ASQ_OPTIMIZED_FATTENING   */

#ifdef  ASQ_OPTIMIZED_FATTENING   /* Asqtad action only, "_fn" executables */
#ifndef FN
BOMB THE COMPILE
#endif
#ifdef DM_DU0
void compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			su3_matrix *link, Real coef, Real coef2) {
#else
void compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			su3_matrix *link, Real coef) {
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
    if staple==NULL then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_field( link, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &staple[i] );
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
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( &(s->link[nu]),&link[i], &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( &staple[i],(su3_matrix *)gen_pt[0][i], 
		      &staple[i] );
      fat1 = &(t_fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 &staple[i], coef, 
				 fat1 );
#ifdef DM_DU0
      fat2 = &(t_dfatlink_du0[4*i+mu]);
      scalar_mult_add_su3_matrix( fat2,
				 &staple[i], coef2, 
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
} /* compute_gen_staple_field */
#endif  /* ASQ_OPTIMIZED_FATTENING   */

/* quark_stuff.c */
