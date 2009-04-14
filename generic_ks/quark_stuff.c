/****** quark_stuff.c  -- ******************/
/* MIMD version 7 */
/* Construct path tables and action coefficients for improved quark actions
 * (Formerly a catch-all file for improved quark action utilities)
 *
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
 * C.D. 10/06 Moved compute_gen_staple to fermion_links_fn.c
 
 * This code combines quark_stuff.c and quark_stuff_tmp.c
 
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
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
#ifndef TADPOLE_IMPROVE
BOMB THE COMPILE
#endif

static Real act_path_coeff_dmdu0[MAX_BASIC_PATHS]; /* coefficient for
						      d(Dslash)/d(u0) */
#endif

/* Array of structures, for each rotation and reflection of each kind of
	path.  */
static Q_path q_paths[MAX_NUM];
#ifdef DM_DU0
static Q_path q_paths_dmdu0[MAX_NUM];
#endif
static int num_q_paths;	/* number of paths in dslash */
static int num_basic_paths;	/* number of paths before rotation/reflection */

static int 
is_path_equal( int *path1, int* path2, int length );
static int 
add_basic_path( int *vec, int length, Real coeff );

/********************************************************************/
/* Make table of paths in action */
/********************************************************************/
int make_path_table(ks_action_paths *ap, ks_action_paths *ap_dmdu0) {

    int i,j;
#ifdef TADPOLE_IMPROVE
    int k;
#endif

    // Don't remake the table if already made

    if(ap->constructed)return 0;
    // node0_printf("MAKING PATH TABLES\n");
  
    /* table of directions, 1 for each kind of path */
    /**int path_ind[MAX_BASIC_PATHS][MAX_LENGTH];**/
    /* table of coefficients in action, for each path */

    node0_printf("%s\n",QUARK_ACTION_DESCRIPTION);
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
	act_path_coeff_dmdu0[j] = this_coeff*(1-path_length_in[j])/u0;
#endif
	i = add_basic_path( path_ind[j], path_length_in[j],
	    this_coeff );
	node0_printf("                    %d      %e     %d\n",
	    j,this_coeff,i);
    }
    ap->num_q_paths = num_q_paths;
    ap->q_paths = q_paths;
    ap->act_path_coeff = act_path_coeff;
#ifdef DM_DU0
    ap_dmdu0->num_q_paths = num_q_paths;
    ap_dmdu0->q_paths = q_paths;
    ap_dmdu0->act_path_coeff = act_path_coeff_dmdu0;
#endif
    ap->constructed = 1;
    return 1;
}

/********************************************************************/
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */
/********************************************************************/
static int 
add_basic_path( int *basic_vec, int length, Real coeff ) {

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
#ifdef DM_DU0
	      q_paths_dmdu0[num_q_paths].length=length;
#endif
	      for(j=0;j<MAX_LENGTH;j++) {
		q_paths[num_q_paths].dir[j]=vec[j];
#ifdef DM_DU0
		q_paths_dmdu0[num_q_paths].dir[j]=vec[j];
#endif
	      }
		/* remember to copy NODIR's, or comparison will get confused */
	      if(ir[0]==0){
		q_paths[num_q_paths].coeff =  coeff;
		q_paths[num_q_paths].forwback =  +1;
#ifdef DM_DU0
		q_paths_dmdu0[num_q_paths].coeff = coeff*(1-length)/u0;
		q_paths_dmdu0[num_q_paths].forwback =  +1;
#endif
	      }
	      else{
		q_paths[num_q_paths].coeff = -coeff;
		q_paths[num_q_paths].forwback = -1;
#ifdef DM_DU0
		q_paths_dmdu0[num_q_paths].coeff = -coeff*(1-length)/u0;
		q_paths_dmdu0[num_q_paths].forwback = -1;
#endif
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
static int 
is_path_equal( int *path1, int* path2, int length ){
   register int i;
   for(i=0;i<length;i++)if(path1[i]!=path2[i])return(0);
   return(1);
}

/********************************************************************/
/* Initialization */
/********************************************************************/
void 
init_path_table(ks_action_paths *ap){
  ap->constructed = 0;
}

/* quark_stuff.c */
