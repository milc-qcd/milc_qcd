/****** ks_action_paths_hisq.c  -- ******************/
/* (Formerly quark_stuff_hisq.c) */
/* MIMD version 7 */
/* Construct path tables and action coefficients for improved quark actions
 * (Formerly a catch-all file for improved quark action utilities)
 *
 * D.T. 1/28/98, starting from gauge_stuff.c
 * K.O. 3/99 Added optimized fattening for Asq actions
 * D.T. 4/99 Combine force calculations for both mass quarks
 * K.O. 4/99 Optimized force for Asq action
 * S.G. 7/01, modified to use Xt_longlink and Xt_fatlink
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
 * A.B. 8/08 Third path table introduced for HISQ to speed up cases
 *           with non-zero epsilon correction to Naik term
 
 * This code combines quark_stuff.c and quark_stuff_tmp.c
 
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
*/


#include <stdlib.h>
#include <stdio.h>
#include "../include/comdefs.h"
#include "../include/complex.h"
#define IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#include "quark_action.h"
#include "../include/ks_action_paths.h"


    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash_site().  Rotations
       and reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at end of
       fermion_forcee_eo_milc.c. */
    /* Include file specifies the basic paths */

static int 
make_component_path_table( char *action_desc, int npaths, int max_paths,
		      int *path_length, Real *coeff,
		      int paths[][MAX_LENGTH], Real *act_coeff, 
		      Q_path *this_q_paths, Real mass, int index_onelink,
		      int index_naik );

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

/* Actual path coefficients     *
 * Without tadpole improvement
 * they are equal to path_coeff   *
 * Otherwise they are  path_coeff*u_0^(L-1) */

static Real act_path_coeff_1[NUM_BASIC_PATHS_1];
static Real act_path_coeff_2[NUM_BASIC_PATHS_2];
static Real act_path_coeff_3[NUM_BASIC_PATHS_3];

/* Array of structures, for each rotation and reflection of each kind of
	path.  */
static Q_path q_paths_1[MAX_NUM_1];
static Q_path q_paths_2[MAX_NUM_2];
static Q_path q_paths_3[MAX_NUM_3];
static int num_q_paths_1;	/* number of paths in dslash */
static int num_q_paths_2;	/* number of paths in dslash */
static int num_q_paths_3;	/* number of paths in dslash */

static int 
is_path_equal( int *path1, int* path2, int length );
static int 
add_basic_path( Q_path *q_paths, int start_path, int *vec, int length, 
		Real coeff, int max_paths );

static void 
build_act_path_coeff(char action_desc[], Real act_coeff[], Real coeff[], 
		     int npaths, int path_length[]){
  int j;

  if(mynode()==0)printf("%s\n",action_desc);
  if(mynode()==0)printf("path coefficients: npath  path_coeff\n");

  for(j=0;j<npaths;j++) {
    Real this_coeff;
    this_coeff = coeff[j];
#ifdef TADPOLE_IMPROVE
    {
      int k;
      for(k=1;k< path_length[j];k++)this_coeff /= u0;
    }
#endif
    act_coeff[j] = this_coeff ;
    if(mynode()==0)printf("                    %d      %e\n", j,this_coeff);
  }
}

/********************************************************************/
/* Get just the path coefficients in the action                     */
/* but don't construct the path tables                              */
/********************************************************************/
void load_act_path_coeff_hisq(ks_action_paths_hisq *ap, int n_naiks,
			      double *eps_naik){
  int i;

  ap->n_naiks = n_naiks;
  for(i = 0; i < n_naiks; i++)
    ap->eps_naik[i] = eps_naik[i];

  build_act_path_coeff(QUARK_ACTION_DESCRIPTION_1, act_path_coeff_1, 
		       path_coeff_1, quark_action_npaths_1, path_length_in_1);
  ap->p1.num_q_paths = 0;
  ap->p1.q_paths = NULL;

  // fat7 stage 
  // convert hisq style coeffs to asqtad style 
  ap->p1.act_path_coeff.one_link     = act_path_coeff_1[0];
  ap->p1.act_path_coeff.three_staple = act_path_coeff_1[1];
  ap->p1.act_path_coeff.five_staple  = act_path_coeff_1[2];
  ap->p1.act_path_coeff.seven_staple = act_path_coeff_1[3];

  // by definition lepage and naik coeffs are 0 for fat7
  ap->p1.act_path_coeff.lepage = 0.;
  ap->p1.act_path_coeff.naik = 0.;

  //  ap->p1.act_path_coeff = act_path_coeff_1;
    
#if ( UNITARIZATION_METHOD==UNITARIZE_NONE )
  if(mynode()==0)printf("Unitarization method = UNITARIZE_NONE\n");
#elif ( UNITARIZATION_METHOD==UNITARIZE_APE )
  if(mynode()==0)printf("UNITARIZE_APE: derivative is not ready for this method\n");
#elif ( UNITARIZATION_METHOD==UNITARIZE_ROOT )
  if(mynode()==0)printf("Unitarization method = UNITARIZE_ROOT\n");
#elif ( UNITARIZATION_METHOD==UNITARIZE_RATIONAL )
  if(mynode()==0)printf("Unitarization method = UNITARIZE_RATIONAL\n");
#elif ( UNITARIZATION_METHOD==UNITARIZE_ANALYTIC )
  if(mynode()==0)printf("Unitarization method = UNITARIZE_ANALYTIC\n");
#elif ( UNITARIZATION_METHOD==UNITARIZE_STOUT )
  if(mynode()==0)printf("Unitarization method = UNITARIZE_STOUT\n");
#else
  if(mynode()==0)printf("Unknown unitarization method\n"); terminate(0);
#endif
  
#if ( UNITARIZATION_GROUP==UNITARIZE_SU3 )
  if(mynode()==0)printf("Unitarizaton group = SU(3)\n");
#elif ( UNITARIZATION_GROUP==UNITARIZE_U3 )
  if(mynode()==0)printf("Unitarizaton group = U(3)\n");
#else
  if(mynode()==0)printf("Unknown unitarization group. Set U(3) or SU(3)\n");
  terminate(0);
#endif

  ap->umethod = UNITARIZATION_METHOD;
  ap->ugroup = UNITARIZATION_GROUP;

  build_act_path_coeff(QUARK_ACTION_DESCRIPTION_2, act_path_coeff_2, 
		       path_coeff_2, quark_action_npaths_2, path_length_in_2);
  ap->p2.num_q_paths = 0;
  ap->p2.q_paths = NULL;

  ap->p2.act_path_coeff.one_link     = act_path_coeff_2[0];
  ap->p2.act_path_coeff.naik         = act_path_coeff_2[1];
  ap->p2.act_path_coeff.three_staple = act_path_coeff_2[2];
  ap->p2.act_path_coeff.five_staple  = act_path_coeff_2[3];
  ap->p2.act_path_coeff.seven_staple = act_path_coeff_2[4];
  ap->p2.act_path_coeff.lepage       = act_path_coeff_2[5];

  // ap->p2.act_path_coeff = act_path_coeff_2;

  build_act_path_coeff(QUARK_ACTION_DESCRIPTION_3, act_path_coeff_3, 
		       path_coeff_3, quark_action_npaths_3, path_length_in_3);
  ap->p3.num_q_paths = 0;
  ap->p3.q_paths = NULL;

  ap->p3.act_path_coeff.one_link     = act_path_coeff_3[0];
  ap->p3.act_path_coeff.naik         = act_path_coeff_3[1];
  ap->p3.act_path_coeff.three_staple = 0.;
  ap->p3.act_path_coeff.five_staple  = 0.;
  ap->p3.act_path_coeff.seven_staple = 0.;
  ap->p3.act_path_coeff.lepage       = 0.;


  // ap->p3.act_path_coeff = act_path_coeff_3;

}

/********************************************************************/
/* Make table of paths in action */
/********************************************************************/
// Returns 0 if the path table did not change. 1 if it did.
int make_path_table_hisq(ks_action_paths_hisq *ap,
			 int n_naiks, double *eps_naik) {

  if(ap->constructed) return 0;

  load_act_path_coeff_hisq(ap, n_naiks, eps_naik);

  if(mynode()==0)printf("MAKING PATH TABLES\n");

  /* Keep a copy of the Naik epsilon values */
  num_q_paths_1 = 
      make_component_path_table( QUARK_ACTION_DESCRIPTION_1, quark_action_npaths_1,
	    MAX_NUM_1, path_length_in_1, path_coeff_1, path_ind_1, 
	    act_path_coeff_1, q_paths_1, 0.0, -1, -1 );
    
  ap->p1.num_q_paths = num_q_paths_1;
  ap->p1.q_paths = q_paths_1;
    
  num_q_paths_2 = 
    make_component_path_table( QUARK_ACTION_DESCRIPTION_2, quark_action_npaths_2,
			  MAX_NUM_2, 
			  path_length_in_2, path_coeff_2, path_ind_2, 
			  act_path_coeff_2, q_paths_2, 0.0, -1, -1 );

  ap->p2.num_q_paths = num_q_paths_2;
  ap->p2.q_paths = q_paths_2;
  //  ap->p2.act_path_coeff = act_path_coeff_2;

  num_q_paths_3 = 
    make_component_path_table( QUARK_ACTION_DESCRIPTION_3, quark_action_npaths_3,
                        MAX_NUM_3, 
                        path_length_in_3, path_coeff_3, path_ind_3, 
                        act_path_coeff_3, q_paths_3, 0.0, -1, -1 );

  ap->p3.num_q_paths = num_q_paths_3;
  ap->p3.q_paths = q_paths_3;
  //  ap->p3.act_path_coeff = act_path_coeff_3;

  ap->constructed = 1;
  return 1;
}

static int 
make_component_path_table( char *action_desc, int npaths, int max_paths,
		      int *path_length, Real *coeff,
		      int paths[][MAX_LENGTH], Real *act_coeff, 
		      Q_path *this_q_paths, Real naik_term_epsilon, 
		      int index_onelink, int index_naik ) 
{

  int i,j;
  int n_basic_paths;	 // number of paths before rotation/reflection
  int n_q_paths; // total number of paths in table
#ifdef TADPOLE_IMPROVE
  int k;
#endif
  
  /* table of directions, 1 for each kind of path */
  /**int paths[npaths][MAX_LENGTH];**/
  /* table of coefficients in action, for each path */
  
  //  if(mynode()==0)printf("%s\n",action_desc);
  n_q_paths = 0;
  n_basic_paths = 0;
  if(MAX_LENGTH > MAX_PATH_LENGTH){
    printf("Path length for this action is too long.  Recompile.\n");
    terminate(1);
  }
  
  /* add rots. and reflects to table, print out the path coefficients */
  //  if(mynode()==0)printf("path coefficients: npath  path_coeff  multiplicity\n");
  for(j=0;j<npaths;j++) {
    Real this_coeff;
//    this_coeff = coeff[j];
//#ifdef TADPOLE_IMPROVE
//    for(k=1;k< path_length[j];k++)this_coeff /= u0;
//#endif
//    act_coeff[j] = this_coeff ;
    this_coeff = act_coeff[j] ;
    //AB THIS SHOULD BE REMOVED
    // Apply mass correction to one-link and Naik coefficients
    if(j == index_onelink){
      ; //this_coeff += onelink_mass_renorm_fact * naik_term_epsilon * naik_term_epsilon;
    }
    if(j == index_naik){
      ; //this_coeff += naik_mass_renorm_fact * naik_term_epsilon * naik_term_epsilon;
    }
    i = add_basic_path( this_q_paths, n_q_paths, paths[j],
			path_length[j], this_coeff, max_paths );
    n_q_paths += i;
    n_basic_paths++;
    //    if(mynode()==0)printf("                    %d      %e     %d\n", j,this_coeff,i);
  }
  return( n_q_paths );
} //make_component_path_table()

/********************************************************************/
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */
/********************************************************************/
static int 
add_basic_path( Q_path *this_q_paths, int path_table_index, 
		int *basic_vec, int length, Real coeff, int max_paths ) {
  // this_q_paths is array of paths we are building
    // path_table_index is starting index when called
    // basic_vec is list of directions in basic path
    // length is length of basic path
    // coeff is coefficient in action

    int perm[8],pp[8],ir[4];
    int j,path_num;
    int vec[MAX_LENGTH];
    int flag;
         //if(mynode()==0)printf("ADD BASIC PATH %d:  ",path_table_index);
	 //printpath( basic_vec, length );

    path_num = 0;  // number of paths made from this basic path so far
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
	    for(j=0;j<path_table_index;j++){
	      flag = is_path_equal( vec, this_q_paths[j].dir, MAX_LENGTH );
	      if(flag==1)break;
	    }
	    if(flag == 0 ){
	      if(path_table_index>=max_paths){
		if(mynode()==0)printf("OOPS: MAX_NUM too small\n");
		exit(0);
	      }
	      this_q_paths[path_table_index].length=length;
	      for(j=0;j<MAX_LENGTH;j++) {
		this_q_paths[path_table_index].dir[j]=vec[j];
	      }
		/* remember to copy NODIR's, or comparison will get confused */
	      if(ir[0]==0){
		this_q_paths[path_table_index].coeff =  coeff;
		this_q_paths[path_table_index].forwback =  +1;
	      }
	      else{
		this_q_paths[path_table_index].coeff = -coeff;
		this_q_paths[path_table_index].forwback = -1;
	      }
	      path_table_index++;
	      path_num++;
	         //if(mynode()==0)printf("ADD PATH %d:  rx=%d ",path_table_index-1,ir[0]);
		 //printpath( vec, length );
	    }

	  } /* end reflection*/
        } /* end permutation if block */
      } /* end permutation */
	//if(mynode()==0)printf("ADD BASIC PATH: added %d entries\n",path_num);
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
/* API */
/********************************************************************/

ks_action_paths_hisq *
create_path_table_hisq(void){
  ks_action_paths_hisq *ap;
  char myname[] = "create_path_table_hisq";

  ap = (ks_action_paths_hisq *)malloc(sizeof(ks_action_paths_hisq));
  if(ap == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  ap->constructed = 0;

  return ap;

}

void destroy_path_table_hisq(ks_action_paths_hisq *ap){
  if(ap == NULL)return;
  free(ap);
}

int get_n_naiks(ks_action_paths_hisq *ap){
  if(ap == NULL)return 0;
  return ap->n_naiks;
}

double *get_eps_naik(ks_action_paths_hisq *ap){
  if(ap == NULL)return NULL;
  return ap->eps_naik;
}

int get_umethod(ks_action_paths_hisq *ap){
  if(ap == NULL)return 0;
  return ap->umethod;
}

int get_ugroup(ks_action_paths_hisq *ap){
  if(ap == NULL)return 0;
  return ap->ugroup;
}

#define MAX_STRING 2048

char *
get_ap_string_hisq(ks_action_paths_hisq *ap){
  static char str[MAX_STRING] = "";
  asqtad_coeffs_t *apc1 = &ap->p1.act_path_coeff;
  asqtad_coeffs_t *apc2 = &ap->p2.act_path_coeff;
  asqtad_coeffs_t *apc3 = &ap->p3.act_path_coeff;
  
  snprintf(str, MAX_STRING, "action.hisq.fat7.one_link %e\naction.hisq.fat7.three_staple %e\naction.hisq.fat7.five_staple %e\naction.hisq.fat7.seven_staple %e\naction.hisq.asqtad.one_link %e\naction.hisq.asqtad.three_staple %e\naction.hisq.asqtad.five_staple %e\naction.hisq.asqtad.seven_staple %e\naction.hisq.asqtad.lepage %e\naction.hisq.asqtad.naik %e\naction.hisq.difference.one_link %e\naction.hisq.difference.naik %e\n",
	   apc1->one_link, apc1->three_staple, apc1->five_staple, apc1->seven_staple, 
	   apc2->one_link, apc2->three_staple, apc2->five_staple, apc2->seven_staple, 
	   apc2->lepage,   apc2->naik, 
	   apc3->one_link, apc3->naik);

  str[MAX_STRING-1] = '\0';
  return str;
}

/* ks_action_paths_hisq.c */
