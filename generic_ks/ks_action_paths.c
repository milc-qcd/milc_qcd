/****** ks_action_paths.c  -- ******************/
/* (Formerly called quark_stuff.c) */
/* (Formerly a catch-all file for improved quark action utilities) */
/* MIMD version 7 */
/* Construct path tables and action coefficients for improved quark actions
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


#include <stdlib.h>
#include <stdio.h>
#include "../include/comdefs.h"
#define IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#include <quark_action.h>
#include "../include/ks_action_paths.h"
#include "lattice.h"   /* Only for u0 */

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash_site().  Rotations
       and reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at end of
       fermion_forcee_eo_milc.c. */
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
#  ifndef TADPOLE_IMPROVE
BOMB THE COMPILE
#  endif
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

static void 
load_asqtad_coeffs_t(asqtad_coeffs_t *apc, Real act_path_coeff[]){
#ifndef ANISOTROPY
    apc->one_link =     act_path_coeff[0];
    apc->naik =	        act_path_coeff[1];
    apc->three_staple = act_path_coeff[2];
    apc->five_staple =  act_path_coeff[3];
    if(MAX_BASIC_PATHS > 4)
      apc->seven_staple = act_path_coeff[4];
    else
      apc->seven_staple = 0.;
    if(MAX_BASIC_PATHS > 5)
      apc->lepage =       act_path_coeff[5];
    else
      apc->lepage = 0.;
#else
#ifndef ABSORB_ANI_XIQ
    // Minimal set of coefficients for a minimal setup without absorbing anisotropy factors. JHW, 2020/10/05
    apc->one_link         = act_path_coeff[ 0];
    apc->naik             = act_path_coeff[ 1];
    apc->three_staple     = act_path_coeff[ 2];
    apc->five_staple      = act_path_coeff[ 3];
    apc->seven_staple     = act_path_coeff[ 4];
    apc->lepage           = act_path_coeff[ 5];
    apc->ani_naik         = act_path_coeff[ 6];
    apc->ani_three_staple = act_path_coeff[ 7];
    apc->ani_five_staple  = act_path_coeff[ 8];
    apc->ani_seven_staple = act_path_coeff[ 9];
    apc->ani_lepage       = act_path_coeff[10];
#else
    // Larger set of coefficients in a complete setup that absorbs the anisotropy factors. JHW, 2020/10/05
    apc->one_link         = act_path_coeff[ 0];
    apc->naik             = act_path_coeff[ 1];
    apc->three_staple     = act_path_coeff[ 2];
    apc->five_staple      = act_path_coeff[ 3];
    apc->lepage           = act_path_coeff[ 4];
    apc->ani_one_link     = act_path_coeff[ 5];
    apc->ani_three_staple = act_path_coeff[ 6];
    apc->ani_five_staple  = act_path_coeff[ 7];
    apc->ani_lepage       = act_path_coeff[ 8];
    apc->ani_seven_staple = act_path_coeff[ 9];
    apc->ani2_three_staple= act_path_coeff[10];
    apc->ani2_five_staple = act_path_coeff[11];
    apc->seven_staple     = act_path_coeff[12]; // Necessarily contains two anisotropic links, JHW 2020/10/04
    apc->ani_naik         = act_path_coeff[13];
    apc->ani4_lepage      = act_path_coeff[14];

#endif
    apc->ani_dir = ani_dir;
#endif
}

/********************************************************************/
/* Make table of paths in action */
/********************************************************************/
int make_path_table(ks_action_paths *ap, ks_action_paths *ap_dmdu0) {

    int i,j;
#ifdef TADPOLE_IMPROVE
    int k;
#endif

    if(ap == NULL){
      printf("make_path_table(%d): Called with NULL ap pointer\n",mynode());
      terminate(1);
    }

#ifdef DM_DU0
    if(ap_dmdu0 == NULL){
      printf("make_path_table(%d): Called with NULL ap_dmdu0 pointer\n",mynode());
      terminate(1);
    }
    if(ap->constructed && ap_dmdu0->constructed)return 0;
#else
    if(ap->constructed)return 0;
#endif

    /* table of directions, 1 for each kind of path */
    /**int path_ind[MAX_BASIC_PATHS][MAX_LENGTH];**/
    /* table of coefficients in action, for each path */

#ifdef ANISOTROPY
     Real u0r1 = u0/ani_u0;
     Real u0r2 = u0r1*u0r1;
     Real u0r4 = u0r2*u0r2;
     node0_printf("Ratio of tadpole factors %f \n",u0r1);
#  ifndef ABSORB_ANI_XIQ
#    define IS_ANI_PATH(j)   ( j >= ANI_NK ) 
#    define IS_ANI_LEPAGE(j) ( j == ANI_LP ) 
#  else
#    define IS_ANI0_PATH(j) ( j <  ANI1_1L ) 
#    define IS_ANI1_PATH(j) ( j >  ANI0_LP && j <  ANI2_3L ) 
#    define IS_ANI2_PATH(j) ( j >  ANI1_LP && j <  ANI3_NK ) 
#    define IS_ANI3_PATH(j) ( j == ANI3_NK ) 
#    define IS_ANI4_PATH(j) ( j == ANI4_LP ) 
#  endif
#endif

    if(mynode()==0)printf("%s\n",QUARK_ACTION_DESCRIPTION);
    num_q_paths = 0;
    num_basic_paths = 0;
    if(MAX_LENGTH > MAX_PATH_LENGTH){
      printf("Path length for this action is too long.  Recompile.\n");
      terminate(1);
    }

    /* add rots. and reflects to table, print out the path coefficients */
    if(mynode()==0)printf("path coefficients: npath  path_coeff  multiplicity\n");
    for(j=0;j<quark_action_npaths;j++) {
	Real this_coeff;
	this_coeff = path_coeff[j];
#ifdef TADPOLE_IMPROVE
	for(k=1;k< path_length_in[j];k++)this_coeff /= u0;
#endif
#  ifndef ANISOTROPY
	act_path_coeff[j] = this_coeff ;
#  else
#    ifndef ABSORB_ANI_XIQ
#ifdef TADPOLE_IMPROVE
        if ( IS_ANI_PATH(j) ) this_coeff *= (IS_ANI_LEPAGE(j)? u0r2 : u0r4 );
        act_path_coeff[j] = this_coeff ; 
#endif
#    else
#ifdef TADPOLE_IMPROVE
	for(k=0;k < path_u0rat_pow[j];k++) this_coeff *= u0r1;
#endif

#  ifdef ONEDIM_ANISO_TEST
        act_path_coeff[j] = this_coeff * ( IS_ANI1_PATH(j) || IS_ANI3_PATH(j) ? ani_xiq : iso_xiq );
#  else
        act_path_coeff[j] = this_coeff * ( IS_ANI1_PATH(j) || IS_ANI3_PATH(j) ? ani_xiq : 1 );
#  endif

#    endif
#  endif

#ifdef DM_DU0
#  ifndef ANISOTROPY
	act_path_coeff_dmdu0[j] = this_coeff*(1-path_length_in[j])/u0;
#  else
#    ifndef ABSORB_ANI_XIQ
	act_path_coeff_dmdu0[j] = this_coeff*(1+(IS_ANI_PATH(j)+IS_ANI_LEPAGE(j))*2-path_length_in[j])/u0;
#    else
        switch ( path_u0_rat_pow[j] ) {
          case 0: 
          case 2: 
            act_path_coeff_dmdu0[j] = act_path_coeff[j] * (1+path_u0_rat_pow[j]-path_length_in[j])/u0 ;
            break;
            break;
          case 4: // No u0 dependence like the isotropic one_link, JHW 2020/10/03
          default:
            act_path_coeff_dmdu0[j] = 0;
            break;
        }
#    endif
#  endif
#endif
	i = add_basic_path( path_ind[j], path_length_in[j], this_coeff );
	if(mynode()==0)printf("                    %d      %e     %d\n",
	    j,this_coeff,i);
    }
    ap->p.num_q_paths = num_q_paths;
    ap->p.q_paths = q_paths;

    load_asqtad_coeffs_t(&ap->p.act_path_coeff, act_path_coeff);

#ifdef DM_DU0
    ap_dmdu0->p.num_q_paths = num_q_paths;
    ap_dmdu0->p.q_paths = q_paths_dmdu0;
    load_asqtad_coeffs_t(&ap_dmdu0->p.act_path_coeff, act_path_coeff_dmdu0);
    ap_dmdu0->constructed = 1;
#endif
    ap->constructed = 1;
#ifdef ANISOTROPY
  char ani_char[]="xyzt";
  ap->ani_dir = ani_dir; 
#  ifndef ABSORB_ANI_XIQ
  node0_printf("Anisotropic action with bare quark anisotropy %.6f in the %c-direction\n", ani_xiq, ani_char[ani_dir]);
  ap->ani_xiq = ani_xiq; 
#  else
  node0_printf("Anisotropic action with bare quark anisotropy %.6f in the %c-direction absorbed into path coefficients\n", ani_xiq, ani_char[ani_dir]);
  ap->ani_xiq = 1; 
#  endif
#  ifdef ONEDIM_ANISO_TEST
  node0_printf("using three isotropic directions with factor %.6f for debugging\n",iso_xiq);
#    ifndef ABSORB_ANI_XIQ
  ap->iso_xiq = iso_xiq; 
#    else
  ap->iso_xiq = 1; 
#    endif
#  endif
#endif
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
#ifdef ANISOTROPY
            /* Distinguish between paths with multiple anisotropic links or without and 
             * extend path table with that criterium. JHW 2020/10/02 */
            int ani_len = 0;
	    for(j=0;j<length;j++) ani_len += ( vec[j] == ani_dir || vec[j] == OPP_DIR(ani_dir) );
#  ifndef ABSORB_ANI_XIQ
            if ( ( ani_len > 1 && num_q_paths < ISO_NUM ) || ( ani_len <= 1 && num_q_paths >= ISO_NUM ) ) continue;
#  else
            /* Did not find any less ugly solution than giving minima and maxima for each combination of 
             * length and length along the anisotropic direction (ani_len). Without separate limits for 
             * lengths would not be able to discriminate between Lepage with one or without anisotropic 
             * links. JHW 2020/10/05 */
            if (  ( ( ani_len == 0 && length == 1 && (                              num_q_paths > ANI0_1L_MAX ) )
                 || ( ani_len == 0 && length == 3 && ( num_q_paths < ANI0_1L_MAX || num_q_paths > ANI0_3L_MAX ) )
                 || ( ani_len == 0 && length == 5 && ( num_q_paths < ANI0_3L_MAX || num_q_paths > ANI0_LP_MAX ) )
                 || ( ani_len == 1 && length == 1 && ( num_q_paths < ANI0_LP_MAX || num_q_paths > ANI1_1L_MAX ) )
                 || ( ani_len == 1 && length == 3 && ( num_q_paths < ANI1_1L_MAX || num_q_paths > ANI1_3L_MAX ) )
                 || ( ani_len == 1 && length == 5 && ( num_q_paths < ANI1_3L_MAX || num_q_paths > ANI1_LP_MAX ) )
                 || ( ani_len == 1 && length == 7 && ( num_q_paths < ANI1_5L_MAX || num_q_paths > ANI1_7L_MAX ) )
                 || ( ani_len == 2 && length == 3 && ( num_q_paths < ANI1_LP_MAX || num_q_paths > ANI2_3L_MAX ) )
                 || ( ani_len == 2 && length == 5 && ( num_q_paths < ANI2_3L_MAX || num_q_paths > ANI2_5L_MAX ) )
                 || ( ani_len == 2 && length == 7 && ( num_q_paths < ANI2_5L_MAX || num_q_paths > ANI2_7L_MAX ) )
                 || ( ani_len == 3 && length == 3 && ( num_q_paths < ANI2_7L_MAX || num_q_paths > ANI3_NK_MAX ) )
                 || ( ani_len == 4 && length == 5 && ( num_q_paths < ANI3_NK_MAX                           ) )
               )  )  continue;
#  endif
#endif
//printf("PATH (%d,%d -- %d,%d): ",path_num,num_q_paths,length,ani_len); for(j=0;j<length;j++) printf(" %d ", vec[j]); printf("\n");
	    /* check if it's a new set: */
	    for(j=0;j<num_q_paths;j++){
	      flag = is_path_equal( vec, q_paths[j].dir, MAX_LENGTH );
	      if(flag==1)break;
	    }
	    if(flag == 0 ){
	      if(num_q_paths>=MAX_NUM){
		if(mynode()==0)printf("OOPS: MAX_NUM too small\n");
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
	      /**if(mynode()==0)printf("ADD PATH %d:  rx=%d ",num_q_paths-1,ir[0]);
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
//void 
//init_path_table(ks_action_paths *ap){
//  ap->constructed = 0;
//}


ks_action_paths *
create_path_table(void){
  ks_action_paths *ap;
  char myname[] = "create_path_table";

  ap = (ks_action_paths *)malloc(sizeof(ks_action_paths));
  if(ap == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  ap->constructed = 0;

  return ap;

}

void destroy_path_table(ks_action_paths *ap){
  if(ap == NULL)return;
  free(ap);
}

#include <string.h>

#define MAX_STRING 2048

char *
get_ap_string(ks_action_paths *ap){
  static char str[MAX_STRING] = "";
#ifdef ASQ_ACTION
  asqtad_coeffs_t *apc = &ap->p.act_path_coeff;

  snprintf(str, MAX_STRING, "action.asqtad.one_link %e\naction.asqtad.three_staple %e\naction.asqtad.five_staple %e\naction.asqtad.seven_staple %e\naction.asqtad.lepage %e\naction.asqtad.naik %e\n",
	   apc->one_link, apc->three_staple, apc->five_staple, apc->seven_staple, 
	   apc->lepage, apc->naik);
#else

  int j;
  char *buf;

  for(j = 0, buf = str; j < quark_action_npaths; j++, buf = str + strlen(str)) {
     snprintf(buf, MAX_STRING-strlen(str), "action.fn.coeff[%d] %e\n", j, act_path_coeff[j]);
  }

#endif

  str[MAX_STRING-1] = '\0';
  return str;
}

/* ks_action_paths.c */
